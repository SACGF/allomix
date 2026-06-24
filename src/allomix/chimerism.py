"""MLE chimerism estimation from informative markers.

Implements maximum-likelihood estimation of donor chimerism fraction using
allele counts at informative SNP markers with known host/donor genotypes.
Based on the mixture genotype likelihood of Crysup & Woerner (2022),
Formula 5, simplified for the case of known contributor genotypes.

This module is the optimisation layer: the grid search, Nelder-Mead refinement,
profile-likelihood CIs, Fisher-information detection limits, and the optional
robust refit. The beta-binomial likelihood and expected-weight model it
optimises live in ``allomix.likelihood``; the result data types it returns live
in ``allomix.results``. Import those names from their own modules, not from here.

The in silico CI-coverage and point-estimate accuracy across depths and noise
conditions are reported in the paper's depth and bias-correction validations.
"""

import math
from dataclasses import replace

import numpy as np
from scipy.optimize import brentq, minimize, minimize_scalar
from scipy.stats import chi2, norm

from allomix.constants import (
    CI_LEVEL,
    DEFAULT_ERROR_RATE,
    PLOIDY,
    ROBUST_K_DEFAULT,
)
from allomix.genotype import InformativeMarker
from allomix.likelihood import (
    PanelCalibration,
    _ll_from_p_alt,
    _MarkerArrays,
    _p_alt_for_f,
    _precompute_marker_arrays,
    _total_ll_vec,
    alt_read_probability,
    expected_weight,
    expected_weight_multi,
    total_log_likelihood_multi_bb,
)
from allomix.results import ChimerismResult, MarkerResult, MultiDonorResult

# Robust-refit tuning. The residual cut itself (ROBUST_K_DEFAULT) lives in
# allomix.constants since it is also the CLI/analysis default. The refit floors
# below keep trimming from gutting sparse panels: "auto" never drops below
# ROBUST_MIN_MARKERS, "force" never below ROBUST_HARD_MIN.
ROBUST_MAX_ITER = 5
ROBUST_MIN_MARKERS = 15
ROBUST_HARD_MIN = 4
ROBUST_MODES = ("off", "auto", "force")

# One-sided robust trim. The symmetric median/MAD cut removes a marker whenever
# its residual is a large outlier in either direction. At low host fraction the
# host-carrying markers sit off the donor-dominated fit *in the host-present
# direction* (their observed VAF is pulled toward the host's own allele dose),
# so a symmetric cut trims the very signal we want and collapses the estimate
# toward the donor-only solution. With this on, a marker whose residual deviates
# toward host presence is never trimmed; only residuals pointing away from host
# presence (genotype miscalls, mapping artifacts, host CNV/LoH in the
# anti-host direction) stay eligible for the cut. This is the intended trade:
# at the limit of detection we would rather keep a few artifacts than discard a
# real low-fraction host signal. See claude/further_improvements.md, Obs 1.
ROBUST_ONE_SIDED = True

# Residual outlier cut for the non-robust per-marker flag, in standard
# deviations: a marker whose residual is more than this many SDs from the mean
# residual is marked excluded (flagged, not refit). Used by both the single-
# and multi-donor estimators below.
OUTLIER_SD_THRESHOLD = 3.0
# "auto" engages the refit only when the first pass finds more residual outliers
# than this trigger: max(ROBUST_TRIGGER_MIN, ceil(ROBUST_TRIGGER_FRAC * n)). At
# k=3.5 a clean panel produces <1 chance outlier on average, so a trigger of ~3
# leaves clean (and the already-validated) samples byte-identical, while
# genuine copy-number / LoH clusters clear the bar. "force" uses a trigger of 1.
ROBUST_TRIGGER_FRAC = 0.03
ROBUST_TRIGGER_MIN = 3


# ---------------------------------------------------------------------------
# MLE estimation
# ---------------------------------------------------------------------------


def _compute_per_marker_results(
    markers: list[InformativeMarker],
    f_mle: float,
    calibration: PanelCalibration,
) -> list[MarkerResult]:
    """Compute per-marker residuals and flag outliers.

    Shared between binomial and beta-binomial estimators.
    """
    per_marker: list[MarkerResult] = []
    residuals: list[float] = []

    for m in markers:
        bias = calibration.bias_for(m)
        w = expected_weight(m.host_gt, m.donor_gts[0], f_mle, bias=bias)
        exp_vaf = 1.0 - w  # ALT VAF = 1 - ref_weight
        obs_vaf = m.admix_ad_alt / m.admix_dp if m.admix_dp > 0 else 0.0
        residual = obs_vaf - exp_vaf
        residuals.append(residual)

        per_marker.append(
            MarkerResult(
                chrom=m.chrom,
                pos=m.pos,
                marker_type=m.marker_type,
                expected_vaf=exp_vaf,
                observed_vaf=obs_vaf,
                residual=residual,
                ad_ref=m.admix_ad_ref,
                ad_alt=m.admix_ad_alt,
                dp=m.admix_dp,
                included=True,
            )
        )

    # Flag outliers (residual > OUTLIER_SD_THRESHOLD SDs)
    if len(residuals) >= 2:
        r_arr = np.array(residuals)
        mean_r = float(np.mean(r_arr))
        sd_r = float(np.std(r_arr, ddof=1))

        if sd_r > 0:
            for i, mr in enumerate(per_marker):
                if abs(residuals[i] - mean_r) > OUTLIER_SD_THRESHOLD * sd_r:
                    per_marker[i] = MarkerResult(
                        chrom=mr.chrom,
                        pos=mr.pos,
                        marker_type=mr.marker_type,
                        expected_vaf=mr.expected_vaf,
                        observed_vaf=mr.observed_vaf,
                        residual=mr.residual,
                        ad_ref=mr.ad_ref,
                        ad_alt=mr.ad_alt,
                        dp=mr.dp,
                        included=False,
                    )

    return per_marker


# One-sided normal quantile at CI_LEVEL (z_0.95 ~= 1.6449), used for EP17-style
# LoB/LoD.
_Z95 = float(norm.ppf(CI_LEVEL))

# Margin used to keep a probability strictly inside the open interval (0, 1),
# so that p * (1 - p) stays positive and the marker variance never collapses to
# zero. This is a safety clamp, not machine epsilon (np.finfo(float).eps).
_PROB_EPS = 1e-9

# Consistency factor converting a median absolute deviation to a standard
# deviation for normally distributed data (1 / norm.ppf(0.75)). Used by the
# robust refit to put the MAD-based residual scale on an SD footing.
_MAD_TO_SD = 1.4826

# Feasible range for the beta-binomial concentration rho during optimisation.
# Values outside this are rejected with a large finite penalty (_INFEASIBLE_PENALTY)
# so Nelder-Mead stays in-bounds without a hard constraint.
_RHO_MIN = 0.5
_RHO_MAX = 50000.0
# Upper rho cap for the initial grid search only. Lower than _RHO_MAX on
# purpose: the grid just seeds Nelder-Mead, which then refines rho up to _RHO_MAX.
_RHO_SEED_MAX = 10000.0
_INFEASIBLE_PENALTY = 1e30

# Per-marker-type overdispersion (issue #33). Below this many informative markers
# a class's rho is not identifiable, so the two-rho mode falls back to shared rho.
MIN_CLASS_MARKERS = 30


def _donor_het_mask(markers: list[InformativeMarker]) -> np.ndarray:
    """Boolean mask, True where the (single) donor genotype is heterozygous.

    Donor-het markers sit at background VAF ~0.5 (symmetric amplification
    scatter); donor-hom markers sit near 0/1 (one-sided host signal). The two
    classes carry different overdispersion, which is what the per-marker-type
    mode fits separately (issue #33).

    Args:
        markers: Informative markers (first donor genotype is used).

    Returns:
        Boolean array, one entry per marker, True at donor-het markers.
    """
    return np.fromiter(
        ((m.donor_gts[0][0] + m.donor_gts[0][1]) == 1 for m in markers),
        dtype=bool,
        count=len(markers),
    )


def _donor_het_mask_multi(markers: list[InformativeMarker], n_donors: int) -> np.ndarray:
    """Multi-donor analogue of ``_donor_het_mask`` (phase 2 stub, issue #33).

    True where the combined donor background ALT balance is intermediate (not
    near 0 or 1), the multi-donor counterpart of donor-het: a marker whose pooled
    donor background sits near VAF 0.5 carries the symmetric overdispersion the
    two-rho mode separates out. To be wired into ``_estimate_multi_donor_core``
    when per-marker-type overdispersion is extended to multi-donor; single-donor
    ships first.

    Args:
        markers: Informative markers.
        n_donors: Number of donors.

    Returns:
        Boolean array, one entry per marker.

    Raises:
        NotImplementedError: Always; multi-donor is a later phase.
    """
    raise NotImplementedError(
        "per-marker-type overdispersion is single-donor only in this phase (issue #33)"
    )


def _profile_rho_at_f(arr: _MarkerArrays, f: float, error_rate: float) -> tuple[float, float]:
    """Max LL over rho in [1, _RHO_MAX] at fixed f for one marker class.

    Mirrors the single-rho profile already used in the grid search and CI, so the
    per-class arithmetic is identical to the shared-rho path's.

    Args:
        arr: Per-marker arrays for one class (from ``_precompute_marker_arrays``).
        f: Donor fraction to hold fixed.
        error_rate: Sequencing error rate.

    Returns:
        ``(ll_max, rho_at_max)``.
    """
    p_alt = _p_alt_for_f(arr, f, error_rate)
    opt = minimize_scalar(
        lambda log_r: -_ll_from_p_alt(arr, p_alt, math.exp(log_r)),
        bounds=(math.log(1.0), math.log(_RHO_MAX)),
        method="bounded",
    )
    return -float(opt.fun), math.exp(float(opt.x))


def _two_rho_profile_ll(
    arr_hom: _MarkerArrays, arr_het: _MarkerArrays, f: float, error_rate: float
) -> float:
    """Total profiled LL at fixed f: ``max_rho_hom LL(hom) + max_rho_het LL(het)``.

    Args:
        arr_hom: Donor-hom class arrays.
        arr_het: Donor-het class arrays.
        f: Donor fraction to hold fixed.
        error_rate: Sequencing error rate.

    Returns:
        Sum of the two per-class profiled log-likelihoods.
    """
    ll_hom, _ = _profile_rho_at_f(arr_hom, f, error_rate)
    ll_het, _ = _profile_rho_at_f(arr_het, f, error_rate)
    return ll_hom + ll_het


def fraction_se(
    markers: list[InformativeMarker],
    f_donor: float,
    error_rate: float = DEFAULT_ERROR_RATE,
    rho: float = float("inf"),
    calibration: PanelCalibration | None = None,
    rho_hom: float | None = None,
    rho_het: float | None = None,
) -> float:
    """Standard error of the donor-fraction estimate at a given fraction.

    Computed from the Fisher information of the beta-binomial model used by
    ``estimate_single_donor_bb``, evaluated at ``f_donor``. Each informative
    marker contributes ``(dp_alt/df)^2 / Var(VAF)``, where the VAF variance
    includes the beta-binomial overdispersion inflation
    ``1 + (n - 1) / (rho + 1)``. Markers whose expected ALT fraction does not
    change with f (host and donor ref-dose equal) carry no information and are
    skipped.

    Args:
        markers: Informative markers with admixture allele counts.
        f_donor: Donor fraction at which to evaluate the SE.
        error_rate: Sequencing error rate.
        rho: Beta-binomial concentration (inf = pure binomial).
        calibration: Optional per-marker bias and error tables (only biases
            are used here).
        rho_hom: Donor-hom class concentration. When both ``rho_hom`` and
            ``rho_het`` are given, each marker uses its class rho instead of the
            scalar ``rho`` (per-marker-type overdispersion, issue #33).
        rho_het: Donor-het class concentration (see ``rho_hom``).

    Returns:
        Standard error of the donor fraction. inf if no marker is informative.
    """
    cal = calibration or PanelCalibration()
    per_class_rho = rho_hom is not None and rho_het is not None
    # P(observe ALT) is linear in the REF weight w, so its slope dp_alt/dw is
    # constant and equals the change across the full weight range w: 0 -> 1.
    dpalt_dw = alt_read_probability(1.0, error_rate) - alt_read_probability(0.0, error_rate)

    info = 0.0
    for m in markers:
        bias = cal.bias_for(m)

        host_ref_dose = PLOIDY - (m.host_gt[0] + m.host_gt[1])
        donor_ref_dose = PLOIDY - (m.donor_gts[0][0] + m.donor_gts[0][1])
        dw_df = (donor_ref_dose - host_ref_dose) / PLOIDY
        if dw_df == 0.0:
            continue

        n = m.admix_ad_ref + m.admix_ad_alt
        if n == 0:
            continue

        w = expected_weight(m.host_gt, m.donor_gts[0], f_donor, bias=bias)
        p_alt = alt_read_probability(w, error_rate)
        p_alt = max(_PROB_EPS, min(1.0 - _PROB_EPS, p_alt))

        if per_class_rho:
            donor_alt_dose = m.donor_gts[0][0] + m.donor_gts[0][1]
            m_rho = rho_het if donor_alt_dose == 1 else rho_hom
        else:
            m_rho = rho
        overdispersion = 1.0 if math.isinf(m_rho) else 1.0 + (n - 1.0) / (m_rho + 1.0)
        var_vaf = p_alt * (1.0 - p_alt) / n * overdispersion
        if var_vaf <= 0.0:
            continue

        dpalt_df = dpalt_dw * dw_df
        info += (dpalt_df * dpalt_df) / var_vaf

    if info <= 0.0:
        return float("inf")
    return 1.0 / math.sqrt(info)


def detection_limit(
    markers: list[InformativeMarker],
    error_rate: float = DEFAULT_ERROR_RATE,
    rho: float = float("inf"),
    calibration: PanelCalibration | None = None,
    rho_hom: float | None = None,
    rho_het: float | None = None,
) -> tuple[float, float]:
    """Per-sample limit of blank and limit of detection (donor fractions).

    Follows the single-replicate Currie / CLSI EP17-A2 construction, using the
    Fisher information of this sample's own markers in place of repeated blank
    and low-level measurements:

        LoB = z * SE(f = 0)
        LoD = LoB + z * SE(f = LoB)

    with ``z`` the one-sided 95% normal quantile. The donor-fraction estimator
    is bounded at 0, so its upper 95th percentile under a true blank is
    ``z * SE(0)`` even though the lower tail piles at 0.

    This is the best achievable sensitivity given the fitted noise model
    (``rho``) and known/corrected biases. It is not a substitute for a
    validated assay LoD from a blank and dilution series, which must come from
    replicated experiments (simulation or wetlab).

    Args:
        markers: Informative markers with admixture allele counts.
        error_rate: Sequencing error rate.
        rho: Beta-binomial concentration from the fit (inf = pure binomial).
        calibration: Optional per-marker bias and error tables (only biases
            are used here).
        rho_hom: Donor-hom class concentration. When both ``rho_hom`` and
            ``rho_het`` are given, each marker uses its class rho instead of the
            scalar ``rho`` (per-marker-type overdispersion, issue #33). Forwarded
            to both ``fraction_se`` calls.
        rho_het: Donor-het class concentration (see ``rho_hom``).

    Returns:
        ``(lob, lod)`` as donor fractions (0.0-1.0). ``(inf, inf)`` if no
        marker is informative.
    """
    se0 = fraction_se(markers, 0.0, error_rate, rho, calibration, rho_hom, rho_het)
    if math.isinf(se0):
        return float("inf"), float("inf")
    lob = _Z95 * se0
    se_lob = fraction_se(markers, lob, error_rate, rho, calibration, rho_hom, rho_het)
    if math.isinf(se_lob):
        se_lob = se0
    lod = lob + _Z95 * se_lob
    return lob, lod


def _estimate_single_donor_bb_core(
    markers: list[InformativeMarker],
    error_rate: float = DEFAULT_ERROR_RATE,
    grid_steps: int = 1001,
    calibration: PanelCalibration | None = None,
    marker_type_overdispersion: bool = True,
) -> ChimerismResult:
    """Single-donor beta-binomial MLE over the given marker set (no robust trim).

    This is the unguarded estimator; ``estimate_single_donor_bb`` wraps it with
    the optional robust refit. Args/returns as in that wrapper, plus
    ``marker_type_overdispersion`` (issue #33): when True (the default), fit a
    separate rho for the donor-hom and donor-het marker classes instead of one
    shared rho, which removes the sub-0.5% MLE floor. Set False to recover the
    legacy shared-rho path (byte-identical to the pre-#33 estimator). When a
    marker class has fewer than ``MIN_CLASS_MARKERS`` markers the per-class rho
    is not identifiable, so the estimator falls back to shared rho for that
    sample and records the reason on the result (diagnostic only, not a warning).
    """
    cal = calibration or PanelCalibration()
    n_informative = len(markers)

    if n_informative == 0:
        return ChimerismResult(
            donor_fraction=0.0,
            donor_fraction_ci=(0.0, 0.0),
            host_fraction=1.0,
            log_likelihood=0.0,
            n_informative=0,
            n_markers_used=0,
            per_marker=[],
            error_rate=error_rate,
        )

    if marker_type_overdispersion:
        het_mask = _donor_het_mask(markers)
        n_het = int(het_mask.sum())
        n_hom = n_informative - n_het
        if n_het >= MIN_CLASS_MARKERS and n_hom >= MIN_CLASS_MARKERS:
            return _estimate_single_donor_two_rho(markers, het_mask, error_rate, grid_steps, cal)
        # Sparse class: its per-class rho is not identifiable, so fall through to
        # the shared-rho path for this sample and record the reason on the result.
        # This is routine for small or hom-dominated panels and is a diagnostic
        # field only, not a QC warning (two-rho is the default, not a request).
        fell_back: str | None = (
            f"a marker class is sparse (hom={n_hom}, het={n_het}, "
            f"min={MIN_CLASS_MARKERS}); used shared rho for this sample"
        )
    else:
        fell_back = None

    # This function is the hot path of every validation sweep in the paper (LoD,
    # relatedness, depth, overdispersion all bottom out here, millions of calls).
    # Per profiling it is ~99% of those builds; simulation and VCF IO are noise.
    # The cost lives in the two scipy searches below (grid rho-profiling and the
    # Nelder-Mead refinement) hammering _ll_from_p_alt. Optimise here, exactly:
    # the f-invariant work (w/bias/p_alt via _p_alt_for_f) is already hoisted out
    # of the rho loops, and the per-call-constant terms are precomputed in
    # _precompute_marker_arrays. Keep any further change bit-identical (this
    # estimator's output is validated against fixtures; ~0.1% drift matters at
    # the low-fraction limit of detection).

    # Precompute the (f, rho)-independent per-marker arrays once and reuse them
    # across the grid search, Nelder-Mead refinement, and profile-likelihood CI.
    arr = _precompute_marker_arrays(markers, cal)

    # Step 1: Grid search over f with rho profiled out at each grid point
    grid = np.linspace(0.0, 1.0, grid_steps)
    best_ll = -math.inf
    best_f = 0.0
    best_rho = 100.0

    # Hoisted out of the grid loop: the rho-profiling bounds are the same at
    # every grid point.
    grid_rho_bounds = (math.log(1.0), math.log(_RHO_SEED_MAX))
    for f in grid:
        # p_alt depends on f only, so compute it once and reuse it across the
        # rho profiling (was recomputed at every rho evaluation).
        p_alt_f = _p_alt_for_f(arr, f, error_rate)
        opt_rho = minimize_scalar(
            lambda log_r, _p=p_alt_f: -_ll_from_p_alt(arr, _p, math.exp(log_r)),
            bounds=grid_rho_bounds,
            method="bounded",
        )
        rho_cand = math.exp(float(opt_rho.x))
        ll_cand = -float(opt_rho.fun)
        if ll_cand > best_ll:
            best_ll = ll_cand
            best_f = float(f)
            best_rho = rho_cand

    # Step 2: Joint Nelder-Mead refinement over (f, log_rho)
    def neg_ll_joint(x):
        f_val, log_rho_val = x
        if f_val < 0.0 or f_val > 1.0:
            return _INFEASIBLE_PENALTY
        rho_val = math.exp(log_rho_val)
        if rho_val < _RHO_MIN or rho_val > _RHO_MAX:
            return _INFEASIBLE_PENALTY
        return -_total_ll_vec(arr, f_val, error_rate, rho_val)

    opt = minimize(
        neg_ll_joint,
        x0=[best_f, math.log(best_rho)],
        method="Nelder-Mead",
        options={"xatol": 1e-6, "fatol": 1e-10, "maxiter": 5000},
    )

    f_mle = max(0.0, min(1.0, float(opt.x[0])))
    rho_mle = math.exp(float(opt.x[1]))

    # Step 3: Profile likelihood CIs for f, profiling out rho at each f.
    # rho upper bound matches the Nelder-Mead constraint (50000) to avoid
    # profile_ll_f(f_mle) < ll_max_joint, which causes brentq sign errors.
    threshold = chi2.ppf(CI_LEVEL, df=1)
    half_threshold = threshold / 2.0

    def profile_ll_f(f_val: float) -> float:
        """Max LL over rho at a given f."""
        p_alt_f = _p_alt_for_f(arr, f_val, error_rate)
        opt_rho = minimize_scalar(
            lambda log_r: -_ll_from_p_alt(arr, p_alt_f, math.exp(log_r)),
            bounds=(math.log(1.0), math.log(_RHO_MAX)),
            method="bounded",
        )
        return -float(opt_rho.fun)

    # Re-derive ll_max from the profile at f_mle so the CI reference is
    # consistent with profile_ll_f across the search interval.
    ll_max = profile_ll_f(f_mle)

    def ci_func(f_val: float) -> float:
        return ll_max - profile_ll_f(f_val) - half_threshold

    # Lower bound
    if f_mle <= 0.0 or ci_func(0.0) <= 0.0:
        f_lo = 0.0
    else:
        f_lo = brentq(ci_func, 0.0, f_mle, xtol=1e-5)

    # Upper bound
    if f_mle >= 1.0 or ci_func(1.0) <= 0.0:
        f_hi = 1.0
    else:
        f_hi = brentq(ci_func, f_mle, 1.0, xtol=1e-5)

    # Step 4: Per-marker residuals
    per_marker = _compute_per_marker_results(markers, f_mle, cal)
    n_markers_used = sum(1 for mr in per_marker if mr.included)

    # Step 5: Per-sample analytical detection limits from the fitted noise model.
    lob, lod = detection_limit(markers, error_rate, rho_mle, cal)

    return ChimerismResult(
        donor_fraction=f_mle,
        donor_fraction_ci=(f_lo, f_hi),
        host_fraction=1.0 - f_mle,
        log_likelihood=ll_max,
        n_informative=n_informative,
        n_markers_used=n_markers_used,
        per_marker=per_marker,
        error_rate=error_rate,
        rho=rho_mle,
        marker_type_overdispersion_fallback=fell_back,
        lob_fraction=lob,
        lod_fraction=lod,
    )


def _estimate_single_donor_two_rho(
    markers: list[InformativeMarker],
    het_mask: np.ndarray,
    error_rate: float,
    grid_steps: int,
    cal: PanelCalibration,
) -> ChimerismResult:
    """Single-donor MLE with a separate rho per marker class (issue #33).

    Fits the donor fraction f jointly with two independent concentration
    parameters, one for the donor-hom markers and one for the donor-het markers,
    profiling each rho out at every f. This down-weights the over-dispersed
    donor-het class (background VAF ~0.5) at low fraction, where its symmetric
    amplification scatter otherwise rectifies into a small positive host
    fraction (the sub-0.5% floor). The structure mirrors the shared-rho core: a
    grid search over f, a Nelder-Mead refinement, a profile-likelihood CI, then
    per-marker residuals over the full marker set and per-class detection limits.

    Args:
        markers: Informative markers (single donor).
        het_mask: Boolean mask from ``_donor_het_mask`` (True at donor-het).
        error_rate: Sequencing error rate.
        grid_steps: Number of f grid points for the initial search.
        cal: Per-marker calibration.

    Returns:
        ChimerismResult with ``rho`` set to the het-class rho and ``rho_hom`` /
        ``rho_het`` populated.
    """
    hom_markers = [m for m, h in zip(markers, het_mask) if not h]
    het_markers = [m for m, h in zip(markers, het_mask) if h]
    arr_hom = _precompute_marker_arrays(hom_markers, cal)
    arr_het = _precompute_marker_arrays(het_markers, cal)

    # Step 1: grid over f, both rhos profiled out at each grid point.
    grid = np.linspace(0.0, 1.0, grid_steps)
    best_ll = -math.inf
    best_f = 0.0
    best_rho_hom = 100.0
    best_rho_het = 100.0
    for f in grid:
        ll_hom, r_hom = _profile_rho_at_f(arr_hom, f, error_rate)
        ll_het, r_het = _profile_rho_at_f(arr_het, f, error_rate)
        ll = ll_hom + ll_het
        if ll > best_ll:
            best_ll = ll
            best_f = float(f)
            best_rho_hom = r_hom
            best_rho_het = r_het

    # Step 2: Nelder-Mead refinement over (f, log_rho_hom, log_rho_het).
    def neg_ll_joint(x):
        f_val, lr_hom, lr_het = x
        if f_val < 0.0 or f_val > 1.0:
            return _INFEASIBLE_PENALTY
        r_hom = math.exp(lr_hom)
        r_het = math.exp(lr_het)
        if not (_RHO_MIN <= r_hom <= _RHO_MAX) or not (_RHO_MIN <= r_het <= _RHO_MAX):
            return _INFEASIBLE_PENALTY
        return -(
            _total_ll_vec(arr_hom, f_val, error_rate, r_hom)
            + _total_ll_vec(arr_het, f_val, error_rate, r_het)
        )

    opt = minimize(
        neg_ll_joint,
        x0=[best_f, math.log(best_rho_hom), math.log(best_rho_het)],
        method="Nelder-Mead",
        options={"xatol": 1e-6, "fatol": 1e-10, "maxiter": 5000},
    )
    f_mle = max(0.0, min(1.0, float(opt.x[0])))
    rho_hom_mle = math.exp(float(opt.x[1]))
    rho_het_mle = math.exp(float(opt.x[2]))

    # Step 3: profile-likelihood CI for f, both rhos profiled out at each f.
    threshold = chi2.ppf(CI_LEVEL, df=1)
    half_threshold = threshold / 2.0

    def profile_ll_f(f_val: float) -> float:
        return _two_rho_profile_ll(arr_hom, arr_het, f_val, error_rate)

    # Re-derive ll_max from the profile at f_mle so the CI reference is
    # consistent with profile_ll_f across the search interval.
    ll_max = profile_ll_f(f_mle)

    def ci_func(f_val: float) -> float:
        return ll_max - profile_ll_f(f_val) - half_threshold

    if f_mle <= 0.0 or ci_func(0.0) <= 0.0:
        f_lo = 0.0
    else:
        f_lo = brentq(ci_func, 0.0, f_mle, xtol=1e-5)

    if f_mle >= 1.0 or ci_func(1.0) <= 0.0:
        f_hi = 1.0
    else:
        f_hi = brentq(ci_func, f_mle, 1.0, xtol=1e-5)

    # Step 4: per-marker residuals over the full marker set at the single f_mle,
    # so per-marker output, residuals and the robust-refit interaction are
    # unchanged in shape; only the f they are evaluated at moves.
    per_marker = _compute_per_marker_results(markers, f_mle, cal)
    n_markers_used = sum(1 for mr in per_marker if mr.included)

    # Step 5: per-sample detection limits under the per-class rhos.
    lob, lod = detection_limit(
        markers,
        error_rate,
        rho=rho_het_mle,
        calibration=cal,
        rho_hom=rho_hom_mle,
        rho_het=rho_het_mle,
    )

    return ChimerismResult(
        donor_fraction=f_mle,
        donor_fraction_ci=(f_lo, f_hi),
        host_fraction=1.0 - f_mle,
        log_likelihood=ll_max,
        n_informative=len(markers),
        n_markers_used=n_markers_used,
        per_marker=per_marker,
        error_rate=error_rate,
        rho=rho_het_mle,  # headline rho = het class (governs the floor and low-f CI)
        rho_hom=rho_hom_mle,
        rho_het=rho_het_mle,
        lob_fraction=lob,
        lod_fraction=lod,
    )


def _marker_key(m: InformativeMarker) -> tuple[str, int, str, str]:
    """Position-based key used to track which markers survive a robust trim."""
    return (m.chrom, m.pos, m.ref, m.alt)


def _robust_trigger(n_markers: int) -> int:
    """Minimum first-pass outliers for "auto" to engage the robust refit."""
    return max(ROBUST_TRIGGER_MIN, math.ceil(ROBUST_TRIGGER_FRAC * n_markers))


def _robust_refit(
    all_markers: list[InformativeMarker],
    core_fn,
    recompute_per_marker_fn,
    robust_k: float,
    min_markers: int,
    min_trigger: int,
):
    """Iteratively drop residual-outlier markers (median/MAD) and refit.

    The built-in 3-SD flag never refits, so a cluster of host copy-number /
    LoH-inconsistent markers biases the estimate and escapes the flag (it sets
    the very SD the flag uses). This trims on a robust median/MAD scale, which a
    minority of contaminating markers cannot inflate, and refits on the
    survivors until the set stabilises.

    Args:
        all_markers: Full informative-marker set.
        core_fn: Callable taking a marker subset and returning a fitted result.
        recompute_per_marker_fn: Callable (markers, result) -> per-marker results
            evaluated at ``result``'s estimate, used to re-flag every original
            marker against the final (survivor) fit.
        robust_k: Residual cut in robust SDs.
        min_markers: Trimming never reduces the surviving set below this.
        min_trigger: Engage the refit only if the first pass finds at least this
            many residual outliers (the gate that keeps clean samples unchanged).

    Returns:
        A result equal to ``core_fn`` on the survivors, but with ``per_marker``
        spanning all original markers (``included`` False for dropped ones),
        ``n_informative`` over all markers, and the robust-accounting fields set.
        When nothing is dropped, returns ``core_fn(all_markers)`` unchanged.
    """
    current = list(all_markers)
    result = core_fn(current)
    for it in range(ROBUST_MAX_ITER):
        if len(current) <= min_markers:
            break
        resids = np.array([mr.residual for mr in result.per_marker], dtype=float)
        med = float(np.median(resids))
        mad = float(np.median(np.abs(resids - med))) * _MAD_TO_SD
        if mad <= 0.0:
            break
        deviation = resids - med
        keep_mask = np.abs(deviation) <= robust_k * mad
        if ROBUST_ONE_SIDED:
            # Protect markers whose residual deviates toward host presence.
            # Increasing the host weight moves a marker's expected ALT VAF toward
            # the host's own ALT dose (host_alt / PLOIDY), so the host-present
            # direction at the current fit is sign(host_alt / PLOIDY - expected_vaf).
            # A residual deviating that way is under-fit host signal, not an
            # artifact, and must not be trimmed. ``current`` and
            # ``result.per_marker`` are in the same order (the core estimator
            # builds per-marker results by iterating the marker list).
            host_alt = np.array([m.host_gt[0] + m.host_gt[1] for m in current], dtype=float)
            exp_vaf = np.array([mr.expected_vaf for mr in result.per_marker], dtype=float)
            host_dir = np.sign(host_alt / PLOIDY - exp_vaf)
            points_to_host = (host_dir != 0.0) & (np.sign(deviation) == host_dir)
            keep_mask = keep_mask | points_to_host
        n_keep = int(keep_mask.sum())
        # Gate: on the first pass, only engage if outliers exceed the trigger,
        # so clean samples (a chance outlier or two) are left untouched.
        if it == 0 and (len(current) - n_keep) < min_trigger:
            break
        if n_keep == len(current) or n_keep < min_markers:
            break
        current = [m for m, k in zip(current, keep_mask) if k]
        result = core_fn(current)

    n_excluded = len(all_markers) - len(current)
    if n_excluded == 0:
        return result  # untouched: identical to the non-robust path

    surviving = {_marker_key(m) for m in current}
    pm_all = recompute_per_marker_fn(all_markers, result)
    pm_final = [
        replace(mr, included=_marker_key(m) in surviving) for m, mr in zip(all_markers, pm_all)
    ]
    return replace(
        result,
        per_marker=pm_final,
        n_informative=len(all_markers),
        n_markers_used=len(current),
        n_robust_excluded=n_excluded,
        robust_drop_fraction=n_excluded / len(all_markers),
    )


def estimate_single_donor_bb(
    markers: list[InformativeMarker],
    error_rate: float = DEFAULT_ERROR_RATE,
    grid_steps: int = 1001,
    calibration: PanelCalibration | None = None,
    robust: str = "off",
    robust_k: float = ROBUST_K_DEFAULT,
    marker_type_overdispersion: bool = True,
) -> ChimerismResult:
    """Estimate single-donor chimerism with beta-binomial likelihood.

    Estimates single-donor chimerism fraction using beta-binomial
    per-marker likelihoods to handle overdispersion. Jointly estimates
    the donor fraction f and concentration parameter rho.

    Args:
        markers: List of informative markers with admixture allele counts.
        error_rate: Sequencing error rate (fallback when a marker is missing
            per-direction rates).
        grid_steps: Number of grid points for initial f search.
        calibration: Optional per-marker bias and error tables. For a marker
            with both per-direction error rates known, the asymmetric
            REF/ALT-only likelihood is used; otherwise the symmetric 4-state
            model with ``error_rate`` is used.
        robust: Robust-refit mode. ``"off"`` (default) is the plain MLE.
            ``"auto"`` iteratively drops median/MAD residual outliers and refits,
            never below ``ROBUST_MIN_MARKERS`` survivors; this protects against
            host copy-number / LoH markers (see ``_robust_refit``). ``"force"``
            is the same but trims down to ``ROBUST_HARD_MIN`` (for experiments).
            On clean data the trim is a no-op and the result is unchanged.
        robust_k: Residual cut in robust SDs for the refit.
        marker_type_overdispersion: Fit a separate beta-binomial concentration
            for the donor-hom and donor-het marker classes instead of one shared
            rho (issue #33). On by default; removes the sub-0.5% MLE floor. Set
            False for the legacy shared-rho path (byte-identical to the pre-#33
            estimator). Falls back to shared rho for a sample when a class has
            fewer than ``MIN_CLASS_MARKERS`` markers.

    Returns:
        ChimerismResult with MLE estimate and beta-binomial CIs.
    """
    if robust not in ROBUST_MODES:
        raise ValueError(f"robust must be one of {ROBUST_MODES}, got {robust!r}")

    cal = calibration or PanelCalibration()

    def core(mk: list[InformativeMarker]) -> ChimerismResult:
        return _estimate_single_donor_bb_core(
            mk, error_rate, grid_steps, cal, marker_type_overdispersion
        )

    if robust == "off" or len(markers) == 0:
        return core(markers)

    min_markers = ROBUST_MIN_MARKERS if robust == "auto" else ROBUST_HARD_MIN
    min_trigger = _robust_trigger(len(markers)) if robust == "auto" else 1
    return _robust_refit(
        markers,
        core,
        lambda mk, res: _compute_per_marker_results(mk, res.donor_fraction, cal),
        robust_k,
        min_markers,
        min_trigger,
    )


# ---------------------------------------------------------------------------
# Multi-donor MLE estimation
# ---------------------------------------------------------------------------


def _estimate_multi_donor_core(
    markers: list[InformativeMarker],
    n_donors: int = 2,
    error_rate: float = DEFAULT_ERROR_RATE,
    grid_steps: int = 101,
    calibration: PanelCalibration | None = None,
) -> MultiDonorResult:
    """Estimate multi-donor chimerism fractions via maximum likelihood.

    Uses a beta-binomial likelihood to handle overdispersion, jointly
    estimating donor fractions and the concentration parameter rho.

    Algorithm:
        1. Triangular grid search over (f1, f2) at fixed rho
        2. Nelder-Mead refinement over (f1, f2, log_rho)
        3. Profile likelihood CI per donor, profiling out other f and rho
        4. Per-marker residuals and outlier flagging

    Args:
        markers: Informative markers (for at least one donor).
        n_donors: Number of donors (currently supports 2).
        error_rate: Sequencing error rate (fallback).
        grid_steps: Grid resolution per dimension.
        calibration: Optional per-marker bias and error tables
            (see ``estimate_single_donor_bb``).

    Returns:
        MultiDonorResult with per-donor fractions and CIs.
    """
    cal = calibration or PanelCalibration()
    if n_donors > 2:
        raise ValueError(
            f"n_donors={n_donors} not supported; "
            "estimate_multi_donor currently supports up to 2 donors"
        )

    n_informative = len(markers)

    if n_informative == 0:
        return MultiDonorResult(
            donor_fractions=[0.0] * n_donors,
            donor_fraction_cis=[(0.0, 0.0)] * n_donors,
            host_fraction=1.0,
            log_likelihood=0.0,
            n_informative=0,
            n_markers_used=0,
            per_marker=[],
            error_rate=error_rate,
            per_donor_n_informative=[0] * n_donors,
        )

    # Step 1: Triangular grid search at fixed rho
    rho_init = 100.0
    best_ll = -math.inf
    best_f = [0.0] * n_donors
    step = 1.0 / (grid_steps - 1)

    for i in range(grid_steps):
        f1 = i * step
        if f1 > 1.0:
            break
        max_f2 = 1.0 - f1
        n_f2_steps = int(max_f2 / step) + 1
        for j in range(n_f2_steps):
            f2 = j * step
            if f1 + f2 > 1.0 + 1e-9:
                break
            ll = total_log_likelihood_multi_bb(
                markers,
                [f1, f2],
                error_rate,
                rho_init,
                cal,
            )
            if ll > best_ll:
                best_ll = ll
                best_f = [f1, f2]

    # Step 2: Nelder-Mead refinement over (f1, f2, log_rho)
    def neg_ll(x):
        f1, f2, log_rho = x
        if f1 < 0 or f2 < 0 or f1 + f2 > 1.0:
            return _INFEASIBLE_PENALTY
        rho = math.exp(log_rho)
        if rho < _RHO_MIN or rho > _RHO_MAX:
            return _INFEASIBLE_PENALTY
        return -total_log_likelihood_multi_bb(
            markers,
            [f1, f2],
            error_rate,
            rho,
            cal,
        )

    opt = minimize(
        neg_ll,
        x0=[best_f[0], best_f[1], math.log(rho_init)],
        method="Nelder-Mead",
        options={"xatol": 1e-5, "fatol": 1e-8, "maxiter": 5000},
    )

    f_mle = [max(0.0, float(x)) for x in opt.x[:2]]
    if sum(f_mle) > 1.0:
        scale = 1.0 / sum(f_mle)
        f_mle = [f * scale for f in f_mle]
    rho_mle = math.exp(float(opt.x[2]))
    ll_max = -float(opt.fun)

    # Step 3: Profile likelihood CIs per donor (profiling out other f and rho)
    cis = _profile_likelihood_cis_multi(
        markers,
        f_mle,
        n_donors,
        error_rate,
        cal,
    )

    # Step 4: Per-marker residuals
    per_marker = _per_marker_results_multi(markers, f_mle, cal)

    # Per-donor informative counts
    per_donor_n_inf = [0] * n_donors
    for m in markers:
        if m.informative_for is not None:
            for d_idx in range(n_donors):
                if d_idx < len(m.informative_for) and m.informative_for[d_idx]:
                    per_donor_n_inf[d_idx] += 1
        else:
            # Legacy: only first donor tracked
            per_donor_n_inf[0] += 1

    n_markers_used = sum(1 for mr in per_marker if mr.included)

    return MultiDonorResult(
        donor_fractions=f_mle,
        donor_fraction_cis=cis,
        host_fraction=1.0 - sum(f_mle),
        log_likelihood=ll_max,
        n_informative=n_informative,
        n_markers_used=n_markers_used,
        per_marker=per_marker,
        error_rate=error_rate,
        per_donor_n_informative=per_donor_n_inf,
        rho=rho_mle,
    )


def estimate_multi_donor(
    markers: list[InformativeMarker],
    n_donors: int = 2,
    error_rate: float = DEFAULT_ERROR_RATE,
    grid_steps: int = 101,
    calibration: PanelCalibration | None = None,
    robust: str = "off",
    robust_k: float = ROBUST_K_DEFAULT,
) -> MultiDonorResult:
    """Estimate multi-donor chimerism with an optional robust refit.

    Wraps the multi-donor beta-binomial MLE with the same median/MAD
    robust-refit logic as ``estimate_single_donor_bb`` (see its ``robust`` and
    ``robust_k`` arguments and ``_robust_refit``). With ``robust="off"`` this is
    the plain estimator.

    Returns:
        MultiDonorResult with per-donor fractions and CIs.
    """
    if robust not in ROBUST_MODES:
        raise ValueError(f"robust must be one of {ROBUST_MODES}, got {robust!r}")

    cal = calibration or PanelCalibration()

    def core(mk: list[InformativeMarker]) -> MultiDonorResult:
        return _estimate_multi_donor_core(mk, n_donors, error_rate, grid_steps, cal)

    if robust == "off" or len(markers) == 0:
        return core(markers)

    min_markers = ROBUST_MIN_MARKERS if robust == "auto" else ROBUST_HARD_MIN
    min_trigger = _robust_trigger(len(markers)) if robust == "auto" else 1
    return _robust_refit(
        markers,
        core,
        lambda mk, res: _per_marker_results_multi(mk, res.donor_fractions, cal),
        robust_k,
        min_markers,
        min_trigger,
    )


def _profile_likelihood_cis_multi(
    markers: list[InformativeMarker],
    f_mle: list[float],
    n_donors: int,
    error_rate: float,
    calibration: PanelCalibration,
) -> list[tuple[float, float]]:
    """Profile likelihood CIs for each donor fraction.

    For donor_i, scan f_i while optimizing f_j (j != i) and rho at each
    point. Uses chi2(df=1) threshold since we profile one parameter at
    a time. Overdispersion is handled by the beta-binomial likelihood
    (via rho profiling) rather than by inflating the threshold.
    """
    threshold = float(chi2.ppf(CI_LEVEL, df=1))
    half_threshold = threshold / 2.0
    cis: list[tuple[float, float]] = []

    for donor_idx in range(n_donors):
        # The estimator is capped at two donors (enforced in
        # _estimate_multi_donor_core), so profiling one donor leaves exactly one
        # "other" donor: the remaining index.
        other_idx = next(j for j in range(n_donors) if j != donor_idx)

        def profile_ll(fi: float, _other=other_idx, _didx=donor_idx) -> float:
            """Max LL over the other donor and rho, with donor_idx fixed at fi."""
            max_fj = max(0.0, 1.0 - fi)
            if max_fj < 1e-9:
                # Only rho to optimise
                fracs = [0.0, 0.0]
                fracs[_didx] = fi
                opt_rho = minimize_scalar(
                    lambda log_r: (
                        -total_log_likelihood_multi_bb(
                            markers,
                            fracs,
                            error_rate,
                            math.exp(log_r),
                            calibration,
                        )
                    ),
                    bounds=(math.log(1.0), math.log(_RHO_MAX)),
                    method="bounded",
                )
                return -float(opt_rho.fun)

            # Optimise over (fj, log_rho) jointly
            def neg_ll_inner(x, _di=_didx):
                fj, log_r = x
                if fj < 0 or fj > max_fj:
                    return _INFEASIBLE_PENALTY
                rho = math.exp(log_r)
                if rho < _RHO_MIN or rho > _RHO_MAX:
                    return _INFEASIBLE_PENALTY
                fracs = [fi, fj] if _di == 0 else [fj, fi]
                return -total_log_likelihood_multi_bb(
                    markers,
                    fracs,
                    error_rate,
                    rho,
                    calibration,
                )

            opt = minimize(
                neg_ll_inner,
                x0=[f_mle[_other], math.log(100.0)],
                method="Nelder-Mead",
                options={"xatol": 1e-5, "fatol": 1e-8, "maxiter": 2000},
            )
            return -float(opt.fun)

        fi_mle = f_mle[donor_idx]
        # Re-derive reference LL from the profile at fi_mle for consistency
        # with profile_ll across the search interval (avoids brentq sign errors).
        pll_at_mle = profile_ll(fi_mle)

        def ci_func(fi: float, _pll=profile_ll, _ref=pll_at_mle) -> float:
            return _ref - _pll(fi) - half_threshold

        # Lower bound
        if fi_mle <= 0.0 or ci_func(0.0) <= 0.0:
            f_lo = 0.0
        else:
            f_lo = brentq(ci_func, 0.0, fi_mle, xtol=1e-5)

        # Upper bound
        if fi_mle >= 1.0 or ci_func(1.0) <= 0.0:
            f_hi = 1.0
        else:
            f_hi = brentq(ci_func, fi_mle, 1.0, xtol=1e-5)

        cis.append((f_lo, f_hi))

    return cis


def _per_marker_results_multi(
    markers: list[InformativeMarker],
    f_mle: list[float],
    calibration: PanelCalibration,
) -> list[MarkerResult]:
    """Compute per-marker residuals for multi-donor model."""
    per_marker: list[MarkerResult] = []
    residuals: list[float] = []

    for m in markers:
        bias = calibration.bias_for(m)
        w = expected_weight_multi(m.host_gt, m.donor_gts, f_mle, bias=bias)
        exp_vaf = 1.0 - w
        obs_vaf = m.admix_ad_alt / m.admix_dp if m.admix_dp > 0 else 0.0
        residual = obs_vaf - exp_vaf
        residuals.append(residual)

        per_marker.append(
            MarkerResult(
                chrom=m.chrom,
                pos=m.pos,
                marker_type=m.marker_type,
                expected_vaf=exp_vaf,
                observed_vaf=obs_vaf,
                residual=residual,
                ad_ref=m.admix_ad_ref,
                ad_alt=m.admix_ad_alt,
                dp=m.admix_dp,
                included=True,
            )
        )

    # Flag outliers (residual > OUTLIER_SD_THRESHOLD SDs)
    if len(residuals) >= 2:
        r_arr = np.array(residuals)
        mean_r = float(np.mean(r_arr))
        sd_r = float(np.std(r_arr, ddof=1))

        if sd_r > 0:
            for i, mr in enumerate(per_marker):
                if abs(residuals[i] - mean_r) > OUTLIER_SD_THRESHOLD * sd_r:
                    per_marker[i] = MarkerResult(
                        chrom=mr.chrom,
                        pos=mr.pos,
                        marker_type=mr.marker_type,
                        expected_vaf=mr.expected_vaf,
                        observed_vaf=mr.observed_vaf,
                        residual=mr.residual,
                        ad_ref=mr.ad_ref,
                        ad_alt=mr.ad_alt,
                        dp=mr.dp,
                        included=False,
                    )

    return per_marker
