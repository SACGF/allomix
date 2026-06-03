"""Core MLE chimerism estimation from informative markers.

Implements maximum-likelihood estimation of donor chimerism fraction using
allele counts at informative SNP markers with known host/donor genotypes.
Based on the mixture genotype likelihood of Crysup & Woerner (2022),
Formula 5, simplified for the case of known contributor genotypes.

Uses a beta-binomial likelihood to handle overdispersion from per-marker
amplification bias and depth variability. A standard binomial model
assumes all variance comes from random sampling; in practice, systematic
effects produce extra-binomial variance that causes binomial CIs to
undercover. The beta-binomial adds a shared concentration parameter rho
that is jointly estimated from the data, naturally widening CIs when
overdispersion is present.

Benchmarking on synthetic data with realistic noise (40 markers, 2000x
depth, bias SD 0.02, depth CV 0.4, 100 replicates x 10 fractions):

    Beta-binomial CI coverage:  88.2%  (vs 79.8% binomial)
    Point estimate MAE:         0.0029 (identical to binomial)
    Mean CI width:              1.3%   (vs 1.1% binomial)

The improvement is largest at low donor fractions (f=1%: 96% vs 86%
coverage; f=2%: 100% vs 84%; f=5%: 88% vs 80%) where accurate CIs
matter most clinically.
"""

import math
from dataclasses import dataclass, replace
from math import lgamma

import numpy as np
from scipy.optimize import brentq, minimize, minimize_scalar
from scipy.special import gammaln
from scipy.stats import chi2, norm

from allomix.detect import HostPresenceResult
from allomix.genotype import InformativeMarker
from allomix.relatedness import AdmixConsistencyResult, RelatednessResult

# ---------------------------------------------------------------------------
# Result types
# ---------------------------------------------------------------------------


@dataclass
class MarkerResult:
    """Per-marker contribution to the chimerism estimate."""

    chrom: str
    pos: int
    marker_type: int
    expected_vaf: float
    observed_vaf: float
    residual: float
    ad_ref: int
    ad_alt: int
    dp: int
    included: bool  # False if outlier-excluded


@dataclass
class ChimerismResult:
    """Result of single-donor chimerism estimation."""

    donor_fraction: float  # MLE point estimate (0.0-1.0)
    donor_fraction_ci: tuple[float, float]  # 95% CI
    host_fraction: float  # 1 - donor_fraction
    log_likelihood: float  # at MLE
    n_informative: int
    n_markers_used: int  # after outlier exclusion
    per_marker: list[MarkerResult]
    error_rate: float
    rho: float = float("inf")  # beta-binomial concentration; inf = no overdispersion
    # Per-sample analytical detection limits (donor fractions, 0.0-1.0), computed
    # from the Fisher information of this sample's own markers. inf = nothing detectable.
    lob_fraction: float = float("inf")  # limit of blank
    lod_fraction: float = float("inf")  # limit of detection
    # Host-presence detector output (see ``allomix.detect``). None when the
    # caller disabled the detector or when there were no donor-homozygous
    # markers to run it on.
    host_presence: HostPresenceResult | None = None
    # Robust-refit accounting (see ``estimate_single_donor_bb`` robust mode).
    # n_robust_excluded is the count of markers dropped as residual outliers and
    # excluded from the final fit; robust_drop_fraction is that count over
    # n_informative. Both are 0 when robust mode is off or nothing was dropped.
    n_robust_excluded: int = 0
    robust_drop_fraction: float = 0.0
    # Identity QC (see ``allomix.relatedness``), attached by ``analyse_sample``.
    # relatedness holds one entry per reference-sample pair (host vs each donor);
    # admix_consistency is the consensus-homozygote swap check. None when the
    # caller did not compute them.
    relatedness: list[RelatednessResult] | None = None
    admix_consistency: AdmixConsistencyResult | None = None


@dataclass
class MultiDonorResult:
    """Result of multi-donor chimerism estimation."""

    donor_fractions: list[float]  # [f_donor1, f_donor2, ...]
    donor_fraction_cis: list[tuple[float, float]]  # [(lo, hi), ...] per donor
    host_fraction: float  # 1 - sum(donor_fractions)
    log_likelihood: float
    n_informative: int
    n_markers_used: int
    per_marker: list[MarkerResult]
    error_rate: float
    per_donor_n_informative: list[int] | None = None  # informative markers per donor
    rho: float = float("inf")  # beta-binomial concentration; inf = no overdispersion
    host_presence: HostPresenceResult | None = None
    n_robust_excluded: int = 0
    robust_drop_fraction: float = 0.0
    # Identity QC; see ChimerismResult above. For multi-donor, relatedness also
    # includes donor-vs-donor pairs.
    relatedness: list[RelatednessResult] | None = None
    admix_consistency: AdmixConsistencyResult | None = None


# Robust-refit tuning. ROBUST_K is the median/MAD residual cut (robust SDs);
# 3.5 leaves clean data essentially untouched (drops <1% of markers by chance)
# while removing copy-number/LoH-inconsistent markers. The refit floors keep
# trimming from gutting sparse panels: "auto" never drops below
# ROBUST_MIN_MARKERS, "force" never below ROBUST_HARD_MIN.
ROBUST_K_DEFAULT = 3.5
ROBUST_MAX_ITER = 5
ROBUST_MIN_MARKERS = 15
ROBUST_HARD_MIN = 4
ROBUST_MODES = ("off", "auto", "force")
# "auto" engages the refit only when the first pass finds more residual outliers
# than this trigger: max(ROBUST_TRIGGER_MIN, ceil(ROBUST_TRIGGER_FRAC * n)). At
# k=3.5 a clean panel produces <1 chance outlier on average, so a trigger of ~3
# leaves clean (and the already-validated) samples byte-identical, while
# genuine copy-number / LoH clusters clear the bar. "force" uses a trigger of 1.
ROBUST_TRIGGER_FRAC = 0.03
ROBUST_TRIGGER_MIN = 3


# ---------------------------------------------------------------------------
# Core likelihood functions
# ---------------------------------------------------------------------------


def expected_weight(
    host_gt: tuple[int, int],
    donor_gt: tuple[int, int],
    f_donor: float,
    bias: float = 0.0,
) -> float:
    """Expected reference allele weight for a given chimerism fraction.

    w = (1 - f) * host_ref_dose / 2 + f * donor_ref_dose / 2

    where ref_dose = 2 - alt_dose.

    When ``bias`` is non-zero, the weight is adjusted to account for
    per-marker amplification bias.  A positive bias means the ALT allele
    is preferentially captured, so the observed REF weight is lower:

        w_corrected = w_true - bias   (clamped to [eps, 1 - eps])

    Args:
        host_gt: Host diploid genotype, e.g. (0, 0), (0, 1), (1, 1).
        donor_gt: Donor diploid genotype.
        f_donor: Donor fraction (0.0 to 1.0).
        bias: Per-marker amplification bias (default 0.0 = no correction).

    Returns:
        Expected reference allele weight (0.0 to 1.0).
    """
    host_ref_dose = 2 - (host_gt[0] + host_gt[1])
    donor_ref_dose = 2 - (donor_gt[0] + donor_gt[1])
    w = (1.0 - f_donor) * host_ref_dose / 2.0 + f_donor * donor_ref_dose / 2.0
    if bias != 0.0:
        # Clamp to avoid 0/1 boundary (log-likelihood needs p > 0)
        w = max(1e-6, min(1.0 - 1e-6, w - bias))
    return w


# DNA has 4 bases, so a sequencing error changes the true base into one of the
# 3 other bases. Assuming errors are spread evenly, a miscall to one specific
# base (e.g. REF read as the ALT allele) has probability error_rate / 3.
_N_OTHER_BASES = 3


def alt_read_probability(w: float, error_rate: float = 0.01) -> float:
    """Probability an observed read is ALT, given expected REF weight ``w``.

    Under the 4-state error model a read is observed as ALT if it comes from a
    true ALT allele and is called correctly (probability ``1 - e``), or from a
    true REF allele miscalled to the ALT base (probability ``e / 3``). The two
    raw probabilities are renormalised so they condition on the read being
    called REF or ALT rather than one of the other two bases.

    Args:
        w: Expected reference allele weight (fraction of REF alleles).
        error_rate: Per-base sequencing error rate ``e``.

    Returns:
        P(observe ALT | w), between 0 and 1.
    """
    e = error_rate
    e_specific = e / _N_OTHER_BASES
    p_alt = (1.0 - w) * (1.0 - e) + w * e_specific
    p_ref = w * (1.0 - e) + (1.0 - w) * e_specific
    return p_alt / (p_ref + p_alt)


def log_likelihood_marker_bb(
    ad_ref: int,
    ad_alt: int,
    w: float,
    error_rate: float = 0.01,
    rho: float = 100.0,
    e_refalt: float | None = None,
    e_altref: float | None = None,
) -> float:
    """Per-marker log-likelihood under a beta-binomial model.

    When ``e_refalt`` and ``e_altref`` are both supplied, uses the asymmetric
    REF/ALT-only error model

        p_alt = w * e_refalt + (1 - w) * (1 - e_altref)

    where the rates are per-direction empirical substitution probabilities
    (typically from ``error_rates.estimate_error_rates``). Either rate may be
    ``None`` individually; in that case the legacy 4-state symmetric model
    with rate ``error_rate`` is used, matching the prior behaviour.

    When rho -> inf this converges to the binomial log-likelihood.

    Args:
        ad_ref: Reference allele read count.
        ad_alt: Alternative allele read count.
        w: Expected reference allele weight.
        error_rate: Symmetric 4-state rate, used only when asymmetric rates
            are not both provided (default 0.01).
        rho: Beta-binomial concentration parameter. Larger = less
            overdispersion. Typical empirical values: 50-500.
        e_refalt: ``P(observe ALT | true REF)`` for this marker.
        e_altref: ``P(observe REF | true ALT)`` for this marker.

    Returns:
        Log-likelihood contribution from this marker.
    """
    if e_refalt is not None and e_altref is not None:
        # Asymmetric REF/ALT-only model. p_alt + p_ref = 1 by construction.
        p_alt = w * e_refalt + (1.0 - w) * (1.0 - e_altref)
    else:
        p_alt = alt_read_probability(w, error_rate)
    p_alt = max(1e-6, min(1.0 - 1e-6, p_alt))

    n = ad_ref + ad_alt
    k = ad_alt

    if n == 0:
        return 0.0

    a = p_alt * rho
    b = (1.0 - p_alt) * rho

    # Clamp to avoid lgamma(0) = inf
    a = max(a, 1e-10)
    b = max(b, 1e-10)

    # log P(k | n, a, b) dropping the constant log C(n,k)
    ll = lgamma(k + a) + lgamma(n - k + b) - lgamma(n + rho) - lgamma(a) - lgamma(b) + lgamma(rho)
    return ll


def total_log_likelihood_bb(
    markers: list[InformativeMarker],
    f_donor: float,
    error_rate: float = 0.01,
    rho: float = 100.0,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
    marker_errors: (
        dict[tuple[str, int, str, str], tuple[float | None, float | None]] | None
    ) = None,
) -> float:
    """Sum of per-marker beta-binomial log-likelihoods.

    Args:
        markers: List of informative markers with admixture allele counts.
        f_donor: Donor fraction to evaluate.
        error_rate: Sequencing error rate (used as the fallback when a marker
            is missing from ``marker_errors`` or only one direction is known).
        rho: Beta-binomial concentration parameter.
        marker_biases: Optional per-marker amplification bias dict.
        marker_errors: Optional per-marker, per-direction empirical error
            rates. Maps ``(chrom, pos, ref, alt)`` to
            ``(e_refalt, e_altref)``. ``None`` in either slot falls through
            to the symmetric ``error_rate`` model for that marker.

    Returns:
        Total log-likelihood.
    """
    if not markers:
        return 0.0
    arr = _precompute_marker_arrays(markers, marker_biases, marker_errors)
    return _total_ll_vec(arr, f_donor, error_rate, rho)


@dataclass(frozen=True)
class _MarkerArrays:
    """Per-marker quantities that do not depend on (f_donor, rho).

    Precomputed once per marker set so the likelihood can be evaluated as a
    single vectorized expression for any (f_donor, rho).
    """

    host_ref_dose: np.ndarray  # 2 - host alt dose, float
    donor_ref_dose: np.ndarray  # 2 - donor[0] alt dose, float
    n: np.ndarray  # ad_ref + ad_alt, float
    k: np.ndarray  # ad_alt, float
    bias: np.ndarray  # per-marker bias, float (0.0 where none)
    bias_mask: np.ndarray  # bias != 0.0, bool
    # Per-marker asymmetric error rates. NaN where the asymmetric model does
    # not apply (either rate missing); the vectorised LL falls back to the
    # symmetric 4-state model at those markers via ``error_mask``.
    e_refalt: np.ndarray
    e_altref: np.ndarray
    error_mask: np.ndarray  # True where both per-direction rates are known


def _precompute_marker_arrays(
    markers: list[InformativeMarker],
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
    marker_errors: (
        dict[tuple[str, int, str, str], tuple[float | None, float | None]] | None
    ) = None,
) -> _MarkerArrays:
    """Build the (f, rho)-independent per-marker arrays for the vectorized LL.

    Uses the first donor genotype (single-donor model), matching
    ``total_log_likelihood_bb``.

    Args:
        markers: List of informative markers with admixture allele counts.
        marker_biases: Optional per-marker amplification bias dict.
        marker_errors: Optional per-marker, per-direction error-rate dict
            (``(e_refalt, e_altref)`` per key, either entry may be ``None``).

    Returns:
        Precomputed per-marker arrays for ``_total_ll_vec``.
    """
    n_markers = len(markers)
    host_ref_dose = np.fromiter(
        (2 - (m.host_gt[0] + m.host_gt[1]) for m in markers),
        dtype=float,
        count=n_markers,
    )
    donor_ref_dose = np.fromiter(
        (2 - (m.donor_gts[0][0] + m.donor_gts[0][1]) for m in markers),
        dtype=float,
        count=n_markers,
    )
    ad_ref = np.fromiter((m.admix_ad_ref for m in markers), dtype=float, count=n_markers)
    ad_alt = np.fromiter((m.admix_ad_alt for m in markers), dtype=float, count=n_markers)
    if marker_biases is not None:
        bias = np.fromiter(
            (marker_biases.get((m.chrom, m.pos, m.ref, m.alt), 0.0) for m in markers),
            dtype=float,
            count=n_markers,
        )
    else:
        bias = np.zeros(n_markers, dtype=float)

    e_refalt = np.full(n_markers, np.nan, dtype=float)
    e_altref = np.full(n_markers, np.nan, dtype=float)
    if marker_errors is not None:
        for i, m in enumerate(markers):
            entry = marker_errors.get((m.chrom, m.pos, m.ref, m.alt))
            if entry is None:
                continue
            e_ra, e_ar = entry
            # Both directions required for the asymmetric path. Storing one
            # only would leave the other side of the likelihood underspecified
            # at hets and at the opposite-homozygous endpoint.
            if e_ra is None or e_ar is None:
                continue
            e_refalt[i] = e_ra
            e_altref[i] = e_ar
    error_mask = ~np.isnan(e_refalt)

    return _MarkerArrays(
        host_ref_dose=host_ref_dose,
        donor_ref_dose=donor_ref_dose,
        n=ad_ref + ad_alt,
        k=ad_alt,
        bias=bias,
        bias_mask=bias != 0.0,
        e_refalt=e_refalt,
        e_altref=e_altref,
        error_mask=error_mask,
    )


def _total_ll_vec(
    arr: _MarkerArrays,
    f_donor: float,
    error_rate: float = 0.01,
    rho: float = 100.0,
) -> float:
    """Vectorized single-donor beta-binomial total log-likelihood.

    Numerically equivalent to ``total_log_likelihood_bb`` (differences only at
    the ``gammaln`` vs ``math.lgamma`` rounding level, ~1e-9). Markers with
    ``n == 0`` contribute 0 automatically.

    Args:
        arr: Per-marker arrays from ``_precompute_marker_arrays``.
        f_donor: Donor fraction to evaluate.
        error_rate: Sequencing error rate.
        rho: Beta-binomial concentration parameter.

    Returns:
        Total log-likelihood.
    """
    # Expected reference-allele weight, with the conditional bias clamp.
    w = (1.0 - f_donor) * arr.host_ref_dose / 2.0 + f_donor * arr.donor_ref_dose / 2.0
    if arr.bias_mask.any():
        w[arr.bias_mask] = np.clip(w[arr.bias_mask] - arr.bias[arr.bias_mask], 1e-6, 1.0 - 1e-6)

    # P(observe ALT | w). Default is the 4-state symmetric model with the
    # global ``error_rate``; per-marker asymmetric rates (where supplied) use
    # the REF/ALT-only form ``p_alt = w * e_refalt + (1 - w) * (1 - e_altref)``.
    e = error_rate
    e_specific = e / _N_OTHER_BASES
    p_alt_raw = (1.0 - w) * (1.0 - e) + w * e_specific
    p_ref_raw = w * (1.0 - e) + (1.0 - w) * e_specific
    p_alt = p_alt_raw / (p_ref_raw + p_alt_raw)
    if arr.error_mask.any():
        # Sub in the asymmetric per-marker rates at the masked positions.
        p_alt_asym = w * arr.e_refalt + (1.0 - w) * (1.0 - arr.e_altref)
        p_alt = np.where(arr.error_mask, p_alt_asym, p_alt)
    p_alt = np.clip(p_alt, 1e-6, 1.0 - 1e-6)

    a = np.maximum(p_alt * rho, 1e-10)
    b = np.maximum((1.0 - p_alt) * rho, 1e-10)
    n, k = arr.n, arr.k

    ll = (
        gammaln(k + a)
        + gammaln(n - k + b)
        - gammaln(n + rho)
        - gammaln(a)
        - gammaln(b)
        + gammaln(rho)
    )
    return float(ll.sum())


def total_log_likelihood_multi_bb(
    markers: list[InformativeMarker],
    donor_fractions: list[float],
    error_rate: float = 0.01,
    rho: float = 100.0,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
    marker_errors: (
        dict[tuple[str, int, str, str], tuple[float | None, float | None]] | None
    ) = None,
) -> float:
    """Sum of per-marker beta-binomial log-likelihoods for multi-donor model.

    Args:
        markers: Informative markers (for at least one donor).
        donor_fractions: [f_donor1, f_donor2, ...].
        error_rate: Sequencing error rate (fallback when a marker's
            per-direction rate is missing).
        rho: Beta-binomial concentration parameter.
        marker_biases: Optional per-marker bias dict.
        marker_errors: Optional per-marker, per-direction error-rate dict
            (see ``total_log_likelihood_bb``).

    Returns:
        Total log-likelihood.
    """
    ll = 0.0
    for m in markers:
        key = (m.chrom, m.pos, m.ref, m.alt)
        bias = marker_biases.get(key, 0.0) if marker_biases is not None else 0.0
        e_ra: float | None = None
        e_ar: float | None = None
        if marker_errors is not None:
            entry = marker_errors.get(key)
            if entry is not None:
                e_ra, e_ar = entry
        w = expected_weight_multi(m.host_gt, m.donor_gts, donor_fractions, bias=bias)
        ll += log_likelihood_marker_bb(
            m.admix_ad_ref, m.admix_ad_alt, w,
            error_rate=error_rate,
            rho=rho,
            e_refalt=e_ra,
            e_altref=e_ar,
        )
    return ll


def expected_weight_multi(
    host_gt: tuple[int, int],
    donor_gts: list[tuple[int, int]],
    donor_fractions: list[float],
    bias: float = 0.0,
) -> float:
    """Expected reference allele weight for multi-donor chimerism.

    w = (1 - f1 - f2) * host_ref_dose/2 + f1 * d1_ref_dose/2 + f2 * d2_ref_dose/2

    Args:
        host_gt: Host diploid genotype.
        donor_gts: List of donor diploid genotypes.
        donor_fractions: List of donor fractions (sum <= 1.0).
        bias: Per-marker amplification bias.

    Returns:
        Expected reference allele weight (0.0 to 1.0).
    """
    host_ref_dose = 2 - (host_gt[0] + host_gt[1])
    f_host = 1.0 - sum(donor_fractions)
    w = f_host * host_ref_dose / 2.0
    for dgt, f in zip(donor_gts, donor_fractions):
        d_ref_dose = 2 - (dgt[0] + dgt[1])
        w += f * d_ref_dose / 2.0
    if bias != 0.0:
        w = max(1e-6, min(1.0 - 1e-6, w - bias))
    return w


# ---------------------------------------------------------------------------
# MLE estimation
# ---------------------------------------------------------------------------


def _compute_per_marker_results(
    markers: list[InformativeMarker],
    f_mle: float,
    marker_biases: dict[tuple[str, int, str, str], float] | None,
) -> list[MarkerResult]:
    """Compute per-marker residuals and flag outliers.

    Shared between binomial and beta-binomial estimators.
    """
    per_marker: list[MarkerResult] = []
    residuals: list[float] = []

    for m in markers:
        bias = 0.0
        if marker_biases is not None:
            bias = marker_biases.get((m.chrom, m.pos, m.ref, m.alt), 0.0)
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

    # Flag outliers (residual > 3 SD)
    if len(residuals) >= 2:
        r_arr = np.array(residuals)
        mean_r = float(np.mean(r_arr))
        sd_r = float(np.std(r_arr, ddof=1))

        if sd_r > 0:
            for i, mr in enumerate(per_marker):
                if abs(residuals[i] - mean_r) > 3.0 * sd_r:
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


# One-sided 95% normal quantile (z_0.95 ~= 1.6449), used for EP17-style LoB/LoD.
_Z95 = float(norm.ppf(0.95))

# Margin used to keep a probability strictly inside the open interval (0, 1),
# so that p * (1 - p) stays positive and the marker variance never collapses to
# zero. This is a safety clamp, not machine epsilon (np.finfo(float).eps).
_PROB_EPS = 1e-9


def fraction_se(
    markers: list[InformativeMarker],
    f_donor: float,
    error_rate: float = 0.01,
    rho: float = float("inf"),
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
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
        marker_biases: Optional per-marker amplification bias dict.

    Returns:
        Standard error of the donor fraction. inf if no marker is informative.
    """
    # P(observe ALT) is linear in the REF weight w, so its slope dp_alt/dw is
    # constant and equals the change across the full weight range w: 0 -> 1.
    dpalt_dw = alt_read_probability(1.0, error_rate) - alt_read_probability(0.0, error_rate)

    info = 0.0
    for m in markers:
        bias = 0.0
        if marker_biases is not None:
            bias = marker_biases.get((m.chrom, m.pos, m.ref, m.alt), 0.0)

        host_ref_dose = 2 - (m.host_gt[0] + m.host_gt[1])
        donor_ref_dose = 2 - (m.donor_gts[0][0] + m.donor_gts[0][1])
        dw_df = (donor_ref_dose - host_ref_dose) / 2.0
        if dw_df == 0.0:
            continue

        n = m.admix_ad_ref + m.admix_ad_alt
        if n == 0:
            continue

        w = expected_weight(m.host_gt, m.donor_gts[0], f_donor, bias=bias)
        p_alt = alt_read_probability(w, error_rate)
        p_alt = max(_PROB_EPS, min(1.0 - _PROB_EPS, p_alt))

        overdispersion = 1.0 if math.isinf(rho) else 1.0 + (n - 1.0) / (rho + 1.0)
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
    error_rate: float = 0.01,
    rho: float = float("inf"),
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
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
        marker_biases: Optional per-marker amplification bias dict.

    Returns:
        ``(lob, lod)`` as donor fractions (0.0-1.0). ``(inf, inf)`` if no
        marker is informative.
    """
    se0 = fraction_se(markers, 0.0, error_rate, rho, marker_biases)
    if math.isinf(se0):
        return float("inf"), float("inf")
    lob = _Z95 * se0
    se_lob = fraction_se(markers, lob, error_rate, rho, marker_biases)
    if math.isinf(se_lob):
        se_lob = se0
    lod = lob + _Z95 * se_lob
    return lob, lod


def _estimate_single_donor_bb_core(
    markers: list[InformativeMarker],
    error_rate: float = 0.01,
    grid_steps: int = 1001,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
    marker_errors: (
        dict[tuple[str, int, str, str], tuple[float | None, float | None]] | None
    ) = None,
) -> ChimerismResult:
    """Single-donor beta-binomial MLE over the given marker set (no robust trim).

    This is the unguarded estimator; ``estimate_single_donor_bb`` wraps it with
    the optional robust refit. Args/returns as in that wrapper.
    """
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

    # Precompute the (f, rho)-independent per-marker arrays once and reuse them
    # across the grid search, Nelder-Mead refinement, and profile-likelihood CI.
    arr = _precompute_marker_arrays(markers, marker_biases, marker_errors)

    # Step 1: Grid search over f with rho profiled out at each grid point
    grid = np.linspace(0.0, 1.0, grid_steps)
    best_ll = -math.inf
    best_f = 0.0
    best_rho = 100.0

    for f in grid:
        # Optimise rho for this f
        opt_rho = minimize_scalar(
            lambda log_r, _f=f: -_total_ll_vec(arr, _f, error_rate, math.exp(log_r)),
            bounds=(math.log(1.0), math.log(10000.0)),
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
            return 1e30
        rho_val = math.exp(log_rho_val)
        if rho_val < 0.5 or rho_val > 50000:
            return 1e30
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
    threshold = chi2.ppf(0.95, df=1)
    half_threshold = threshold / 2.0

    def profile_ll_f(f_val: float) -> float:
        """Max LL over rho at a given f."""
        opt_rho = minimize_scalar(
            lambda log_r: -_total_ll_vec(arr, f_val, error_rate, math.exp(log_r)),
            bounds=(math.log(1.0), math.log(50000.0)),
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
    per_marker = _compute_per_marker_results(markers, f_mle, marker_biases)
    n_markers_used = sum(1 for mr in per_marker if mr.included)

    # Step 5: Per-sample analytical detection limits from the fitted noise model.
    lob, lod = detection_limit(markers, error_rate, rho_mle, marker_biases)

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
        mad = float(np.median(np.abs(resids - med))) * 1.4826
        if mad <= 0.0:
            break
        keep_mask = np.abs(resids - med) <= robust_k * mad
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
        replace(mr, included=_marker_key(m) in surviving)
        for m, mr in zip(all_markers, pm_all)
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
    error_rate: float = 0.01,
    grid_steps: int = 1001,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
    marker_errors: (
        dict[tuple[str, int, str, str], tuple[float | None, float | None]] | None
    ) = None,
    robust: str = "off",
    robust_k: float = ROBUST_K_DEFAULT,
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
        marker_biases: Optional per-marker bias dict.
        marker_errors: Optional per-marker, per-direction empirical error
            rates (``(e_refalt, e_altref)`` per key, either may be ``None``).
            When present and both directions are known, the asymmetric
            REF/ALT-only likelihood is used at that marker; otherwise the
            symmetric 4-state model with ``error_rate`` is used.
        robust: Robust-refit mode. ``"off"`` (default) is the plain MLE.
            ``"auto"`` iteratively drops median/MAD residual outliers and refits,
            never below ``ROBUST_MIN_MARKERS`` survivors; this protects against
            host copy-number / LoH markers (see ``_robust_refit``). ``"force"``
            is the same but trims down to ``ROBUST_HARD_MIN`` (for experiments).
            On clean data the trim is a no-op and the result is unchanged.
        robust_k: Residual cut in robust SDs for the refit.

    Returns:
        ChimerismResult with MLE estimate and beta-binomial CIs.
    """
    if robust not in ROBUST_MODES:
        raise ValueError(f"robust must be one of {ROBUST_MODES}, got {robust!r}")

    def core(mk: list[InformativeMarker]) -> ChimerismResult:
        return _estimate_single_donor_bb_core(
            mk, error_rate, grid_steps, marker_biases, marker_errors
        )

    if robust == "off" or len(markers) == 0:
        return core(markers)

    min_markers = ROBUST_MIN_MARKERS if robust == "auto" else ROBUST_HARD_MIN
    min_trigger = _robust_trigger(len(markers)) if robust == "auto" else 1
    return _robust_refit(
        markers,
        core,
        lambda mk, res: _compute_per_marker_results(mk, res.donor_fraction, marker_biases),
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
    error_rate: float = 0.01,
    grid_steps: int = 101,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
    marker_errors: (
        dict[tuple[str, int, str, str], tuple[float | None, float | None]] | None
    ) = None,
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
        marker_biases: Optional per-marker bias dict.
        marker_errors: Optional per-marker, per-direction empirical error
            rates (see ``estimate_single_donor_bb``).

    Returns:
        MultiDonorResult with per-donor fractions and CIs.
    """
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
                markers, [f1, f2], error_rate, rho_init, marker_biases,
                marker_errors,
            )
            if ll > best_ll:
                best_ll = ll
                best_f = [f1, f2]

    # Step 2: Nelder-Mead refinement over (f1, f2, log_rho)
    def neg_ll(x):
        f1, f2, log_rho = x
        if f1 < 0 or f2 < 0 or f1 + f2 > 1.0:
            return 1e30
        rho = math.exp(log_rho)
        if rho < 0.5 or rho > 50000:
            return 1e30
        return -total_log_likelihood_multi_bb(
            markers, [f1, f2], error_rate, rho, marker_biases, marker_errors,
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
        markers, f_mle, n_donors, error_rate, marker_biases, marker_errors,
    )

    # Step 4: Per-marker residuals
    per_marker = _per_marker_results_multi(markers, f_mle, marker_biases)

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
    error_rate: float = 0.01,
    grid_steps: int = 101,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
    marker_errors: (
        dict[tuple[str, int, str, str], tuple[float | None, float | None]] | None
    ) = None,
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

    def core(mk: list[InformativeMarker]) -> MultiDonorResult:
        return _estimate_multi_donor_core(
            mk, n_donors, error_rate, grid_steps, marker_biases, marker_errors
        )

    if robust == "off" or len(markers) == 0:
        return core(markers)

    min_markers = ROBUST_MIN_MARKERS if robust == "auto" else ROBUST_HARD_MIN
    min_trigger = _robust_trigger(len(markers)) if robust == "auto" else 1
    return _robust_refit(
        markers,
        core,
        lambda mk, res: _per_marker_results_multi(mk, res.donor_fractions, marker_biases),
        robust_k,
        min_markers,
        min_trigger,
    )


def _profile_likelihood_cis_multi(
    markers: list[InformativeMarker],
    f_mle: list[float],
    n_donors: int,
    error_rate: float,
    marker_biases: dict[tuple[str, int, str, str], float] | None,
    marker_errors: (
        dict[tuple[str, int, str, str], tuple[float | None, float | None]] | None
    ) = None,
) -> list[tuple[float, float]]:
    """Profile likelihood CIs for each donor fraction.

    For donor_i, scan f_i while optimizing f_j (j != i) and rho at each
    point. Uses chi2(df=1) threshold since we profile one parameter at
    a time. Overdispersion is handled by the beta-binomial likelihood
    (via rho profiling) rather than by inflating the threshold.
    """
    threshold = float(chi2.ppf(0.95, df=1))
    half_threshold = threshold / 2.0
    cis: list[tuple[float, float]] = []

    for donor_idx in range(n_donors):
        other_idx = 1 - donor_idx  # works for 2 donors

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
                            markers, fracs, error_rate, math.exp(log_r),
                            marker_biases, marker_errors,
                        )
                    ),
                    bounds=(math.log(1.0), math.log(50000.0)),
                    method="bounded",
                )
                return -float(opt_rho.fun)

            # Optimise over (fj, log_rho) jointly
            def neg_ll_inner(x, _di=_didx):
                fj, log_r = x
                if fj < 0 or fj > max_fj:
                    return 1e30
                rho = math.exp(log_r)
                if rho < 0.5 or rho > 50000:
                    return 1e30
                fracs = [fi, fj] if _di == 0 else [fj, fi]
                return -total_log_likelihood_multi_bb(
                    markers, fracs, error_rate, rho, marker_biases,
                    marker_errors,
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
    marker_biases: dict[tuple[str, int, str, str], float] | None,
) -> list[MarkerResult]:
    """Compute per-marker residuals for multi-donor model."""
    per_marker: list[MarkerResult] = []
    residuals: list[float] = []

    for m in markers:
        bias = 0.0
        if marker_biases is not None:
            bias = marker_biases.get((m.chrom, m.pos, m.ref, m.alt), 0.0)
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

    # Flag outliers (residual > 3 SD)
    if len(residuals) >= 2:
        r_arr = np.array(residuals)
        mean_r = float(np.mean(r_arr))
        sd_r = float(np.std(r_arr, ddof=1))

        if sd_r > 0:
            for i, mr in enumerate(per_marker):
                if abs(residuals[i] - mean_r) > 3.0 * sd_r:
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
