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
overdispersion is present, with the largest improvement in CI coverage at
low donor fractions where accurate CIs matter most clinically. Per-marker
amplification bias (het-site VAF deviation, SD ~0.02 empirically) is
corrected multiplicatively in logit space (see ``apply_bias``). The in
silico CI-coverage and point-estimate accuracy across depths and noise
conditions are reported in the paper's depth and bias-correction
validations.
"""

import math
from dataclasses import dataclass, field, replace
from math import lgamma

import numpy as np
from scipy.optimize import brentq, minimize, minimize_scalar
from scipy.special import expit, gammaln, logit
from scipy.stats import chi2, norm

from allomix.constants import (
    CI_LEVEL,
    DEFAULT_ERROR_RATE,
    N_OTHER_BASES,
    ROBUST_K_DEFAULT,
)
from allomix.contamination import ContaminationResult
from allomix.detect import HostPresenceResult
from allomix.error_rates import MarkerErrorRates
from allomix.genotype import InformativeMarker, MarkerKey
from allomix.relatedness import AdmixConsistencyResult, RelatednessResult
from allomix.runmeta import RunUnitInfo


@dataclass(frozen=True)
class PanelCalibration:
    """Per-marker calibration applied during chimerism estimation.

    Bundles the two optional per-marker tables the estimators consume:

    - ``biases``: amplification bias per marker (see ``allomix.bias``). A
      positive value means the ALT allele is preferentially captured.
    - ``errors``: per-direction empirical substitution rates per marker (see
      ``allomix.error_rates``). Used for the asymmetric REF/ALT-only likelihood
      where both directions are known.

    Both default to empty, so an uncalibrated run is ``PanelCalibration()``.
    Markers absent from a table fall through to the uncorrected weight (bias 0)
    or the symmetric global ``error_rate`` (errors).
    """

    biases: dict[MarkerKey, float] = field(default_factory=dict)
    errors: dict[MarkerKey, MarkerErrorRates] = field(default_factory=dict)

    def __post_init__(self) -> None:
        # Treat an explicit ``None`` (a caller's "no table" sentinel) the same as
        # an omitted argument, so callers can pass an optional dict straight
        # through without a local ``or {}`` guard.
        if self.biases is None:
            object.__setattr__(self, "biases", {})
        if self.errors is None:
            object.__setattr__(self, "errors", {})

    def bias_for(self, m: InformativeMarker) -> float:
        """Amplification bias for a marker, or 0.0 if it is not in the table."""
        return self.biases.get((m.chrom, m.pos, m.ref, m.alt), 0.0)

    def error_for(self, m: InformativeMarker) -> MarkerErrorRates | None:
        """Per-direction rates for a marker, or ``None`` if not in the table."""
        return self.errors.get((m.chrom, m.pos, m.ref, m.alt))


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
    # In-data contamination estimate (see ``allomix.contamination``), attached by
    # ``analyse_sample``. None when the caller did not compute it.
    contamination: ContaminationResult | None = None
    # Sequencing run-unit metadata for this admix sample, read from the admix VCF
    # header (see ``allomix.runmeta``). None when the VCF carried no run metadata.
    run_unit: RunUnitInfo | None = None


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
    contamination: ContaminationResult | None = None
    run_unit: RunUnitInfo | None = None


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
# Core likelihood functions
# ---------------------------------------------------------------------------


# Clamp keeping expected weights off the 0/1 boundary (log-likelihood needs p > 0).
W_EPS = 1e-6


def apply_bias(w: float | np.ndarray, bias: float | np.ndarray) -> float | np.ndarray:
    """Correct an expected REF weight for per-marker amplification bias.

    ``bias`` is the het-site convention used throughout: ``median(observed het
    ALT VAF) - 0.5`` (positive = ALT preferentially captured). It is estimated
    at heterozygous sites, where the expected REF weight is 0.5. Applying it as
    a flat additive shift (``w - bias``) is only valid near 0.5; at informative
    markers whose expected VAF is near 0 or 1 (the norm at low chimerism) it
    overcorrects badly (issue #20). Instead apply it multiplicatively, in logit
    space, where it is valid at any expected VAF:

        w_corrected = expit(logit(w) - logit(0.5 + bias))

    At a het site (w = 0.5) this reduces to ``0.5 - bias`` (ALT VAF 0.5 + bias),
    matching the estimate; at an extreme expected VAF it is only a small
    multiplicative nudge rather than a large additive jump. Works for scalars
    and numpy arrays. The result is clamped to ``[W_EPS, 1 - W_EPS]``.

    Args:
        w: Expected reference allele weight(s) before correction.
        bias: Per-marker bias(es) in het-site VAF units.

    Returns:
        Bias-corrected reference allele weight(s).
    """
    p = np.clip(0.5 + bias, W_EPS, 1.0 - W_EPS)  # observed het ALT-favouring as a probability
    w_clamped = np.clip(w, W_EPS, 1.0 - W_EPS)
    return np.clip(expit(logit(w_clamped) - logit(p)), W_EPS, 1.0 - W_EPS)


def inject_bias(alt_vaf: float | np.ndarray, bias: float | np.ndarray) -> float | np.ndarray:
    """Shift a true ALT VAF by a het-site bias, the simulator-side counterpart
    of ``apply_bias``.

    Used by the simulator to inject per-marker amplification bias the same way
    the estimator corrects for it, so simulation and estimation stay
    self-consistent: at the true parameters the injected ALT VAF equals the
    estimator's expected biased ALT VAF, ``1 - apply_bias(w_true, bias)``.
    Equivalent to ``expit(logit(alt_vaf) + logit(0.5 + bias))``; at a het site
    (true VAF 0.5) the observed VAF becomes ``0.5 + bias``. (Note this is not an
    algebraic inverse of ``apply_bias``: both shift in the same direction, so
    composing them would double-correct.)

    Args:
        alt_vaf: True ALT-allele VAF(s) before bias.
        bias: Per-marker bias(es) in het-site VAF units.

    Returns:
        Biased ALT VAF(s), clamped to ``[W_EPS, 1 - W_EPS]``.
    """
    return 1.0 - apply_bias(1.0 - np.asarray(alt_vaf, dtype=float), bias)


def expected_weight(
    host_gt: tuple[int, int],
    donor_gt: tuple[int, int],
    f_donor: float,
    bias: float = 0.0,
) -> float:
    """Expected reference allele weight for a given chimerism fraction.

    w = (1 - f) * host_ref_dose / 2 + f * donor_ref_dose / 2

    where ref_dose = 2 - alt_dose.

    When ``bias`` is non-zero, the weight is corrected for per-marker
    amplification bias in logit space (see ``apply_bias``).

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
        w = float(apply_bias(w, bias))
    return w


def alt_read_probability(w: float, error_rate: float = DEFAULT_ERROR_RATE) -> float:
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
    e_specific = e / N_OTHER_BASES
    p_alt = (1.0 - w) * (1.0 - e) + w * e_specific
    p_ref = w * (1.0 - e) + (1.0 - w) * e_specific
    return p_alt / (p_ref + p_alt)


def log_likelihood_marker_bb(
    ad_ref: int,
    ad_alt: int,
    w: float,
    error_rate: float = DEFAULT_ERROR_RATE,
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
    error_rate: float = DEFAULT_ERROR_RATE,
    rho: float = 100.0,
    calibration: PanelCalibration | None = None,
) -> float:
    """Sum of per-marker beta-binomial log-likelihoods.

    Args:
        markers: List of informative markers with admixture allele counts.
        f_donor: Donor fraction to evaluate.
        error_rate: Sequencing error rate (used as the fallback when a marker
            is missing from ``calibration.errors`` or only one direction is
            known).
        rho: Beta-binomial concentration parameter.
        calibration: Optional per-marker bias and error tables. A marker's
            per-direction rates are used (asymmetric model) only when both are
            known; otherwise the symmetric ``error_rate`` model applies.

    Returns:
        Total log-likelihood.
    """
    if not markers:
        return 0.0
    arr = _precompute_marker_arrays(markers, calibration or PanelCalibration())
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
    # (f, rho)-independent precomputations, hoisted out of the per-evaluation
    # hot path. ``has_bias``/``has_error`` cache the per-call ``.any()`` checks;
    # ``logit_bias_masked`` is ``logit(clip(0.5 + bias))`` at the biased markers
    # (the only bias-dependent term in ``apply_bias``), so the likelihood loop
    # never recomputes it. All exact: same values the per-call code produced.
    has_bias: bool
    has_error: bool
    logit_bias_masked: np.ndarray


def _precompute_marker_arrays(
    markers: list[InformativeMarker],
    calibration: PanelCalibration,
) -> _MarkerArrays:
    """Build the (f, rho)-independent per-marker arrays for the vectorized LL.

    Uses the first donor genotype (single-donor model), matching
    ``total_log_likelihood_bb``.

    Args:
        markers: List of informative markers with admixture allele counts.
        calibration: Per-marker bias and error tables (empty tables = no
            correction).

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
    if calibration.biases:
        bias = np.fromiter(
            (calibration.bias_for(m) for m in markers),
            dtype=float,
            count=n_markers,
        )
    else:
        bias = np.zeros(n_markers, dtype=float)

    e_refalt = np.full(n_markers, np.nan, dtype=float)
    e_altref = np.full(n_markers, np.nan, dtype=float)
    if calibration.errors:
        for i, m in enumerate(markers):
            entry = calibration.error_for(m)
            if entry is None:
                continue
            # Both directions required for the asymmetric path. Storing one
            # only would leave the other side of the likelihood underspecified
            # at hets and at the opposite-homozygous endpoint.
            if entry.e_refalt is None or entry.e_altref is None:
                continue
            e_refalt[i] = entry.e_refalt
            e_altref[i] = entry.e_altref
    error_mask = ~np.isnan(e_refalt)

    bias_mask = bias != 0.0
    # logit(clip(0.5 + bias)) at the biased markers, in array order, so it lines
    # up with ``w[bias_mask]``. This is exactly the ``logit(p)`` term inside
    # ``apply_bias``; lifting it here keeps the per-evaluation work to the
    # w-dependent part only.
    logit_bias_masked = logit(np.clip(0.5 + bias[bias_mask], W_EPS, 1.0 - W_EPS))

    return _MarkerArrays(
        host_ref_dose=host_ref_dose,
        donor_ref_dose=donor_ref_dose,
        n=ad_ref + ad_alt,
        k=ad_alt,
        bias=bias,
        bias_mask=bias_mask,
        e_refalt=e_refalt,
        e_altref=e_altref,
        error_mask=error_mask,
        has_bias=bool(bias_mask.any()),
        has_error=bool(error_mask.any()),
        logit_bias_masked=logit_bias_masked,
    )


def _total_ll_vec(
    arr: _MarkerArrays,
    f_donor: float,
    error_rate: float = DEFAULT_ERROR_RATE,
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
    return _ll_from_p_alt(arr, _p_alt_for_f(arr, f_donor, error_rate), rho)


def _p_alt_for_f(
    arr: _MarkerArrays,
    f_donor: float,
    error_rate: float = DEFAULT_ERROR_RATE,
) -> np.ndarray:
    """Per-marker P(observe ALT) at donor fraction ``f_donor`` (rho-independent).

    This is the f-dependent half of ``_total_ll_vec``. Splitting it out lets the
    rho profiling (which holds f fixed and sweeps rho) compute this once per f
    instead of once per (f, rho). The arithmetic is identical to the inline
    version, so the returned array is bit-for-bit the same.
    """
    # Expected reference-allele weight, with the conditional logit-space bias
    # correction (see apply_bias) applied only at markers that have an estimate.
    w = (1.0 - f_donor) * arr.host_ref_dose / 2.0 + f_donor * arr.donor_ref_dose / 2.0
    if arr.has_bias:
        # Inlined apply_bias with the bias-only ``logit(p)`` term precomputed.
        wm = np.clip(w[arr.bias_mask], W_EPS, 1.0 - W_EPS)
        w[arr.bias_mask] = np.clip(
            expit(logit(wm) - arr.logit_bias_masked), W_EPS, 1.0 - W_EPS
        )

    # P(observe ALT | w). Default is the 4-state symmetric model with the
    # global ``error_rate``; per-marker asymmetric rates (where supplied) use
    # the REF/ALT-only form ``p_alt = w * e_refalt + (1 - w) * (1 - e_altref)``.
    e = error_rate
    e_specific = e / N_OTHER_BASES
    p_alt_raw = (1.0 - w) * (1.0 - e) + w * e_specific
    p_ref_raw = w * (1.0 - e) + (1.0 - w) * e_specific
    p_alt = p_alt_raw / (p_ref_raw + p_alt_raw)
    if arr.has_error:
        # Sub in the asymmetric per-marker rates at the masked positions.
        p_alt_asym = w * arr.e_refalt + (1.0 - w) * (1.0 - arr.e_altref)
        p_alt = np.where(arr.error_mask, p_alt_asym, p_alt)
    return np.clip(p_alt, 1e-6, 1.0 - 1e-6)


def _ll_from_p_alt(arr: _MarkerArrays, p_alt: np.ndarray, rho: float) -> float:
    """Total beta-binomial log-likelihood from a precomputed ``p_alt`` and rho.

    The rho-dependent half of ``_total_ll_vec``; identical arithmetic to the
    inline version.
    """
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
    error_rate: float = DEFAULT_ERROR_RATE,
    rho: float = 100.0,
    calibration: PanelCalibration | None = None,
) -> float:
    """Sum of per-marker beta-binomial log-likelihoods for multi-donor model.

    Args:
        markers: Informative markers (for at least one donor).
        donor_fractions: [f_donor1, f_donor2, ...].
        error_rate: Sequencing error rate (fallback when a marker's
            per-direction rate is missing).
        rho: Beta-binomial concentration parameter.
        calibration: Optional per-marker bias and error tables
            (see ``total_log_likelihood_bb``).

    Returns:
        Total log-likelihood.
    """
    cal = calibration or PanelCalibration()
    ll = 0.0
    for m in markers:
        bias = cal.bias_for(m)
        entry = cal.error_for(m)
        e_ra = entry.e_refalt if entry is not None else None
        e_ar = entry.e_altref if entry is not None else None
        w = expected_weight_multi(m.host_gt, m.donor_gts, donor_fractions, bias=bias)
        ll += log_likelihood_marker_bb(
            m.admix_ad_ref,
            m.admix_ad_alt,
            w,
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
        w = float(apply_bias(w, bias))
    return w


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
_INFEASIBLE_PENALTY = 1e30


def fraction_se(
    markers: list[InformativeMarker],
    f_donor: float,
    error_rate: float = DEFAULT_ERROR_RATE,
    rho: float = float("inf"),
    calibration: PanelCalibration | None = None,
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

    Returns:
        Standard error of the donor fraction. inf if no marker is informative.
    """
    cal = calibration or PanelCalibration()
    # P(observe ALT) is linear in the REF weight w, so its slope dp_alt/dw is
    # constant and equals the change across the full weight range w: 0 -> 1.
    dpalt_dw = alt_read_probability(1.0, error_rate) - alt_read_probability(0.0, error_rate)

    info = 0.0
    for m in markers:
        bias = cal.bias_for(m)

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
    error_rate: float = DEFAULT_ERROR_RATE,
    rho: float = float("inf"),
    calibration: PanelCalibration | None = None,
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

    Returns:
        ``(lob, lod)`` as donor fractions (0.0-1.0). ``(inf, inf)`` if no
        marker is informative.
    """
    se0 = fraction_se(markers, 0.0, error_rate, rho, calibration)
    if math.isinf(se0):
        return float("inf"), float("inf")
    lob = _Z95 * se0
    se_lob = fraction_se(markers, lob, error_rate, rho, calibration)
    if math.isinf(se_lob):
        se_lob = se0
    lod = lob + _Z95 * se_lob
    return lob, lod


def _estimate_single_donor_bb_core(
    markers: list[InformativeMarker],
    error_rate: float = DEFAULT_ERROR_RATE,
    grid_steps: int = 1001,
    calibration: PanelCalibration | None = None,
) -> ChimerismResult:
    """Single-donor beta-binomial MLE over the given marker set (no robust trim).

    This is the unguarded estimator; ``estimate_single_donor_bb`` wraps it with
    the optional robust refit. Args/returns as in that wrapper.
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

    for f in grid:
        # p_alt depends on f only, so compute it once and reuse it across the
        # rho profiling (was recomputed at every rho evaluation).
        p_alt_f = _p_alt_for_f(arr, f, error_rate)
        opt_rho = minimize_scalar(
            lambda log_r, _p=p_alt_f: -_ll_from_p_alt(arr, _p, math.exp(log_r)),
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
        lob_fraction=lob,
        lod_fraction=lod,
    )


# ---------------------------------------------------------------------------
# Fast vectorized grid estimator (opt-in, single donor)
# ---------------------------------------------------------------------------
#
# The exact estimator above profiles rho with ``minimize_scalar`` at every grid
# f and then runs a joint Nelder-Mead refine, so each call makes hundreds of
# scalar scipy optimisations. For the LoD sweep we only need ``donor_fraction``,
# and that is recovered to well within 1e-4 by maximising the same beta-binomial
# log-likelihood on a 2-D (f, rho) grid, then a single 1-D bounded f-search
# bracketed around the grid argmax (with rho profiled out). This path replaces
# the per-f scipy calls of the grid search with one vectorized array pass and a
# short local refine.
#
# Accuracy is validated against the exact estimator in
# ``scripts/validate_grid_estimator.py``: across the LoD parameter space the
# fraction agrees to < 1e-4 and the resulting LoD-summary percentages to well
# under 0.01 pp. The exact estimator remains the default; this is opt-in only.


@dataclass
class GridChimerismResult:
    """Lightweight result of the fast grid single-donor estimator.

    Carries only what the LoD sweep consumes (``donor_fraction``) plus a few
    fields so callers that read ``ChimerismResult`` attributes still work. CI is
    a coarse profile-likelihood bracket on the grid; ``detection_limit`` style
    fields are not computed.
    """

    donor_fraction: float
    donor_fraction_ci: tuple[float, float]
    host_fraction: float
    log_likelihood: float
    n_informative: int
    n_markers_used: int
    error_rate: float
    rho: float = float("inf")


def _p_alt_grid(
    arr: _MarkerArrays,
    f_grid: np.ndarray,
    error_rate: float = DEFAULT_ERROR_RATE,
) -> np.ndarray:
    """Per-marker P(observe ALT) for a whole f-grid at once.

    Vectorized form of ``_p_alt_for_f`` over an ``(n_f,)`` fraction grid; returns
    an ``(n_f, M)`` array whose row i is exactly ``_p_alt_for_f(arr, f_grid[i])``.

    Args:
        arr: Per-marker arrays from ``_precompute_marker_arrays``.
        f_grid: Donor fractions to evaluate, shape ``(n_f,)``.
        error_rate: Sequencing error rate.

    Returns:
        ``(n_f, M)`` array of P(observe ALT).
    """
    f = f_grid[:, None]  # (n_f, 1)
    host = arr.host_ref_dose[None, :]  # (1, M)
    donor = arr.donor_ref_dose[None, :]
    w = (1.0 - f) * host / 2.0 + f * donor / 2.0  # (n_f, M)

    if arr.has_bias:
        # Inlined apply_bias with the bias-only logit(p) term precomputed.
        bm = arr.bias_mask
        wm = np.clip(w[:, bm], W_EPS, 1.0 - W_EPS)
        w[:, bm] = np.clip(
            expit(logit(wm) - arr.logit_bias_masked[None, :]), W_EPS, 1.0 - W_EPS
        )

    e = error_rate
    e_specific = e / N_OTHER_BASES
    p_alt_raw = (1.0 - w) * (1.0 - e) + w * e_specific
    p_ref_raw = w * (1.0 - e) + (1.0 - w) * e_specific
    p_alt = p_alt_raw / (p_ref_raw + p_alt_raw)
    if arr.has_error:
        p_alt_asym = w * arr.e_refalt[None, :] + (1.0 - w) * (1.0 - arr.e_altref[None, :])
        p_alt = np.where(arr.error_mask[None, :], p_alt_asym, p_alt)
    return np.clip(p_alt, 1e-6, 1.0 - 1e-6)


def _ll_grid_over_rho(
    arr: _MarkerArrays,
    p_alt: np.ndarray,
    rho_grid: np.ndarray,
) -> np.ndarray:
    """Total log-likelihood on the full ``(n_f, n_rho)`` grid.

    ``p_alt`` is the ``(n_f, M)`` array from ``_p_alt_grid``. Loops over the
    rho-grid (cheap: ``n_rho`` is small) so the largest temporary is one
    ``(n_f, M)`` block per rho, keeping memory bounded. Each rho column equals
    summing the per-marker beta-binomial terms over the marker axis, i.e. the
    vectorized counterpart of ``_ll_from_p_alt`` evaluated for every f at once.

    Returns:
        ``(n_f, n_rho)`` array of total log-likelihoods.
    """
    n, k = arr.n[None, :], arr.k[None, :]
    out = np.empty((p_alt.shape[0], rho_grid.shape[0]), dtype=float)
    for j, rho in enumerate(rho_grid):
        a = np.maximum(p_alt * rho, 1e-10)
        b = np.maximum((1.0 - p_alt) * rho, 1e-10)
        ll = (
            gammaln(k + a)
            + gammaln(n - k + b)
            - gammaln(n + rho)
            - gammaln(a)
            - gammaln(b)
            + math.lgamma(rho)
        )
        out[:, j] = ll.sum(axis=1)
    return out


# rho range for the fast grid and its profile. The lower bound matches the exact
# estimator's grid rho-profiling; the upper bound is _RHO_MAX (50000), the same
# ceiling the exact estimator's Nelder-Mead refine allows. Capping the fast
# profile at the lower 10000 (the exact estimator's *grid* bound) leaves a ~0.01
# pp bias at high-concentration samples where the optimum rho sits above 10000;
# matching _RHO_MAX removes it.
_GRID_RHO_LO = 1.0


def estimate_single_donor_bb_grid(
    markers: list[InformativeMarker],
    error_rate: float = DEFAULT_ERROR_RATE,
    calibration: PanelCalibration | None = None,
    n_f: int = 201,
    n_rho: int = 32,
    refine: bool = True,
) -> GridChimerismResult:
    """Fast approximate single-donor MLE via a vectorized (f, rho) grid.

    Opt-in alternative to ``estimate_single_donor_bb`` for parameter sweeps
    where only ``donor_fraction`` is needed. Builds an f-grid on ``[0, 1]`` and a
    log-spaced rho-grid on ``[1, _RHO_MAX]``, evaluates the beta-binomial total
    log-likelihood on the full grid in a few vectorized array passes, and takes
    the grid argmax. When ``refine`` (the default) it then polishes the donor
    fraction with a 1-D bounded search bracketed to +/- 2 f-grid steps around the
    argmax, profiling rho out at each candidate f (the same profile the exact
    estimator uses). Because the total log-likelihood is unimodal in f, this
    local solve lands on the exact estimator's fraction to < 1e-4 across the LoD
    space (see ``scripts/validate_grid_estimator.py``). The 1-D profiled refine
    is both faster and tighter than a joint Nelder-Mead polish, which can stall
    in the very flat rho direction.

    The exact estimator (``estimate_single_donor_bb``) stays the default and is
    untouched; this path is selected explicitly by the caller.

    Performance and accuracy (measured on the LoD sweep parameter space, see
    ``scripts/validate_grid_estimator.py``):

      - About 6.5x faster per call than the exact estimator (~30 ms vs ~191 ms;
        4-7x depending on panel size, with the grid build dominating and the
        refine ~5 ms). The win comes from replacing the exact estimator's
        Python-level f-grid loop (a bounded rho profile per f point) and the
        joint Nelder-Mead refinement with a few vectorized array passes.
      - Donor-fraction agreement with the exact estimator: median 1e-6 pp,
        worst case 0.0115 pp over fractions <= 5% (the whole LoD-sweep regime).
        The only deviations above 0.01 pp are at f=0.5 (outside the LoD range),
        and there the grid finds the strictly higher likelihood, i.e. it is the
        more accurate of the two, not the reverse.
      - End-to-end effect on the reported LoD: across the full 60-cell
        ``lod_summary`` grid (relatedness x depth x panel size), the per-cell
        ``lod_pct`` matches the exact estimator to a median of 0.0000 pp and a
        maximum of 0.0011 pp, comfortably under a 0.01 pp tolerance.

    So this path is appropriate for fast iteration on the LoD sweeps; run the
    exact estimator (the default) for the final publication figures.

    Args:
        markers: Informative markers with admixture allele counts.
        error_rate: Sequencing error rate (fallback when a marker lacks
            per-direction rates).
        calibration: Optional per-marker bias and error tables.
        n_f: Number of f-grid points on ``[0, 1]``.
        n_rho: Number of log-spaced rho-grid points on ``[1, _RHO_MAX]``.
        refine: Run the 1-D profiled local polish from the grid argmax
            (default True).

    Returns:
        ``GridChimerismResult`` with the donor-fraction estimate and a coarse CI.
    """
    cal = calibration or PanelCalibration()
    n_informative = len(markers)
    if n_informative == 0:
        return GridChimerismResult(
            donor_fraction=0.0,
            donor_fraction_ci=(0.0, 0.0),
            host_fraction=1.0,
            log_likelihood=0.0,
            n_informative=0,
            n_markers_used=0,
            error_rate=error_rate,
        )

    arr = _precompute_marker_arrays(markers, cal)

    f_grid = np.linspace(0.0, 1.0, n_f)
    rho_grid = np.exp(np.linspace(math.log(_GRID_RHO_LO), math.log(_RHO_MAX), n_rho))

    p_alt = _p_alt_grid(arr, f_grid, error_rate)  # (n_f, M)
    ll = _ll_grid_over_rho(arr, p_alt, rho_grid)  # (n_f, n_rho)

    flat = int(np.argmax(ll))
    fi, ri = np.unravel_index(flat, ll.shape)
    best_f = float(f_grid[fi])
    best_rho = float(rho_grid[ri])
    best_ll = float(ll[fi, ri])

    def profile_ll_f(f_val: float) -> float:
        """Max LL over rho at a fixed f (rho profiled out, as in the exact MLE)."""
        p = _p_alt_for_f(arr, f_val, error_rate)
        opt_rho = minimize_scalar(
            lambda log_r, _p=p: -_ll_from_p_alt(arr, _p, math.exp(log_r)),
            bounds=(math.log(_GRID_RHO_LO), math.log(_RHO_MAX)),
            method="bounded",
        )
        return -float(opt_rho.fun)

    if refine:
        # The LL is unimodal in f, so the optimum lies within one grid step of the
        # argmax; bracket +/- 2 steps for safety and solve f with rho profiled out.
        step = 1.0 / (n_f - 1)
        lo = max(0.0, best_f - 2.0 * step)
        hi = min(1.0, best_f + 2.0 * step)
        if hi > lo:
            opt = minimize_scalar(
                lambda f: -profile_ll_f(f),
                bounds=(lo, hi),
                method="bounded",
                options={"xatol": 1e-8},
            )
            f_ref = max(0.0, min(1.0, float(opt.x)))
            ll_ref = -float(opt.fun)
            if ll_ref >= best_ll:
                best_f = f_ref
                best_ll = ll_ref
                # Profile rho at the refined f for the reported concentration.
                p = _p_alt_for_f(arr, best_f, error_rate)
                opt_rho = minimize_scalar(
                    lambda log_r: -_ll_from_p_alt(arr, p, math.exp(log_r)),
                    bounds=(math.log(_GRID_RHO_LO), math.log(_RHO_MAX)),
                    method="bounded",
                )
                best_rho = math.exp(float(opt_rho.x))

    # Coarse profile-likelihood CI from the grid: profile rho out of the (f, rho)
    # grid (max over the rho axis at each f), then bracket where the profile drops
    # by chi2(0.95, df=1)/2 from the maximum. Cheap and good enough for a sweep.
    prof = ll.max(axis=1)
    half_threshold = float(chi2.ppf(CI_LEVEL, df=1)) / 2.0
    above = prof >= (best_ll - half_threshold)
    idx = np.nonzero(above)[0]
    if idx.size:
        f_lo = float(f_grid[idx[0]])
        f_hi = float(f_grid[idx[-1]])
    else:
        f_lo = f_hi = best_f

    return GridChimerismResult(
        donor_fraction=best_f,
        donor_fraction_ci=(f_lo, f_hi),
        host_fraction=1.0 - best_f,
        log_likelihood=best_ll,
        n_informative=n_informative,
        n_markers_used=n_informative,
        error_rate=error_rate,
        rho=best_rho,
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
            # the host's own ALT dose (host_alt / 2), so the host-present
            # direction at the current fit is sign(host_alt / 2 - expected_vaf).
            # A residual deviating that way is under-fit host signal, not an
            # artifact, and must not be trimmed. ``current`` and
            # ``result.per_marker`` are in the same order (the core estimator
            # builds per-marker results by iterating the marker list).
            host_alt = np.array([m.host_gt[0] + m.host_gt[1] for m in current], dtype=float)
            exp_vaf = np.array([mr.expected_vaf for mr in result.per_marker], dtype=float)
            host_dir = np.sign(host_alt / 2.0 - exp_vaf)
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

    Returns:
        ChimerismResult with MLE estimate and beta-binomial CIs.
    """
    if robust not in ROBUST_MODES:
        raise ValueError(f"robust must be one of {ROBUST_MODES}, got {robust!r}")

    cal = calibration or PanelCalibration()

    def core(mk: list[InformativeMarker]) -> ChimerismResult:
        return _estimate_single_donor_bb_core(mk, error_rate, grid_steps, cal)

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
