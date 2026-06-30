"""Likelihood and weight model for chimerism estimation (no optimisation).

The beta-binomial per-marker likelihood and the expected-weight model it is
built on; ``allomix.chimerism`` imports these to build the MLE estimators.

The beta-binomial handles overdispersion from per-marker amplification bias and
depth variability. A binomial model assumes all variance is random sampling, but
systematic effects add extra-binomial variance that makes binomial CIs
undercover. The shared concentration parameter rho, jointly estimated from the
data, widens CIs when overdispersion is present, with the largest coverage gain
at low donor fractions where accurate CIs matter most clinically. Per-marker
amplification bias (het-site VAF deviation, SD ~0.02 empirically) is corrected
multiplicatively in logit space (see ``apply_bias``).
"""

from dataclasses import dataclass, field
from math import lgamma

import numpy as np
from scipy.special import expit, gammaln, logit

from allomix.constants import DEFAULT_ERROR_RATE, N_OTHER_BASES, PLOIDY
from allomix.contamination_table import ContaminationCorrection
from allomix.error_rates import MarkerErrorRates
from allomix.genotype import InformativeMarker, MarkerKey


@dataclass(frozen=True)
class PanelCalibration:
    """Per-marker calibration applied during chimerism estimation.

    Bundles the optional per-marker tables the estimators consume:

    - ``biases``: amplification bias per marker (see ``allomix.bias``). Positive
      means the ALT allele is preferentially captured.
    - ``errors``: per-direction empirical substitution rates (see
      ``allomix.error_rates``). Used for the asymmetric REF/ALT-only likelihood
      where both directions are known.
    - ``contamination_correction``: per-marker co-pooled contamination correction
      (Step 30, see ``allomix.contamination_table``). Applied by
      ``estimate_single_donor_bb`` before the MLE. None (the default) and a
      gated-out table both leave estimation byte-identical.

    Tables default to empty, so an uncalibrated run is ``PanelCalibration()``.
    Markers absent from a table fall through to the uncorrected weight (bias 0)
    or the symmetric global ``error_rate``.
    """

    biases: dict[MarkerKey, float] = field(default_factory=dict)
    errors: dict[MarkerKey, MarkerErrorRates] = field(default_factory=dict)
    contamination_correction: ContaminationCorrection | None = None

    def __post_init__(self) -> None:
        # Treat an explicit ``None`` the same as an omitted argument, so callers
        # can pass an optional dict straight through without an ``or {}`` guard.
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
    """
    p = np.clip(0.5 + bias, W_EPS, 1.0 - W_EPS)  # observed het ALT-favouring as a probability
    w_clamped = np.clip(w, W_EPS, 1.0 - W_EPS)
    return np.clip(expit(logit(w_clamped) - logit(p)), W_EPS, 1.0 - W_EPS)


def inject_bias(alt_vaf: float | np.ndarray, bias: float | np.ndarray) -> float | np.ndarray:
    """Shift a true ALT VAF by a het-site bias, the simulator-side counterpart
    of ``apply_bias``.

    The simulator injects bias the same way the estimator corrects it, so the
    two stay self-consistent: at the true parameters the injected ALT VAF equals
    the estimator's expected biased ALT VAF, ``1 - apply_bias(w_true, bias)``.
    Equivalent to ``expit(logit(alt_vaf) + logit(0.5 + bias))``; at a het site
    (true VAF 0.5) the observed VAF becomes ``0.5 + bias``. Not an algebraic
    inverse of ``apply_bias``: both shift in the same direction, so composing
    them would double-correct.
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

    where ref_dose = 2 - alt_dose. The single-donor case of
    ``expected_weight_multi``, delegating to it to keep the dose math in one
    place (the two are algebraically identical, bias branch included).
    """
    return expected_weight_multi(host_gt, [donor_gt], [f_donor], bias=bias)


def alt_read_probability(w: float, error_rate: float = DEFAULT_ERROR_RATE) -> float:
    """Probability an observed read is ALT, given expected REF weight ``w``.

    Under the 4-state error model a read is observed as ALT if it comes from a
    true ALT allele and is called correctly (probability ``1 - e``), or from a
    true REF allele miscalled to the ALT base (probability ``e / 3``). The two
    raw probabilities are renormalised so they condition on the read being
    called REF or ALT rather than one of the other two bases.
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
        error_rate: Symmetric 4-state rate, used only when asymmetric rates are
            not both provided.
        rho: Beta-binomial concentration. Larger = less overdispersion. Typical
            empirical values: 50-500.
        e_refalt: ``P(observe ALT | true REF)`` for this marker.
        e_altref: ``P(observe REF | true ALT)`` for this marker.
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
        error_rate: Fallback when a marker is missing from ``calibration.errors``
            or only one direction is known.
        calibration: Per-marker bias and error tables. A marker's per-direction
            rates are used (asymmetric model) only when both are known; otherwise
            the symmetric ``error_rate`` model applies.
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
    # (f, rho)-independent precomputations hoisted out of the hot path.
    # ``has_bias``/``has_error`` cache the per-call ``.any()`` checks;
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
    """
    n_markers = len(markers)
    host_ref_dose = np.fromiter(
        (PLOIDY - (m.host_gt[0] + m.host_gt[1]) for m in markers),
        dtype=float,
        count=n_markers,
    )
    donor_ref_dose = np.fromiter(
        (PLOIDY - (m.donor_gts[0][0] + m.donor_gts[0][1]) for m in markers),
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
            # Both directions required: one alone leaves the likelihood
            # underspecified at hets and at the opposite-homozygous endpoint.
            if entry.e_refalt is None or entry.e_altref is None:
                continue
            e_refalt[i] = entry.e_refalt
            e_altref[i] = entry.e_altref
    error_mask = ~np.isnan(e_refalt)

    bias_mask = bias != 0.0
    # The ``logit(p)`` term inside ``apply_bias`` at the biased markers (array
    # order, lines up with ``w[bias_mask]``), lifted here so the per-evaluation
    # loop does only the w-dependent part.
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
    w = (1.0 - f_donor) * arr.host_ref_dose / PLOIDY + f_donor * arr.donor_ref_dose / PLOIDY
    if arr.has_bias:
        # Inlined apply_bias with the bias-only ``logit(p)`` term precomputed.
        wm = np.clip(w[arr.bias_mask], W_EPS, 1.0 - W_EPS)
        w[arr.bias_mask] = np.clip(expit(logit(wm) - arr.logit_bias_masked), W_EPS, 1.0 - W_EPS)

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
    """Sum of per-marker beta-binomial log-likelihoods for the multi-donor model.

    Args:
        error_rate: Fallback when a marker's per-direction rate is missing.
    """
    cal = calibration or PanelCalibration()
    ll = 0.0
    for m in markers:
        bias = cal.bias_for(m)
        entry = cal.error_for(m)
        err_kwargs = {}
        if entry is not None:
            err_kwargs = {"e_refalt": entry.e_refalt, "e_altref": entry.e_altref}
        w = expected_weight_multi(m.host_gt, m.donor_gts, donor_fractions, bias=bias)
        ll += log_likelihood_marker_bb(
            m.admix_ad_ref,
            m.admix_ad_alt,
            w,
            error_rate=error_rate,
            rho=rho,
            **err_kwargs,
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

    ``donor_fractions`` must sum to <= 1.0. Returns the expected reference allele
    weight (0.0 to 1.0).
    """
    host_ref_dose = PLOIDY - (host_gt[0] + host_gt[1])
    f_host = 1.0 - sum(donor_fractions)
    w = f_host * host_ref_dose / PLOIDY
    for dgt, f in zip(donor_gts, donor_fractions):
        d_ref_dose = PLOIDY - (dgt[0] + dgt[1])
        w += f * d_ref_dose / PLOIDY
    if bias != 0.0:
        w = float(apply_bias(w, bias))
    return w
