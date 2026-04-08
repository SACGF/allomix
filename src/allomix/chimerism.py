"""Core MLE chimerism estimation from informative markers.

Implements maximum-likelihood estimation of donor chimerism fraction using
allele counts at informative SNP markers with known host/donor genotypes.
Based on the Demixtify likelihood model (Formula 5) simplified for known
genotypes.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq, minimize, minimize_scalar
from scipy.stats import chi2

from allomix.genotype import InformativeMarker

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


def log_likelihood_marker(
    ad_ref: int,
    ad_alt: int,
    w: float,
    error_rate: float = 0.01,
) -> float:
    """Per-marker log-likelihood using Demixtify Formula 5 with known genotypes.

    LL = n_ref * log(w*(1-e) + (1-w)*e/3) + n_alt * log((1-w)*(1-e) + w*e/3)

    The error rate prevents log(0) when w is exactly 0 or 1.

    Args:
        ad_ref: Reference allele read count.
        ad_alt: Alternative allele read count.
        w: Expected reference allele weight.
        error_rate: Sequencing error rate (default 0.01).

    Returns:
        Log-likelihood contribution from this marker.
    """
    e = error_rate
    p_ref = w * (1.0 - e) + (1.0 - w) * e / 3.0
    p_alt = (1.0 - w) * (1.0 - e) + w * e / 3.0

    # Clamp to avoid log(0) in degenerate cases
    p_ref = max(p_ref, 1e-300)
    p_alt = max(p_alt, 1e-300)

    ll = 0.0
    if ad_ref > 0:
        ll += ad_ref * math.log(p_ref)
    if ad_alt > 0:
        ll += ad_alt * math.log(p_alt)
    return ll


def total_log_likelihood(
    markers: list[InformativeMarker],
    f_donor: float,
    error_rate: float = 0.01,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
) -> float:
    """Sum of per-marker log-likelihoods across all informative markers.

    Args:
        markers: List of informative markers with admixture allele counts.
        f_donor: Donor fraction to evaluate.
        error_rate: Sequencing error rate.
        marker_biases: Optional dict mapping (chrom, pos, ref, alt) to per-marker
            amplification bias. When provided, the expected weight at each marker
            is adjusted to account for systematic capture bias.

    Returns:
        Total log-likelihood.
    """
    ll = 0.0
    for m in markers:
        bias = 0.0
        if marker_biases is not None:
            bias = marker_biases.get((m.chrom, m.pos, m.ref, m.alt), 0.0)
        w = expected_weight(m.host_gt, m.donor_gts[0], f_donor, bias=bias)
        ll += log_likelihood_marker(m.admix_ad_ref, m.admix_ad_alt, w, error_rate)
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


def total_log_likelihood_multi(
    markers: list[InformativeMarker],
    donor_fractions: list[float],
    error_rate: float = 0.01,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
) -> float:
    """Total log-likelihood for multi-donor model.

    A marker contributes if it is informative for any donor. The expected
    weight uses all donor genotypes simultaneously.

    Args:
        markers: Informative markers (for at least one donor).
        donor_fractions: [f_donor1, f_donor2, ...].
        error_rate: Sequencing error rate.
        marker_biases: Optional per-marker bias dict.

    Returns:
        Total log-likelihood.
    """
    ll = 0.0
    for m in markers:
        bias = 0.0
        if marker_biases is not None:
            bias = marker_biases.get((m.chrom, m.pos, m.ref, m.alt), 0.0)
        w = expected_weight_multi(m.host_gt, m.donor_gts, donor_fractions, bias=bias)
        ll += log_likelihood_marker(m.admix_ad_ref, m.admix_ad_alt, w, error_rate)
    return ll


def _compute_overdispersion(
    markers: list[InformativeMarker],
    f_donor: float,
    error_rate: float,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
) -> float:
    """Compute Pearson overdispersion factor for CI calibration.

    Compares observed per-marker residual variance to expected binomial
    variance. Values > 1 indicate model misspecification (e.g. uncorrected
    per-marker amplification bias) and are used to inflate the profile
    likelihood CI threshold proportionally.

    Args:
        markers: Informative markers with admixture allele counts.
        f_donor: MLE donor fraction.
        error_rate: Sequencing error rate.
        marker_biases: Optional per-marker bias corrections.

    Returns:
        Pearson overdispersion factor (phi).
    """
    if len(markers) <= 1:
        return 1.0

    pearson_sum = 0.0
    n_valid = 0

    for m in markers:
        dp = m.admix_ad_ref + m.admix_ad_alt
        if dp == 0:
            continue

        bias = 0.0
        if marker_biases is not None:
            bias = marker_biases.get((m.chrom, m.pos, m.ref, m.alt), 0.0)
        w = expected_weight(m.host_gt, m.donor_gts[0], f_donor, bias=bias)

        e = error_rate
        p_alt = (1.0 - w) * (1.0 - e) + w * e / 3.0
        var_vaf = p_alt * (1.0 - p_alt) / dp
        if var_vaf < 1e-12:
            continue

        observed_alt_frac = m.admix_ad_alt / dp
        residual = observed_alt_frac - p_alt
        pearson_sum += residual * residual / var_vaf
        n_valid += 1

    if n_valid <= 1:
        return 1.0

    return pearson_sum / (n_valid - 1)


def _compute_overdispersion_multi(
    markers: list[InformativeMarker],
    donor_fractions: list[float],
    error_rate: float,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
) -> float:
    """Compute Pearson overdispersion factor for multi-donor model."""
    if len(markers) <= 1:
        return 1.0

    pearson_sum = 0.0
    n_valid = 0

    for m in markers:
        dp = m.admix_ad_ref + m.admix_ad_alt
        if dp == 0:
            continue

        bias = 0.0
        if marker_biases is not None:
            bias = marker_biases.get((m.chrom, m.pos, m.ref, m.alt), 0.0)
        w = expected_weight_multi(m.host_gt, m.donor_gts, donor_fractions, bias=bias)

        e = error_rate
        p_alt = (1.0 - w) * (1.0 - e) + w * e / 3.0
        var_vaf = p_alt * (1.0 - p_alt) / dp
        if var_vaf < 1e-12:
            continue

        observed_alt_frac = m.admix_ad_alt / dp
        residual = observed_alt_frac - p_alt
        pearson_sum += residual * residual / var_vaf
        n_valid += 1

    if n_valid <= 1:
        return 1.0

    return pearson_sum / (n_valid - 1)


# ---------------------------------------------------------------------------
# MLE estimation
# ---------------------------------------------------------------------------


def estimate_error_rate(markers: list[InformativeMarker]) -> float:
    """Estimate sequencing error rate from marker data.

    Currently returns the default value of 0.01. Future versions may estimate
    empirically from non-informative marker data.

    Args:
        markers: List of informative markers (unused in v1).

    Returns:
        Error rate estimate.
    """
    return 0.01


def estimate_single_donor(
    markers: list[InformativeMarker],
    error_rate: float = 0.01,
    grid_steps: int = 1001,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
) -> ChimerismResult:
    """Estimate single-donor chimerism fraction via maximum likelihood.

    Algorithm:
        1. Grid search over f in [0, 1] at 1/grid_steps resolution
        2. Brent refinement via scipy.optimize.minimize_scalar in +/-1% window
        3. Profile likelihood CI using chi-squared threshold (df=1, alpha=0.05)
        4. Per-marker residuals and outlier flagging (>3 SD)

    Args:
        markers: List of informative markers with admixture allele counts.
        error_rate: Sequencing error rate.
        grid_steps: Number of grid points for initial search.
        marker_biases: Optional dict mapping (chrom, pos, ref, alt) to per-marker
            amplification bias for likelihood correction.

    Returns:
        ChimerismResult with MLE estimate, CI, and per-marker details.
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

    # Step 1: Grid search
    grid = np.linspace(0.0, 1.0, grid_steps)
    ll_values = np.array(
        [total_log_likelihood(markers, f, error_rate, marker_biases) for f in grid]
    )
    best_idx = int(np.argmax(ll_values))
    f_grid = float(grid[best_idx])

    # Step 2: Brent refinement in +/-1% window around grid max
    lo = max(0.0, f_grid - 0.01)
    hi = min(1.0, f_grid + 0.01)

    result = minimize_scalar(
        lambda f: -total_log_likelihood(markers, f, error_rate, marker_biases),
        bounds=(lo, hi),
        method="bounded",
    )
    f_mle = float(result.x)
    ll_max = -float(result.fun)

    # Compute overdispersion for robust CIs
    phi = max(1.0, _compute_overdispersion(markers, f_mle, error_rate, marker_biases))

    # Step 3: Profile likelihood CI via root-finding
    # Find f where: ll_max - ll(f) - half_threshold = 0
    # phi > 1 inflates CIs to account for model misspecification (e.g. marker bias)
    threshold = chi2.ppf(0.95, df=1) * phi
    half_threshold = threshold / 2.0

    def ci_func(f: float) -> float:
        """Zero-crossing where LL drops below CI threshold."""
        return ll_max - total_log_likelihood(markers, f, error_rate, marker_biases) - half_threshold

    # Lower bound: find root on [0, f_mle], or 0 if LL never drops enough
    if f_mle <= 0.0 or ci_func(0.0) <= 0.0:
        f_lo = 0.0
    else:
        f_lo = brentq(ci_func, 0.0, f_mle, xtol=1e-5)

    # Upper bound: find root on [f_mle, 1], or 1 if LL never drops enough
    if f_mle >= 1.0 or ci_func(1.0) <= 0.0:
        f_hi = 1.0
    else:
        f_hi = brentq(ci_func, f_mle, 1.0, xtol=1e-5)

    # Step 4: Per-marker results
    per_marker: list[MarkerResult] = []
    residuals: list[float] = []

    for m in markers:
        bias = 0.0
        if marker_biases is not None:
            bias = marker_biases.get((m.chrom, m.pos, m.ref, m.alt), 0.0)
        w = expected_weight(m.host_gt, m.donor_gts[0], f_mle, bias=bias)
        expected_vaf = 1.0 - w  # ALT VAF = 1 - ref_weight
        observed_vaf = m.admix_ad_alt / m.admix_dp if m.admix_dp > 0 else 0.0
        residual = observed_vaf - expected_vaf
        residuals.append(residual)

        per_marker.append(
            MarkerResult(
                chrom=m.chrom,
                pos=m.pos,
                marker_type=m.marker_type,
                expected_vaf=expected_vaf,
                observed_vaf=observed_vaf,
                residual=residual,
                ad_ref=m.admix_ad_ref,
                ad_alt=m.admix_ad_alt,
                dp=m.admix_dp,
                included=True,  # will be updated below
            )
        )

    # Flag outliers (residual > 3 SD) but do NOT exclude them in v1
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

    n_markers_used = sum(1 for mr in per_marker if mr.included)

    return ChimerismResult(
        donor_fraction=f_mle,
        donor_fraction_ci=(f_lo, f_hi),
        host_fraction=1.0 - f_mle,
        log_likelihood=ll_max,
        n_informative=n_informative,
        n_markers_used=n_markers_used,
        per_marker=per_marker,
        error_rate=error_rate,
    )


# ---------------------------------------------------------------------------
# Multi-donor MLE estimation
# ---------------------------------------------------------------------------


def estimate_multi_donor(
    markers: list[InformativeMarker],
    n_donors: int = 2,
    error_rate: float = 0.01,
    grid_steps: int = 101,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
) -> MultiDonorResult:
    """Estimate multi-donor chimerism fractions via maximum likelihood.

    Algorithm:
        1. Triangular grid search over (f1, f2) with f1 + f2 <= 1
        2. Nelder-Mead refinement from grid maximum
        3. Profile likelihood CI per donor (chi-squared df=1 per donor)
        4. Per-marker residuals and outlier flagging

    Args:
        markers: Informative markers (for at least one donor).
        n_donors: Number of donors (currently supports 2).
        error_rate: Sequencing error rate.
        grid_steps: Grid resolution per dimension.
        marker_biases: Optional per-marker bias dict.

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

    # Step 1: Triangular grid search
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
            ll = total_log_likelihood_multi(markers, [f1, f2], error_rate, marker_biases)
            if ll > best_ll:
                best_ll = ll
                best_f = [f1, f2]

    # Step 2: Nelder-Mead refinement
    def neg_ll(x):
        f1, f2 = x
        if f1 < 0 or f2 < 0 or f1 + f2 > 1.0:
            return 1e30
        return -total_log_likelihood_multi(markers, [f1, f2], error_rate, marker_biases)

    opt = minimize(
        neg_ll,
        x0=best_f,
        method="Nelder-Mead",
        options={"xatol": 1e-5, "fatol": 1e-8, "maxiter": 2000},
    )

    f_mle = [max(0.0, float(x)) for x in opt.x]
    if sum(f_mle) > 1.0:
        scale = 1.0 / sum(f_mle)
        f_mle = [f * scale for f in f_mle]
    ll_max = -float(opt.fun)

    # Compute overdispersion for robust CIs
    phi = max(1.0, _compute_overdispersion_multi(markers, f_mle, error_rate, marker_biases))

    # Step 3: Profile likelihood CIs per donor
    cis = _profile_likelihood_cis_multi(
        markers, f_mle, ll_max, n_donors, error_rate, marker_biases, phi
    )

    # Step 4: Per-marker residuals
    per_marker = _per_marker_results_multi(markers, f_mle, error_rate, marker_biases)

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
    )


def _profile_likelihood_cis_multi(
    markers: list[InformativeMarker],
    f_mle: list[float],
    ll_max: float,
    n_donors: int,
    error_rate: float,
    marker_biases: dict[tuple[str, int, str, str], float] | None,
    phi: float = 1.0,
) -> list[tuple[float, float]]:
    """Profile likelihood CIs for each donor fraction.

    For donor_i, scan f_i while optimizing f_j (j != i) at each point.
    Uses chi2(df=1) threshold since we profile one parameter at a time,
    inflated by phi to account for overdispersion.
    """
    threshold = float(chi2.ppf(0.95, df=1)) * phi
    half_threshold = threshold / 2.0
    cis: list[tuple[float, float]] = []

    for donor_idx in range(n_donors):
        other_idx = 1 - donor_idx  # works for 2 donors

        def profile_ll(fi: float, _other=other_idx) -> float:
            """Max LL over the other donor, with donor_idx fixed at fi."""
            max_fj = max(0.0, 1.0 - fi)
            if max_fj < 1e-9:
                fracs = [0.0, 0.0]
                fracs[donor_idx] = fi
                return total_log_likelihood_multi(markers, fracs, error_rate, marker_biases)
            res = minimize_scalar(
                lambda fj: (
                    -total_log_likelihood_multi(
                        markers,
                        [fi, fj] if donor_idx == 0 else [fj, fi],
                        error_rate,
                        marker_biases,
                    )
                ),
                bounds=(0.0, max_fj),
                method="bounded",
            )
            return -float(res.fun)

        def ci_func(fi: float) -> float:
            return ll_max - profile_ll(fi) - half_threshold

        fi_mle = f_mle[donor_idx]

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
    error_rate: float,
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
