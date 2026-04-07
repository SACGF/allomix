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
from scipy.optimize import minimize_scalar
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

    # Step 3: Profile likelihood CI
    # Threshold: 2 * (LL_max - LL(f)) = chi2.ppf(0.95, df=1)
    threshold = chi2.ppf(0.95, df=1)  # ~3.84
    half_threshold = threshold / 2.0

    # Scan left from MLE
    f_lo = f_mle
    step = 0.001
    while f_lo > 0.0:
        f_test = max(0.0, f_lo - step)
        ll_test = total_log_likelihood(markers, f_test, error_rate, marker_biases)
        if (ll_max - ll_test) > half_threshold:
            # Interpolate between f_test and f_lo for better precision
            f_lo = f_test
            break
        f_lo = f_test
        if f_test == 0.0:
            break

    # Scan right from MLE
    f_hi = f_mle
    while f_hi < 1.0:
        f_test = min(1.0, f_hi + step)
        ll_test = total_log_likelihood(markers, f_test, error_rate, marker_biases)
        if (ll_max - ll_test) > half_threshold:
            f_hi = f_test
            break
        f_hi = f_test
        if f_test == 1.0:
            break

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
        mean_r = sum(residuals) / len(residuals)
        var_r = sum((r - mean_r) ** 2 for r in residuals) / (len(residuals) - 1)
        sd_r = math.sqrt(var_r) if var_r > 0 else 0.0

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
