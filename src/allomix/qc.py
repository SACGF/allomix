"""Quality control assessment for chimerism results.

Evaluates marker counts, sequencing depth, confidence intervals,
and goodness-of-fit to flag potential issues in chimerism estimates.
"""

from __future__ import annotations

import statistics
from dataclasses import dataclass, field

from scipy.stats import chi2

try:
    from allomix.chimerism import ChimerismResult, MarkerResult
except ImportError:

    @dataclass
    class MarkerResult:
        """Per-marker result from chimerism estimation."""

        chrom: str
        pos: int
        marker_type: int
        expected_vaf: float
        observed_vaf: float
        residual: float
        ad_ref: int
        ad_alt: int
        dp: int
        included: bool

    @dataclass
    class ChimerismResult:
        """Result of single-donor chimerism estimation."""

        donor_fraction: float
        donor_fraction_ci: tuple[float, float]
        host_fraction: float
        log_likelihood: float
        n_informative: int
        n_markers_used: int
        per_marker: list[MarkerResult]
        error_rate: float


from allomix.genotype import MarkerGenotypes


@dataclass
class QCReport:
    """Quality control assessment of a chimerism result.

    Attributes:
        n_total_markers: Total markers in input VCFs.
        n_shared_markers: Markers present across all samples.
        n_informative: Markers with differing host/donor genotypes.
        n_used: Markers with included=True in the chimerism fit.
        n_excluded_depth: Markers excluded for low depth.
        n_excluded_quality: Markers excluded for quality issues.
        n_excluded_outlier: Markers excluded as outliers.
        mean_depth: Mean sequencing depth across informative markers.
        median_depth: Median sequencing depth across informative markers.
        min_depth: Minimum sequencing depth across informative markers.
        goodness_of_fit_pval: Chi-squared goodness-of-fit p-value, or None.
        warnings: List of warning messages.
        pass_: Overall QC pass/fail status.
        per_donor_n_informative: Per-donor informative marker counts (multi-donor).
    """

    n_total_markers: int
    n_shared_markers: int
    n_informative: int
    n_used: int
    n_excluded_depth: int
    n_excluded_quality: int
    n_excluded_outlier: int
    mean_depth: float
    median_depth: float
    min_depth: int
    goodness_of_fit_pval: float | None
    warnings: list[str] = field(default_factory=list)
    pass_: bool = True
    per_donor_n_informative: list[int] | None = None


def _compute_gof_pval(per_marker: list[MarkerResult]) -> float | None:
    """Compute chi-squared goodness-of-fit p-value from per-marker residuals.

    Uses the sum of squared Pearson residuals as the test statistic
    with degrees of freedom equal to the number of included markers minus 1.

    Args:
        per_marker: List of per-marker results from chimerism estimation.

    Returns:
        p-value from chi-squared survival function, or None if fewer than
        2 included markers.
    """
    included = [m for m in per_marker if m.included]
    if len(included) < 2:
        return None

    chi_sq = sum(m.residual**2 for m in included)
    df = len(included) - 1
    pval: float = chi2.sf(chi_sq, df)
    return pval


def assess_quality(
    result: ChimerismResult,
    genotypes: MarkerGenotypes,
    min_informative: int = 3,
) -> QCReport:
    """Assess quality of a chimerism result and produce a QC report.

    Checks marker counts, sequencing depth, confidence interval width,
    and goodness-of-fit. Sets pass/fail status and collects warnings.
    Handles both ChimerismResult (single-donor) and MultiDonorResult.

    Args:
        result: The chimerism estimation result to evaluate.
        genotypes: Marker genotype classification data.
        min_informative: Minimum number of informative markers required.

    Returns:
        QCReport with metrics, warnings, and pass/fail status.
    """
    warnings: list[str] = []
    pass_ = True

    # Marker counts
    n_informative = result.n_informative
    n_used = sum(1 for m in result.per_marker if m.included)
    n_excluded_outlier = sum(1 for m in result.per_marker if not m.included)

    # Depth statistics from all per-marker results
    depths = [m.dp for m in result.per_marker]
    if depths:
        mean_depth = statistics.mean(depths)
        median_depth = statistics.median(depths)
        min_depth = min(depths)
    else:
        mean_depth = 0.0
        median_depth = 0.0
        min_depth = 0

    # Goodness of fit
    gof_pval = _compute_gof_pval(result.per_marker)

    # --- QC checks ---

    # Insufficient informative markers
    if n_informative < min_informative:
        pass_ = False
        warnings.append(f"Insufficient informative markers: {n_informative} < {min_informative}")

    # Low mean depth
    if mean_depth < 100:
        warnings.append(f"Low mean depth: {mean_depth:.0f}x < 100x")

    # CI width check — handle single-donor and multi-donor
    per_donor_n_inf = None
    if hasattr(result, "donor_fraction_cis"):
        # Multi-donor result
        for i, (ci_lo, ci_hi) in enumerate(result.donor_fraction_cis):
            ci_width = (ci_hi - ci_lo) * 100
            if ci_width > 20:
                warnings.append(f"Wide CI for donor {i + 1}: {ci_width:.1f}% > 20%")
        per_donor_n_inf = getattr(result, "per_donor_n_informative", None)
        if per_donor_n_inf:
            for i, n_inf in enumerate(per_donor_n_inf):
                if n_inf < min_informative:
                    warnings.append(
                        f"Few informative markers for donor {i + 1}: {n_inf} < {min_informative}"
                    )
    else:
        # Single-donor result
        ci_lo, ci_hi = result.donor_fraction_ci
        ci_width = (ci_hi - ci_lo) * 100
        if ci_width > 20:
            warnings.append(f"Wide confidence interval: {ci_width:.1f}% > 20%")

    # Poor goodness of fit
    if gof_pval is not None and gof_pval < 0.01:
        warnings.append(
            "Poor model fit (goodness-of-fit p < 0.01) — "
            "possible genotyping error, CNV, or sample issue"
        )

    return QCReport(
        n_total_markers=genotypes.n_total,
        n_shared_markers=genotypes.n_shared,
        n_informative=n_informative,
        n_used=n_used,
        n_excluded_depth=genotypes.n_filtered,
        n_excluded_quality=0,
        n_excluded_outlier=n_excluded_outlier,
        mean_depth=mean_depth,
        median_depth=median_depth,
        min_depth=min_depth,
        goodness_of_fit_pval=gof_pval,
        warnings=warnings,
        pass_=pass_,
        per_donor_n_informative=per_donor_n_inf,
    )
