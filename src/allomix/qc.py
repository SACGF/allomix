"""Quality control assessment for chimerism results.

Evaluates marker counts, sequencing depth, confidence intervals,
and goodness-of-fit to flag potential issues in chimerism estimates.
"""

from __future__ import annotations

import math
import statistics
from dataclasses import dataclass, field

from scipy.stats import chi2

from allomix.chimerism import ChimerismResult, MarkerResult
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


def _error_adjusted_p_alt(expected_vaf: float, error_rate: float) -> float:
    """Error-model-adjusted expected ALT fraction.

    Mirrors the 4-state error model used by ``log_likelihood_marker_bb``,
    so that GoF compares the observed VAF to the same quantity the
    likelihood does. The raw ``expected_vaf`` stored on ``MarkerResult``
    is ``1 - w`` without sequencing-error adjustment; at homozygous
    markers this is 0 or 1 exactly, clamping the variance floor and
    producing spurious chi-sq blow-ups against observations that differ
    by a typical sequencing error (~1%).
    """
    e = error_rate
    w = 1.0 - expected_vaf  # reference weight
    p_alt = (1.0 - w) * (1.0 - e) + w * e / 3.0
    p_ref = w * (1.0 - e) + (1.0 - w) * e / 3.0
    return p_alt / (p_ref + p_alt)


def _compute_gof_pval(
    per_marker: list[MarkerResult],
    rho: float = float("inf"),
    n_fitted_params: int = 2,
    error_rate: float = 0.0,
) -> float | None:
    """Compute chi-squared goodness-of-fit p-value from per-marker residuals.

    Uses the sum of squared Pearson residuals standardised by the
    beta-binomial variance at each marker:

        Var(k/n) = p(1-p) * (n + rho) / (n * (rho + 1))

    As rho -> inf this collapses to the binomial Pearson chi-squared.

    When ``error_rate > 0``, the expected VAF used for the variance
    floor is adjusted via the 4-state error model (``p_alt`` from
    ``log_likelihood_marker_bb``) so that at f near 0 or 1 the floor is
    the actual error rate, not ``1 - 1e-6``. Without this, a typical
    ~1% sequencing-error residual against a saturated ``expected_vaf``
    of 0 or 1 produces a spurious chi-sq blow-up. The stored
    ``m.residual`` is used as-is so synthetic test fixtures that set a
    residual independently of the observed VAF keep working.

    Args:
        per_marker: Per-marker results from chimerism estimation.
        rho: Fitted beta-binomial concentration parameter. Pass
            ``math.inf`` for binomial scaling.
        n_fitted_params: Number of parameters jointly estimated with the
            fit, used to set degrees of freedom. Single-donor BB: 2
            (f, rho). Multi-donor BB with k donors: k + 1.
        error_rate: Sequencing error rate used by the likelihood. Pass
            0.0 to use the raw ``expected_vaf`` for the variance floor.

    Returns:
        p-value from chi-squared survival function, or None if there are
        not enough included markers (<= n_fitted_params).
    """
    included = [m for m in per_marker if m.included]
    if len(included) <= n_fitted_params:
        return None

    chi_sq = 0.0
    for m in included:
        n = m.dp
        if n <= 0:
            continue
        if error_rate > 0:
            ev_raw = _error_adjusted_p_alt(m.expected_vaf, error_rate)
        else:
            ev_raw = m.expected_vaf
        ev = max(1e-6, min(1.0 - 1e-6, ev_raw))
        if math.isinf(rho):
            var_vaf = ev * (1.0 - ev) / n
        else:
            var_vaf = ev * (1.0 - ev) * (n + rho) / (n * (rho + 1.0))
        if var_vaf <= 0:
            continue
        chi_sq += m.residual**2 / var_vaf

    df = len(included) - n_fitted_params
    if df <= 0:
        return None
    pval: float = chi2.sf(chi_sq, df)
    return pval


def _marker_loss_diagnosis(g: MarkerGenotypes, n_informative: int) -> str:
    """Explain which input starved the informative-marker set.

    Walks the per-input marker counts recorded by ``classify_markers`` and names
    the dominant bottleneck: a low-coverage admixture sample, sparse host or donor
    genotyping, or a per-marker filter (depth, GQ, PASS). Returns "" when no
    counts were recorded (``g.marker_counts`` is None, e.g. a hand-built
    MarkerGenotypes).

    Args:
        g: Marker genotype classification, with ``g.marker_counts`` populated.
        n_informative: Number of informative markers that survived.

    Returns:
        A one-line diagnostic, or "" if no count data is available.
    """
    mc = g.marker_counts
    if mc is None or mc.n_admix == 0:
        return ""

    donor_str = "/".join(str(x) for x in mc.n_donor_markers) if mc.n_donor_markers else "?"
    summary = (
        f"counts: host {mc.n_host}, donor {donor_str}, admix {mc.n_admix}; "
        f"{g.n_shared} shared; {n_informative} informative"
    )

    # 1. Admixture sample itself sparse (many no-calls / low input).
    geno_counts = [mc.n_host, *mc.n_donor_markers]
    if geno_counts and mc.n_admix < 0.5 * min(geno_counts):
        return (
            f"Few informative markers: admixture sample has only {mc.n_admix} genotyped "
            f"markers vs host {mc.n_host}/donor {donor_str} (low coverage or many no-calls) "
            f"[{summary}]."
        )

    # 2. A genotyping input covers few of the admixture markers.
    cover = [("host", mc.n_admix_in_host)]
    for i, c in enumerate(mc.n_admix_in_donor):
        cover.append(("donor" if len(mc.n_admix_in_donor) == 1 else f"donor{i + 1}", c))
    name, cov = min(cover, key=lambda kv: kv[1])
    if cov < 0.5 * mc.n_admix:
        return (
            f"Few informative markers: {name} genotyping covers only {cov}/{mc.n_admix} "
            f"admixture markers (genotyping likely failed) [{summary}]."
        )

    # 3. Sharing is fine; blame the dominant per-marker filter.
    drops = [
        ("low admixture depth (DP<min)", mc.n_drop_admix_dp),
        ("host GQ below threshold", mc.n_drop_gq_host),
        ("donor GQ below threshold", mc.n_drop_gq_donor),
        ("non-PASS calls", mc.n_drop_pass),
    ]
    reason, n = max(drops, key=lambda kv: kv[1])
    if n > 0:
        return (
            f"Few informative markers: {n}/{g.n_shared} shared markers dropped at "
            f"{reason} [{summary}]."
        )

    return f"Few informative markers [{summary}]."


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
    rho = getattr(result, "rho", float("inf"))
    if hasattr(result, "donor_fractions"):
        n_fitted = len(result.donor_fractions) + 1  # k donors + rho
    else:
        n_fitted = 2  # f + rho
    gof_pval = _compute_gof_pval(
        result.per_marker,
        rho=rho,
        n_fitted_params=n_fitted,
        error_rate=result.error_rate,
    )

    # --- QC checks ---

    # Insufficient informative markers
    if n_informative < min_informative:
        pass_ = False
        warnings.append(f"Insufficient informative markers: {n_informative} < {min_informative}")
        diagnosis = _marker_loss_diagnosis(genotypes, n_informative)
        if diagnosis:
            warnings.append(diagnosis)

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
