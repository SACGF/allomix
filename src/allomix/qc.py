"""Quality control assessment for chimerism results.

Evaluates marker counts, sequencing depth, confidence intervals,
and goodness-of-fit to flag potential issues in chimerism estimates.
"""

import math
import statistics
from dataclasses import dataclass, field

from scipy.stats import chi2

from allomix.chimerism import ChimerismResult, MarkerResult
from allomix.constants import N_OTHER_BASES
from allomix.detect import HostPresenceResult
from allomix.genotype import MarkerGenotypes
from allomix.relatedness import (
    DEGREE_IDENTICAL,
    MIN_CONSENSUS,
    AdmixConsistencyResult,
    RelatednessResult,
    evaluate_expected,
)

# Thresholds for the optional REVIEW warning when the host-presence detector
# fires significantly but the global MLE does not echo the signal. Tunable
# here rather than buried in ``assess_quality`` so an operator can audit them
# without reading the call site.
HOST_PRESENCE_REVIEW_P = 0.01
# Treat the MLE as "below the detector estimate" when it's < 1/3 of f_h_hat.
# A factor of 3 absorbs sampling noise; tighter ratios produce a lot of
# noise warnings at borderline cells.
HOST_PRESENCE_RATIO_GAP = 3.0
# When the robust refit excludes more than this fraction of informative markers,
# promote the result to REVIEW: a large exclusion points at host copy-number /
# LoH (or a genotyping problem), and the robust refit itself is unreliable once
# the aberrant markers are no longer a clear minority.
ROBUST_REVIEW_FRACTION = 0.15
# The admixture carries alleles in neither host nor donor: when the
# consensus-homozygote swap test is significant at this level (and rests on at
# least MIN_CONSENSUS markers) the sample is promoted to REVIEW.
SWAP_REVIEW_P = 1e-3

# Sample-level QC warning thresholds (soft warnings, not the per-marker filters).
LOW_MEAN_DEPTH_WARN = 100  # warn when mean admixture depth is below this
WIDE_CI_WARN_PCT = 20  # warn when a donor-fraction CI spans more than this (%)
GOF_REVIEW_P = 0.01  # goodness-of-fit p below this promotes the result to REVIEW

# Clamp keeping the expected VAF strictly inside (0, 1) so the GoF variance never
# collapses to zero at homozygous markers (see _error_adjusted_p_alt rationale).
_VAF_EPS = 1e-6


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
        status: Overall QC status, one of "PASS", "REVIEW", or "FAIL". "FAIL"
            means the result is unusable (e.g. too few informative markers).
            "REVIEW" means it was computed but a reliability check failed (poor
            model fit or wide CI), so it needs manual interpretation rather than
            being trusted or discarded automatically.
        per_donor_n_informative: Per-donor informative marker counts (multi-donor).
        relatedness: Estimated relatedness per reference-sample pair, or None.
        admix_consistency: Consensus-homozygote swap check result, or None.
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
    status: str = "PASS"
    per_donor_n_informative: list[int] | None = None
    relatedness: list[RelatednessResult] | None = None
    admix_consistency: AdmixConsistencyResult | None = None

    @property
    def pass_(self) -> bool:
        """True unless the result is a hard FAIL (i.e. PASS or REVIEW)."""
        return self.status != "FAIL"


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
    e_specific = e / N_OTHER_BASES  # miscall to one specific base
    p_alt = (1.0 - w) * (1.0 - e) + w * e_specific
    p_ref = w * (1.0 - e) + (1.0 - w) * e_specific
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
        ev = max(_VAF_EPS, min(1.0 - _VAF_EPS, ev_raw))
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


def _mle_host_estimate(result: ChimerismResult) -> float:
    """Extract the MLE's host-fraction estimate.

    Single-donor results report ``donor_fraction`` directly; multi-donor
    results carry an explicit ``host_fraction`` field. Either way the host
    fraction is what we need to compare against the dedicated detector.
    """
    if hasattr(result, "host_fraction"):
        return float(result.host_fraction)
    return 1.0 - float(result.donor_fraction)


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
    expected_relatedness: list[str] | None = None,
    relatedness_tolerance: int = 1,
) -> QCReport:
    """Assess quality of a chimerism result and produce a QC report.

    Checks marker counts, sequencing depth, confidence interval width,
    goodness-of-fit, and sample identity (relatedness against a declared
    expectation, and an admixture-vs-(host+donor) swap test). Sets a three-state
    status (PASS / REVIEW / FAIL) and collects warnings. Too few informative
    markers is a FAIL (unusable); a relatedness declaration that crosses the
    related/unrelated boundary is a FAIL (sample swap or mislabel); a poor model
    fit, a wide CI, a 2-level relatedness mismatch, or a significant admix swap
    test is a REVIEW (computed, but flagged for manual interpretation). Handles
    both ChimerismResult and MultiDonorResult.

    Identity inputs are read off ``result`` (``result.relatedness`` and
    ``result.admix_consistency``, attached by ``analyse_sample``) via getattr,
    so callers that skip them still get a valid report.

    Args:
        result: The chimerism estimation result to evaluate.
        genotypes: Marker genotype classification data.
        min_informative: Minimum number of informative markers required.
        expected_relatedness: Declared host-vs-donor relationships, aligned with
            the leading host-vs-donor entries of ``result.relatedness`` (one per
            donor, in donor order). "NA"/None entries are skipped.
        relatedness_tolerance: Allowed degree distance for a relatedness PASS.

    Returns:
        QCReport with metrics, warnings, and pass/fail status.
    """
    warnings: list[str] = []
    status = "PASS"
    wide_ci = False
    poor_gof = False

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
        status = "FAIL"
        warnings.append(f"Insufficient informative markers: {n_informative} < {min_informative}")
        diagnosis = _marker_loss_diagnosis(genotypes, n_informative)
        if diagnosis:
            warnings.append(diagnosis)

    # Low mean depth
    if mean_depth < LOW_MEAN_DEPTH_WARN:
        warnings.append(f"Low mean depth: {mean_depth:.0f}x < 100x")

    # CI width check — handle single-donor and multi-donor
    per_donor_n_inf = None
    if hasattr(result, "donor_fraction_cis"):
        # Multi-donor result
        for i, (ci_lo, ci_hi) in enumerate(result.donor_fraction_cis):
            ci_width = (ci_hi - ci_lo) * 100
            if ci_width > WIDE_CI_WARN_PCT:
                wide_ci = True
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
        if ci_width > WIDE_CI_WARN_PCT:
            wide_ci = True
            warnings.append(f"Wide confidence interval: {ci_width:.1f}% > 20%")

    # Poor goodness of fit
    if gof_pval is not None and gof_pval < GOF_REVIEW_P:
        poor_gof = True
        warnings.append(
            "Poor model fit (goodness-of-fit p < 0.01) — "
            "possible genotyping error, CNV, or sample issue"
        )

    # Host-presence vs MLE disagreement: a significant presence test that the
    # global MLE host-fraction estimate does not echo is the clinically
    # interesting "low-level host signal below the MLE's resolution" case.
    # Soft warning only — v1 does not promote qc.status because the operating
    # characteristics on real samples are still being mapped.
    hp: HostPresenceResult | None = getattr(result, "host_presence", None)
    if hp is not None and hp.n_markers > 0 and hp.lrt_pval < HOST_PRESENCE_REVIEW_P:
        mle_host = _mle_host_estimate(result)
        lob = getattr(result, "lob_fraction", float("inf"))
        gap_ratio = hp.f_host_mle / HOST_PRESENCE_RATIO_GAP
        below_lob = math.isfinite(lob) and mle_host < lob
        below_ratio = mle_host < gap_ratio
        if below_lob or below_ratio:
            warnings.append(
                "Low-level host signal detected below the fraction estimate's "
                f"resolution (host_present_p={hp.lrt_pval:.2e}, "
                f"f_host_est={hp.f_host_mle:.4%}, mle_host={mle_host:.4%})"
            )

    # Robust refit exclusions: a few dropped markers is routine, but a large
    # fraction signals host copy-number / LoH (or genotyping error) and means the
    # trimmed fit should not be trusted blindly.
    high_robust_drop = False
    n_robust = getattr(result, "n_robust_excluded", 0)
    drop_frac = getattr(result, "robust_drop_fraction", 0.0)
    if n_robust > 0:
        warnings.append(
            f"Robust refit excluded {n_robust} marker(s) ({drop_frac:.0%}) as "
            "residual outliers (possible host copy-number/LoH or genotyping error)"
        )
        if drop_frac > ROBUST_REVIEW_FRACTION:
            high_robust_drop = True

    # Sample-identity QC. Two checks on the reference samples plus the
    # admixture-vs-(host+donor) swap test.
    identity_review = False
    relatedness: list[RelatednessResult] | None = getattr(result, "relatedness", None)
    duplicate_pairs: set[str] = set()
    if relatedness:
        # Duplicate / sample reuse: two reference samples that should be distinct
        # individuals reading as the same genome. Checked unconditionally (no
        # declaration needed) since it is intrinsically an error, and it makes
        # host-vs-donor chimerism meaningless. Hard FAIL with a clear message.
        for rel in relatedness:
            if rel.degree == DEGREE_IDENTICAL:
                status = "FAIL"
                duplicate_pairs.add(rel.pair)
                coef = "" if rel.coefficient is None else f"r={rel.coefficient:.2f}, "
                warnings.append(
                    f"Identical reference samples: {rel.pair} estimate as the "
                    f"same genome ({coef}{rel.n_sites} markers) — sample "
                    "reuse/mislabel, or an identical-twin (syngeneic) donor; "
                    "either way genotype-based chimerism cannot be measured"
                )

    # Relatedness vs a declared expectation. A close relationship declared but
    # detected unrelated (a likely random swap) is a hard FAIL; a milder
    # mismatch folds into the REVIEW block below. Identical pairs already
    # reported as duplicates above are skipped to avoid a redundant message.
    if relatedness and expected_relatedness:
        # Declarations are per donor; align them to the host-vs-donor pairs,
        # which lead `relatedness` in donor order. strict=True turns a
        # count mismatch (one declaration too few/many) into an error rather
        # than silently leaving a donor unchecked.
        host_pairs = [r for r in relatedness if r.a_name == "host"]
        for rel, declared in zip(host_pairs, expected_relatedness, strict=True):
            if rel.pair in duplicate_pairs:
                continue
            verdict = evaluate_expected(rel, declared, tolerance=relatedness_tolerance)
            if verdict is None:
                continue
            warnings.append(verdict.message)
            if verdict.status == "FAIL":
                status = "FAIL"
            elif verdict.status == "REVIEW":
                identity_review = True

    ac: AdmixConsistencyResult | None = getattr(result, "admix_consistency", None)
    if ac is not None and ac.n_consensus_hom >= MIN_CONSENSUS and ac.swap_pval < SWAP_REVIEW_P:
        identity_review = True
        warnings.append(
            f"Possible sample swap: admixture carries alleles in neither host "
            f"nor donor at {ac.n_discordant}/{ac.n_consensus_hom} "
            f"consensus-homozygous markers (swap p={ac.swap_pval:.2e})"
        )

    # A computed-but-questionable result (poor fit, imprecise, heavily trimmed,
    # or a softer identity flag) is flagged for review rather than passed
    # silently or failed.
    if status != "FAIL" and (poor_gof or wide_ci or high_robust_drop or identity_review):
        status = "REVIEW"

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
        status=status,
        per_donor_n_informative=per_donor_n_inf,
        relatedness=relatedness,
        admix_consistency=ac,
    )
