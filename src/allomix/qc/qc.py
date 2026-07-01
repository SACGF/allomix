"""Quality control assessment for chimerism results.

Evaluates marker counts, sequencing depth, confidence intervals,
and goodness-of-fit to flag potential issues in chimerism estimates.
"""

import math
import statistics
from dataclasses import dataclass, field

from scipy.stats import chi2

from allomix.constants import N_OTHER_BASES
from allomix.genotype import MarkerGenotypes
from allomix.qc.host_presence import HostPresenceResult
from allomix.qc.relatedness import (
    MIN_CONSENSUS,
    AdmixConsistencyResult,
    Relatedness,
    RelatednessResult,
    SharedHetBalanceResult,
    evaluate_expected,
)
from allomix.qc.runmeta import RunUnitInfo
from allomix.qc.sample_contamination import ContaminationResult
from allomix.results import ChimerismResult, MarkerResult

# Warn when the host-presence detector fires but the global MLE does not echo it.
HOST_PRESENCE_REVIEW_P = 0.01
# MLE counts as "below the detector" under 1/3 of f_h_hat; the factor of 3
# absorbs sampling noise without flooding borderline cells with warnings.
HOST_PRESENCE_RATIO_GAP = 3.0
# Above this robust-trim fraction, promote to REVIEW: a large exclusion points at
# host CNV/LoH (or genotyping error) and the trimmed fit is itself unreliable.
ROBUST_REVIEW_FRACTION = 0.15
# Consensus-homozygote swap test significant at this level (>= MIN_CONSENSUS
# markers) promotes to REVIEW.
SWAP_REVIEW_P = 1e-3

# In-data contamination: third-party signal at consensus-hom markers (see
# ``allomix.qc.sample_contamination``), the low-level floor the gross swap test misses.
# Distinct from the index-hopping provenance flag. At panel depth the pooled test
# is significant for any real excess, so the magnitude thresholds below gate.
CONTAMINATION_P = 0.01
# Warn at/above this floor (significant). Below it is benign vs the <1% target;
# tunable per lab.
CONTAMINATION_WARN_FRACTION = 0.002  # 0.2%
# Above this the floor biases low-level host detection -> REVIEW.
CONTAMINATION_REVIEW_FRACTION = 0.01  # 1%

# Sample-level soft warnings (not per-marker filters).
LOW_MEAN_DEPTH_WARN = 100  # mean admixture depth below this
WIDE_CI_WARN_PCT = 20  # donor-fraction CI wider than this (%)
GOF_REVIEW_P = 0.01  # goodness-of-fit p below this -> REVIEW

# Coverage uniformity: fraction of the sample's per-marker depths that must clear
# UNIFORMITY_DEPTH_FRACTION * mean depth. A lopsided depth distribution (partial
# capture failure, a few markers carrying most reads) can pass the mean-depth check
# while most markers are starved; warn when too few markers clear the bar. Healthy
# panel data keeps >= 0.88 of markers above the bar (SRP434573), so the floor sits
# below that. Warning-only (does not change status).
UNIFORMITY_DEPTH_FRACTION = 0.2
UNIFORMITY_MIN_PASS_FRACTION = 0.8

# Shared-het allele balance (see ``allomix.qc.relatedness.shared_het_balance``). Below
# MIN_SHARED_HET markers the imbalanced fraction is too noisy to act on. At/above
# SHARED_HET_REVIEW_FRACTION imbalanced -> REVIEW (contamination, CNV/allelic
# imbalance, or a sample mix-up). Healthy panel data sits <= 0.07 imbalanced at ~170
# sites (SRP434573), so the gate keeps ~2x margin over sampling noise.
MIN_SHARED_HET = 20
SHARED_HET_REVIEW_FRACTION = 0.15

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
        goodness_of_fit_pval_pretrim: Goodness-of-fit p-value computed on the
            full marker set before any robust trim, or None. Equal to
            goodness_of_fit_pval when no markers were trimmed. A fit that trimmed
            away its outliers can report a clean post-trim GoF, so the REVIEW
            gate uses the worse of the two.
        warnings: List of warning messages.
        status: Overall QC status, one of "PASS", "REVIEW", or "FAIL". "FAIL"
            means the result is unusable (e.g. too few informative markers).
            "REVIEW" means it was computed but a reliability check failed (poor
            model fit or wide CI), so it needs manual interpretation rather than
            being trusted or discarded automatically.
        coverage_uniformity: Fraction of per-marker depths above
            ``UNIFORMITY_DEPTH_FRACTION`` of the sample mean. 1.0 when there are no
            markers. A low value means a lopsided depth distribution.
        per_donor_n_informative: Per-donor informative marker counts (multi-donor).
        relatedness: Estimated relatedness per reference-sample pair, or None.
        admix_consistency: Consensus-homozygote swap check result, or None.
        shared_het_balance: Consensus-het allele-balance check result, or None.
        contamination: In-data third-party contamination estimate, or None.
        run_unit: Sequencing run-unit metadata for the sample, or None.
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
    goodness_of_fit_pval_pretrim: float | None = None
    warnings: list[str] = field(default_factory=list)
    status: str = "PASS"
    coverage_uniformity: float = 1.0
    per_donor_n_informative: list[int] | None = None
    relatedness: list[RelatednessResult] | None = None
    admix_consistency: AdmixConsistencyResult | None = None
    shared_het_balance: SharedHetBalanceResult | None = None
    contamination: ContaminationResult | None = None
    run_unit: RunUnitInfo | None = None

    @property
    def pass_(self) -> bool:
        """True unless the result is a hard FAIL (i.e. PASS or REVIEW)."""
        return self.status != "FAIL"


def _error_adjusted_p_alt(expected_vaf: float, error_rate: float) -> float:
    """Error-model-adjusted expected ALT fraction.

    Mirrors the 4-state error model used by ``log_likelihood_marker_bb`` so GoF
    compares the observed VAF to the same quantity the likelihood does. The raw
    ``expected_vaf`` is 0 or 1 at homozygous markers, which clamps the variance
    floor and produces spurious chi-sq blow-ups against ~1% sequencing-error
    observations.
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
    pretrim: bool = False,
) -> float | None:
    """Compute chi-squared goodness-of-fit p-value from per-marker residuals.

    Sum of squared Pearson residuals standardised by the beta-binomial variance:

        Var(k/n) = p(1-p) * (n + rho) / (n * (rho + 1))

    As rho -> inf this collapses to the binomial Pearson chi-squared. When
    ``error_rate > 0`` the variance-floor expected VAF is adjusted via the 4-state
    error model so the floor at f near 0 or 1 is the real error rate, not
    ``1 - 1e-6`` (which would blow up chi-sq against ~1% error residuals). The
    stored ``m.residual`` is used as-is so synthetic fixtures keep working.

    Args:
        rho: Fitted BB concentration. Pass ``math.inf`` for binomial scaling.
        n_fitted_params: Parameters jointly estimated, setting the degrees of
            freedom. Single-donor BB: 2 (f, rho). Multi-donor with k donors: k + 1.
        error_rate: Pass 0.0 to use the raw ``expected_vaf`` for the variance floor.
        pretrim: When True, evaluate over every marker regardless of its
            ``included`` flag, so a robust trim cannot hide a poor fit by dropping
            its own outliers. When False (default), only ``included`` markers count.

    Returns:
        p-value, or None when there are too few markers (<= n_fitted_params).
    """
    included = per_marker if pretrim else [m for m in per_marker if m.included]
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


def _coverage_uniformity(depths: list[int], depth_fraction: float) -> float:
    """Fraction of ``depths`` exceeding ``depth_fraction`` of the mean depth.

    An evenness metric complementary to mean/median/min: a sample can hold an
    acceptable mean while a lopsided distribution (partial capture failure) starves
    most markers. Returns 1.0 for an empty list or a non-positive mean (nothing to
    flag).
    """
    if not depths:
        return 1.0
    mean = statistics.mean(depths)
    if mean <= 0:
        return 1.0
    cutoff = depth_fraction * mean
    return sum(1 for d in depths if d > cutoff) / len(depths)


def _mle_host_estimate(result: ChimerismResult) -> float:
    """Extract the MLE's host-fraction estimate.

    Multi-donor results carry an explicit ``host_fraction``; single-donor results
    give it as ``1 - donor_fraction``.
    """
    if hasattr(result, "host_fraction"):
        return float(result.host_fraction)
    return 1.0 - float(result.donor_fraction)


def _marker_loss_diagnosis(g: MarkerGenotypes, n_informative: int) -> str:
    """Explain which input starved the informative-marker set.

    Walks the per-input marker counts from ``classify_markers`` and names the
    dominant bottleneck: a low-coverage admixture sample, sparse host or donor
    genotyping, or a per-marker filter (depth, GQ, PASS). Returns "" when no counts
    were recorded (``g.marker_counts`` is None, e.g. a hand-built MarkerGenotypes).
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
        ("low admixture depth (DP<min)", mc.n_dropped_low_admix_dp),
        ("host GQ below threshold", mc.n_dropped_low_gq_host),
        ("donor GQ below threshold", mc.n_dropped_low_gq_donor),
        ("non-PASS calls", mc.n_dropped_failed_filter),
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
    expected_relatedness: list[Relatedness | None] | None = None,
    relatedness_tolerance: int = 1,
) -> QCReport:
    """Assess quality of a chimerism result and produce a QC report.

    Checks marker counts, depth, CI width, goodness-of-fit, and sample identity
    (relatedness against a declared expectation, and an admixture-vs-(host+donor)
    swap test), setting a three-state status and collecting warnings. FAIL means
    unusable (too few informative markers, or a relatedness declaration that
    crosses the related/unrelated boundary). REVIEW means computed but flagged for
    manual interpretation (poor fit, wide CI, 2-level relatedness mismatch, or a
    significant swap test). Handles both ChimerismResult and MultiDonorResult.

    Identity inputs are read off ``result`` (``relatedness``, ``admix_consistency``,
    attached by ``analyse_sample``) via getattr, so callers that skip them still
    get a valid report.

    Args:
        expected_relatedness: Declared host-vs-donor relationships as ``Relatedness``
            members, aligned with the leading host-vs-donor entries of
            ``result.relatedness`` (one per donor, in donor order). None entries (no
            expectation) are skipped.
        relatedness_tolerance: Allowed degree distance for a relatedness PASS.
    """
    warnings: list[str] = []
    status = "PASS"
    wide_ci = False
    poor_gof = False

    n_informative = result.n_informative
    n_used = sum(1 for m in result.per_marker if m.included)
    n_excluded_outlier = sum(1 for m in result.per_marker if not m.included)

    depths = [m.dp for m in result.per_marker]
    if depths:
        mean_depth = statistics.mean(depths)
        median_depth = statistics.median(depths)
        min_depth = min(depths)
    else:
        mean_depth = 0.0
        median_depth = 0.0
        min_depth = 0

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
    # Also evaluate the full pre-trim set: the robust refit can trim the worst
    # markers (CNV/LoH, miscalls, or low-fraction host signal) and leave a clean
    # post-trim GoF on a poorly-fitting sample. Recompute only when trimmed.
    if any(not m.included for m in result.per_marker):
        gof_pval_pretrim = _compute_gof_pval(
            result.per_marker,
            rho=rho,
            n_fitted_params=n_fitted,
            error_rate=result.error_rate,
            pretrim=True,
        )
    else:
        gof_pval_pretrim = gof_pval

    if n_informative < min_informative:
        status = "FAIL"
        warnings.append(f"Insufficient informative markers: {n_informative} < {min_informative}")
        diagnosis = _marker_loss_diagnosis(genotypes, n_informative)
        if diagnosis:
            warnings.append(diagnosis)

    if mean_depth < LOW_MEAN_DEPTH_WARN:
        warnings.append(f"Low mean depth: {mean_depth:.0f}x < 100x")

    # Coverage uniformity: a lopsided depth distribution can pass the mean-depth
    # check while most markers are starved (partial capture failure). Warning-only.
    coverage_uniformity = _coverage_uniformity(depths, UNIFORMITY_DEPTH_FRACTION)
    if depths and coverage_uniformity < UNIFORMITY_MIN_PASS_FRACTION:
        warnings.append(
            f"Uneven coverage: only {coverage_uniformity:.0%} of markers exceed "
            f"{UNIFORMITY_DEPTH_FRACTION:.0%} of the mean depth "
            f"(< {UNIFORMITY_MIN_PASS_FRACTION:.0%}); possible partial capture failure"
        )

    per_donor_n_inf = None
    if hasattr(result, "donor_fraction_cis"):
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
        ci_lo, ci_hi = result.donor_fraction_ci
        ci_width = (ci_hi - ci_lo) * 100
        if ci_width > WIDE_CI_WARN_PCT:
            wide_ci = True
            warnings.append(f"Wide confidence interval: {ci_width:.1f}% > 20%")

    # Gate on the worse of the post-trim and pre-trim fits so a sample cannot pass
    # by trimming away its inconvenient markers.
    gof_candidates = [p for p in (gof_pval, gof_pval_pretrim) if p is not None]
    gof_for_review = min(gof_candidates) if gof_candidates else None
    if gof_for_review is not None and gof_for_review < GOF_REVIEW_P:
        poor_gof = True
        if (
            gof_pval_pretrim is not None
            and gof_pval_pretrim < GOF_REVIEW_P
            and (gof_pval is None or gof_pval >= GOF_REVIEW_P)
        ):
            # Post-trim fit looks fine; only the full set fails (the trim removed
            # the misfitting markers).
            warnings.append(
                "Poor model fit on the full marker set (pre-trim "
                f"goodness-of-fit p={gof_pval_pretrim:.1e} < 0.01) that the "
                "robust trim masks — possible genotyping error, CNV, or, at low "
                "host fraction, trimmed host signal"
            )
        else:
            warnings.append(
                "Poor model fit (goodness-of-fit p < 0.01) — "
                "possible genotyping error, CNV, or sample issue"
            )

    # Host-presence vs MLE disagreement: a significant presence test the MLE does
    # not echo is the clinically interesting "host below the MLE's resolution"
    # case. Soft warning only; v1 does not promote status (real-sample operating
    # characteristics still being mapped).
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

    # A two-rho -> shared-rho fallback (``marker_type_overdispersion_fallback``) is
    # deliberately NOT warned: a sparse class is routine and shared-rho is the
    # validated baseline, so it stays a diagnostic field only.

    # A few exclusions are routine; a large fraction signals host CNV/LoH (or
    # genotyping error) and an untrustworthy trimmed fit.
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

    # Sample-identity QC: reference-sample checks plus the admix swap test.
    identity_review = False
    relatedness: list[RelatednessResult] | None = getattr(result, "relatedness", None)
    duplicate_pairs: set[str] = set()
    if relatedness:
        # Duplicate/reuse: two reference samples reading as the same genome.
        # Intrinsically an error (chimerism is meaningless), so a hard FAIL
        # checked unconditionally.
        for rel in relatedness:
            if rel.degree == Relatedness.IDENTICAL:
                status = "FAIL"
                duplicate_pairs.add(rel.pair)
                coef = "" if rel.coefficient is None else f"r={rel.coefficient:.2f}, "
                warnings.append(
                    f"Identical reference samples: {rel.pair} estimate as the "
                    f"same genome ({coef}{rel.n_sites} markers) — sample "
                    "reuse/mislabel, or an identical-twin (syngeneic) donor; "
                    "either way genotype-based chimerism cannot be measured"
                )

    # Relatedness vs declared expectation. Close-declared-but-unrelated (likely
    # swap) is a hard FAIL; milder mismatches go to REVIEW. Skip duplicate pairs.
    if relatedness and expected_relatedness:
        # Declarations are per donor; align to the host-vs-donor pairs that lead
        # `relatedness` in donor order. strict=True makes a count mismatch an error.
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

    # In-data contamination floor: low-level third-party signal the swap test
    # misses. Magnitude-gated warning; REVIEW when large enough to bias host
    # detection.
    contamination_review = False
    contamination: ContaminationResult | None = getattr(result, "contamination", None)
    if (
        contamination is not None
        and contamination.n_markers > 0
        and contamination.p_value < CONTAMINATION_P
        and contamination.contamination_fraction >= CONTAMINATION_WARN_FRACTION
    ):
        warnings.append(
            f"Contamination: third-party signal at {contamination.contamination_fraction:.4%} "
            f"above the sequencing-error floor at "
            f"{contamination.n_markers} consensus-homozygous markers "
            f"(contamination_p={contamination.p_value:.2e}); possible index hopping or "
            "cross-contamination, limits low-fraction host detection"
        )
        if contamination.contamination_fraction >= CONTAMINATION_REVIEW_FRACTION:
            contamination_review = True

    # Shared-het allele balance: at sites het in host and every donor the admix VAF
    # sits near 0.5 whatever the mixing fraction, so a raised imbalanced fraction is
    # orthogonal signal (contamination, CNV/allelic imbalance, or a sample mix-up)
    # the consensus-hom checks miss. REVIEW above the gate; needs enough sites.
    shared_het_review = False
    shb: SharedHetBalanceResult | None = getattr(result, "shared_het_balance", None)
    if (
        shb is not None
        and shb.n_shared_het >= MIN_SHARED_HET
        and shb.imbalanced_fraction >= SHARED_HET_REVIEW_FRACTION
    ):
        shared_het_review = True
        warnings.append(
            f"Shared-het allele imbalance: {shb.n_imbalanced}/{shb.n_shared_het} "
            f"markers het in all parties fall outside the "
            f"{shb.band:.0%}:{1 - shb.band:.0%} balance band "
            f"(imbalanced={shb.imbalanced_fraction:.0%}, pooled VAF "
            f"{shb.pooled_vaf:.3f}); possible contamination, CNV/allelic imbalance, "
            "or sample mix-up"
        )

    # Index-hopping provenance: shares a sequencing run unit with the host, so
    # hopped host reads could leak in. Soft warning only (a risk, not a defect;
    # the contamination estimate above measures whether it bit). Silent when the
    # admix VCF carried no run metadata.
    run_unit: RunUnitInfo | None = getattr(result, "run_unit", None)
    if run_unit is not None and run_unit.shares_run_with_host:
        warnings.append(
            f"Index-hopping risk: shares sequencing run unit "
            f"{run_unit.run_unit} with the host, so barcode hopping could leak "
            "host reads into this sample; cross-check the contamination estimate"
        )

    # A computed-but-questionable result (poor fit, imprecise, heavily trimmed, or
    # a softer identity flag) is flagged for review rather than passed or failed.
    if status != "FAIL" and (
        poor_gof
        or wide_ci
        or high_robust_drop
        or identity_review
        or contamination_review
        or shared_het_review
    ):
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
        goodness_of_fit_pval_pretrim=gof_pval_pretrim,
        warnings=warnings,
        status=status,
        coverage_uniformity=coverage_uniformity,
        per_donor_n_informative=per_donor_n_inf,
        relatedness=relatedness,
        admix_consistency=ac,
        shared_het_balance=shb,
        contamination=contamination,
        run_unit=run_unit,
    )
