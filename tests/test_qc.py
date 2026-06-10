"""Tests for allomix.qc — quality control assessment."""

import math
import random
import re

import pytest

import allomix
from allomix.chimerism import estimate_single_donor_bb
from allomix.contamination import ContaminationResult
from allomix.genotype import InformativeMarker, MarkerCounts, MarkerGenotypes
from allomix.qc import (
    ChimerismResult,
    MarkerResult,
    _compute_gof_pval,
    _marker_loss_diagnosis,
    assess_quality,
)
from allomix.relatedness import (
    DEGREE_FIRST,
    DEGREE_IDENTICAL,
    DEGREE_UNRELATED,
    AdmixConsistencyResult,
    RelatednessResult,
)
from allomix.runmeta import RunUnitInfo

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_marker_result(
    chrom: str = "chr1",
    pos: int = 100,
    marker_type: int = 0,
    expected_vaf: float = 0.10,
    observed_vaf: float = 0.10,
    residual: float = 0.0,
    ad_ref: int = 900,
    ad_alt: int = 100,
    dp: int = 1000,
    included: bool = True,
) -> MarkerResult:
    return MarkerResult(
        chrom=chrom,
        pos=pos,
        marker_type=marker_type,
        expected_vaf=expected_vaf,
        observed_vaf=observed_vaf,
        residual=residual,
        ad_ref=ad_ref,
        ad_alt=ad_alt,
        dp=dp,
        included=included,
    )


def _make_chimerism_result(
    donor_fraction: float = 0.10,
    ci: tuple[float, float] = (0.08, 0.12),
    n_informative: int = 10,
    per_marker: list[MarkerResult] | None = None,
    rho: float = float("inf"),
) -> ChimerismResult:
    if per_marker is None:
        per_marker = [_make_marker_result(pos=i * 100, dp=1000) for i in range(n_informative)]
    n_used = sum(1 for m in per_marker if m.included)
    return ChimerismResult(
        donor_fraction=donor_fraction,
        donor_fraction_ci=ci,
        host_fraction=1.0 - donor_fraction,
        log_likelihood=-100.0,
        n_informative=n_informative,
        n_markers_used=n_used,
        per_marker=per_marker,
        error_rate=0.01,
        rho=rho,
    )


def _make_genotypes(
    n_total: int = 76,
    n_shared: int = 50,
    n_filtered: int = 5,
) -> MarkerGenotypes:
    return MarkerGenotypes(
        informative=[],
        non_informative=[],
        n_total=n_total,
        n_shared=n_shared,
        n_filtered=n_filtered,
        sample_name="test_sample",
    )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestInsufficientMarkers:
    """n_informative < min_informative should cause pass_=False."""

    def test_too_few_markers_fails(self):
        result = _make_chimerism_result(
            n_informative=2,
            per_marker=[
                _make_marker_result(pos=100),
                _make_marker_result(pos=200),
            ],
        )
        genotypes = _make_genotypes()
        qc = assess_quality(result, genotypes, min_informative=3)
        assert qc.pass_ is False

    def test_warning_message_present(self):
        result = _make_chimerism_result(
            n_informative=1,
            per_marker=[
                _make_marker_result(pos=100),
            ],
        )
        genotypes = _make_genotypes()
        qc = assess_quality(result, genotypes, min_informative=3)
        assert any("Insufficient" in w for w in qc.warnings)

    def test_exact_threshold_passes(self):
        result = _make_chimerism_result(
            n_informative=3, per_marker=[_make_marker_result(pos=i * 100) for i in range(3)]
        )
        genotypes = _make_genotypes()
        qc = assess_quality(result, genotypes, min_informative=3)
        assert qc.pass_ is True


class TestContaminationFlag:
    """QC reads result.contamination and warns / REVIEWs by magnitude."""

    @staticmethod
    def _contamination(frac: float, p: float = 1e-12, n: int = 100) -> ContaminationResult:
        return ContaminationResult(
            n_markers=n,
            contamination_fraction=frac,
            median_minor_frac=frac,
            error_floor=0.0,
            floor_empirical=True,
            pooled_minor_frac=frac,
            n_minor_reads=int(frac * n * 2000),
            total_depth=n * 2000,
            p_value=p,
            n_excluded_high=0,
            used_per_site_error=False,
            error_rate_source="global-fallback",
        )

    def test_clean_no_warning(self):
        result = _make_chimerism_result()
        result.contamination = self._contamination(0.0, p=1.0)
        qc = assess_quality(result, _make_genotypes())
        assert not any("Contamination" in w for w in qc.warnings)
        assert qc.status == "PASS"

    def test_moderate_contamination_warns_not_review(self):
        result = _make_chimerism_result()
        result.contamination = self._contamination(0.004)  # 0.4%: warn, below REVIEW
        qc = assess_quality(result, _make_genotypes())
        assert any("Contamination" in w for w in qc.warnings)
        assert qc.status == "PASS"

    def test_high_contamination_promotes_review(self):
        result = _make_chimerism_result()
        result.contamination = self._contamination(0.015)  # 1.5%: REVIEW
        qc = assess_quality(result, _make_genotypes())
        assert any("Contamination" in w for w in qc.warnings)
        assert qc.status == "REVIEW"

    def test_significant_but_tiny_not_flagged(self):
        # p is significant (high depth) but the magnitude is below the warn gate.
        result = _make_chimerism_result()
        result.contamination = self._contamination(0.0005, p=1e-9)
        qc = assess_quality(result, _make_genotypes())
        assert not any("Contamination" in w for w in qc.warnings)

    def test_no_markers_not_flagged(self):
        result = _make_chimerism_result()
        result.contamination = self._contamination(0.0, p=1.0, n=0)
        qc = assess_quality(result, _make_genotypes())
        assert not any("Contamination" in w for w in qc.warnings)


class TestIndexHoppingFlag:
    """QC reads result.run_unit and warns (softly) when it shares the host run."""

    def test_shared_run_warns_without_review(self):
        result = _make_chimerism_result()
        result.run_unit = RunUnitInfo(
            run_unit="FC1:1", source="RG:PU", shares_run_with_host=True
        )
        qc = assess_quality(result, _make_genotypes())
        assert any("Index-hopping" in w and "FC1:1" in w for w in qc.warnings)
        assert qc.status == "PASS"  # soft warning only

    def test_not_shared_no_warning(self):
        result = _make_chimerism_result()
        result.run_unit = RunUnitInfo(
            run_unit="FC2:2", source="RG:PU", shares_run_with_host=False
        )
        qc = assess_quality(result, _make_genotypes())
        assert not any("Index-hopping" in w for w in qc.warnings)

    def test_undetermined_no_warning(self):
        result = _make_chimerism_result()
        result.run_unit = RunUnitInfo(
            run_unit="FC3:3", source="RG:PU", shares_run_with_host=None
        )
        qc = assess_quality(result, _make_genotypes())
        assert not any("Index-hopping" in w for w in qc.warnings)

    def test_absent_metadata_no_warning(self):
        result = _make_chimerism_result()  # run_unit defaults to None
        qc = assess_quality(result, _make_genotypes())
        assert not any("Index-hopping" in w for w in qc.warnings)
        assert qc.run_unit is None


class TestMarkerLossDiagnosis:
    """_marker_loss_diagnosis should name the input that starved the markers."""

    def test_no_funnel_data_returns_empty(self):
        # Hand-built genotypes with default (zero) funnel fields.
        assert _marker_loss_diagnosis(_make_genotypes(), 1) == ""

    def test_sparse_donor_genotyping(self):
        g = MarkerGenotypes(
            informative=[],
            non_informative=[],
            n_total=64,
            n_shared=8,
            n_filtered=0,
            marker_counts=MarkerCounts(
                n_host=64, n_donor_markers=[9], n_admix=64, n_admix_in_host=64, n_admix_in_donor=[9]
            ),
        )
        msg = _marker_loss_diagnosis(g, 1)
        assert "donor" in msg and "9/64" in msg

    def test_low_admixture_depth(self):
        g = MarkerGenotypes(
            informative=[],
            non_informative=[],
            n_total=64,
            n_shared=60,
            n_filtered=58,
            marker_counts=MarkerCounts(
                n_host=64,
                n_donor_markers=[62],
                n_admix=64,
                n_admix_in_host=64,
                n_admix_in_donor=[62],
                n_drop_admix_dp=58,
            ),
        )
        msg = _marker_loss_diagnosis(g, 2)
        assert "admixture depth" in msg

    def test_sparse_admixture_sample(self):
        g = MarkerGenotypes(
            informative=[],
            non_informative=[],
            n_total=5,
            n_shared=5,
            n_filtered=0,
            marker_counts=MarkerCounts(
                n_host=70, n_donor_markers=[68], n_admix=5, n_admix_in_host=5, n_admix_in_donor=[5]
            ),
        )
        msg = _marker_loss_diagnosis(g, 1)
        assert "admixture sample has only 5" in msg

    def test_diagnosis_appears_in_warnings(self):
        result = _make_chimerism_result(n_informative=1, per_marker=[_make_marker_result()])
        g = MarkerGenotypes(
            informative=[],
            non_informative=[],
            n_total=64,
            n_shared=8,
            n_filtered=0,
            marker_counts=MarkerCounts(
                n_host=64, n_donor_markers=[9], n_admix=64, n_admix_in_host=64, n_admix_in_donor=[9]
            ),
        )
        qc = assess_quality(result, g, min_informative=3)
        assert any("donor genotyping covers only" in w for w in qc.warnings)


class TestLowDepth:
    """Mean depth < 100 should produce a warning."""

    def test_low_depth_warning(self):
        markers = [_make_marker_result(pos=i * 100, dp=50, ad_ref=45, ad_alt=5) for i in range(5)]
        result = _make_chimerism_result(n_informative=5, per_marker=markers)
        genotypes = _make_genotypes()
        qc = assess_quality(result, genotypes)
        assert any("depth" in w.lower() for w in qc.warnings)

    def test_adequate_depth_no_warning(self):
        markers = [_make_marker_result(pos=i * 100, dp=500) for i in range(5)]
        result = _make_chimerism_result(n_informative=5, per_marker=markers)
        genotypes = _make_genotypes()
        qc = assess_quality(result, genotypes)
        assert not any("depth" in w.lower() for w in qc.warnings)


class TestWideCi:
    """CI width > 20% should produce a warning."""

    def test_wide_ci_warning(self):
        result = _make_chimerism_result(ci=(0.01, 0.35))
        genotypes = _make_genotypes()
        qc = assess_quality(result, genotypes)
        assert any("confidence interval" in w.lower() for w in qc.warnings)
        assert qc.status == "REVIEW"

    def test_narrow_ci_no_warning(self):
        result = _make_chimerism_result(ci=(0.09, 0.11))
        genotypes = _make_genotypes()
        qc = assess_quality(result, genotypes)
        assert not any("confidence interval" in w.lower() for w in qc.warnings)
        assert qc.status == "PASS"


class TestGoodData:
    """Well-formed data should produce pass_=True with no warnings."""

    def test_good_data_passes(self):
        markers = [
            _make_marker_result(
                pos=i * 100,
                dp=1500,
                ad_ref=1350,
                ad_alt=150,
                residual=0.005,
            )
            for i in range(20)
        ]
        result = _make_chimerism_result(
            donor_fraction=0.10,
            ci=(0.08, 0.12),
            n_informative=20,
            per_marker=markers,
        )
        genotypes = _make_genotypes()
        qc = assess_quality(result, genotypes)
        assert qc.pass_ is True
        assert qc.status == "PASS"
        assert len(qc.warnings) == 0

    def test_counts_are_correct(self):
        markers = [_make_marker_result(pos=i * 100, dp=1000) for i in range(10)]
        # Mark two as excluded
        markers[0].included = False
        markers[1].included = False
        result = _make_chimerism_result(
            n_informative=10,
            per_marker=markers,
        )
        genotypes = _make_genotypes(n_total=76, n_shared=50, n_filtered=5)
        qc = assess_quality(result, genotypes)
        assert qc.n_total_markers == 76
        assert qc.n_shared_markers == 50
        assert qc.n_informative == 10
        assert qc.n_used == 8
        assert qc.n_excluded_outlier == 2
        assert qc.n_excluded_depth == 5


class TestGoodnessOfFit:
    """Bad residuals should trigger a GOF warning."""

    def test_bad_residuals_warning(self):
        markers = [
            _make_marker_result(
                pos=i * 100,
                dp=1000,
                residual=5.0,
                included=True,
            )
            for i in range(10)
        ]
        result = _make_chimerism_result(
            n_informative=10,
            per_marker=markers,
        )
        genotypes = _make_genotypes()
        qc = assess_quality(result, genotypes)
        assert qc.goodness_of_fit_pval is not None
        assert qc.goodness_of_fit_pval < 0.01
        assert any("model fit" in w.lower() for w in qc.warnings)
        assert qc.status == "REVIEW"
        assert qc.pass_ is True  # REVIEW is not a hard fail

    def test_small_residuals_no_warning(self):
        markers = [
            _make_marker_result(
                pos=i * 100,
                dp=1000,
                residual=0.005,
                included=True,
            )
            for i in range(10)
        ]
        result = _make_chimerism_result(
            n_informative=10,
            per_marker=markers,
        )
        genotypes = _make_genotypes()
        qc = assess_quality(result, genotypes)
        assert qc.goodness_of_fit_pval is not None
        assert qc.goodness_of_fit_pval >= 0.01
        assert not any("model fit" in w.lower() for w in qc.warnings)

    def test_gof_none_with_too_few_included(self):
        markers = [
            _make_marker_result(pos=100, included=True),
        ]
        result = _make_chimerism_result(
            n_informative=1,
            per_marker=markers,
        )
        genotypes = _make_genotypes()
        qc = assess_quality(result, genotypes, min_informative=1)
        assert qc.goodness_of_fit_pval is None

    def test_overdispersed_residuals_not_flagged_under_bb(self):
        """With a small rho capturing overdispersion, residuals that fail
        under binomial variance should pass under beta-binomial variance."""
        markers = [
            _make_marker_result(pos=i * 100, dp=1000, residual=0.03, included=True)
            for i in range(20)
        ]
        result = _make_chimerism_result(
            n_informative=20,
            per_marker=markers,
            rho=5.0,
        )
        genotypes = _make_genotypes()
        qc = assess_quality(result, genotypes)
        assert qc.goodness_of_fit_pval is not None
        assert qc.goodness_of_fit_pval > 0.01
        assert not any("model fit" in w.lower() for w in qc.warnings)

    def test_same_residuals_fail_under_binomial(self):
        """Sanity: those same residuals fail under the binomial limit."""
        markers = [
            _make_marker_result(pos=i * 100, dp=1000, residual=0.03, included=True)
            for i in range(20)
        ]
        result = _make_chimerism_result(
            n_informative=20,
            per_marker=markers,
            rho=float("inf"),
        )
        genotypes = _make_genotypes()
        qc = assess_quality(result, genotypes)
        assert qc.goodness_of_fit_pval is not None
        assert qc.goodness_of_fit_pval < 0.01


# ---------------------------------------------------------------------------
# Pearson GoF statistic tests
# ---------------------------------------------------------------------------


def _make_informative_marker(
    host_gt,
    donor_gt,
    ad_ref,
    ad_alt,
    marker_type=0,
    chrom="chr1",
    pos=100,
):
    return InformativeMarker(
        chrom=chrom,
        pos=pos,
        ref="A",
        alt="T",
        host_gt=host_gt,
        donor_gts=[donor_gt],
        marker_type=marker_type,
        admix_ad_ref=ad_ref,
        admix_ad_alt=ad_alt,
        admix_dp=ad_ref + ad_alt,
    )


class TestPearsonGoF:
    """Pearson chi-squared GoF should detect model–data mismatch."""

    def test_detects_systematic_vaf_shift(self):
        """A 5% systematic shift at dp=2000 across 20 markers should be detected."""
        markers = [
            _make_marker_result(
                pos=i * 100,
                dp=2000,
                expected_vaf=0.10,
                observed_vaf=0.15,
                residual=0.05,
                ad_ref=1700,
                ad_alt=300,
                included=True,
            )
            for i in range(20)
        ]
        pval = _compute_gof_pval(markers)
        assert pval is not None
        assert pval < 0.01

    def test_detects_single_outlier_marker(self):
        """One marker with residual=0.40 among 19 good ones should be detected."""
        good = [
            _make_marker_result(pos=i * 100, dp=2000, residual=0.001, included=True)
            for i in range(19)
        ]
        bad = [
            _make_marker_result(
                pos=1900,
                dp=2000,
                expected_vaf=0.10,
                observed_vaf=0.50,
                residual=0.40,
                included=True,
            )
        ]
        pval = _compute_gof_pval(good + bad)
        assert pval is not None
        assert pval < 0.01

    def test_calibration_under_null(self):
        """Under the null, p should not always be ~1.0."""
        rng = random.Random(42)
        markers = []
        for i in range(30):
            exp_vaf = 0.10
            dp = 2000
            sd = math.sqrt(exp_vaf * (1 - exp_vaf) / dp)
            obs_vaf = max(0.0, min(1.0, rng.gauss(exp_vaf, sd)))
            markers.append(
                _make_marker_result(
                    pos=i * 100,
                    dp=dp,
                    expected_vaf=exp_vaf,
                    observed_vaf=obs_vaf,
                    residual=obs_vaf - exp_vaf,
                    included=True,
                )
            )
        pval = _compute_gof_pval(markers)
        assert pval is not None
        assert pval < 0.999, f"GoF p-value is {pval:.6f} — not calibrated"


class TestGoFEndToEnd:
    """GoF should catch a swapped genotype through the full estimation pipeline."""

    def test_catches_wrong_genotype(self):
        """A single swapped genotype should be detected end-to-end, either
        by the 3-SD outlier rule (excluding it from the fit) or by a GoF
        failure. Under the beta-binomial likelihood a lone outlier is
        normally caught by the outlier rule; GoF only flags it when the
        outlier rule fails to exclude it."""
        rng = random.Random(42)
        true_f = 0.20
        dp = 2000
        markers = []

        for i in range(19):
            alt_count = sum(1 for _ in range(dp) if rng.random() < true_f)
            markers.append(
                _make_informative_marker(
                    host_gt=(0, 0),
                    donor_gt=(1, 1),
                    ad_ref=dp - alt_count,
                    ad_alt=alt_count,
                    marker_type=0,
                    chrom=f"chr{i + 1}",
                    pos=1000 * (i + 1),
                )
            )

        # One marker with swapped host/donor genotypes
        alt_count = sum(1 for _ in range(dp) if rng.random() < true_f)
        markers.append(
            _make_informative_marker(
                host_gt=(1, 1),
                donor_gt=(0, 0),
                ad_ref=dp - alt_count,
                ad_alt=alt_count,
                marker_type=1,
                chrom="chr20",
                pos=20000,
            )
        )

        result = estimate_single_donor_bb(markers)
        genotypes = MarkerGenotypes(
            informative=[],
            non_informative=[],
            n_total=20,
            n_shared=20,
            n_filtered=0,
        )
        qc = assess_quality(result, genotypes)
        outlier_excluded = result.n_markers_used < result.n_informative
        gof_fail = qc.goodness_of_fit_pval is not None and qc.goodness_of_fit_pval < 0.01
        assert outlier_excluded or gof_fail, (
            f"Swapped genotype not detected: n_used={result.n_markers_used}/"
            f"{result.n_informative}, gof_pval={qc.goodness_of_fit_pval}"
        )


# ---------------------------------------------------------------------------
# Version consistency
# ---------------------------------------------------------------------------


class TestVersionConsistency:
    def test_versions_match(self):
        """pyproject.toml version should match __init__.py __version__."""
        with open("pyproject.toml", encoding="utf-8") as f:
            content = f.read()
        match = re.search(r'^version\s*=\s*"([^"]+)"', content, re.MULTILINE)
        assert match is not None, "Could not find version in pyproject.toml"
        assert match.group(1) == allomix.__version__


class TestRobustExclusionReview:
    """A large robust-refit exclusion fraction should warn and flag REVIEW."""

    def test_high_drop_promotes_review(self):
        per_marker = [_make_marker_result(pos=i * 100) for i in range(20)]
        result = _make_chimerism_result(n_informative=20, per_marker=per_marker)
        result.n_robust_excluded = 5
        result.robust_drop_fraction = 0.25  # > ROBUST_REVIEW_FRACTION (0.15)
        qc = assess_quality(result, _make_genotypes())
        assert qc.status == "REVIEW"
        assert any("Robust refit excluded" in w for w in qc.warnings)

    def test_small_drop_warns_but_passes(self):
        per_marker = [_make_marker_result(pos=i * 100) for i in range(40)]
        result = _make_chimerism_result(n_informative=40, per_marker=per_marker)
        result.n_robust_excluded = 1
        result.robust_drop_fraction = 0.025  # below review threshold
        qc = assess_quality(result, _make_genotypes())
        assert qc.status == "PASS"
        assert any("Robust refit excluded" in w for w in qc.warnings)

    def test_no_robust_no_warning(self):
        result = _make_chimerism_result(n_informative=30)
        qc = assess_quality(result, _make_genotypes())
        assert not any("Robust refit" in w for w in qc.warnings)


def _rel(degree: int, coefficient: float) -> RelatednessResult:
    """A RelatednessResult with a fixed degree, for QC-integration tests."""
    from allomix.relatedness import DEGREE_LABELS

    return RelatednessResult(
        a_name="host",
        b_name="donor",
        coefficient=coefficient,
        ci_low=coefficient - 0.05,
        ci_high=coefficient + 0.05,
        confidence="high",
        relationship=DEGREE_LABELS[degree],
        degree=degree,
        n_sites=80,
        het_a=60,
        het_b=60,
        shared_hets=40,
        ibs0=2,
    )


class TestRelatednessQC:
    """Declared-vs-detected relatedness drives the QC status."""

    def test_boundary_mismatch_fails(self):
        result = _make_chimerism_result(n_informative=30)
        result.relatedness = [_rel(DEGREE_UNRELATED, 0.03)]
        qc = assess_quality(result, _make_genotypes(), expected_relatedness=["first-degree"])
        assert qc.status == "FAIL"
        assert any("relatedness check FAIL" in w for w in qc.warnings)

    def test_within_tolerance_passes(self):
        result = _make_chimerism_result(n_informative=30)
        result.relatedness = [_rel(DEGREE_FIRST, 0.5)]
        qc = assess_quality(result, _make_genotypes(), expected_relatedness=["first-degree"])
        assert qc.status == "PASS"

    def test_no_expectation_no_verdict(self):
        result = _make_chimerism_result(n_informative=30)
        result.relatedness = [_rel(DEGREE_FIRST, 0.5)]
        qc = assess_quality(result, _make_genotypes(), expected_relatedness=["NA"])
        assert qc.status == "PASS"
        assert not any("relatedness check" in w for w in qc.warnings)
        # The estimate is still carried on the report even without a verdict.
        assert qc.relatedness is not None

    def test_relatedness_attached_to_report(self):
        result = _make_chimerism_result(n_informative=30)
        result.relatedness = [_rel(DEGREE_FIRST, 0.5)]
        qc = assess_quality(result, _make_genotypes())
        assert qc.relatedness == result.relatedness

    def test_count_mismatch_raises(self):
        # More declarations than host-vs-donor pairs: strict zip errors rather
        # than silently leaving a donor unchecked.
        result = _make_chimerism_result(n_informative=30)
        result.relatedness = [_rel(DEGREE_FIRST, 0.5)]  # one host-vs-donor pair
        with pytest.raises(ValueError):
            assess_quality(
                result,
                _make_genotypes(),
                expected_relatedness=["first-degree", "unrelated"],
            )

    def test_gaining_relatedness_reviews_not_fails(self):
        # Declared unrelated but a non-identical relationship detected: REVIEW,
        # not FAIL (not a random-swap signature).
        result = _make_chimerism_result(n_informative=30)
        result.relatedness = [_rel(DEGREE_FIRST, 0.5)]
        qc = assess_quality(result, _make_genotypes(), expected_relatedness=["unrelated"])
        assert qc.status == "REVIEW"
        assert not any("FAIL" in w for w in qc.warnings)

    def test_duplicate_fails_without_declaration(self):
        # An identical reference pair is an unconditional FAIL with a clear
        # message, even when no expected relationship is declared.
        result = _make_chimerism_result(n_informative=30)
        result.relatedness = [_rel(DEGREE_IDENTICAL, 1.0)]
        qc = assess_quality(result, _make_genotypes())
        assert qc.status == "FAIL"
        assert any("Identical reference samples" in w for w in qc.warnings)

    def test_duplicate_reported_once_with_declaration(self):
        # When identical and a declaration is given, only the duplicate message
        # is emitted (the declared-verdict message is skipped).
        result = _make_chimerism_result(n_informative=30)
        result.relatedness = [_rel(DEGREE_IDENTICAL, 1.0)]
        qc = assess_quality(result, _make_genotypes(), expected_relatedness=["unrelated"])
        assert qc.status == "FAIL"
        assert any("Identical reference samples" in w for w in qc.warnings)
        assert not any("relatedness check" in w for w in qc.warnings)


class TestAdmixSwapQC:
    """A significant consensus-homozygote swap test promotes REVIEW."""

    def test_swap_promotes_review(self):
        result = _make_chimerism_result(n_informative=30)
        result.admix_consistency = AdmixConsistencyResult(
            n_consensus_hom=40, n_discordant=40, discordant_fraction=1.0, swap_pval=1e-30
        )
        qc = assess_quality(result, _make_genotypes())
        assert qc.status == "REVIEW"
        assert any("Possible sample swap" in w for w in qc.warnings)

    def test_clean_admix_no_warning(self):
        result = _make_chimerism_result(n_informative=30)
        result.admix_consistency = AdmixConsistencyResult(
            n_consensus_hom=40, n_discordant=0, discordant_fraction=0.0, swap_pval=1.0
        )
        qc = assess_quality(result, _make_genotypes())
        assert qc.status == "PASS"
        assert not any("sample swap" in w for w in qc.warnings)

    def test_too_few_consensus_sites_no_action(self):
        result = _make_chimerism_result(n_informative=30)
        # Significant but below MIN_CONSENSUS -> not acted on.
        result.admix_consistency = AdmixConsistencyResult(
            n_consensus_hom=5, n_discordant=5, discordant_fraction=1.0, swap_pval=1e-10
        )
        qc = assess_quality(result, _make_genotypes())
        assert qc.status == "PASS"
