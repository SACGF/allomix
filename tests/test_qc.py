"""Tests for allomix.qc — quality control assessment."""

from __future__ import annotations

import math
import random
import re

from allomix.chimerism import estimate_single_donor_bb
from allomix.genotype import InformativeMarker, MarkerGenotypes
from allomix.qc import ChimerismResult, MarkerResult, _compute_gof_pval, assess_quality

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

    def test_narrow_ci_no_warning(self):
        result = _make_chimerism_result(ci=(0.09, 0.11))
        genotypes = _make_genotypes()
        qc = assess_quality(result, genotypes)
        assert not any("confidence interval" in w.lower() for w in qc.warnings)


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
                pos=i * 100, dp=2000, expected_vaf=0.10, observed_vaf=0.15,
                residual=0.05, ad_ref=1700, ad_alt=300, included=True,
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
                pos=1900, dp=2000, expected_vaf=0.10, observed_vaf=0.50,
                residual=0.40, included=True,
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
                    pos=i * 100, dp=dp, expected_vaf=exp_vaf,
                    observed_vaf=obs_vaf, residual=obs_vaf - exp_vaf, included=True,
                )
            )
        pval = _compute_gof_pval(markers)
        assert pval is not None
        assert pval < 0.999, f"GoF p-value is {pval:.6f} — not calibrated"


class TestGoFEndToEnd:
    """GoF should catch a swapped genotype through the full estimation pipeline."""

    def test_catches_wrong_genotype(self):
        rng = random.Random(42)
        true_f = 0.20
        dp = 2000
        markers = []

        for i in range(19):
            alt_count = sum(1 for _ in range(dp) if rng.random() < true_f)
            markers.append(
                _make_informative_marker(
                    host_gt=(0, 0), donor_gt=(1, 1),
                    ad_ref=dp - alt_count, ad_alt=alt_count,
                    marker_type=0, chrom=f"chr{i + 1}", pos=1000 * (i + 1),
                )
            )

        # One marker with swapped host/donor genotypes
        alt_count = sum(1 for _ in range(dp) if rng.random() < true_f)
        markers.append(
            _make_informative_marker(
                host_gt=(1, 1), donor_gt=(0, 0),
                ad_ref=dp - alt_count, ad_alt=alt_count,
                marker_type=1, chrom="chr20", pos=20000,
            )
        )

        result = estimate_single_donor_bb(markers)
        genotypes = MarkerGenotypes(
            informative=[], non_informative=[], n_total=20, n_shared=20, n_filtered=0,
        )
        qc = assess_quality(result, genotypes)
        assert qc.goodness_of_fit_pval is not None
        assert qc.goodness_of_fit_pval < 0.01


# ---------------------------------------------------------------------------
# Version consistency
# ---------------------------------------------------------------------------


class TestVersionConsistency:
    def test_versions_match(self):
        """pyproject.toml version should match __init__.py __version__."""
        import allomix

        with open("pyproject.toml") as f:
            content = f.read()
        match = re.search(r'^version\s*=\s*"([^"]+)"', content, re.MULTILINE)
        assert match is not None, "Could not find version in pyproject.toml"
        assert match.group(1) == allomix.__version__
