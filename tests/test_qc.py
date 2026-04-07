"""Tests for allomix.qc — quality control assessment."""

from __future__ import annotations

from allomix.genotype import MarkerGenotypes
from allomix.qc import ChimerismResult, MarkerResult, assess_quality

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
        per_marker = [
            _make_marker_result(pos=i * 100, dp=1000)
            for i in range(n_informative)
        ]
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
        result = _make_chimerism_result(n_informative=2, per_marker=[
            _make_marker_result(pos=100),
            _make_marker_result(pos=200),
        ])
        genotypes = _make_genotypes()
        qc = assess_quality(result, genotypes, min_informative=3)
        assert qc.pass_ is False

    def test_warning_message_present(self):
        result = _make_chimerism_result(n_informative=1, per_marker=[
            _make_marker_result(pos=100),
        ])
        genotypes = _make_genotypes()
        qc = assess_quality(result, genotypes, min_informative=3)
        assert any("Insufficient" in w for w in qc.warnings)

    def test_exact_threshold_passes(self):
        result = _make_chimerism_result(n_informative=3, per_marker=[
            _make_marker_result(pos=i * 100) for i in range(3)
        ])
        genotypes = _make_genotypes()
        qc = assess_quality(result, genotypes, min_informative=3)
        assert qc.pass_ is True


class TestLowDepth:
    """Mean depth < 100 should produce a warning."""

    def test_low_depth_warning(self):
        markers = [
            _make_marker_result(pos=i * 100, dp=50, ad_ref=45, ad_alt=5)
            for i in range(5)
        ]
        result = _make_chimerism_result(n_informative=5, per_marker=markers)
        genotypes = _make_genotypes()
        qc = assess_quality(result, genotypes)
        assert any("depth" in w.lower() for w in qc.warnings)

    def test_adequate_depth_no_warning(self):
        markers = [
            _make_marker_result(pos=i * 100, dp=500)
            for i in range(5)
        ]
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
                residual=0.1,
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
        markers = [
            _make_marker_result(pos=i * 100, dp=1000)
            for i in range(10)
        ]
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
                residual=0.1,
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
