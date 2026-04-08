"""Tests for allomix.chimerism — MLE chimerism estimation."""

from __future__ import annotations

import math
import random

import pytest

from allomix.chimerism import (
    ChimerismResult,
    MarkerResult,
    estimate_error_rate,
    estimate_single_donor,
    expected_weight,
    log_likelihood_marker,
    total_log_likelihood,
)
from allomix.genotype import InformativeMarker

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_marker(
    host_gt: tuple[int, int],
    donor_gt: tuple[int, int],
    ad_ref: int,
    ad_alt: int,
    marker_type: int = 0,
    chrom: str = "chr1",
    pos: int = 100,
) -> InformativeMarker:
    """Create an InformativeMarker for testing."""
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


def _make_markers_for_fraction(
    f_donor: float,
    n_markers: int,
    dp: int,
    seed: int = 42,
) -> list[InformativeMarker]:
    """Create type-0 markers (host 0/0, donor 1/1) with simulated allele counts.

    For type-0 markers: expected ALT VAF = f_donor.
    """
    rng = random.Random(seed)
    markers = []
    for i in range(n_markers):
        # Expected REF weight = 1 - f_donor, so ALT count ~ Binomial(dp, f_donor)
        alt_count = sum(1 for _ in range(dp) if rng.random() < f_donor)
        ref_count = dp - alt_count
        markers.append(
            _make_marker(
                host_gt=(0, 0),
                donor_gt=(1, 1),
                ad_ref=ref_count,
                ad_alt=alt_count,
                marker_type=0,
                chrom=f"chr{i + 1}",
                pos=1000 * (i + 1),
            )
        )
    return markers


# ---------------------------------------------------------------------------
# Tests: expected_weight
# ---------------------------------------------------------------------------


class TestExpectedWeight:
    """Test expected_weight for all 9 GT combinations at f=0, 0.5, 1.0.

    expected_weight computes the REF allele weight:
        w = (1-f) * host_ref_dose/2 + f * donor_ref_dose/2
    where ref_dose = 2 - alt_dose.
    """

    @pytest.mark.parametrize(
        "host_gt, donor_gt, f, expected_w",
        [
            # host 0/0 (ref_dose=2), donor 0/0 (ref_dose=2) -> always 1.0
            ((0, 0), (0, 0), 0.0, 1.0),
            ((0, 0), (0, 0), 0.5, 1.0),
            ((0, 0), (0, 0), 1.0, 1.0),
            # host 0/0 (ref_dose=2), donor 0/1 (ref_dose=1)
            # w = (1-f)*1.0 + f*0.5
            ((0, 0), (0, 1), 0.0, 1.0),
            ((0, 0), (0, 1), 0.5, 0.75),
            ((0, 0), (0, 1), 1.0, 0.5),
            # host 0/0 (ref_dose=2), donor 1/1 (ref_dose=0)
            # w = (1-f)*1.0 + f*0.0 = 1-f
            ((0, 0), (1, 1), 0.0, 1.0),
            ((0, 0), (1, 1), 0.5, 0.5),
            ((0, 0), (1, 1), 1.0, 0.0),
            # host 0/1 (ref_dose=1), donor 0/0 (ref_dose=2)
            # w = (1-f)*0.5 + f*1.0
            ((0, 1), (0, 0), 0.0, 0.5),
            ((0, 1), (0, 0), 0.5, 0.75),
            ((0, 1), (0, 0), 1.0, 1.0),
            # host 0/1, donor 0/1 -> always 0.5
            ((0, 1), (0, 1), 0.0, 0.5),
            ((0, 1), (0, 1), 0.5, 0.5),
            ((0, 1), (0, 1), 1.0, 0.5),
            # host 0/1 (ref_dose=1), donor 1/1 (ref_dose=0)
            # w = (1-f)*0.5 + f*0.0 = 0.5*(1-f)
            ((0, 1), (1, 1), 0.0, 0.5),
            ((0, 1), (1, 1), 0.5, 0.25),
            ((0, 1), (1, 1), 1.0, 0.0),
            # host 1/1 (ref_dose=0), donor 0/0 (ref_dose=2)
            # w = (1-f)*0.0 + f*1.0 = f
            ((1, 1), (0, 0), 0.0, 0.0),
            ((1, 1), (0, 0), 0.5, 0.5),
            ((1, 1), (0, 0), 1.0, 1.0),
            # host 1/1 (ref_dose=0), donor 0/1 (ref_dose=1)
            # w = (1-f)*0.0 + f*0.5 = 0.5*f
            ((1, 1), (0, 1), 0.0, 0.0),
            ((1, 1), (0, 1), 0.5, 0.25),
            ((1, 1), (0, 1), 1.0, 0.5),
            # host 1/1, donor 1/1 -> always 0.0
            ((1, 1), (1, 1), 0.0, 0.0),
            ((1, 1), (1, 1), 0.5, 0.0),
            ((1, 1), (1, 1), 1.0, 0.0),
        ],
    )
    def test_expected_weight(
        self,
        host_gt: tuple[int, int],
        donor_gt: tuple[int, int],
        f: float,
        expected_w: float,
    ) -> None:
        result = expected_weight(host_gt, donor_gt, f)
        assert result == pytest.approx(expected_w, abs=1e-10)


# ---------------------------------------------------------------------------
# Tests: log_likelihood_marker
# ---------------------------------------------------------------------------


class TestLogLikelihoodMarker:
    """Test log_likelihood_marker with hand-checked values."""

    def test_perfect_ref_only(self) -> None:
        """All reads are REF and w=1.0 -> highest LL."""
        ll = log_likelihood_marker(ad_ref=1000, ad_alt=0, w=1.0, error_rate=0.01)
        # p_ref = 1.0*(1-0.01) + 0.0*0.01/3 = 0.99
        # LL = 1000 * log(0.99)
        expected = 1000 * math.log(0.99)
        assert ll == pytest.approx(expected, rel=1e-10)

    def test_perfect_alt_only(self) -> None:
        """All reads are ALT and w=0.0 -> highest LL for ALT."""
        ll = log_likelihood_marker(ad_ref=0, ad_alt=1000, w=0.0, error_rate=0.01)
        # p_alt = 1.0*(1-0.01) + 0.0*0.01/3 = 0.99
        # LL = 1000 * log(0.99)
        expected = 1000 * math.log(0.99)
        assert ll == pytest.approx(expected, rel=1e-10)

    def test_balanced_at_half(self) -> None:
        """Equal REF/ALT reads at w=0.5."""
        ll = log_likelihood_marker(ad_ref=500, ad_alt=500, w=0.5, error_rate=0.01)
        # p_ref = 0.5*(0.99) + 0.5*(0.01/3) = 0.495 + 0.001667 = 0.496667
        # p_alt = 0.5*(0.99) + 0.5*(0.01/3) = 0.496667
        # LL = 500*log(0.496667) + 500*log(0.496667) = 1000*log(0.496667)
        p = 0.5 * 0.99 + 0.5 * 0.01 / 3.0
        expected = 1000 * math.log(p)
        assert ll == pytest.approx(expected, rel=1e-10)

    def test_mismatch_penalty(self) -> None:
        """Alt reads when w=1.0 should have lower LL than ref reads when w=1.0."""
        ll_good = log_likelihood_marker(ad_ref=1000, ad_alt=0, w=1.0, error_rate=0.01)
        ll_bad = log_likelihood_marker(ad_ref=0, ad_alt=1000, w=1.0, error_rate=0.01)
        assert ll_good > ll_bad

    def test_zero_reads(self) -> None:
        """No reads should give LL = 0."""
        ll = log_likelihood_marker(ad_ref=0, ad_alt=0, w=0.5, error_rate=0.01)
        assert ll == 0.0

    def test_error_rate_prevents_log_zero(self) -> None:
        """Even with w=0, REF reads still get a small probability from error."""
        ll = log_likelihood_marker(ad_ref=100, ad_alt=0, w=0.0, error_rate=0.01)
        # p_ref = 0.0*(0.99) + 1.0*(0.01/3) = 0.003333
        expected = 100 * math.log(0.01 / 3.0)
        assert ll == pytest.approx(expected, rel=1e-10)
        assert math.isfinite(ll)

    def test_known_hand_calculation(self) -> None:
        """Verify with a fully worked example.

        w=0.8, e=0.02, ad_ref=800, ad_alt=200.
        p_ref = 0.8*0.98 + 0.2*0.02/3 = 0.784 + 0.001333 = 0.785333
        p_alt = 0.2*0.98 + 0.8*0.02/3 = 0.196 + 0.005333 = 0.201333
        LL = 800*log(0.785333) + 200*log(0.201333)
        """
        e = 0.02
        w = 0.8
        p_ref = w * (1 - e) + (1 - w) * e / 3
        p_alt = (1 - w) * (1 - e) + w * e / 3
        expected = 800 * math.log(p_ref) + 200 * math.log(p_alt)
        ll = log_likelihood_marker(ad_ref=800, ad_alt=200, w=w, error_rate=e)
        assert ll == pytest.approx(expected, rel=1e-10)


# ---------------------------------------------------------------------------
# Tests: total_log_likelihood
# ---------------------------------------------------------------------------


class TestTotalLogLikelihood:
    """Test total_log_likelihood consistency."""

    def test_sum_matches_individual(self) -> None:
        """Total LL should equal sum of individual marker LLs."""
        markers = [
            _make_marker((0, 0), (1, 1), 900, 100, marker_type=0, chrom="chr1", pos=100),
            _make_marker((0, 0), (1, 1), 800, 200, marker_type=0, chrom="chr2", pos=200),
            _make_marker((1, 1), (0, 0), 100, 900, marker_type=1, chrom="chr3", pos=300),
        ]

        f = 0.10
        e = 0.01
        total = total_log_likelihood(markers, f, e)

        individual_sum = 0.0
        for m in markers:
            w = expected_weight(m.host_gt, m.donor_gts[0], f)
            individual_sum += log_likelihood_marker(m.admix_ad_ref, m.admix_ad_alt, w, e)

        assert total == pytest.approx(individual_sum, rel=1e-10)

    def test_empty_markers(self) -> None:
        """Empty marker list should give LL = 0."""
        assert total_log_likelihood([], 0.5) == 0.0

    def test_maximum_at_true_fraction(self) -> None:
        """LL should be highest near the true donor fraction."""
        markers = _make_markers_for_fraction(0.20, n_markers=20, dp=2000, seed=99)

        ll_at_true = total_log_likelihood(markers, 0.20)
        ll_at_wrong_lo = total_log_likelihood(markers, 0.05)
        ll_at_wrong_hi = total_log_likelihood(markers, 0.50)

        assert ll_at_true > ll_at_wrong_lo
        assert ll_at_true > ll_at_wrong_hi


# ---------------------------------------------------------------------------
# Tests: estimate_single_donor
# ---------------------------------------------------------------------------


class TestEstimateSingleDonor:
    """Test MLE estimation with synthetic data."""

    def test_ten_percent_chimerism(self) -> None:
        """10 type-0 markers at f=0.10, dp=2000 -> estimate near 0.10."""
        markers = _make_markers_for_fraction(0.10, n_markers=10, dp=2000, seed=42)
        result = estimate_single_donor(markers)

        assert result.donor_fraction == pytest.approx(0.10, abs=0.02)
        assert result.host_fraction == pytest.approx(1.0 - result.donor_fraction, rel=1e-10)
        assert result.n_informative == 10
        assert result.n_markers_used <= 10
        assert len(result.per_marker) == 10

    def test_pure_host(self) -> None:
        """f=0.0 (pure host) -> estimate near 0.0."""
        markers = _make_markers_for_fraction(0.0, n_markers=10, dp=2000, seed=42)
        result = estimate_single_donor(markers)

        assert result.donor_fraction == pytest.approx(0.0, abs=0.005)
        assert result.host_fraction == pytest.approx(1.0, abs=0.005)

    def test_pure_donor(self) -> None:
        """f=1.0 (pure donor) -> estimate near 1.0."""
        markers = _make_markers_for_fraction(1.0, n_markers=10, dp=2000, seed=42)
        result = estimate_single_donor(markers)

        assert result.donor_fraction == pytest.approx(1.0, abs=0.005)
        assert result.host_fraction == pytest.approx(0.0, abs=0.005)

    def test_fifty_percent(self) -> None:
        """f=0.50 -> estimate near 0.50."""
        markers = _make_markers_for_fraction(0.50, n_markers=20, dp=2000, seed=77)
        result = estimate_single_donor(markers)

        assert result.donor_fraction == pytest.approx(0.50, abs=0.03)

    def test_ci_contains_true_value(self) -> None:
        """95% CI should contain the true fraction in most seeds.

        With 30 markers at dp=2000, the CI is wide enough that this should
        hold reliably for a well-behaved seed.
        """
        true_f = 0.15
        markers = _make_markers_for_fraction(true_f, n_markers=30, dp=2000, seed=123)
        result = estimate_single_donor(markers)

        lo, hi = result.donor_fraction_ci
        assert lo <= true_f <= hi, f"CI [{lo:.4f}, {hi:.4f}] does not contain true f={true_f}"

    def test_ci_narrows_with_depth(self) -> None:
        """CI should be narrower at higher depth."""
        true_f = 0.15

        markers_low = _make_markers_for_fraction(true_f, n_markers=10, dp=500, seed=42)
        result_low = estimate_single_donor(markers_low)
        ci_width_low = result_low.donor_fraction_ci[1] - result_low.donor_fraction_ci[0]

        markers_high = _make_markers_for_fraction(true_f, n_markers=10, dp=5000, seed=42)
        result_high = estimate_single_donor(markers_high)
        ci_width_high = result_high.donor_fraction_ci[1] - result_high.donor_fraction_ci[0]

        assert ci_width_high < ci_width_low, (
            f"High-depth CI width ({ci_width_high:.4f}) should be smaller "
            f"than low-depth CI width ({ci_width_low:.4f})"
        )

    def test_ci_narrows_with_more_markers(self) -> None:
        """CI should be narrower with more informative markers."""
        true_f = 0.15

        markers_few = _make_markers_for_fraction(true_f, n_markers=5, dp=2000, seed=42)
        result_few = estimate_single_donor(markers_few)
        ci_width_few = result_few.donor_fraction_ci[1] - result_few.donor_fraction_ci[0]

        markers_many = _make_markers_for_fraction(true_f, n_markers=40, dp=2000, seed=42)
        result_many = estimate_single_donor(markers_many)
        ci_width_many = result_many.donor_fraction_ci[1] - result_many.donor_fraction_ci[0]

        assert ci_width_many < ci_width_few, (
            f"Many-marker CI width ({ci_width_many:.4f}) should be smaller "
            f"than few-marker CI width ({ci_width_few:.4f})"
        )

    def test_per_marker_results(self) -> None:
        """Per-marker results should have correct structure and reasonable values."""
        markers = _make_markers_for_fraction(0.10, n_markers=5, dp=2000, seed=42)
        result = estimate_single_donor(markers)

        for mr in result.per_marker:
            assert isinstance(mr, MarkerResult)
            assert 0.0 <= mr.expected_vaf <= 1.0
            assert 0.0 <= mr.observed_vaf <= 1.0
            assert mr.dp == mr.ad_ref + mr.ad_alt
            assert mr.dp > 0

    def test_log_likelihood_is_finite(self) -> None:
        """MLE log-likelihood should be finite and negative."""
        markers = _make_markers_for_fraction(0.10, n_markers=10, dp=2000, seed=42)
        result = estimate_single_donor(markers)

        assert math.isfinite(result.log_likelihood)
        assert result.log_likelihood < 0.0

    def test_empty_markers(self) -> None:
        """Empty marker list should return zero-fraction result."""
        result = estimate_single_donor([])

        assert result.donor_fraction == 0.0
        assert result.host_fraction == 1.0
        assert result.n_informative == 0
        assert result.n_markers_used == 0
        assert len(result.per_marker) == 0

    def test_single_marker(self) -> None:
        """Should work with a single marker (edge case)."""
        marker = _make_marker((0, 0), (1, 1), ad_ref=1800, ad_alt=200, marker_type=0)
        result = estimate_single_donor([marker])

        assert result.n_informative == 1
        assert result.donor_fraction == pytest.approx(0.10, abs=0.02)

    def test_all_same_marker_type(self) -> None:
        """Should work when all markers are the same type."""
        # All type-1: host 1/1, donor 0/0 at f=0.20
        # Expected ALT VAF = 1 - f = 0.80, so ad_alt should be ~80%
        rng = random.Random(42)
        markers = []
        for i in range(10):
            ad_alt = sum(1 for _ in range(2000) if rng.random() < 0.80)
            ad_ref = 2000 - ad_alt
            markers.append(
                _make_marker(
                    host_gt=(1, 1),
                    donor_gt=(0, 0),
                    ad_ref=ad_ref,
                    ad_alt=ad_alt,
                    marker_type=1,
                    chrom=f"chr{i + 1}",
                    pos=1000 * (i + 1),
                )
            )

        result = estimate_single_donor(markers)
        assert result.donor_fraction == pytest.approx(0.20, abs=0.02)

    def test_mixed_marker_types(self) -> None:
        """Should work with a mix of type-0 and type-1 markers."""
        rng = random.Random(42)
        true_f = 0.15
        markers = []

        # Type-0: host 0/0, donor 1/1 -> expected ALT VAF = f
        for i in range(5):
            alt_count = sum(1 for _ in range(2000) if rng.random() < true_f)
            markers.append(
                _make_marker(
                    host_gt=(0, 0),
                    donor_gt=(1, 1),
                    ad_ref=2000 - alt_count,
                    ad_alt=alt_count,
                    marker_type=0,
                    chrom=f"chr{i + 1}",
                    pos=1000 * (i + 1),
                )
            )

        # Type-1: host 1/1, donor 0/0 -> expected ALT VAF = 1-f
        for i in range(5):
            alt_count = sum(1 for _ in range(2000) if rng.random() < (1.0 - true_f))
            markers.append(
                _make_marker(
                    host_gt=(1, 1),
                    donor_gt=(0, 0),
                    ad_ref=2000 - alt_count,
                    ad_alt=alt_count,
                    marker_type=1,
                    chrom=f"chr{i + 6}",
                    pos=1000 * (i + 6),
                )
            )

        result = estimate_single_donor(markers)
        assert result.donor_fraction == pytest.approx(true_f, abs=0.03)

    def test_het_markers(self) -> None:
        """Test with partially informative markers (host het, donor hom)."""
        rng = random.Random(42)
        true_f = 0.20
        markers = []

        # Type-20: host 0/0, donor 0/1 -> ref_weight = (1-f)*1.0 + f*0.5
        # Expected ALT VAF = 1 - ref_weight = f * 0.5
        for i in range(10):
            expected_alt_vaf = true_f * 0.5
            alt_count = sum(1 for _ in range(2000) if rng.random() < expected_alt_vaf)
            markers.append(
                _make_marker(
                    host_gt=(0, 0),
                    donor_gt=(0, 1),
                    ad_ref=2000 - alt_count,
                    ad_alt=alt_count,
                    marker_type=20,
                    chrom=f"chr{i + 1}",
                    pos=1000 * (i + 1),
                )
            )

        result = estimate_single_donor(markers)
        # Less precise due to partial informativeness
        assert result.donor_fraction == pytest.approx(true_f, abs=0.05)

    def test_result_type(self) -> None:
        """Verify the return type is ChimerismResult."""
        markers = _make_markers_for_fraction(0.10, n_markers=5, dp=1000, seed=42)
        result = estimate_single_donor(markers)
        assert isinstance(result, ChimerismResult)
        assert isinstance(result.donor_fraction_ci, tuple)
        assert len(result.donor_fraction_ci) == 2


# ---------------------------------------------------------------------------
# Tests: estimate_error_rate
# ---------------------------------------------------------------------------


class TestEstimateErrorRate:
    def test_returns_default(self) -> None:
        """Should return 0.01 in v1."""
        markers = _make_markers_for_fraction(0.10, n_markers=5, dp=1000, seed=42)
        assert estimate_error_rate(markers) == 0.01

    def test_empty_markers(self) -> None:
        assert estimate_error_rate([]) == 0.01


# ---------------------------------------------------------------------------
# Tests: end-to-end with simulate module
# ---------------------------------------------------------------------------


class TestEndToEndWithSimulate:
    """Integration tests using allomix.simulate for data generation."""

    def test_roundtrip_multiple_fractions(self) -> None:
        """Generate data at several fractions and verify MLE recovers them."""
        from allomix.simulate import expected_vaf, sample_allele_counts

        fractions = [0.0, 0.05, 0.10, 0.25, 0.50, 0.75, 1.0]
        dp = 3000

        # Define a set of host/donor GT pairs (all type-0 for simplicity)
        gt_pairs = [((0, 0), (1, 1))] * 10 + [((1, 1), (0, 0))] * 10

        for true_f in fractions:
            rng = random.Random(int(true_f * 1000) + 42)
            markers = []

            for i, (h_gt, d_gt) in enumerate(gt_pairs):
                vaf = expected_vaf(h_gt, d_gt, true_f)
                ref_count, alt_count = sample_allele_counts(vaf, dp, rng)
                from allomix.genotype import marker_type as get_mtype

                mtype = get_mtype(h_gt, d_gt)
                markers.append(
                    InformativeMarker(
                        chrom=f"chr{i + 1}",
                        pos=1000 * (i + 1),
                        ref="A",
                        alt="T",
                        host_gt=h_gt,
                        donor_gts=[d_gt],
                        marker_type=mtype if mtype is not None else 0,
                        admix_ad_ref=ref_count,
                        admix_ad_alt=alt_count,
                        admix_dp=ref_count + alt_count,
                    )
                )

            result = estimate_single_donor(markers)

            # Point estimate should be close to truth
            assert result.donor_fraction == pytest.approx(true_f, abs=0.03), (
                f"At true_f={true_f}: estimated {result.donor_fraction}"
            )

            # CI should be in the neighbourhood of the true value.
            # Due to sampling noise the true value may occasionally fall
            # just outside the 95% CI, so we use a generous tolerance.
            lo, hi = result.donor_fraction_ci
            assert lo <= true_f + 0.01 and hi >= true_f - 0.01, (
                f"At true_f={true_f}: CI [{lo:.4f}, {hi:.4f}]"
            )


class TestProfileLikelihoodCIPrecision:
    """CI width should not be inflated by step-size artifacts."""

    def test_ci_precision_at_high_depth(self):
        """With 50 type-0 markers at dp=5000, CI width should be tight.

        Overdispersion adjustment may widen slightly beyond pure profile
        likelihood (phi can be ~1.2 even for ideal binomial data), so we
        allow up to 0.008 instead of the raw ~0.004.
        """
        rng = random.Random(42)
        true_f = 0.50
        markers = []
        for i in range(50):
            alt_count = sum(1 for _ in range(5000) if rng.random() < true_f)
            markers.append(
                InformativeMarker(
                    chrom=f"chr{i + 1}",
                    pos=1000 * (i + 1),
                    ref="A",
                    alt="T",
                    host_gt=(0, 0),
                    donor_gts=[(1, 1)],
                    marker_type=0,
                    admix_ad_ref=5000 - alt_count,
                    admix_ad_alt=alt_count,
                    admix_dp=5000,
                )
            )
        result = estimate_single_donor(markers)
        ci_width = result.donor_fraction_ci[1] - result.donor_fraction_ci[0]
        assert ci_width < 0.008, f"CI width {ci_width:.6f} unexpectedly wide"


class TestRoundtripWithErrorRate:
    """Round-trip simulate -> estimate with matched error model (e=0.01)."""

    def test_roundtrip_with_error_rate(self) -> None:
        from allomix.simulate import expected_vaf, sample_allele_counts

        error_rate = 0.01
        fractions = [0.0, 0.10, 0.50, 0.90, 1.0]
        dp = 5000
        gt_pairs = [((0, 0), (1, 1))] * 15 + [((1, 1), (0, 0))] * 15

        for true_f in fractions:
            rng = random.Random(int(true_f * 1000) + 99)
            markers = []

            for i, (h_gt, d_gt) in enumerate(gt_pairs):
                vaf = expected_vaf(h_gt, d_gt, true_f)
                ref_count, alt_count = sample_allele_counts(
                    vaf, dp, rng, error_rate=error_rate,
                )
                from allomix.genotype import marker_type as get_mtype

                mtype = get_mtype(h_gt, d_gt)
                markers.append(
                    InformativeMarker(
                        chrom=f"chr{i + 1}",
                        pos=1000 * (i + 1),
                        ref="A",
                        alt="T",
                        host_gt=h_gt,
                        donor_gts=[d_gt],
                        marker_type=mtype if mtype is not None else 0,
                        admix_ad_ref=ref_count,
                        admix_ad_alt=alt_count,
                        admix_dp=ref_count + alt_count,
                    )
                )

            result = estimate_single_donor(markers, error_rate=error_rate)
            assert result.donor_fraction == pytest.approx(true_f, abs=0.03), (
                f"At true_f={true_f}: estimated {result.donor_fraction}"
            )


class TestErrorModelMatchValidation:
    """Validate that simulator + estimator produce unbiased results.

    Generates data with error_rate=0.01 (where the old symmetric model had
    measurable bias) and checks that the MLE is unbiased. Under the old
    symmetric error model, the simulator generated ~1% ALT reads at hom-ref
    sites while the estimator expected ~0.33%, producing a ~0.6% upward bias
    at f=0. With matched models, the bias should be negligible.
    """

    def _run_replicates(
        self, true_f: float, n_replicates: int = 20, seed_offset: int = 0,
    ) -> list[float]:
        from allomix.simulate import expected_vaf, sample_allele_counts

        error_rate = 0.01
        dp = 2000
        gt_pairs = [((0, 0), (1, 1))] * 20 + [((1, 1), (0, 0))] * 20

        estimates = []
        for rep in range(n_replicates):
            rng = random.Random(rep * 1000 + seed_offset)
            markers = []
            for i, (h_gt, d_gt) in enumerate(gt_pairs):
                vaf = expected_vaf(h_gt, d_gt, true_f)
                ref_count, alt_count = sample_allele_counts(
                    vaf, dp, rng, error_rate=error_rate,
                )
                from allomix.genotype import marker_type as get_mtype

                mtype = get_mtype(h_gt, d_gt)
                markers.append(
                    InformativeMarker(
                        chrom=f"chr{i + 1}",
                        pos=1000 * (i + 1),
                        ref="A",
                        alt="T",
                        host_gt=h_gt,
                        donor_gts=[d_gt],
                        marker_type=mtype if mtype is not None else 0,
                        admix_ad_ref=ref_count,
                        admix_ad_alt=alt_count,
                        admix_dp=ref_count + alt_count,
                    )
                )
            result = estimate_single_donor(markers, error_rate=error_rate)
            estimates.append(result.donor_fraction)
        return estimates

    def test_no_bias_at_zero_fraction(self) -> None:
        """At f=0 with error_rate=0.01, estimate should be near 0, not ~0.6%."""
        estimates = self._run_replicates(0.0, n_replicates=20, seed_offset=0)
        mean_estimate = sum(estimates) / len(estimates)
        # With matched models, mean should be near 0 (within sampling noise).
        # Under the old mismatch, this would be ~0.006 (0.6%).
        assert mean_estimate < 0.003, (
            f"Mean estimate at f=0 is {mean_estimate:.4f}; "
            "possible error model mismatch (old symmetric model gave ~0.006)"
        )

    def test_no_bias_at_full_fraction(self) -> None:
        """At f=1.0 with error_rate=0.01, estimate should be near 1.0."""
        estimates = self._run_replicates(1.0, n_replicates=20, seed_offset=500)
        mean_estimate = sum(estimates) / len(estimates)
        assert mean_estimate > 0.997, (
            f"Mean estimate at f=1 is {mean_estimate:.4f}; "
            "possible error model mismatch (old symmetric model gave ~0.994)"
        )
