"""Tests for allomix.chimerism — MLE chimerism estimation."""

import math
import random

import numpy as np
import pytest

from allomix import chimerism
from allomix.chimerism import (
    ChimerismResult,
    MarkerResult,
    PanelCalibration,
    _precompute_marker_arrays,
    _total_ll_vec,
    detection_limit,
    estimate_multi_donor,
    estimate_single_donor_bb,
    expected_weight,
    fraction_se,
    log_likelihood_marker_bb,
)
from allomix.genotype import InformativeMarker
from allomix.simulate import (
    expected_vaf,
    generate_marker_biases_realistic,
    sample_allele_counts,
    sample_marker_depths,
)

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


def _make_markers_overdispersed(
    f_donor: float,
    n_markers: int,
    dp: int,
    seed: int = 42,
    bias_sd: float = 0.02,
    depth_cv: float = 0.4,
) -> list[InformativeMarker]:
    """Create markers with realistic overdispersion sources.

    Unlike _make_markers_for_fraction (pure binomial), this adds:
    - Per-marker amplification bias drawn from a heavy-tailed distribution
    - Per-marker depth variability (log-normal with given CV)

    These are the noise sources that cause the binomial CI to fail.
    """
    rng = random.Random(seed)
    biases = generate_marker_biases_realistic(n_markers, rng, sd=bias_sd * 0.6)
    depths = sample_marker_depths(n_markers, dp, depth_cv, rng)

    markers = []
    for i in range(n_markers):
        # Alternate type-0 and type-1
        if i % 2 == 0:
            h_gt, d_gt, mtype = (0, 0), (1, 1), 0
        else:
            h_gt, d_gt, mtype = (1, 1), (0, 0), 1

        vaf = expected_vaf(h_gt, d_gt, f_donor) + biases[i]
        vaf = max(0.0, min(1.0, vaf))
        ref_count, alt_count = sample_allele_counts(vaf, depths[i], rng)

        markers.append(
            _make_marker(
                host_gt=h_gt,
                donor_gt=d_gt,
                ad_ref=ref_count,
                ad_alt=alt_count,
                marker_type=mtype,
                chrom=f"chr{(i % 22) + 1}",
                pos=1_000_000 + i * 100_000,
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
# Tests: beta-binomial likelihood
# ---------------------------------------------------------------------------


class TestBetaBinomialLikelihood:
    """Test beta-binomial likelihood functions."""

    def test_lower_rho_gives_flatter_likelihood(self) -> None:
        """Lower rho (more overdispersion) should flatten the likelihood.

        The LL difference between the best and a wrong w value should
        be smaller when rho is small (flat surface) vs large (sharp peak).
        """
        ad_ref, ad_alt, e = 700, 300, 0.01

        # Sharp (high rho, binomial-like)
        diff_sharp = log_likelihood_marker_bb(
            ad_ref, ad_alt, 0.7, e, rho=10000
        ) - log_likelihood_marker_bb(ad_ref, ad_alt, 0.5, e, rho=10000)
        # Flat (low rho, overdispersed)
        diff_flat = log_likelihood_marker_bb(
            ad_ref, ad_alt, 0.7, e, rho=10
        ) - log_likelihood_marker_bb(ad_ref, ad_alt, 0.5, e, rho=10)

        assert abs(diff_flat) < abs(diff_sharp), (
            "Lower rho should produce a flatter likelihood surface"
        )

    def test_zero_reads_returns_zero(self) -> None:
        assert log_likelihood_marker_bb(0, 0, 0.5, 0.01, 100.0) == 0.0

    def test_ll_negative_for_nonzero_reads(self) -> None:
        """LL should be negative for any non-zero read counts."""
        ll = log_likelihood_marker_bb(700, 300, 0.7, 0.01, 100.0)
        assert math.isfinite(ll)
        assert ll < 0.0

    def test_ll_maximised_near_true_w(self) -> None:
        """LL should be higher when w matches the data."""
        ad_ref, ad_alt, e = 700, 300, 0.01
        rho = 100.0

        ll_good = log_likelihood_marker_bb(ad_ref, ad_alt, 0.7, e, rho)
        ll_bad = log_likelihood_marker_bb(ad_ref, ad_alt, 0.3, e, rho)
        assert ll_good > ll_bad


# ---------------------------------------------------------------------------
# Tests: estimate_single_donor_bb
# ---------------------------------------------------------------------------


class TestEstimateSingleDonorBB:
    """Test MLE estimation with beta-binomial model."""

    def test_ten_percent_chimerism(self) -> None:
        """10 type-0 markers at f=0.10, dp=2000 -> estimate near 0.10."""
        markers = _make_markers_for_fraction(0.10, n_markers=10, dp=2000, seed=42)
        result = estimate_single_donor_bb(markers)

        assert result.donor_fraction == pytest.approx(0.10, abs=0.02)
        assert result.host_fraction == pytest.approx(1.0 - result.donor_fraction, rel=1e-10)
        assert result.n_informative == 10
        assert result.n_markers_used <= 10
        assert len(result.per_marker) == 10

    def test_pure_host(self) -> None:
        """f=0.0 (pure host) -> estimate near 0.0."""
        markers = _make_markers_for_fraction(0.0, n_markers=10, dp=2000, seed=42)
        result = estimate_single_donor_bb(markers)

        assert result.donor_fraction == pytest.approx(0.0, abs=0.005)

    def test_pure_donor(self) -> None:
        """f=1.0 (pure donor) -> estimate near 1.0."""
        markers = _make_markers_for_fraction(1.0, n_markers=10, dp=2000, seed=42)
        result = estimate_single_donor_bb(markers)

        assert result.donor_fraction == pytest.approx(1.0, abs=0.005)

    def test_fifty_percent(self) -> None:
        """f=0.50 -> estimate near 0.50."""
        markers = _make_markers_for_fraction(0.50, n_markers=20, dp=2000, seed=77)
        result = estimate_single_donor_bb(markers)

        assert result.donor_fraction == pytest.approx(0.50, abs=0.03)

    def test_ci_contains_true_value(self) -> None:
        """95% CI should contain the true fraction for a well-behaved seed."""
        true_f = 0.15
        markers = _make_markers_for_fraction(true_f, n_markers=30, dp=2000, seed=123)
        result = estimate_single_donor_bb(markers)

        lo, hi = result.donor_fraction_ci
        assert lo <= true_f <= hi, f"CI [{lo:.4f}, {hi:.4f}] does not contain true f={true_f}"

    def test_ci_narrows_with_depth(self) -> None:
        """CI should be narrower at higher depth."""
        true_f = 0.15

        markers_low = _make_markers_for_fraction(true_f, n_markers=10, dp=500, seed=42)
        result_low = estimate_single_donor_bb(markers_low)
        ci_width_low = result_low.donor_fraction_ci[1] - result_low.donor_fraction_ci[0]

        markers_high = _make_markers_for_fraction(true_f, n_markers=10, dp=5000, seed=42)
        result_high = estimate_single_donor_bb(markers_high)
        ci_width_high = result_high.donor_fraction_ci[1] - result_high.donor_fraction_ci[0]

        assert ci_width_high < ci_width_low

    def test_per_marker_results(self) -> None:
        """Per-marker results should have correct structure and reasonable values."""
        markers = _make_markers_for_fraction(0.10, n_markers=5, dp=2000, seed=42)
        result = estimate_single_donor_bb(markers)

        for mr in result.per_marker:
            assert isinstance(mr, MarkerResult)
            assert 0.0 <= mr.expected_vaf <= 1.0
            assert 0.0 <= mr.observed_vaf <= 1.0
            assert mr.dp == mr.ad_ref + mr.ad_alt
            assert mr.dp > 0

    def test_log_likelihood_is_finite(self) -> None:
        """MLE log-likelihood should be finite and negative."""
        markers = _make_markers_for_fraction(0.10, n_markers=10, dp=2000, seed=42)
        result = estimate_single_donor_bb(markers)

        assert math.isfinite(result.log_likelihood)
        assert result.log_likelihood < 0.0

    def test_empty_markers(self) -> None:
        """Empty marker list should return zero-fraction result."""
        result = estimate_single_donor_bb([])

        assert result.donor_fraction == 0.0
        assert result.host_fraction == 1.0
        assert result.n_informative == 0
        assert len(result.per_marker) == 0

    def test_single_marker(self) -> None:
        """Should work with a single marker (edge case)."""
        marker = _make_marker((0, 0), (1, 1), ad_ref=1800, ad_alt=200, marker_type=0)
        result = estimate_single_donor_bb([marker])

        assert result.n_informative == 1
        assert result.donor_fraction == pytest.approx(0.10, abs=0.02)

    def test_result_type(self) -> None:
        """Verify the return type is ChimerismResult."""
        markers = _make_markers_for_fraction(0.10, n_markers=5, dp=1000, seed=42)
        result = estimate_single_donor_bb(markers)
        assert isinstance(result, ChimerismResult)
        assert isinstance(result.donor_fraction_ci, tuple)
        assert len(result.donor_fraction_ci) == 2


# ---------------------------------------------------------------------------
# Tests: BB estimator on overdispersed data
# ---------------------------------------------------------------------------


class TestBetaBinomialEstimator:
    """Smoke tests for estimate_single_donor_bb on overdispersed data.

    Monte Carlo CI coverage validation is in scripts/benchmark_ci_models.py.
    These tests verify basic correctness without the multi-minute runtime.
    """

    def test_bb_on_overdispersed_data(self) -> None:
        """BB estimator should produce a reasonable estimate on overdispersed data."""
        markers = _make_markers_overdispersed(0.20, n_markers=30, dp=2000, seed=42)
        res = estimate_single_donor_bb(markers)

        assert res.donor_fraction == pytest.approx(0.20, abs=0.05)
        lo, hi = res.donor_fraction_ci
        assert lo < hi
        assert hi - lo < 0.10  # CI should not be absurdly wide

    def test_bb_ci_wider_on_overdispersed_than_clean(self) -> None:
        """BB CIs should be wider on overdispersed data than clean data.

        This tests that rho adapts to the overdispersion level.
        """
        clean_widths = []
        od_widths = []
        for seed in range(5):
            clean = _make_markers_for_fraction(0.20, n_markers=20, dp=2000, seed=seed + 200)
            od = _make_markers_overdispersed(0.20, n_markers=20, dp=2000, seed=seed + 200)
            res_clean = estimate_single_donor_bb(clean)
            res_od = estimate_single_donor_bb(od)
            clean_widths.append(res_clean.donor_fraction_ci[1] - res_clean.donor_fraction_ci[0])
            od_widths.append(res_od.donor_fraction_ci[1] - res_od.donor_fraction_ci[0])

        assert sum(od_widths) / len(od_widths) > sum(clean_widths) / len(clean_widths), (
            "BB CIs should be wider on overdispersed data (rho should adapt)"
        )


# ---------------------------------------------------------------------------
# Tests: per-sample analytical detection limit (LoB / LoD)
# ---------------------------------------------------------------------------


class TestDetectionLimit:
    """Tests for fraction_se and detection_limit."""

    def test_blank_sample_has_finite_positive_lod(self) -> None:
        """A pure-host sample still yields a finite, positive LoB/LoD."""
        markers = _make_markers_for_fraction(0.0, n_markers=40, dp=2000, seed=7)
        lob, lod = detection_limit(markers)
        assert 0.0 < lob < lod < 1.0
        assert math.isfinite(lod)

    def test_lod_reported_on_result(self) -> None:
        """estimate_single_donor_bb populates lob_fraction/lod_fraction."""
        markers = _make_markers_for_fraction(0.0, n_markers=40, dp=2000, seed=7)
        result = estimate_single_donor_bb(markers)
        assert math.isfinite(result.lod_fraction)
        assert 0.0 < result.lob_fraction <= result.lod_fraction

    def test_lod_decreases_with_more_markers(self) -> None:
        """More informative markers should lower the detection limit."""
        few = detection_limit(_make_markers_for_fraction(0.0, n_markers=10, dp=1000, seed=1))[1]
        many = detection_limit(_make_markers_for_fraction(0.0, n_markers=100, dp=1000, seed=1))[1]
        assert many < few

    def test_lod_decreases_with_depth(self) -> None:
        """Higher depth should lower the detection limit."""
        shallow = detection_limit(_make_markers_for_fraction(0.0, n_markers=40, dp=250, seed=2))[1]
        deep = detection_limit(_make_markers_for_fraction(0.0, n_markers=40, dp=4000, seed=2))[1]
        assert deep < shallow

    def test_overdispersion_widens_lod(self) -> None:
        """Smaller rho (more overdispersion) should raise the detection limit."""
        markers = _make_markers_for_fraction(0.0, n_markers=40, dp=2000, seed=3)
        lod_clean = detection_limit(markers, rho=float("inf"))[1]
        lod_od = detection_limit(markers, rho=50.0)[1]
        assert lod_od > lod_clean

    def test_empty_markers_infinite_lod(self) -> None:
        """No markers means nothing is detectable."""
        lob, lod = detection_limit([])
        assert math.isinf(lob)
        assert math.isinf(lod)
        # And the estimator surfaces inf on its result.
        result = estimate_single_donor_bb([])
        assert math.isinf(result.lod_fraction)

    def test_se_positive_and_finite(self) -> None:
        """fraction_se returns a positive finite SE for informative markers."""
        markers = _make_markers_for_fraction(0.0, n_markers=30, dp=1000, seed=9)
        se = fraction_se(markers, 0.0)
        assert math.isfinite(se)
        assert se > 0.0


# ---------------------------------------------------------------------------
# Vectorized likelihood numeric-equivalence tests
# ---------------------------------------------------------------------------


def _scalar_total(markers, f, err, rho, biases):
    """Independent re-derivation of the scalar reference total log-likelihood.

    Does not call total_log_likelihood_bb, which now routes through the
    vectorized path, so this stays an independent check.
    """
    ll = 0.0
    for m in markers:
        bias = biases.get((m.chrom, m.pos, m.ref, m.alt), 0.0) if biases else 0.0
        w = expected_weight(m.host_gt, m.donor_gts[0], f, bias=bias)
        ll += log_likelihood_marker_bb(m.admix_ad_ref, m.admix_ad_alt, w, err, rho)
    return ll


def _mk(chrom, pos, hgt, dgt, ref_n, alt_n):
    return InformativeMarker(
        chrom=chrom,
        pos=pos,
        ref="A",
        alt="G",
        host_gt=hgt,
        donor_gts=[dgt],
        marker_type=0,
        admix_ad_ref=ref_n,
        admix_ad_alt=alt_n,
        admix_dp=ref_n + alt_n,
    )


def test_vectorized_ll_matches_scalar() -> None:
    """Vectorized total LL matches the scalar reference across an (f, rho) grid."""
    markers = [
        _mk("chr1", 100, (0, 0), (0, 1), 480, 20),
        _mk("chr1", 200, (0, 1), (1, 1), 250, 250),
        _mk("chr1", 300, (1, 1), (0, 0), 15, 985),
        _mk("chr1", 400, (0, 0), (1, 1), 0, 0),  # dropout: n==0 -> 0 contribution
    ]
    biases = {("chr1", 200, "A", "G"): 0.03}  # one biased marker
    err = 0.01
    arr = _precompute_marker_arrays(markers, PanelCalibration(biases=biases))
    max_diff = 0.0
    for f in np.linspace(0.0, 1.0, 51):
        for rho in (1.0, 10.0, 50.0, 200.0, 1000.0, 50000.0):
            s = _scalar_total(markers, float(f), err, rho, biases)
            v = _total_ll_vec(arr, float(f), err, rho)
            max_diff = max(max_diff, abs(s - v))
    assert max_diff < 1e-6


def test_vectorized_ll_no_biases_path() -> None:
    """Vectorized LL matches scalar on the no-bias code path."""
    markers = [_mk("chr1", 100, (0, 0), (0, 1), 480, 20)]
    arr = _precompute_marker_arrays(markers, PanelCalibration())
    assert (
        abs(
            _total_ll_vec(arr, 0.05, 0.01, 200.0)
            - _scalar_total(markers, 0.05, 0.01, 200.0, None)
        )
        < 1e-9
    )


# ---------------------------------------------------------------------------
# Tests: robust refit (host CNV / LoH mitigation)
# ---------------------------------------------------------------------------


def _make_contaminated(f_donor, n_clean, n_bad, dp=1000, seed=0):
    """Clean type-0 markers at VAF=f_donor plus n_bad gross outliers (VAF~0.9)."""
    rng = random.Random(seed)
    markers = []
    for i in range(n_clean):
        alt = sum(1 for _ in range(dp) if rng.random() < f_donor)
        markers.append(_make_marker((0, 0), (1, 1), dp - alt, alt, 0, f"chr{i+1}", 1000))
    for j in range(n_bad):
        alt = sum(1 for _ in range(dp) if rng.random() < 0.9)  # CNV/LoH-like outlier
        markers.append(_make_marker((0, 0), (1, 1), dp - alt, alt, 0, f"chr{j+1}", 9000))
    return markers


class TestRobustRefit:
    def test_invalid_mode_raises(self) -> None:
        markers = _make_markers_for_fraction(0.3, 30, 1000)
        with pytest.raises(ValueError, match="robust must be one of"):
            estimate_single_donor_bb(markers, robust="bogus")

    def test_off_is_default(self) -> None:
        markers = _make_markers_for_fraction(0.3, 30, 1000)
        res = estimate_single_donor_bb(markers)
        assert res.n_robust_excluded == 0
        assert res.robust_drop_fraction == 0.0

    def test_clean_data_unchanged(self) -> None:
        """On clean data, robust auto returns the same estimate (no exclusions)."""
        markers = _make_markers_for_fraction(0.3, 60, 1000, seed=7)
        std = estimate_single_donor_bb(markers, robust="off")
        rob = estimate_single_donor_bb(markers, robust="auto")
        assert rob.n_robust_excluded == 0
        assert rob.donor_fraction == pytest.approx(std.donor_fraction, abs=1e-9)

    def test_recovers_from_contamination(self) -> None:
        markers = _make_contaminated(0.3, 52, 8, seed=1)
        std = estimate_single_donor_bb(markers, robust="off")
        rob = estimate_single_donor_bb(markers, robust="auto")
        # Standard is pulled off; robust recovers close to truth and drops the bad ones.
        assert abs(rob.donor_fraction - 0.3) < abs(std.donor_fraction - 0.3)
        assert rob.donor_fraction == pytest.approx(0.3, abs=0.02)
        assert rob.n_robust_excluded == 8
        assert rob.robust_drop_fraction == pytest.approx(8 / 60, abs=1e-6)
        # All original markers are still reported; the 8 outliers are excluded.
        assert rob.n_informative == 60
        assert rob.n_markers_used == 52
        excluded = [m for m in rob.per_marker if not m.included]
        assert len(excluded) == 8
        assert all(m.pos == 9000 for m in excluded)

    def test_min_marker_floor(self) -> None:
        """'auto' will not trim a tiny panel below the floor."""
        markers = _make_contaminated(0.3, 6, 4, seed=2)  # 10 markers, < floor
        rob = estimate_single_donor_bb(markers, robust="auto")
        assert rob.n_robust_excluded == 0  # floor protects the small panel

    def test_multi_donor_robust_runs(self) -> None:
        markers = _make_markers_for_fraction(0.3, 30, 1000)
        res = estimate_multi_donor(markers, n_donors=2, robust="auto")
        assert res.n_robust_excluded == 0
        with pytest.raises(ValueError, match="robust must be one of"):
            estimate_multi_donor(markers, n_donors=2, robust="bogus")

    def test_one_sided_keeps_host_direction_marker_in_fit(self) -> None:
        """A marker deviating toward host presence is kept in the fit by the
        one-sided trim but dropped from the fit by the symmetric trim.

        host(1,1)/donor(0,0) markers have expected ALT VAF = f_host, so a marker
        with ALT VAF well above f_host carries excess host signal (the
        host-present direction). The symmetric MAD cut sees a large residual and
        trims it from the fit (pulling the host estimate down); the one-sided cut
        must keep it, so the host estimate is not biased low.
        """
        f_host = 0.05
        rng = random.Random(3)
        dp = 1000
        markers = [
            _make_marker((1, 1), (0, 0), dp - alt, alt, 0, f"chr{i + 1}", 1000)
            for i, alt in enumerate(
                sum(1 for _ in range(dp) if rng.random() < f_host) for _ in range(30)
            )
        ]
        # One marker carrying far more host ALT than the fit expects.
        big_alt = int(dp * (f_host + 0.15))
        markers.append(_make_marker((1, 1), (0, 0), dp - big_alt, big_alt, 0, "chrZ", 9000))

        saved = chimerism.ROBUST_ONE_SIDED
        try:
            chimerism.ROBUST_ONE_SIDED = False
            sym = estimate_single_donor_bb(markers, robust="force")
            chimerism.ROBUST_ONE_SIDED = True
            one = estimate_single_donor_bb(markers, robust="force")
        finally:
            chimerism.ROBUST_ONE_SIDED = saved

        # Symmetric trims the host-signal marker from the fit; one-sided keeps it.
        assert sym.n_robust_excluded >= 1
        assert one.n_robust_excluded == 0
        # Keeping it leaves the host estimate higher (less biased low).
        assert (1.0 - one.donor_fraction) > (1.0 - sym.donor_fraction)

    def test_one_sided_still_trims_anti_host_outlier(self) -> None:
        """The one-sided trim does not protect outliers pointing away from host.

        The contamination fixture's bad markers (host(0,0)/donor(1,1) at VAF~0.9)
        deviate away from the host-present direction, so one-sided trims them
        exactly like the symmetric refit and recovers the true fraction.
        """
        markers = _make_contaminated(0.3, 52, 8, seed=1)
        saved = chimerism.ROBUST_ONE_SIDED
        try:
            chimerism.ROBUST_ONE_SIDED = True
            rob = estimate_single_donor_bb(markers, robust="auto")
        finally:
            chimerism.ROBUST_ONE_SIDED = saved
        assert rob.n_robust_excluded == 8
        assert rob.donor_fraction == pytest.approx(0.3, abs=0.02)
