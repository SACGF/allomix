"""Tests for allomix.calibration.bias — per-marker amplification bias estimation and correction."""

import random

import pytest

from allomix.calibration.bias import (
    MarkerBias,
    biases_to_simple_dict,
    estimate_biases,
    estimate_biases_both_het,
    load_bias_table,
    save_bias_table,
)
from allomix.estimate.chimerism import estimate_single_donor_bb
from allomix.estimate.likelihood import (
    PanelCalibration,
    expected_weight,
    inject_bias,
    total_log_likelihood_bb,
)
from allomix.genotype import InformativeMarker, MarkerData

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_het_marker(
    chrom: str = "chr1",
    pos: int = 100,
    ref: str = "A",
    alt: str = "T",
    ad_ref: int = 500,
    ad_alt: int = 500,
) -> MarkerData:
    """Create a heterozygous MarkerData for bias estimation tests."""
    return MarkerData(
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        gt=(0, 1),
        ad_ref=ad_ref,
        ad_alt=ad_alt,
        dp=ad_ref + ad_alt,
    )


def _make_informative_marker(
    host_gt: tuple[int, int],
    donor_gt: tuple[int, int],
    ad_ref: int,
    ad_alt: int,
    marker_type: int = 0,
    chrom: str = "chr1",
    pos: int = 100,
) -> InformativeMarker:
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


# ---------------------------------------------------------------------------
# estimate_biases
# ---------------------------------------------------------------------------


class TestEstimateBiases:
    def test_no_het_markers_returns_empty(self):
        """Hom-only samples produce no bias estimates."""
        m = MarkerData(
            chrom="chr1",
            pos=100,
            ref="A",
            alt="T",
            gt=(0, 0),
            ad_ref=1000,
            ad_alt=0,
            dp=1000,
        )
        biases = estimate_biases([[m]])
        assert biases == {}

    def test_single_het_unbiased(self):
        """A het marker at VAF=0.5 has bias=0."""
        m = _make_het_marker(ad_ref=500, ad_alt=500)
        biases = estimate_biases([[m]])
        key = ("chr1", 100, "A", "T")
        assert key in biases
        assert biases[key].bias == 0.0
        assert biases[key].n_het == 1

    def test_single_het_positive_bias(self):
        """ALT-favoured marker has positive bias."""
        # VAF = 520/1000 = 0.52, bias = 0.52 - 0.5 = 0.02
        m = _make_het_marker(ad_ref=480, ad_alt=520)
        biases = estimate_biases([[m]])
        key = ("chr1", 100, "A", "T")
        assert abs(biases[key].bias - 0.02) < 1e-9

    def test_single_het_negative_bias(self):
        """REF-favoured marker has negative bias."""
        # VAF = 480/1000 = 0.48, bias = -0.02
        m = _make_het_marker(ad_ref=520, ad_alt=480)
        biases = estimate_biases([[m]])
        key = ("chr1", 100, "A", "T")
        assert abs(biases[key].bias - (-0.02)) < 1e-9

    def test_multiple_samples_uses_median(self):
        """Median across multiple het observations."""
        # Three samples: VAFs of 0.48, 0.50, 0.54
        # Deviations: -0.02, 0.0, 0.04 -> median = 0.0
        markers_lists = [
            [_make_het_marker(ad_ref=520, ad_alt=480)],
            [_make_het_marker(ad_ref=500, ad_alt=500)],
            [_make_het_marker(ad_ref=460, ad_alt=540)],
        ]
        biases = estimate_biases(markers_lists)
        key = ("chr1", 100, "A", "T")
        assert biases[key].n_het == 3
        assert abs(biases[key].bias - 0.0) < 1e-9

    def test_min_het_filter(self):
        """Markers with fewer than min_het observations are excluded."""
        m = _make_het_marker()
        biases = estimate_biases([[m]], min_het=2)
        assert biases == {}

    def test_multiple_markers(self):
        """Different markers get independent bias estimates."""
        m1 = _make_het_marker(chrom="chr1", pos=100, ad_ref=490, ad_alt=510)
        m2 = _make_het_marker(chrom="chr2", pos=200, ad_ref=530, ad_alt=470)
        biases = estimate_biases([[m1, m2]])
        assert len(biases) == 2
        assert biases[("chr1", 100, "A", "T")].bias == pytest.approx(0.01)
        assert biases[("chr2", 200, "A", "T")].bias == pytest.approx(-0.03)

    def test_non_het_markers_ignored(self):
        """Hom-ref and hom-alt markers are not used for bias estimation."""
        het = _make_het_marker(pos=100, ad_ref=490, ad_alt=510)
        hom_ref = MarkerData(
            chrom="chr1",
            pos=200,
            ref="A",
            alt="T",
            gt=(0, 0),
            ad_ref=1000,
            ad_alt=5,
            dp=1005,
        )
        hom_alt = MarkerData(
            chrom="chr1",
            pos=300,
            ref="A",
            alt="T",
            gt=(1, 1),
            ad_ref=5,
            ad_alt=1000,
            dp=1005,
        )
        biases = estimate_biases([[het, hom_ref, hom_alt]])
        assert len(biases) == 1
        assert ("chr1", 100, "A", "T") in biases


class TestEstimateBiasesBothHet:
    """Estimate bias from admix samples at host+donor both-het markers (issue #11)."""

    def _gt_marker(self, pos, gt, ad_ref, ad_alt, chrom="chr1"):
        return MarkerData(
            chrom=chrom, pos=pos, ref="A", alt="T", gt=gt,
            ad_ref=ad_ref, ad_alt=ad_alt, dp=ad_ref + ad_alt,
        )

    def test_only_both_het_markers_used(self):
        """Only markers where host and every donor are het contribute."""
        host = [
            self._gt_marker(100, (0, 1), 500, 500),  # both-het
            self._gt_marker(200, (0, 0), 1000, 0),   # host hom -> excluded
            self._gt_marker(300, (0, 1), 500, 500),  # donor hom -> excluded
        ]
        donors = [[
            self._gt_marker(100, (0, 1), 500, 500),
            self._gt_marker(200, (0, 1), 500, 500),
            self._gt_marker(300, (1, 1), 0, 1000),
        ]]
        # admix VAF 0.53 at the both-het marker -> bias +0.03
        admix = [[
            self._gt_marker(100, (0, 1), 470, 530),
            self._gt_marker(200, (0, 1), 480, 520),
            self._gt_marker(300, (0, 1), 480, 520),
        ]]
        biases = estimate_biases_both_het(host, donors, admix)
        assert list(biases.keys()) == [("chr1", 100, "A", "T")]
        assert biases[("chr1", 100, "A", "T")].bias == pytest.approx(0.03)

    def test_multiple_admix_samples_median(self):
        """Each admix sample contributes one observation; bias is their median."""
        host = [self._gt_marker(100, (0, 1), 500, 500)]
        donors = [[self._gt_marker(100, (0, 1), 500, 500)]]
        admix = [
            [self._gt_marker(100, (0, 1), 490, 510)],  # +0.01
            [self._gt_marker(100, (0, 1), 480, 520)],  # +0.02
            [self._gt_marker(100, (0, 1), 460, 540)],  # +0.04
        ]
        biases = estimate_biases_both_het(host, donors, admix)
        mb = biases[("chr1", 100, "A", "T")]
        assert mb.n_het == 3
        assert mb.bias == pytest.approx(0.02)  # median of 0.01, 0.02, 0.04

    def test_min_het_filter(self):
        """Markers below min_het are dropped."""
        host = [self._gt_marker(100, (0, 1), 500, 500)]
        donors = [[self._gt_marker(100, (0, 1), 500, 500)]]
        admix = [[self._gt_marker(100, (0, 1), 480, 520)]]
        assert estimate_biases_both_het(host, donors, admix, min_het=2) == {}

    def test_low_depth_excluded(self):
        """Observations below min_dp are not counted."""
        host = [self._gt_marker(100, (0, 1), 500, 500)]
        donors = [[self._gt_marker(100, (0, 1), 500, 500)]]
        admix = [[self._gt_marker(100, (0, 1), 5, 5)]]
        assert estimate_biases_both_het(host, donors, admix, min_dp=100) == {}


# ---------------------------------------------------------------------------
# Bias table I/O
# ---------------------------------------------------------------------------


class TestBiasTableIO:
    def test_roundtrip(self, tmp_path):
        """Save and load produces the same biases."""
        biases = {
            ("chr1", 100, "A", "T"): MarkerBias("chr1", 100, "A", "T", 0.015, 5),
            ("chr2", 200, "C", "G"): MarkerBias("chr2", 200, "C", "G", -0.008, 3),
        }
        path = tmp_path / "bias.tsv"
        save_bias_table(biases, path)
        loaded = load_bias_table(path)

        assert len(loaded) == 2
        assert loaded[("chr1", 100, "A", "T")] == pytest.approx(0.015)
        assert loaded[("chr2", 200, "C", "G")] == pytest.approx(-0.008)

    def test_biases_to_simple_dict(self):
        """Convert MarkerBias dict to simple float dict."""
        biases = {
            ("chr1", 100, "A", "T"): MarkerBias("chr1", 100, "A", "T", 0.02, 5),
        }
        simple = biases_to_simple_dict(biases)
        assert simple == {("chr1", 100, "A", "T"): 0.02}


# ---------------------------------------------------------------------------
# expected_weight with bias
# ---------------------------------------------------------------------------


class TestExpectedWeightWithBias:
    def test_no_bias_unchanged(self):
        """With bias=0, expected_weight is unchanged."""
        w = expected_weight((0, 0), (1, 1), 0.3, bias=0.0)
        # host ref_dose=2, donor ref_dose=0 -> w = 0.7*1.0 + 0.3*0.0 = 0.7
        assert w == pytest.approx(0.7)

    def test_positive_bias_decreases_ref_weight(self):
        """Positive bias (ALT-favoured) decreases expected ref weight.

        Correction is multiplicative in logit/odds space (issue #20), not a flat
        0.02 shift: w' = odds(0.7) / (odds(0.7) + odds(0.52)) = (7/3)/(7/3 + 13/12)
        = 28/41.
        """
        w_no_bias = expected_weight((0, 0), (1, 1), 0.3, bias=0.0)
        w_biased = expected_weight((0, 0), (1, 1), 0.3, bias=0.02)
        assert w_biased < w_no_bias
        assert w_biased == pytest.approx(28 / 41)

    def test_negative_bias_increases_ref_weight(self):
        """Negative bias (REF-favoured) increases expected ref weight.

        w' = odds(0.7) / (odds(0.7) + odds(0.48)) = (7/3)/(7/3 + 12/13) = 91/127.
        """
        w_biased = expected_weight((0, 0), (1, 1), 0.3, bias=-0.02)
        assert w_biased > expected_weight((0, 0), (1, 1), 0.3, bias=0.0)
        assert w_biased == pytest.approx(91 / 127)

    def test_het_site_matches_estimate(self):
        """At a het site (w=0.5) the corrected ALT VAF is exactly 0.5 + bias."""
        for bias in (-0.05, -0.02, 0.0, 0.02, 0.05):
            w = expected_weight((0, 1), (1, 1), 0.0, bias=bias)  # host het, f=0 -> w_true=0.5
            assert (1.0 - w) == pytest.approx(0.5 + bias, abs=1e-9)

    def test_injection_matches_model_expectation(self):
        """Simulator injection equals the estimator's expected biased weight.

        At the true parameters the observed (injected) ALT VAF must equal what
        the estimator's bias-corrected model expects, so simulation and
        estimation stay self-consistent.
        """
        from allomix.estimate.likelihood import apply_bias

        for w_true in (0.05, 0.3, 0.5, 0.7, 0.95):
            for bias in (-0.04, 0.0, 0.03):
                alt_observed = float(inject_bias(1.0 - w_true, bias))
                expected_alt = 1.0 - float(apply_bias(w_true, bias))
                assert alt_observed == pytest.approx(expected_alt, abs=1e-9)

    def test_extreme_vaf_not_overcorrected(self):
        """At extreme expected VAF the logit shift is small, unlike additive."""
        # host 0/0, donor 1/1, f=0.03 -> expected ALT VAF 0.03 (w=0.97)
        w = expected_weight((0, 0), (1, 1), 0.03, bias=0.02)
        alt = 1.0 - w
        assert alt == pytest.approx(0.0324, abs=1e-3)  # additive would give 0.05
        assert alt < 0.04

    def test_bias_clamped_near_zero(self):
        """Bias doesn't push expected weight below epsilon."""
        # w_true = 0.0 (host 1/1, donor 1/1 wouldn't be informative, but test the math)
        # With no bias: w=0, bias=0.05 -> w - bias = -0.05 -> clamped to 1e-6
        w = expected_weight((1, 1), (1, 1), 0.0, bias=0.05)
        assert w == pytest.approx(1e-6)

    def test_bias_clamped_near_one(self):
        """Bias doesn't push expected weight above 1-epsilon."""
        w = expected_weight((0, 0), (0, 0), 0.0, bias=-0.05)
        assert w == pytest.approx(1.0 - 1e-6)


# ---------------------------------------------------------------------------
# total_log_likelihood with bias
# ---------------------------------------------------------------------------


class TestTotalLogLikelihoodWithBias:
    def test_bias_shifts_likelihood(self):
        """Bias correction changes the total log-likelihood."""
        markers = [
            _make_informative_marker((0, 0), (1, 1), 700, 300, marker_type=0),
        ]
        ll_no_bias = total_log_likelihood_bb(markers, 0.3, 0.01)
        biases = {("chr1", 100, "A", "T"): 0.02}
        ll_biased = total_log_likelihood_bb(
            markers, 0.3, 0.01, calibration=PanelCalibration(biases=biases)
        )
        assert ll_no_bias != ll_biased

    def test_correct_bias_improves_likelihood(self):
        """When data has a known bias, correcting for it should improve LL at the true f."""
        # Simulate a marker with bias +0.02: at f=0.3, true ALT VAF = 0.3
        # With bias: observed ALT VAF ≈ 0.32 -> ref=680, alt=320
        markers = [
            _make_informative_marker((0, 0), (1, 1), 680, 320, marker_type=0),
        ]
        # Without bias correction, the MLE will be pulled toward 0.32
        ll_uncorrected = total_log_likelihood_bb(markers, 0.30, 0.01)
        # With bias correction (0.02), the model expects VAF 0.32 at f=0.30
        biases = {("chr1", 100, "A", "T"): 0.02}
        ll_corrected = total_log_likelihood_bb(
            markers, 0.30, 0.01, calibration=PanelCalibration(biases=biases)
        )
        assert ll_corrected > ll_uncorrected


# ---------------------------------------------------------------------------
# estimate_single_donor with bias
# ---------------------------------------------------------------------------


class TestEstimateSingleDonorWithBias:
    def _make_biased_markers(self, true_f, biases_map, depth=2000):
        """Create markers where observed counts reflect a systematic bias."""
        rng = random.Random(42)
        markers = []
        # Create type 0 markers (host 0/0, donor 1/1)
        for pos, bias in biases_map.items():
            # True ALT VAF at this marker = true_f. Inject bias the same
            # (multiplicative, logit-space) way the estimator corrects for it.
            obs_vaf = float(inject_bias(true_f, bias))
            # Sample allele counts
            ad_alt = rng.binomialvariate(depth, obs_vaf)
            ad_ref = depth - ad_alt
            markers.append(
                _make_informative_marker(
                    (0, 0),
                    (1, 1),
                    ad_ref,
                    ad_alt,
                    marker_type=0,
                    chrom="chr1",
                    pos=pos,
                )
            )
        return markers

    def test_bias_correction_improves_accuracy(self):
        """With known biases, correction reduces estimation error."""
        true_f = 0.30
        # 20 markers with systematic biases
        biases_map = {i * 1000: 0.02 * ((-1) ** i) for i in range(20)}
        markers = self._make_biased_markers(true_f, biases_map)

        # Without bias correction
        result_no_bias = estimate_single_donor_bb(markers, error_rate=0.01)

        # With bias correction
        marker_biases = {("chr1", pos, "A", "T"): b for pos, b in biases_map.items()}
        result_biased = estimate_single_donor_bb(
            markers,
            error_rate=0.01,
            calibration=PanelCalibration(biases=marker_biases),
        )

        # Bias-corrected estimate should be closer to truth
        err_no_bias = abs(result_no_bias.donor_fraction - true_f)
        err_biased = abs(result_biased.donor_fraction - true_f)
        assert err_biased <= err_no_bias + 0.001  # allow tiny float tolerance

    def test_no_bias_dict_same_as_none(self):
        """Passing an empty bias dict should give same result as None."""
        markers = [
            _make_informative_marker((0, 0), (1, 1), 700, 300, marker_type=0),
        ]
        result_none = estimate_single_donor_bb(markers, error_rate=0.01)
        result_empty = estimate_single_donor_bb(
            markers,
            error_rate=0.01,
            calibration=PanelCalibration(biases={}),
        )
        assert result_none.donor_fraction == pytest.approx(
            result_empty.donor_fraction,
            abs=1e-6,
        )

    def test_bias_correction_all_same_direction(self):
        """When all biases push in the same direction, correction is essential."""
        true_f = 0.20
        # All markers biased +0.03 (ALT favoured)
        biases_map = {i * 1000: 0.03 for i in range(30)}
        markers = self._make_biased_markers(true_f, biases_map)

        result_no_bias = estimate_single_donor_bb(markers, error_rate=0.01)
        marker_biases = {("chr1", pos, "A", "T"): b for pos, b in biases_map.items()}
        result_corrected = estimate_single_donor_bb(
            markers,
            error_rate=0.01,
            calibration=PanelCalibration(biases=marker_biases),
        )

        # Without correction, estimate is biased high. The logit-space shift at
        # VAF 0.2 for bias +0.03 is ~0.0199, so the uncorrected pull is ~0.018.
        assert result_no_bias.donor_fraction > true_f + 0.01
        # With correction, estimate should be close to truth and clearly better.
        assert abs(result_corrected.donor_fraction - true_f) < 0.01
        assert abs(result_corrected.donor_fraction - true_f) < abs(
            result_no_bias.donor_fraction - true_f
        )

    def test_ci_coverage_improves_with_bias_correction(self):
        """Bias correction should help CIs cover the truth more often."""
        rng = random.Random(123)

        covers_no_bias = 0
        covers_corrected = 0
        n_trials = 20

        for _ in range(n_trials):
            true_f = 0.25
            # Random biases per marker
            biases_map = {i * 1000: rng.gauss(0, 0.02) for i in range(40)}
            markers = self._make_biased_markers(true_f, biases_map)

            result_no = estimate_single_donor_bb(markers, error_rate=0.01)
            if result_no.donor_fraction_ci[0] <= true_f <= result_no.donor_fraction_ci[1]:
                covers_no_bias += 1

            mb = {("chr1", p, "A", "T"): b for p, b in biases_map.items()}
            result_yes = estimate_single_donor_bb(
                markers,
                error_rate=0.01,
                calibration=PanelCalibration(biases=mb),
            )
            if result_yes.donor_fraction_ci[0] <= true_f <= result_yes.donor_fraction_ci[1]:
                covers_corrected += 1

        # Corrected should have better or equal CI coverage
        assert covers_corrected >= covers_no_bias
