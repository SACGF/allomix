"""Tests for allomix.error_rates — per-site empirical error rate estimation."""

from __future__ import annotations

import random
import tempfile
from pathlib import Path

import pytest

from allomix.chimerism import (
    estimate_single_donor_bb,
    log_likelihood_marker_bb,
    total_log_likelihood_bb,
)
from allomix.error_rates import (
    DEFAULT_ERROR_FLOOR,
    MarkerError,
    errors_to_simple_dict,
    estimate_error_rates,
    load_error_table,
    save_error_table,
)
from allomix.genotype import InformativeMarker, MarkerData

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _homref(pos: int, depth: int, n_alt: int, chrom: str = "chr1") -> MarkerData:
    """A hom-ref MarkerData with the given depth and observed ALT count."""
    return MarkerData(
        chrom=chrom, pos=pos, ref="A", alt="T",
        gt=(0, 0), ad_ref=depth - n_alt, ad_alt=n_alt, dp=depth,
    )


def _homalt(pos: int, depth: int, n_ref: int, chrom: str = "chr1") -> MarkerData:
    """A hom-alt MarkerData with the given depth and observed REF count."""
    return MarkerData(
        chrom=chrom, pos=pos, ref="A", alt="T",
        gt=(1, 1), ad_ref=n_ref, ad_alt=depth - n_ref, dp=depth,
    )


def _het(pos: int, depth: int, n_alt: int, chrom: str = "chr1") -> MarkerData:
    return MarkerData(
        chrom=chrom, pos=pos, ref="A", alt="T",
        gt=(0, 1), ad_ref=depth - n_alt, ad_alt=n_alt, dp=depth,
    )


def _informative(
    pos: int,
    host_gt: tuple[int, int],
    donor_gt: tuple[int, int],
    ad_ref: int,
    ad_alt: int,
    chrom: str = "chr1",
) -> InformativeMarker:
    return InformativeMarker(
        chrom=chrom, pos=pos, ref="A", alt="T",
        host_gt=host_gt, donor_gts=[donor_gt],
        marker_type=0,
        admix_ad_ref=ad_ref, admix_ad_alt=ad_alt,
        admix_dp=ad_ref + ad_alt,
    )


# ---------------------------------------------------------------------------
# estimate_error_rates
# ---------------------------------------------------------------------------


class TestEstimateErrorRates:
    def test_no_homozygous_calls_returns_empty(self) -> None:
        """Het-only training input has no observations of either direction."""
        samples = [[_het(100, 1000, 500)] for _ in range(20)]
        out = estimate_error_rates(samples)
        assert out == {}

    def test_clean_homref_recovers_input_rate(self) -> None:
        """Synthetic data drawn binomially at a known rate is recovered.

        Generates per-sample hom-ref observations with a fixed per-direction
        error rate, then checks the pooled estimate matches within a few
        binomial standard errors.
        """
        rng = random.Random(0)
        true_rate = 3.3e-3  # ~ 0.01/3
        depth = 2000
        n_samples = 50
        samples = []
        for _ in range(n_samples):
            n_alt = sum(1 for _ in range(depth) if rng.random() < true_rate)
            samples.append([_homref(100, depth, n_alt)])
        out = estimate_error_rates(samples, min_reads=1000)
        key = ("chr1", 100, "A", "T")
        assert key in out
        est = out[key].e_refalt
        assert est is not None
        # Total reads = 50 * 2000 = 1e5. SE = sqrt(p(1-p)/n) ~ 1.3e-4. Allow 5x.
        assert abs(est - true_rate) < 7e-4

    def test_max_vaf_homref_drops_miscalled_het(self) -> None:
        """A hom-ref observation with vaf > max_vaf_homref is excluded entirely."""
        # Clean baseline + one obviously-miscalled-het sample at the same site.
        clean = [[_homref(100, 2000, 5)] for _ in range(10)]  # 5/2000 ~ 0.0025
        miscalled_het = [[_homref(100, 2000, 800)]]  # vaf 0.4 -- a het, not hom-ref
        out = estimate_error_rates(clean + miscalled_het, min_reads=1000)
        key = ("chr1", 100, "A", "T")
        est = out[key].e_refalt
        # Without the filter the rate would be inflated by the 800 ALT-read
        # outlier: (50 + 800) / (11 * 2000) ~ 0.039. With the filter it stays
        # near the clean baseline 50 / (10 * 2000) = 0.0025.
        assert est is not None and est < 0.005

    def test_min_vaf_homalt_drops_miscalled_het(self) -> None:
        """Same as above for the hom-alt direction."""
        clean = [[_homalt(100, 2000, 5)] for _ in range(10)]
        miscalled_het = [[_homalt(100, 2000, 800)]]  # vaf 0.6 != hom-alt
        out = estimate_error_rates(clean + miscalled_het, min_reads=1000)
        key = ("chr1", 100, "A", "T")
        est = out[key].e_altref
        assert est is not None and est < 0.005

    def test_min_reads_filter_per_direction(self) -> None:
        """A site below min_reads in one direction returns None for that direction."""
        # Site A: enough hom-ref reads (>1000), no hom-alt
        # Site B: enough hom-alt reads (>1000), no hom-ref
        samples = [
            [_homref(100, 2000, 6), _homalt(200, 500, 1)],  # 500 < 1000 for site B
        ]
        out = estimate_error_rates(samples, min_reads=1000)
        a = out[("chr1", 100, "A", "T")]
        assert a.e_refalt is not None
        assert a.e_altref is None  # no hom-alt obs
        assert ("chr1", 200, "A", "T") not in out  # both directions below floor

    def test_no_usable_observations_omits_site(self) -> None:
        """A site with only-het observations is omitted from the output."""
        samples = [[_het(100, 1000, 500)] for _ in range(50)]
        out = estimate_error_rates(samples, min_reads=1)
        assert out == {}

    def test_pooling_weights_by_depth(self) -> None:
        """Two samples at the same site pool proportionally to depth.

        Sample with depth 100, 1 ALT (rate 0.01) and sample with depth 10000,
        10 ALT (rate 0.001) should pool to (1+10)/(100+10000) ~ 0.00109, much
        closer to the deep sample's rate than the shallow one's.
        """
        samples = [
            [_homref(100, 100, 1)],
            [_homref(100, 10000, 10)],
        ]
        out = estimate_error_rates(samples, min_reads=100)
        est = out[("chr1", 100, "A", "T")].e_refalt
        assert est is not None
        # Read-pooled MLE: 11 / 10100 ~ 0.001089
        assert abs(est - 11.0 / 10100.0) < 1e-9


# ---------------------------------------------------------------------------
# save / load round-trip
# ---------------------------------------------------------------------------


class TestSaveLoadRoundtrip:
    def test_roundtrip_with_na(self) -> None:
        """NA entries round-trip to None; floor is applied at load."""
        errors = {
            ("chr1", 100, "A", "T"): MarkerError(
                "chr1", 100, "A", "T",
                e_refalt=1e-3, e_altref=None,
                n_reads_homref=5000, n_reads_homalt=0,
            ),
            ("chr1", 200, "C", "G"): MarkerError(
                "chr1", 200, "C", "G",
                e_refalt=0.0, e_altref=2e-4,  # zero rate hits the floor
                n_reads_homref=5000, n_reads_homalt=5000,
            ),
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "errors.tsv"
            save_error_table(errors, path)
            loaded = load_error_table(path)

        assert loaded[("chr1", 100, "A", "T")] == (1e-3, None)
        # Floor applied to the 0.0 entry.
        e_ra, e_ar = loaded[("chr1", 200, "C", "G")]
        assert e_ra == DEFAULT_ERROR_FLOOR
        assert e_ar == 2e-4

    def test_loader_disable_floor(self) -> None:
        """Setting error_floor=0 disables the floor."""
        errors = {
            ("chr1", 100, "A", "T"): MarkerError(
                "chr1", 100, "A", "T",
                e_refalt=0.0, e_altref=0.0,
                n_reads_homref=5000, n_reads_homalt=5000,
            ),
        }
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "errors.tsv"
            save_error_table(errors, path)
            loaded = load_error_table(path, error_floor=0.0)
        assert loaded[("chr1", 100, "A", "T")] == (0.0, 0.0)


# ---------------------------------------------------------------------------
# errors_to_simple_dict
# ---------------------------------------------------------------------------


class TestSimpleDict:
    def test_floor_applied_in_helper(self) -> None:
        errors = {
            ("chr1", 100, "A", "T"): MarkerError(
                "chr1", 100, "A", "T",
                e_refalt=1e-7,  # below default floor
                e_altref=1e-3,
                n_reads_homref=5000, n_reads_homalt=5000,
            ),
        }
        simple = errors_to_simple_dict(errors)
        assert simple[("chr1", 100, "A", "T")] == (DEFAULT_ERROR_FLOOR, 1e-3)


# ---------------------------------------------------------------------------
# Likelihood integration
# ---------------------------------------------------------------------------


class TestLikelihoodIntegration:
    def test_symmetric_per_marker_matches_global(self) -> None:
        """Per-marker (e/3, e/3) should approximately match the legacy 4-state
        path. The renormalisation in the 4-state model adds ~O(e^2) so we use
        a loose tolerance.
        """
        m = _informative(100, host_gt=(0, 0), donor_gt=(1, 1),
                         ad_ref=900, ad_alt=100)
        ll_symmetric = log_likelihood_marker_bb(
            m.admix_ad_ref, m.admix_ad_alt, w=0.5, error_rate=0.01, rho=100.0,
        )
        ll_asymmetric = log_likelihood_marker_bb(
            m.admix_ad_ref, m.admix_ad_alt, w=0.5, error_rate=0.01, rho=100.0,
            e_refalt=0.01 / 3.0, e_altref=0.01 / 3.0,
        )
        # The two are not identical (the 4-state path renormalises onto
        # {REF, ALT}); they agree at the O(e^2) level.
        assert abs(ll_symmetric - ll_asymmetric) < 0.1

    def test_asymmetric_shifts_likelihood(self) -> None:
        """Switching to a much larger e_refalt at a hom-ref-like marker (w=1)
        increases the likelihood of seeing ALT reads.
        """
        # ad_alt = 10 at depth 1000 is unlikely under e_refalt = 1e-4 but
        # easy under e_refalt = 1e-2.
        ll_low = log_likelihood_marker_bb(
            ad_ref=990, ad_alt=10, w=1.0, rho=100.0,
            e_refalt=1e-4, e_altref=1e-4,
        )
        ll_high = log_likelihood_marker_bb(
            ad_ref=990, ad_alt=10, w=1.0, rho=100.0,
            e_refalt=1e-2, e_altref=1e-4,
        )
        assert ll_high > ll_low

    def test_marker_errors_missing_falls_back(self) -> None:
        """A marker absent from marker_errors uses the global 4-state path
        (i.e. matches the same call with marker_errors=None).
        """
        m = _informative(100, host_gt=(0, 1), donor_gt=(1, 1),
                         ad_ref=400, ad_alt=600)
        ll_with_empty = total_log_likelihood_bb(
            [m], f_donor=0.5, error_rate=0.01, rho=100.0, marker_errors={},
        )
        ll_no_table = total_log_likelihood_bb(
            [m], f_donor=0.5, error_rate=0.01, rho=100.0,
        )
        assert ll_with_empty == pytest.approx(ll_no_table)

    def test_marker_errors_one_direction_missing_falls_back(self) -> None:
        """If only one per-direction rate is known for a marker, the marker
        falls back to the symmetric path (cannot specify p_alt fully without
        both).
        """
        m = _informative(100, host_gt=(0, 1), donor_gt=(1, 1),
                         ad_ref=400, ad_alt=600)
        partial = {("chr1", 100, "A", "T"): (1e-3, None)}
        ll_partial = total_log_likelihood_bb(
            [m], f_donor=0.5, error_rate=0.01, rho=100.0,
            marker_errors=partial,
        )
        ll_baseline = total_log_likelihood_bb(
            [m], f_donor=0.5, error_rate=0.01, rho=100.0,
        )
        assert ll_partial == pytest.approx(ll_baseline)

    def test_estimator_default_unchanged(self) -> None:
        """estimate_single_donor_bb with default marker_errors=None matches
        the pre-Step-14 behaviour (regression guard).
        """
        rng = random.Random(7)
        # Build a small informative set at f_donor = 0.10.
        f = 0.10
        markers = []
        for i in range(40):
            host_gt = (0, 0)
            donor_gt = (1, 1)
            depth = 1000
            true_alt_prob = (1 - f) * 0 + f * 1.0  # w = 1 - f
            n_alt = sum(1 for _ in range(depth) if rng.random() < true_alt_prob)
            markers.append(
                _informative(i, host_gt, donor_gt, depth - n_alt, n_alt)
            )
        res_default = estimate_single_donor_bb(markers, error_rate=0.01)
        res_with_empty = estimate_single_donor_bb(
            markers, error_rate=0.01, marker_errors={},
        )
        # Empty error-table dict = no asymmetric markers = identical fit.
        assert res_default.donor_fraction == pytest.approx(
            res_with_empty.donor_fraction
        )
        assert res_default.log_likelihood == pytest.approx(
            res_with_empty.log_likelihood, abs=1e-8,
        )
