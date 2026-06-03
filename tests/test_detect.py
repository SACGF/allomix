"""Tests for allomix.detect — host-presence detection at donor-hom markers.

The detector is calibrated by construction when the simulator's symmetric
error rate ``e`` is mirrored by ``e_i = e/3`` on the detector side; these
tests exploit that by building ``InformativeMarker`` lists directly with
binomially-sampled donor-absent counts. The realistic-bias / overdispersion
path is exercised by the prototype in ``paper/scripts/run_presence_lod.py``.
"""

import math
import random

import numpy as np
import pytest

from allomix.chimerism import estimate_single_donor_bb
from allomix.detect import (
    HostPresenceResult,
    host_presence_test,
    select_donor_hom_markers,
)
from allomix.genotype import InformativeMarker, marker_type

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _imarker(
    host_gt: tuple[int, int],
    donor_gts: list[tuple[int, int]],
    ad_ref: int,
    ad_alt: int,
    chrom: str = "chr1",
    pos: int = 100,
    ref: str = "A",
    alt: str = "G",
) -> InformativeMarker:
    """Build an InformativeMarker with explicit allele counts.

    Constructs directly (no simulator) so the per-marker bias / depth
    variability noise sources documented in the plan's prototype-results
    section don't contaminate the calibration check.
    """
    mtypes = [marker_type(host_gt, d) for d in donor_gts]
    mt = mtypes[0] if mtypes and mtypes[0] is not None else (
        next((m for m in mtypes if m is not None), 0)
    )
    return InformativeMarker(
        chrom=chrom,
        pos=pos,
        ref=ref,
        alt=alt,
        host_gt=host_gt,
        donor_gts=list(donor_gts),
        marker_type=mt,
        admix_ad_ref=ad_ref,
        admix_ad_alt=ad_alt,
        admix_dp=ad_ref + ad_alt,
        marker_types=mtypes,
        informative_for=[m is not None for m in mtypes],
    )


def _pure_donor_panel(
    n_per_type: int,
    depth: int,
    e: float,
    rng: random.Random,
) -> list[InformativeMarker]:
    """A donor-homozygous panel under H0 (no host).

    For each Vynck type 0/1/10/11 we generate ``n_per_type`` markers and draw
    the donor-absent count from Binomial(depth, e/3). Other reads are
    assigned to the donor-present allele (or one of two alleles for host-het
    types — irrelevant to the detector).
    """
    out: list[InformativeMarker] = []
    pos = 100
    for _ in range(n_per_type):
        # Type 0: host 0/0, donor 1/1 -> donor-absent = REF
        y = sum(1 for _ in range(depth) if rng.random() < e / 3.0)
        out.append(_imarker((0, 0), [(1, 1)], ad_ref=y, ad_alt=depth - y, pos=pos))
        pos += 100
        # Type 1: host 1/1, donor 0/0 -> donor-absent = ALT
        y = sum(1 for _ in range(depth) if rng.random() < e / 3.0)
        out.append(_imarker((1, 1), [(0, 0)], ad_ref=depth - y, ad_alt=y, pos=pos))
        pos += 100
        # Type 10: host 0/1, donor 0/0 -> donor-absent = ALT
        y = sum(1 for _ in range(depth) if rng.random() < e / 3.0)
        out.append(_imarker((0, 1), [(0, 0)], ad_ref=depth - y, ad_alt=y, pos=pos))
        pos += 100
        # Type 11: host 0/1, donor 1/1 -> donor-absent = REF
        y = sum(1 for _ in range(depth) if rng.random() < e / 3.0)
        out.append(_imarker((0, 1), [(1, 1)], ad_ref=y, ad_alt=depth - y, pos=pos))
        pos += 100
    return out


def _spiked_panel(
    n_per_type: int,
    depth: int,
    e: float,
    f_h: float,
    rng: random.Random,
) -> list[InformativeMarker]:
    """Donor-homozygous panel with a known host fraction ``f_h`` injected.

    Donor-absent reads are drawn from ``Binomial(depth, q)`` with
    ``q = e/3 + (h/2) * f_h`` so the LRT MLE should track ``f_h`` and the
    pooled-Poisson detector should fire well below alpha at clinically
    plausible depths.
    """
    out: list[InformativeMarker] = []
    pos = 100
    for _ in range(n_per_type):
        for host_gt, donor_gt, h in (
            ((0, 0), (1, 1), 2),
            ((1, 1), (0, 0), 2),
            ((0, 1), (0, 0), 1),
            ((0, 1), (1, 1), 1),
        ):
            q = e / 3.0 + (h / 2.0) * f_h
            y = sum(1 for _ in range(depth) if rng.random() < q)
            donor_absent_is_alt = donor_gt == (0, 0)
            if donor_absent_is_alt:
                out.append(
                    _imarker(host_gt, [donor_gt], ad_ref=depth - y, ad_alt=y, pos=pos)
                )
            else:
                out.append(
                    _imarker(host_gt, [donor_gt], ad_ref=y, ad_alt=depth - y, pos=pos)
                )
            pos += 100
    return out


# ---------------------------------------------------------------------------
# Marker selection
# ---------------------------------------------------------------------------


class TestSelectDonorHomMarkers:
    def test_includes_types_0_1_10_11(self):
        markers = [
            _imarker((0, 0), [(1, 1)], 5, 995, pos=100),  # type 0
            _imarker((1, 1), [(0, 0)], 995, 5, pos=200),  # type 1
            _imarker((0, 1), [(0, 0)], 990, 10, pos=300),  # type 10
            _imarker((0, 1), [(1, 1)], 10, 990, pos=400),  # type 11
        ]
        rows = select_donor_hom_markers(markers)
        assert len(rows) == 4
        # Host doses: type 0/1 -> 2, types 10/11 -> 1.
        h_by_pos = {r.key[1]: r.h for r in rows}
        assert h_by_pos == {100: 2, 200: 2, 300: 1, 400: 1}

    def test_excludes_donor_het_types(self):
        # Type 20: host 0/0, donor 0/1 -> donor het, no donor-absent allele.
        # Type 21: host 1/1, donor 0/1.
        markers = [
            _imarker((0, 0), [(0, 1)], 990, 10, pos=100),  # type 20
            _imarker((1, 1), [(0, 1)], 10, 990, pos=200),  # type 21
        ]
        assert select_donor_hom_markers(markers) == []

    def test_multi_donor_requires_absence_from_all(self):
        # First donor 0/0 (would be type 1 with host 1/1), but second donor
        # is 0/1 and carries the ALT allele — disqualifies the marker.
        bad = _imarker(
            host_gt=(1, 1),
            donor_gts=[(0, 0), (0, 1)],
            ad_ref=990,
            ad_alt=10,
        )
        # Both donors 0/0: usable.
        good = _imarker(
            host_gt=(1, 1),
            donor_gts=[(0, 0), (0, 0)],
            ad_ref=990,
            ad_alt=10,
            pos=200,
        )
        rows = select_donor_hom_markers([bad, good])
        assert [r.key[1] for r in rows] == [200]

    def test_multi_donor_requires_same_homozygous_allele(self):
        # Donor 1 is 0/0, donor 2 is 1/1: between them they carry both
        # alleles, so no donor-absent allele exists.
        marker = _imarker(
            host_gt=(0, 1),
            donor_gts=[(0, 0), (1, 1)],
            ad_ref=500,
            ad_alt=500,
        )
        assert select_donor_hom_markers([marker]) == []

    def test_direction_matches_donor_homozygote(self):
        # Donor hom-ref -> direction ref->alt; donor hom-alt -> alt->ref.
        markers = [
            _imarker((1, 1), [(0, 0)], 995, 5, pos=100),
            _imarker((0, 0), [(1, 1)], 5, 995, pos=200),
        ]
        rows = select_donor_hom_markers(markers)
        by_pos = {r.key[1]: r.direction for r in rows}
        assert by_pos == {100: "ref->alt", 200: "alt->ref"}


# ---------------------------------------------------------------------------
# Calibration (false-positive rate on pure donor)
# ---------------------------------------------------------------------------


class TestFalsePositiveRate:
    def test_pure_donor_returns_non_significant(self):
        """One pure-donor panel should not produce a small p-value."""
        rng = random.Random(0)
        markers = _pure_donor_panel(n_per_type=20, depth=2000, e=1e-3, rng=rng)
        res = host_presence_test(markers, error_rate=1e-3)
        assert res.n_markers == 80
        assert res.lrt_pval > 0.05
        assert res.poisson_pval > 0.05

    def test_calibrated_across_replicates(self):
        """FP rate at alpha=0.05 should be in the Wald 95% band for n=200.

        Mirrors the prototype's gate-1 check (LRT FP rate ~ 0.05 across cells)
        but at a single cell with enough replicates to keep the test fast and
        the band wide enough to absorb sampling noise.
        """
        rng = random.Random(1234)
        n_reps = 200
        n_fp = 0
        for _ in range(n_reps):
            markers = _pure_donor_panel(n_per_type=20, depth=1000, e=1e-3, rng=rng)
            res = host_presence_test(markers, error_rate=1e-3)
            if res.lrt_pval < 0.05:
                n_fp += 1
        fp_rate = n_fp / n_reps
        # Wald band for p=0.05, n=200: roughly [0.020, 0.080]. Use a slightly
        # wider band to keep the test resilient to a tail draw.
        assert 0.015 <= fp_rate <= 0.10, f"FP rate {fp_rate:.3f} out of band"


# ---------------------------------------------------------------------------
# Power (spiked positives)
# ---------------------------------------------------------------------------


class TestPower:
    def test_spiked_positive_is_significant(self):
        rng = random.Random(42)
        # 80 markers, 2000x, f_h=1e-3 is well above the LoD this cell can
        # support; expect a tiny p-value.
        markers = _spiked_panel(
            n_per_type=20, depth=2000, e=3e-4, f_h=1e-3, rng=rng,
        )
        res = host_presence_test(markers, error_rate=3e-4)
        assert res.lrt_pval < 1e-3
        assert res.poisson_pval < 1e-3
        # MLE should land in roughly the right neighbourhood.
        assert res.f_host_mle == pytest.approx(1e-3, abs=5e-4)
        # CI should contain the truth.
        lo, hi = res.f_host_ci
        assert lo <= 1e-3 <= hi

    def test_mle_tracks_one_minus_donor_fraction(self):
        """LRT MLE should agree with 1 - donor_fraction from the BB estimator.

        Cross-check between the dedicated detector and the global MLE on the
        same admixture sample — when the host is present at a modest fraction
        both should land near the truth, with the detector typically tighter
        because it ignores noisy het markers.
        """
        rng = random.Random(7)
        markers = _spiked_panel(
            n_per_type=15, depth=2000, e=1e-3, f_h=5e-3, rng=rng,
        )
        res = host_presence_test(markers, error_rate=1e-3)
        mle = estimate_single_donor_bb(markers, error_rate=1e-3, grid_steps=201)
        mle_host = 1.0 - mle.donor_fraction
        # Both should land near the spiked 5e-3; tolerate sampling jitter.
        assert res.f_host_mle == pytest.approx(5e-3, abs=2e-3)
        assert mle_host == pytest.approx(5e-3, abs=3e-3)
        # The two should agree with each other to within sampling noise.
        assert res.f_host_mle == pytest.approx(mle_host, abs=3e-3)


# ---------------------------------------------------------------------------
# Per-site error table effect
# ---------------------------------------------------------------------------


class TestErrorTable:
    def test_per_site_rate_shifts_pvalue_in_expected_direction(self):
        """Inflating the per-site rate makes a given count look more like
        background and pushes the p-value up. Deflating it does the opposite.

        We hold the data fixed and only change the error table the detector
        sees; the comparison isolates the per-site lookup.
        """
        rng = random.Random(99)
        markers = _spiked_panel(
            n_per_type=10, depth=1000, e=1e-3, f_h=5e-4, rng=rng,
        )
        baseline = host_presence_test(markers, error_rate=1e-3)

        # Build a table that tells the detector the background is 10x higher
        # everywhere (per direction). The pooled count Y now sits closer to
        # Lam, so the p-value rises.
        # Same set of keys covered in both directions, defensively.
        big = 1e-2
        table_high = {
            (m.chrom, m.pos, m.ref, m.alt): (big, big) for m in markers
        }
        inflated = host_presence_test(markers, marker_errors=table_high, error_rate=1e-3)
        assert inflated.lrt_pval >= baseline.lrt_pval
        assert inflated.used_per_site_error
        assert inflated.error_rate_source == "per-site"

        # And a table that says background is 10x lower drives the p-value
        # downward (the same count looks more anomalous).
        small = 1e-4
        table_low = {
            (m.chrom, m.pos, m.ref, m.alt): (small, small) for m in markers
        }
        deflated = host_presence_test(markers, marker_errors=table_low, error_rate=1e-3)
        assert deflated.lrt_pval <= baseline.lrt_pval

    def test_missing_sites_fall_back_to_global(self):
        """Sites missing from the table use error_rate / 3 and the source
        flag should reflect that. Drop half the sites from the table; expect
        ``mixed``.
        """
        rng = random.Random(123)
        markers = _pure_donor_panel(n_per_type=10, depth=1000, e=1e-3, rng=rng)
        # Cover only half of the donor-hom markers in the table.
        rows = select_donor_hom_markers(markers)
        half = rows[: len(rows) // 2]
        table = {r.key: (1e-3, 1e-3) for r in half}
        res = host_presence_test(markers, marker_errors=table, error_rate=1e-3)
        assert res.used_per_site_error
        assert res.error_rate_source == "mixed"

    def test_no_table_marks_global_fallback(self):
        rng = random.Random(55)
        markers = _pure_donor_panel(n_per_type=5, depth=1000, e=1e-3, rng=rng)
        res = host_presence_test(markers, error_rate=1e-3)
        assert not res.used_per_site_error
        assert res.error_rate_source == "global-fallback"

    def test_per_direction_missing_falls_back(self):
        """A site present in the table but with ``None`` for the relevant
        direction falls back to global (and counts as fallback for the
        source flag).
        """
        rng = random.Random(8)
        markers = _pure_donor_panel(n_per_type=5, depth=1000, e=1e-3, rng=rng)
        rows = select_donor_hom_markers(markers)
        # Only fill the *wrong* direction per row; the detector should still
        # fall back to the global rate.
        table: dict[tuple[str, int, str, str], tuple[float | None, float | None]] = {}
        for r in rows:
            if r.direction == "ref->alt":
                table[r.key] = (None, 1e-3)  # only alt->ref filled
            else:
                table[r.key] = (1e-3, None)
        res = host_presence_test(markers, marker_errors=table, error_rate=1e-3)
        assert res.error_rate_source == "global-fallback"
        assert not res.used_per_site_error


# ---------------------------------------------------------------------------
# Degenerate inputs
# ---------------------------------------------------------------------------


class TestDegenerate:
    def test_empty_input_returns_neutral_result(self):
        res = host_presence_test([], error_rate=1e-3)
        assert isinstance(res, HostPresenceResult)
        assert res.n_markers == 0
        assert res.poisson_pval == 1.0
        assert res.lrt_pval == 1.0
        assert res.f_host_mle == 0.0
        assert res.f_host_ci == (0.0, 0.0)
        assert res.error_rate_source == "none"

    def test_all_het_donor_markers_excluded(self):
        markers = [
            _imarker((0, 0), [(0, 1)], 990, 10),  # type 20
            _imarker((1, 1), [(0, 1)], 10, 990),  # type 21
        ]
        res = host_presence_test(markers, error_rate=1e-3)
        assert res.n_markers == 0
        assert res.error_rate_source == "none"

    def test_zero_donor_absent_reads_gives_pval_one(self):
        """No stray reads at all -> Y=0 -> Poisson p=1, LRT collapses to 0."""
        markers = [
            _imarker((0, 0), [(1, 1)], 0, 1000, pos=100),
            _imarker((1, 1), [(0, 0)], 1000, 0, pos=200),
        ]
        res = host_presence_test(markers, error_rate=1e-3)
        assert res.poisson_pval == 1.0
        assert res.lrt_pval == 1.0
        assert res.f_host_mle == 0.0


# ---------------------------------------------------------------------------
# Sanity: numpy / math availability under the dataclass machinery
# ---------------------------------------------------------------------------


def test_module_smoke():
    """Importing detect and round-tripping a tiny call shouldn't blow up."""
    assert math.isfinite(np.float64(1.0))
    res = host_presence_test([], error_rate=1e-3)
    assert res.n_markers == 0
