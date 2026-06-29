"""Tests for allomix.marker_contamination — Step 30 per-marker contamination."""

import math

from allomix.chimerism import estimate_single_donor_bb
from allomix.genotype import InformativeMarker, MarkerData
from allomix.likelihood import PanelCalibration
from allomix.marker_contamination import (
    ContaminationCorrection,
    apply_contamination_correction,
    estimate_carrier_counts,
    estimate_contamination_table,
    load_contamination_table,
    save_contamination_table,
)


def _im(
    pos: int,
    marker_type: int,
    ad_ref: int,
    ad_alt: int,
    host_gt: tuple[int, int],
    donor_gt: tuple[int, int],
    chrom: str = "chr1",
) -> InformativeMarker:
    """An informative marker with admix counts for correction tests."""
    return InformativeMarker(
        chrom=chrom,
        pos=pos,
        ref="A",
        alt="G",
        host_gt=host_gt,
        donor_gts=[donor_gt],
        marker_type=marker_type,
        admix_ad_ref=ad_ref,
        admix_ad_alt=ad_alt,
        admix_dp=ad_ref + ad_alt,
    )


def _md(
    pos: int,
    gt: tuple[int, int],
    ad_ref: int = 0,
    ad_alt: int = 0,
    chrom: str = "chr1",
) -> MarkerData:
    return MarkerData(
        chrom=chrom,
        pos=pos,
        ref="A",
        alt="G",
        gt=gt,
        ad_ref=ad_ref,
        ad_alt=ad_alt,
        dp=ad_ref + ad_alt,
    )


class TestCarrierDose:
    """The dose excludes the host's and donor's own carriage and caps."""

    def test_type0_excludes_host(self):
        # type 0: host 0/0 (carries REF), donor 1/1 (no REF). 4 cohort REF
        # carriers -> co-pooled = 4 - 1 (host) = 3.
        m = _im(100, 0, 1000, 0, (0, 0), (1, 1))
        corr = ContaminationCorrection(carriers={("chr1", 100, "A", "G"): (4, 2)})
        assert corr.carrier_dose(m, host_allele=0) == 3

    def test_type1_excludes_host(self):
        # type 1: host 1/1 (carries ALT), donor 0/0 (no ALT). 5 ALT carriers ->
        # co-pooled = 5 - 1 = 4.
        m = _im(100, 1, 0, 1000, (1, 1), (0, 0))
        corr = ContaminationCorrection(carriers={("chr1", 100, "A", "G"): (3, 5)})
        assert corr.carrier_dose(m, host_allele=1) == 4

    def test_cap(self):
        m = _im(100, 0, 1000, 0, (0, 0), (1, 1))
        corr = ContaminationCorrection(carriers={("chr1", 100, "A", "G"): (20, 0)}, dose_cap=5)
        assert corr.carrier_dose(m, host_allele=0) == 5

    def test_missing_marker_is_zero(self):
        m = _im(100, 0, 1000, 0, (0, 0), (1, 1))
        corr = ContaminationCorrection(carriers={})
        assert corr.carrier_dose(m, host_allele=0) == 0


class TestApplyNoOp:
    """The correction is a no-op (same list object) unless gated and positive."""

    def _markers(self):
        return [_im(100, 0, 1000, 5, (0, 0), (1, 1))]

    def test_none_correction(self):
        mk = self._markers()
        assert apply_contamination_correction(mk, None) is mk

    def test_not_gated(self):
        mk = self._markers()
        corr = ContaminationCorrection(
            carriers={("chr1", 100, "A", "G"): (4, 0)}, slope=0.001, gated=False
        )
        assert apply_contamination_correction(mk, corr) is mk

    def test_zero_slope(self):
        mk = self._markers()
        corr = ContaminationCorrection(
            carriers={("chr1", 100, "A", "G"): (4, 0)}, slope=0.0, gated=True
        )
        assert apply_contamination_correction(mk, corr) is mk


class TestCorrectionMath:
    """Subtracts slope * dose * depth from the host allele and rebuilds depth."""

    def test_type0_subtracts_from_ref(self):
        # dose = 4 carriers - 1 host = 3; slope 0.001; depth 2000 -> subtract 6.
        m = _im(100, 0, ad_ref=1990, ad_alt=10, host_gt=(0, 0), donor_gt=(1, 1))
        corr = ContaminationCorrection(
            carriers={("chr1", 100, "A", "G"): (4, 0)}, slope=0.001, gated=True
        )
        out = apply_contamination_correction([m], corr)[0]
        assert out.admix_ad_ref == 1990 - round(0.001 * 3 * 2000)  # 1984
        assert out.admix_ad_alt == 10  # untouched
        assert out.admix_dp == out.admix_ad_ref + out.admix_ad_alt

    def test_type1_subtracts_from_alt(self):
        m = _im(100, 1, ad_ref=10, ad_alt=1990, host_gt=(1, 1), donor_gt=(0, 0))
        corr = ContaminationCorrection(
            carriers={("chr1", 100, "A", "G"): (0, 4)}, slope=0.001, gated=True
        )
        out = apply_contamination_correction([m], corr)[0]
        assert out.admix_ad_alt == 1990 - round(0.001 * 3 * 2000)
        assert out.admix_ad_ref == 10
        assert out.admix_dp == out.admix_ad_ref + out.admix_ad_alt

    def test_donor_het_untouched(self):
        # marker_type 10 (donor het) is not a donor-hom class -> passes through.
        m = _im(100, 10, ad_ref=500, ad_alt=500, host_gt=(0, 1), donor_gt=(0, 0))
        corr = ContaminationCorrection(
            carriers={("chr1", 100, "A", "G"): (4, 4)}, slope=0.01, gated=True
        )
        assert apply_contamination_correction([m], corr)[0] is m

    def test_never_below_zero(self):
        m = _im(100, 0, ad_ref=5, ad_alt=1995, host_gt=(0, 0), donor_gt=(1, 1))
        corr = ContaminationCorrection(
            carriers={("chr1", 100, "A", "G"): (6, 0)}, slope=0.5, gated=True
        )
        out = apply_contamination_correction([m], corr)[0]
        assert out.admix_ad_ref == 0


class TestEstimatorByteIdentical:
    """Default estimation is unchanged: no correction == no contamination field."""

    def _markers(self):
        # A small informative set with a low host signal and a couple of carriers.
        mk = []
        for i in range(40):
            # donor-hom type 0, mostly donor (host fraction ~ 1%).
            mk.append(_im(100 + i, 0, ad_ref=20, ad_alt=1980, host_gt=(0, 0), donor_gt=(1, 1)))
        return mk

    def test_none_matches_default(self):
        mk = self._markers()
        base = estimate_single_donor_bb(mk, calibration=PanelCalibration())
        explicit_none = estimate_single_donor_bb(
            mk, calibration=PanelCalibration(contamination_correction=None)
        )
        assert base.host_fraction == explicit_none.host_fraction

    def test_gated_out_matches_default(self):
        mk = self._markers()
        base = estimate_single_donor_bb(mk, calibration=PanelCalibration())
        corr = ContaminationCorrection(
            carriers={(m.chrom, m.pos, m.ref, m.alt): (4, 0) for m in mk},
            slope=0.01,
            gated=False,  # clean flowcell: no-op even with a slope present
        )
        gated_out = estimate_single_donor_bb(
            mk, calibration=PanelCalibration(contamination_correction=corr)
        )
        assert base.host_fraction == gated_out.host_fraction


class TestEstimateCarrierCounts:
    def test_counts_ref_and_alt_carriers(self):
        # Three cohort individuals: 0/0, 0/1, 1/1 at one marker.
        cohort = [
            [_md(100, (0, 0))],
            [_md(100, (0, 1))],
            [_md(100, (1, 1))],
        ]
        counts = estimate_carrier_counts(cohort)
        # REF carried by 0/0 and 0/1 -> 2; ALT carried by 0/1 and 1/1 -> 2.
        assert counts[("chr1", 100, "A", "G")] == (2, 2)


def _consensus_admix(pos, gt, minor_frac, dp=2000):
    """Admix MarkerData at a consensus-hom site with a given minor fraction."""
    minor = round(minor_frac * dp)
    if gt == (0, 0):
        return _md(pos, gt, ad_ref=dp - minor, ad_alt=minor)
    return _md(pos, gt, ad_ref=minor, ad_alt=dp - minor)


class TestGate:
    """Clean flowcell (flat consensus) gates out; a dose response gates in."""

    def _cohort(self, n_carriers_per_marker):
        """Cohort where marker i has ``n`` carriers of the ALT (minor) allele."""
        cohort = []
        # Build n carriers as separate individuals carrying ALT (0/1), the rest 0/0.
        max_n = max(n_carriers_per_marker)
        for indiv in range(max_n):
            markers = []
            for i, n in enumerate(n_carriers_per_marker):
                gt = (0, 1) if indiv < n else (0, 0)
                markers.append(_md(100 + i, gt))
            cohort.append(markers)
        return cohort

    def test_clean_flat_gates_out(self):
        # 60 consensus-hom markers, minor fraction flat at the error floor
        # regardless of carrier count -> no dose response -> not gated.
        n_car = [(i % 5) for i in range(60)]
        cohort = self._cohort([max(1, n) for n in n_car])
        host = [_md(100 + i, (0, 0)) for i in range(60)]
        donor = [_md(100 + i, (0, 0)) for i in range(60)]
        admix = [_consensus_admix(100 + i, (0, 0), 0.002) for i in range(60)]
        corr = estimate_contamination_table(host, [donor], [admix], cohort, min_dp=1)
        assert not corr.gated
        assert corr.slope == 0.0

    def test_dose_response_gates_in(self):
        # Minor fraction rises with carrier count -> significant positive slope.
        n_car = [(i % 5) + 1 for i in range(120)]
        cohort = self._cohort(n_car)
        host = [_md(100 + i, (0, 0)) for i in range(120)]
        donor = [_md(100 + i, (0, 0)) for i in range(120)]
        admix = [_consensus_admix(100 + i, (0, 0), 0.001 + 0.002 * n_car[i]) for i in range(120)]
        corr = estimate_contamination_table(host, [donor], [admix], cohort, min_dp=1)
        assert corr.gated
        assert corr.gate_p_value < 0.05
        assert corr.gate_slope > 0.0


class TestSaveLoadRoundTrip:
    def test_round_trip(self, tmp_path):
        corr = ContaminationCorrection(
            carriers={
                ("chr1", 100, "A", "G"): (4, 2),
                ("chr2", 250, "C", "T"): (0, 7),
            },
            slope=0.00123,
            dose_cap=5,
            gated=True,
            gate_slope=0.0009,
            gate_p_value=1.2e-8,
            n_consensus=300,
            n_informative=120,
        )
        path = tmp_path / "contam.tsv"
        save_contamination_table(corr, path)
        loaded = load_contamination_table(path)
        assert loaded.carriers == corr.carriers
        assert math.isclose(loaded.slope, corr.slope, rel_tol=1e-6)
        assert loaded.dose_cap == corr.dose_cap
        assert loaded.gated is True
        assert math.isclose(loaded.gate_slope, corr.gate_slope, rel_tol=1e-6)
        assert loaded.n_consensus == 300
        assert loaded.n_informative == 120

    def test_gated_out_round_trip(self, tmp_path):
        corr = ContaminationCorrection(carriers={("chr1", 1, "A", "G"): (1, 1)})
        path = tmp_path / "contam.tsv"
        save_contamination_table(corr, path)
        loaded = load_contamination_table(path)
        assert loaded.gated is False
        assert loaded.slope == 0.0
