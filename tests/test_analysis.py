"""Tests for allomix.analysis wiring: the contamination LoD floor and sex inference.

The contamination-floor-into-LoD logic is validated here against a synthetic
contamination scalar, which isolates the flooring rule from the in-data
contamination estimator (tested separately in ``test_sample_contamination.py``).
"""

import math

from allomix.analysis import _floor_detection_limits, analyse_sample
from allomix.genotype import ContigPolicy, MarkerData
from allomix.qc.sex import PairStatus, Sex
from allomix.results import ChimerismResult, MultiDonorResult


def _single_result(lob: float, lod: float) -> ChimerismResult:
    """A minimal single-donor result carrying only the LoB/LoD under test."""
    return ChimerismResult(
        donor_fraction=0.99,
        donor_fraction_ci=(0.98, 1.0),
        host_fraction=0.01,
        log_likelihood=0.0,
        n_informative=50,
        n_markers_used=50,
        per_marker=[],
        error_rate=1e-3,
        lob_fraction=lob,
        lod_fraction=lod,
    )


class TestFloorDetectionLimits:
    def test_floor_raises_both_when_above(self):
        """A floor above the analytical LoD replaces both limits."""
        res = _single_result(lob=0.0008, lod=0.0015)
        _floor_detection_limits(res, 0.002)
        assert res.lob_fraction == 0.002
        assert res.lod_fraction == 0.002

    def test_floor_below_lod_leaves_lod(self):
        """A floor between LoB and LoD raises only LoB; LoD stays and remains the
        larger of the two.
        """
        res = _single_result(lob=0.0008, lod=0.005)
        _floor_detection_limits(res, 0.002)
        assert res.lob_fraction == 0.002
        assert res.lod_fraction == 0.005

    def test_floor_below_both_is_noop(self):
        res = _single_result(lob=0.003, lod=0.006)
        _floor_detection_limits(res, 0.001)
        assert res.lob_fraction == 0.003
        assert res.lod_fraction == 0.006

    def test_zero_floor_is_noop(self):
        res = _single_result(lob=0.003, lod=0.006)
        _floor_detection_limits(res, 0.0)
        assert res.lob_fraction == 0.003
        assert res.lod_fraction == 0.006

    def test_infinite_lod_stays_infinite_above_floor(self):
        """An uninformative sample (inf LoD) stays inf; a finite floor cannot
        lower it.
        """
        res = _single_result(lob=float("inf"), lod=float("inf"))
        _floor_detection_limits(res, 0.002)
        assert math.isinf(res.lob_fraction)
        assert math.isinf(res.lod_fraction)

    def test_multidonor_result_unchanged(self):
        """Multi-donor results carry no LoB/LoD fields, so flooring is a no-op
        and must not raise.
        """
        res = MultiDonorResult(
            donor_fractions=[0.5, 0.49],
            donor_fraction_cis=[(0.45, 0.55), (0.44, 0.54)],
            host_fraction=0.01,
            log_likelihood=0.0,
            n_informative=50,
            n_markers_used=50,
            per_marker=[],
            error_rate=1e-3,
        )
        _floor_detection_limits(res, 0.002)  # must not raise
        assert not hasattr(res, "lob_fraction")


# ---------------------------------------------------------------------------
# Sex inference wiring (#50)
# ---------------------------------------------------------------------------


def _ref(chrom: str, pos: int, gt: tuple[int, int]) -> MarkerData:
    return MarkerData(chrom, pos, "A", "G", gt, 50, 50, 100, 99)


def _admix(chrom: str, pos: int, alt_frac: float, dp: int = 1000) -> MarkerData:
    alt = round(dp * alt_frac)
    return MarkerData(chrom, pos, "A", "G", (0, 1), dp - alt, alt, dp, None)


def _synthetic_trio(n_auto: int = 30, host_x_het: int = 12, donor_x_het: int = 0):
    """Autosomal hom-ref host / hom-alt donor markers plus 25 non-PAR chrX sites.

    The admixture is 10% donor on the autosomes. chrX host hets are at the
    first ``host_x_het`` sites, donor hets at the first ``donor_x_het``. Host
    and donor also carry ``n_auto`` heterozygous chr2 sites absent from the
    admixture (so they never enter classification): the female model's het
    probability comes from the sample's own autosomal het rate, which a real
    sample has at about 0.5 and an all-homozygous fixture would not.
    """
    host, donor, admix = [], [], []
    for i in range(n_auto):
        pos = 1000 * (i + 1)
        host.append(_ref("chr1", pos, (0, 0)))
        donor.append(_ref("chr1", pos, (1, 1)))
        admix.append(_admix("chr1", pos, 0.10))
        host.append(_ref("chr2", pos, (0, 1)))
        donor.append(_ref("chr2", pos, (0, 1)))
    for i in range(25):
        pos = 10_000_000 + 1000 * i
        host.append(_ref("chrX", pos, (0, 1) if i < host_x_het else (1, 1)))
        donor.append(_ref("chrX", pos, (0, 1) if i < donor_x_het else (1, 1)))
        admix.append(_admix("chrX", pos, 0.45))
    return host, donor, admix


class TestSexWiring:
    def test_sex_attached_and_chrx_excluded(self):
        host, donor, admix = _synthetic_trio()
        a = analyse_sample(host, [donor], admix, min_dp=0, min_gq=0, error_rate=0.01)
        sex = a.result.sex
        assert sex is not None
        assert sex.host.sex is Sex.FEMALE
        assert sex.donors[0].sex is Sex.MALE
        assert sex.pair is PairStatus.MISMATCHED
        # Phase 1: chrX markers never reach the estimate, even when informative.
        assert all(m.chrom == "chr1" for m in a.genotypes.informative)
        assert a.genotypes.n_informative_sex_chrom_excluded == 12
        assert a.genotypes.n_chrx_used == 0
        assert a.qc.sex is sex
        assert a.qc.n_informative_sex_chrom_excluded == 12
        assert a.qc.status != "FAIL"

    def test_declared_sex_resolves_autosome_only_panel(self):
        host, donor, admix = _synthetic_trio()
        host = [m for m in host if m.chrom == "chr1"]
        donor = [m for m in donor if m.chrom == "chr1"]
        admix = [m for m in admix if m.chrom == "chr1"]
        a = analyse_sample(
            host,
            [donor],
            admix,
            min_dp=0,
            min_gq=0,
            error_rate=0.01,
            declared_host_sex=Sex.FEMALE,
            declared_donor_sexes=[None],
        )
        sex = a.result.sex
        assert sex.host.sex is Sex.UNAVAILABLE
        assert sex.host.effective is Sex.FEMALE
        assert sex.host.source == "declared"
        assert sex.donors[0].source == "unavailable"
        assert sex.pair is PairStatus.UNKNOWN

    def test_conflict_fails_qc(self):
        host, donor, admix = _synthetic_trio()
        a = analyse_sample(
            host,
            [donor],
            admix,
            min_dp=0,
            min_gq=0,
            error_rate=0.01,
            declared_host_sex=Sex.MALE,
        )
        assert a.result.sex.host.source == "conflict"
        assert a.qc.status == "FAIL"
        assert any("Sex mismatch for host" in w for w in a.qc.warnings)

    def test_all_primary_policy_admits_chrx(self):
        host, donor, admix = _synthetic_trio()
        a = analyse_sample(
            host,
            [donor],
            admix,
            min_dp=0,
            min_gq=0,
            error_rate=0.01,
            contig_policy=ContigPolicy.ALL_PRIMARY,
        )
        assert a.genotypes.n_chrx_used == 12
        assert a.genotypes.n_informative_sex_chrom_excluded == 0
