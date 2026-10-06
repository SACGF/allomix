"""Tests for allomix.analysis wiring: the contamination LoD floor and sex inference.

The contamination-floor-into-LoD logic is validated here against a synthetic
contamination scalar, which isolates the flooring rule from the in-data
contamination estimator (tested separately in ``test_sample_contamination.py``).
"""

import math

import pytest

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
        # Mismatched pair: chrX markers never reach the estimate, even when informative.
        assert all(m.chrom == "chr1" for m in a.genotypes.informative)
        assert a.genotypes.n_informative_sex_chrom_excluded == 12
        assert a.genotypes.n_chrx_used == 0
        assert a.qc.sex is sex
        assert a.qc.n_informative_sex_chrom_excluded == 12
        assert a.qc.status != "FAIL"

    def test_host_chry_depth_resolves_autosome_only_panel(self):
        """``host_region_depth`` reaches sex inference as the chrY secondary signal."""
        host, donor, admix = _synthetic_trio()
        host = [m for m in host if m.chrom == "chr1"]
        donor = [m for m in donor if m.chrom == "chr1"]
        admix = [m for m in admix if m.chrom == "chr1"]
        depth = [("chr1", 1_000_000 + i, 1000) for i in range(10)] + [("chrY", 2_787_394, 900)]
        a = analyse_sample(
            host, [donor], admix, min_dp=0, min_gq=0, error_rate=0.01, host_region_depth=depth
        )
        sex = a.result.sex
        assert sex.host.sex is Sex.MALE
        assert sex.host.chry_resolved
        assert sex.donors[0].sex is Sex.UNAVAILABLE
        assert sex.pair is PairStatus.UNKNOWN

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


def _sex_matched_trio(host_sex: str, donor_sex: str, f: float = 0.10):
    """30 autosomal type-0 markers plus 25 non-PAR chrX markers for a sex-matched pair.

    chrX genotypes (sites i = 0..24): FF host het at i < 12, donor het at
    6 <= i < 18, hom-alt elsewhere, so 12 chrX sites are informative. MM host
    hom-alt at i < 12 then hom-ref, donor hom-ref at i < 20 then hom-alt, with
    one spurious host het at i = 24 (the hemizygous-het genotyping error), so
    16 hom/hom contrasts are informative and one site is dropped. Each party
    also carries het chr2 sites absent from the admixture, for the autosomal
    het rate the female model uses. Admixture VAFs follow the diploid dosage
    model at ``f`` donor, which is exact for a sex-matched pair.
    """
    host, donor, admix = [], [], []
    for i in range(30):
        pos = 1000 * (i + 1)
        host.append(_ref("chr1", pos, (0, 0)))
        donor.append(_ref("chr1", pos, (1, 1)))
        admix.append(_admix("chr1", pos, f))
        host.append(_ref("chr2", pos, (0, 1)))
        donor.append(_ref("chr2", pos, (0, 1)))
    for i in range(25):
        pos = 10_000_000 + 1000 * i
        if host_sex == "F":
            h = (0, 1) if i < 12 else (1, 1)
            d = (0, 1) if 6 <= i < 18 else (1, 1)
        else:
            h = (1, 1) if i < 12 else (0, 0)
            d = (0, 0) if i < 20 else (1, 1)
            if i == 24:
                h = (0, 1)
        host.append(_ref("chrX", pos, h))
        donor.append(_ref("chrX", pos, d))
        admix.append(_admix("chrX", pos, ((1 - f) * sum(h) + f * sum(d)) / 2))
    return host, donor, admix


class TestSexAwareRoutingWiring:
    """assess_sex runs first and its pair status drives classify_markers (#46)."""

    def test_matched_female_uses_chrx(self):
        host, donor, admix = _sex_matched_trio("F", "F")
        a = analyse_sample(host, [donor], admix, min_dp=0, min_gq=0, error_rate=0.01)
        assert a.result.sex.pair is PairStatus.MATCHED_FEMALE
        assert a.genotypes.n_chrx_used == 12
        assert a.genotypes.n_informative_sex_chrom_excluded == 0
        assert a.genotypes.n_chrx_male_het_dropped == 0
        assert sum(m.chrom == "chrX" for m in a.genotypes.informative) == 12
        assert a.qc.n_chrx_used == 12
        assert a.result.donor_fraction == pytest.approx(0.10, abs=0.01)

    def test_matched_male_hom_hom_used_het_dropped(self):
        host, donor, admix = _sex_matched_trio("M", "M")
        a = analyse_sample(host, [donor], admix, min_dp=0, min_gq=0, error_rate=0.01)
        assert a.result.sex.host.sex is Sex.MALE
        assert a.result.sex.donors[0].sex is Sex.MALE
        assert a.result.sex.pair is PairStatus.MATCHED_MALE
        assert a.genotypes.n_chrx_used == 16
        assert a.genotypes.n_chrx_male_het_dropped == 1
        assert a.qc.n_chrx_male_het_dropped == 1
        chrx = [m for m in a.genotypes.informative if m.chrom == "chrX"]
        assert all(m.host_gt[0] == m.host_gt[1] for m in chrx)
        assert a.result.donor_fraction == pytest.approx(0.10, abs=0.01)

    def test_autosomes_only_policy_overrides_matched_pair(self):
        host, donor, admix = _sex_matched_trio("F", "F")
        a = analyse_sample(
            host,
            [donor],
            admix,
            min_dp=0,
            min_gq=0,
            error_rate=0.01,
            contig_policy=ContigPolicy.AUTOSOMES_ONLY,
        )
        assert a.result.sex.pair is PairStatus.MATCHED_FEMALE
        assert a.genotypes.n_chrx_used == 0
        assert a.genotypes.n_informative_sex_chrom_excluded == 12

    def test_declared_sex_resolves_pair_for_routing(self):
        """A declared donor sex turns an unknown pair into a matched one."""
        host, donor, admix = _sex_matched_trio("F", "F")
        donor_auto = [m for m in donor if m.chrom != "chrX"]  # donor chrX ungenotyped
        a = analyse_sample(host, [donor_auto], admix, min_dp=0, min_gq=0, error_rate=0.01)
        assert a.result.sex.pair is PairStatus.UNKNOWN
        assert a.genotypes.n_chrx_used == 0
        b = analyse_sample(
            host,
            [donor_auto],
            admix,
            min_dp=0,
            min_gq=0,
            error_rate=0.01,
            declared_donor_sexes=[Sex.FEMALE],
        )
        assert b.result.sex.pair is PairStatus.MATCHED_FEMALE
        # No shared chrX markers (donor has none), so nothing is routed, but the
        # pair status now permits it.
        assert b.genotypes.n_chrx_used == 0
        assert b.genotypes.n_informative_sex_chrom_excluded == 0
