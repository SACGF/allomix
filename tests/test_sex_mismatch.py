"""Tests for allomix.qc.sex_mismatch: the sex-chromosome cross-check (#48)."""

import random
from pathlib import Path

import pytest

from allomix.analysis import analyse_sample
from allomix.constants import SEXCHROM_CONCORDANCE_PP, SEXCHROM_MIN_MARKERS
from allomix.estimate.likelihood import expected_weight
from allomix.genotype import InformativeMarker, MarkerData, parse_vcf
from allomix.qc.sex import Sex, SexInference, SexResult, assess_sex, pair_status
from allomix.qc.sex_mismatch import (
    CHRX_BASIS,
    CHRY_BASIS,
    SexMismatchResult,
    chry_depth_ratio,
    chry_male_fraction,
    expected_ref_fraction,
    fit_chrx_fraction,
    is_concordant,
    ref_fraction,
    region_depths_from_vcf,
    sex_mismatch_check,
    usable_chrx_markers,
)
from allomix.simulate import build_joint_vcf_from_genotype_dicts, write_joint_vcf

SEXCHROM_DIR = Path(__file__).resolve().parent / "test_data" / "sexchrom"
MF_VCF = SEXCHROM_DIR / "joint_MF.vcf"

# ---------------------------------------------------------------------------
# Expected-fraction formula
# ---------------------------------------------------------------------------


class TestExpectedRefFraction:
    def test_male_minor_at_female_hom_site(self):
        """Male (cn 1) hom-alt donor at a female hom-ref host site: ALT = f / (2 - f)."""
        for f in (0.01, 0.1, 0.5):
            w = expected_ref_fraction(f, ref_fraction((0, 0)), ref_fraction((1, 1)), 2, 1)
            assert 1.0 - w == pytest.approx(f / (2.0 - f))

    def test_female_hom_alt_minor_at_male_hom_ref(self):
        """Female (cn 2) hom-alt donor at a male hom-ref host site: ALT = 2f / (1 + f)."""
        for f in (0.01, 0.1, 0.5):
            w = expected_ref_fraction(f, ref_fraction((0, 0)), ref_fraction((1, 1)), 1, 2)
            assert 1.0 - w == pytest.approx(2.0 * f / (1.0 + f))

    @pytest.mark.parametrize("cn", [1, 2])
    @pytest.mark.parametrize(
        "host_gt,donor_gt", [((0, 0), (1, 1)), ((0, 1), (1, 1)), ((0, 0), (0, 1)), ((1, 1), (0, 0))]
    )
    def test_equal_copy_number_reduces_to_diploid(self, cn, host_gt, donor_gt):
        """FF and MM pairs: the copy numbers cancel and the diploid weight comes back."""
        for f in (0.0, 0.05, 0.3, 1.0):
            w = expected_ref_fraction(f, ref_fraction(host_gt), ref_fraction(donor_gt), cn, cn)
            assert w == pytest.approx(expected_weight(host_gt, donor_gt, f))

    def test_endpoints_are_pure_parties(self):
        assert expected_ref_fraction(0.0, 1.0, 0.0, 2, 1) == 1.0
        assert expected_ref_fraction(1.0, 1.0, 0.0, 2, 1) == 0.0

    def test_ref_fraction_values(self):
        assert ref_fraction((0, 0)) == 1.0
        assert ref_fraction((0, 1)) == 0.5
        assert ref_fraction((1, 1)) == 0.0


# ---------------------------------------------------------------------------
# chrX estimator on synthetic counts
# ---------------------------------------------------------------------------


def _chrx_marker(
    i: int, host_gt: tuple[int, int], donor_gt: tuple[int, int], ad_ref: int, ad_alt: int
) -> InformativeMarker:
    return InformativeMarker(
        chrom="chrX",
        pos=10_000_000 + i * 1_000_000,
        ref="A",
        alt="G",
        host_gt=host_gt,
        donor_gts=[donor_gt],
        marker_type=1,
        admix_ad_ref=ad_ref,
        admix_ad_alt=ad_alt,
        admix_dp=ad_ref + ad_alt,
    )


def _synthetic_chrx(f_true: float, cn_host: int, cn_donor: int, n: int, seed: int, depth: int):
    """Binomial counts at the copy-number-weighted expectation (4-state error 1%)."""
    rng = random.Random(seed)
    gts = [((0, 0), (1, 1)), ((1, 1), (0, 0)), ((0, 0), (0, 1)), ((1, 1), (0, 1))]
    markers = []
    for i in range(n):
        hg, dg = gts[i % len(gts)]
        # A hemizygous party must be homozygous; swap a het onto the diploid side.
        if cn_host == 1 and hg[0] != hg[1]:
            hg, dg = dg, hg
        if cn_donor == 1 and dg[0] != dg[1]:
            hg, dg = dg, hg
        w = expected_ref_fraction(f_true, ref_fraction(hg), ref_fraction(dg), cn_host, cn_donor)
        p_alt = (1.0 - w) * 0.99 + w * (0.01 / 3)
        p_ref = w * 0.99 + (1.0 - w) * (0.01 / 3)
        p = p_alt / (p_alt + p_ref)
        k = sum(1 for _ in range(depth) if rng.random() < p)
        markers.append(_chrx_marker(i, hg, dg, depth - k, k))
    return markers


class TestFitChrxFraction:
    @pytest.mark.parametrize("cn_host,cn_donor", [(2, 1), (1, 2)])
    @pytest.mark.parametrize("f_true", [0.02, 0.10, 0.90])
    def test_recovers_known_fraction(self, cn_host, cn_donor, f_true):
        markers = _synthetic_chrx(f_true, cn_host, cn_donor, n=20, seed=3, depth=2000)
        f, (lo, hi), rho = fit_chrx_fraction(markers, cn_host, cn_donor)
        assert abs(f - f_true) < 0.01
        assert lo <= f <= hi
        assert hi - lo < 0.05
        assert rho > 0

    def test_diploid_model_would_be_biased(self):
        """Fitting the same male-minor counts with cn 2/2 misses the truth: the
        copy-number model is what makes the estimate correct."""
        markers = _synthetic_chrx(0.10, 2, 1, n=20, seed=5, depth=2000)
        f_cn, _, _ = fit_chrx_fraction(markers, 2, 1)
        f_diploid, _, _ = fit_chrx_fraction(markers, 2, 2)
        assert abs(f_cn - 0.10) < 0.01
        assert abs(f_diploid - 0.10) > 0.02

    def test_empty_raises(self):
        with pytest.raises(ValueError):
            fit_chrx_fraction([], 2, 1)


# ---------------------------------------------------------------------------
# Marker selection
# ---------------------------------------------------------------------------


def _md(chrom: str, pos: int, gt: tuple[int, int], ad: tuple[int, int] = (500, 500)) -> MarkerData:
    return MarkerData(
        chrom=chrom, pos=pos, ref="A", alt="G", gt=gt, ad_ref=ad[0], ad_alt=ad[1], dp=1000, gq=99
    )


class TestUsableChrxMarkers:
    def test_drops_male_het_and_par_and_autosomes(self):
        host = [
            _md("chrX", 20_000_000, (0, 0)),  # usable
            _md("chrX", 21_000_000, (0, 1)),  # male het -> dropped
            _md("chrX", 1_000_000, (0, 0)),  # PAR -> never considered
            _md("chr1", 5_000_000, (0, 0)),  # autosome -> not chrX
            _md("chrX", 22_000_000, (1, 1)),  # host == donor -> not informative
        ]
        donor = [
            _md("chrX", 20_000_000, (1, 1)),
            _md("chrX", 21_000_000, (1, 1)),
            _md("chrX", 1_000_000, (1, 1)),
            _md("chr1", 5_000_000, (1, 1)),
            _md("chrX", 22_000_000, (1, 1)),
        ]
        admix = [_md(m.chrom, m.pos, (0, 1)) for m in host]
        kept, dropped = usable_chrx_markers(host, donor, admix, "host", min_dp=0, min_gq=0)
        assert [m.pos for m in kept] == [20_000_000]
        assert dropped == 1

    def test_female_het_is_kept(self):
        host = [_md("chrX", 20_000_000, (0, 1))]
        donor = [_md("chrX", 20_000_000, (1, 1))]
        admix = [_md("chrX", 20_000_000, (0, 1))]
        kept, dropped = usable_chrx_markers(host, donor, admix, "donor", min_dp=0, min_gq=0)
        assert len(kept) == 1 and dropped == 0


# ---------------------------------------------------------------------------
# chrY depth ratio (toy tables)
# ---------------------------------------------------------------------------


def _depth_table(y_depths: list[int], auto_depth: int = 1000, n_auto: int = 10):
    rows = [("chr1", 1_000_000 + i, auto_depth) for i in range(n_auto)]
    rows += [("chrY", 2_787_394 + i * 100_000, d) for i, d in enumerate(y_depths)]
    rows.append(("chrY", 1_000_000, 10_000))  # PAR1: must be ignored
    return rows


class TestChrYDepth:
    def test_ratio_is_max_over_median(self):
        r, n = chry_depth_ratio(_depth_table([100, 400, 250]))
        assert r == pytest.approx(0.4)
        assert n == 3

    def test_ratio_none_without_chry(self):
        r, n = chry_depth_ratio([("chr1", 1, 100), ("chr2", 2, 200)])
        assert r is None and n == 0

    def test_male_fraction_from_ratio(self):
        ref = _depth_table([800, 1000, 900])  # R_ref = 1.0
        admix = _depth_table([80, 100, 90])  # R_admix = 0.1
        est = chry_male_fraction(admix, ref)
        assert est is not None
        assert est.frac_male == pytest.approx(0.1)
        assert est.n_sites == 3
        assert not est.ci_unreliable
        assert est.ci[0] <= 0.1 <= est.ci[1]

    def test_clipped_to_unit_interval(self):
        est = chry_male_fraction(_depth_table([2000, 2000, 2000]), _depth_table([1000] * 3))
        assert est.frac_male == 1.0

    def test_fewer_than_three_sites_gives_placeholder_ci(self):
        est = chry_male_fraction(_depth_table([100, 100]), _depth_table([1000, 1000]))
        assert est is not None
        assert est.ci == (0.0, 1.0)
        assert est.ci_unreliable

    def test_none_when_reference_has_no_chry_depth(self):
        assert chry_male_fraction(_depth_table([100]), _depth_table([0])) is None

    def test_region_depths_from_vcf(self):
        depths = region_depths_from_vcf(str(MF_VCF), "HOST")
        assert depths and all(len(t) == 3 for t in depths)
        assert all(dp > 0 for _, _, dp in depths)
        with pytest.raises(ValueError):
            region_depths_from_vcf(str(MF_VCF), "NOBODY")


# ---------------------------------------------------------------------------
# Concordance
# ---------------------------------------------------------------------------


class TestConcordance:
    def test_overlapping_cis_concordant(self):
        assert is_concordant(0.10, (0.06, 0.14), 0.15, (0.13, 0.17))

    def test_within_tolerance_without_overlap(self):
        tol = SEXCHROM_CONCORDANCE_PP / 100.0
        assert is_concordant(0.10, (0.099, 0.101), 0.10 + tol * 0.9, (0.115, 0.12))

    def test_discordant(self):
        assert not is_concordant(0.10, (0.09, 0.11), 0.15, (0.14, 0.16))


# ---------------------------------------------------------------------------
# Top-level check
# ---------------------------------------------------------------------------


def _inf(sex: Sex) -> SexInference:
    return SexInference(
        sex=sex,
        n_x_sites=20,
        n_x_het=10 if sex is Sex.FEMALE else 0,
        x_het_rate=0.5 if sex is Sex.FEMALE else 0.0,
        log10_lr=8.0 if sex is Sex.FEMALE else -4.0,
        chry_rel_depth=None,
        n_y_sites=0,
        declared=None,
        source="inferred",
    )


def _sex_result(host: Sex, donors: list[Sex]) -> SexResult:
    return SexResult(
        host=_inf(host), donors=[_inf(d) for d in donors], pair=pair_status(host, donors)
    )


def _mf_inputs():
    host = parse_vcf(MF_VCF, sample="HOST", min_gq=0, gt_ad_consistency=True)
    donor = parse_vcf(MF_VCF, sample="DONOR", min_gq=0, gt_ad_consistency=True)
    admix = parse_vcf(MF_VCF, sample="ADMIX_F0.80", min_dp=0)
    return host, donor, admix


class TestSexMismatchCheck:
    def test_matched_pair_returns_none(self):
        host, donor, admix = _mf_inputs()
        sex = _sex_result(Sex.FEMALE, [Sex.FEMALE])
        assert (
            sex_mismatch_check(
                host,
                [donor],
                admix,
                sex,
                mle_frac_donor=0.8,
                mle_ci=(0.79, 0.81),
                min_dp=0,
                min_gq=0,
            )
            is None
        )

    def test_multi_donor_returns_none(self):
        host, donor, admix = _mf_inputs()
        sex = _sex_result(Sex.MALE, [Sex.FEMALE, Sex.FEMALE])
        assert sex.pair.value == "mismatched"
        assert (
            sex_mismatch_check(
                host,
                [donor, donor],
                admix,
                sex,
                mle_frac_donor=0.8,
                mle_ci=(0.79, 0.81),
                min_dp=0,
                min_gq=0,
            )
            is None
        )

    def test_no_sex_returns_none(self):
        host, donor, admix = _mf_inputs()
        assert (
            sex_mismatch_check(
                host,
                [donor],
                admix,
                None,
                mle_frac_donor=0.8,
                mle_ci=(0.79, 0.81),
                min_dp=0,
                min_gq=0,
            )
            is None
        )

    def test_chrx_basis_on_mf_fixture(self):
        host, donor, admix = _mf_inputs()
        sex = _sex_result(Sex.MALE, [Sex.FEMALE])
        r = sex_mismatch_check(
            host,
            [donor],
            admix,
            sex,
            mle_frac_donor=0.80,
            mle_ci=(0.793, 0.807),
            min_dp=0,
            min_gq=0,
        )
        assert isinstance(r, SexMismatchResult)
        assert r.male_party == "host"
        assert r.basis == CHRX_BASIS
        assert r.n == r.chrx_n >= SEXCHROM_MIN_MARKERS
        assert abs(r.frac_donor - 0.80) < 0.02
        assert r.ci_low <= r.frac_donor <= r.ci_high
        assert r.concordant is True
        assert r.chry_frac_donor is None and r.chry_n_sites == 0

    def test_discordant_against_a_wrong_mle(self):
        host, donor, admix = _mf_inputs()
        sex = _sex_result(Sex.MALE, [Sex.FEMALE])
        r = sex_mismatch_check(
            host,
            [donor],
            admix,
            sex,
            mle_frac_donor=0.95,
            mle_ci=(0.945, 0.955),
            min_dp=0,
            min_gq=0,
        )
        assert r.concordant is False

    def test_too_few_markers_gives_no_basis(self):
        host, donor, admix = _mf_inputs()
        # Keep only four chrX markers (below SEXCHROM_MIN_MARKERS).
        chrx = [m for m in host if m.chrom == "chrX"][: SEXCHROM_MIN_MARKERS - 1]
        keep = {(m.chrom, m.pos) for m in chrx} | {
            (m.chrom, m.pos) for m in host if m.chrom != "chrX"
        }
        host2 = [m for m in host if (m.chrom, m.pos) in keep]
        sex = _sex_result(Sex.MALE, [Sex.FEMALE])
        r = sex_mismatch_check(
            host2, [donor], admix, sex, mle_frac_donor=0.8, mle_ci=(0.79, 0.81), min_dp=0, min_gq=0
        )
        assert r is not None
        assert r.basis is None and r.frac_donor is None and r.concordant is None
        assert r.chrx_n < SEXCHROM_MIN_MARKERS

    def test_chry_depth_takes_priority_and_maps_to_donor_fraction(self):
        host, donor, admix = _mf_inputs()
        sex = _sex_result(Sex.MALE, [Sex.FEMALE])  # host is male: male fraction = host fraction
        ref = _depth_table([1000, 900, 950])
        adm = _depth_table([200, 180, 190])  # male (host) fraction 0.2 -> donor 0.8
        r = sex_mismatch_check(
            host,
            [donor],
            admix,
            sex,
            mle_frac_donor=0.80,
            mle_ci=(0.79, 0.81),
            min_dp=0,
            min_gq=0,
            admix_region_depth=adm,
            male_ref_region_depth=ref,
        )
        assert r.basis == CHRY_BASIS
        assert r.frac_donor == pytest.approx(0.8)
        assert r.n == 3
        assert r.chrx_frac_donor is not None  # both readouts kept
        assert r.concordant is True

    def test_chry_with_too_few_sites_falls_back_to_chrx(self):
        host, donor, admix = _mf_inputs()
        sex = _sex_result(Sex.MALE, [Sex.FEMALE])
        r = sex_mismatch_check(
            host,
            [donor],
            admix,
            sex,
            mle_frac_donor=0.80,
            mle_ci=(0.79, 0.81),
            min_dp=0,
            min_gq=0,
            admix_region_depth=_depth_table([200, 180]),
            male_ref_region_depth=_depth_table([1000, 900]),
        )
        assert r.basis == CHRX_BASIS
        assert r.chry_frac_donor is not None and r.chry_ci_unreliable


# ---------------------------------------------------------------------------
# Through analyse_sample and the simulator
# ---------------------------------------------------------------------------


class TestAnalyseSample:
    def test_mf_fixture_attaches_result(self):
        host, donor, admix = _mf_inputs()
        a = analyse_sample(host, [donor], admix, min_dp=0, min_gq=0, error_rate=0.01)
        sm = a.result.sex_mismatch
        assert sm is not None and sm.basis == CHRX_BASIS
        assert sm.mle_frac_donor == a.result.donor_fraction
        assert sm.concordant is True
        assert a.qc.status != "FAIL"

    def test_ff_fixture_has_none(self):
        v = SEXCHROM_DIR / "joint_FF.vcf"
        host = parse_vcf(v, sample="HOST", min_gq=0)
        donor = parse_vcf(v, sample="DONOR", min_gq=0)
        admix = parse_vcf(v, sample="ADMIX_F0.80", min_dp=0)
        a = analyse_sample(host, [donor], admix, min_dp=0, min_gq=0, error_rate=0.01)
        assert a.result.sex_mismatch is None

    def test_female_minor_direction_in_silico(self, tmp_path):
        """Host F into donor M at 10% host: the double-contrast direction."""
        from allomix.simulate import generate_related_genotypes, generate_sex_chrom_genotypes

        rng = random.Random(21)
        autosomal = generate_related_genotypes(40, "unrelated", rng)
        for i, m in enumerate(autosomal):
            m["chrom"] = f"chr{(i % 22) + 1}"
            m["pos"] = 1_000_000 + i * 100_000
        chrx = generate_sex_chrom_genotypes(20, "F", "M", rng, n_par=0)
        joint = build_joint_vcf_from_genotype_dicts(
            [*autosomal, *chrx], [0.90], ["ADMIX"], target_depth=1500, seed=7
        )
        path = tmp_path / "fm.vcf"
        write_joint_vcf(joint, path)
        host = parse_vcf(path, sample="HOST", min_gq=0)
        donor = parse_vcf(path, sample="DONOR", min_gq=0)
        admix = parse_vcf(path, sample="ADMIX", min_dp=0)
        a = analyse_sample(
            host,
            [donor],
            admix,
            min_dp=0,
            min_gq=0,
            error_rate=0.01,
            declared_host_sex=Sex.FEMALE,
            declared_donor_sexes=[Sex.MALE],
        )
        sm = a.result.sex_mismatch
        assert sm is not None and sm.basis == CHRX_BASIS
        assert sm.male_party == "donor"
        assert abs(sm.frac_donor - 0.90) < 0.02
        assert sm.concordant is True
        assert assess_sex(host, [donor], min_dp=20, min_gq=0).pair.value in (
            "mismatched",
            "unknown",
        )
