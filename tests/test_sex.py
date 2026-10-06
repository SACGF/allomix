"""Tests for allomix.qc.sex: chrX-heterozygosity sex inference and pair status."""

import pytest

from allomix.constants import (
    CHRY_REL_DEPTH_FEMALE_WEAK,
    CHRY_REL_DEPTH_MALE,
    CHRY_REL_DEPTH_MALE_WEAK,
    MIN_X_SITES,
    SEX_LR_THRESHOLD,
)
from allomix.genotype import MarkerData
from allomix.qc.sex import (
    PairStatus,
    Sex,
    assess_sex,
    grade_chry_depth,
    infer_sex,
    pair_status,
    parse_declared_sex,
)
from allomix.sex_types import PairStatus as PairStatusLeaf
from allomix.sex_types import Sex as SexLeaf

# Non-PAR chrX on both builds (PAR1 ends at 2.78 Mb, PAR2 starts at 154.9 Mb).
X_NONPAR_START = 10_000_000
# Inside PAR1 on both builds.
X_PAR1_POS = 1_000_000


def _m(chrom: str, pos: int, gt: tuple[int, int], dp: int = 100, gq: int | None = 99) -> MarkerData:
    half = dp // 2
    return MarkerData(chrom, pos, "A", "G", gt, half, dp - half, dp, gq)


def _autosomes(n: int = 40, n_het: int = 20, **kw) -> list[MarkerData]:
    """``n`` chr1 markers, the first ``n_het`` heterozygous."""
    return [_m("chr1", 1000 * (i + 1), (0, 1) if i < n_het else (0, 0), **kw) for i in range(n)]


def _x(n: int, n_het: int, pos0: int = X_NONPAR_START, **kw) -> list[MarkerData]:
    """``n`` chrX markers from ``pos0`` upward, the first ``n_het`` heterozygous."""
    return [_m("chrX", pos0 + 1000 * i, (0, 1) if i < n_het else (1, 1), **kw) for i in range(n)]


class TestInferSex:
    def test_clean_female(self):
        inf = infer_sex(_autosomes() + _x(25, 12), min_dp=20, min_gq=20)
        assert inf.sex is Sex.FEMALE
        assert inf.n_x_sites == 25
        assert inf.n_x_het == 12
        assert inf.x_het_rate == pytest.approx(12 / 25)
        assert inf.log10_lr is not None and inf.log10_lr >= SEX_LR_THRESHOLD
        assert inf.source == "inferred"
        assert inf.effective is Sex.FEMALE
        assert inf.chry_rel_depth is None

    def test_clean_male(self):
        inf = infer_sex(_autosomes() + _x(25, 0), min_dp=20, min_gq=20)
        assert inf.sex is Sex.MALE
        assert inf.n_x_het == 0
        assert inf.log10_lr is not None and inf.log10_lr <= -SEX_LR_THRESHOLD
        assert inf.effective is Sex.MALE

    @pytest.mark.parametrize("n_spurious", [1, 2])
    def test_male_with_spurious_hets(self, n_spurious):
        """A couple of diploid-caller het errors on a hemizygous X stay MALE."""
        inf = infer_sex(_autosomes() + _x(25, n_spurious), min_dp=20, min_gq=20)
        assert inf.sex is Sex.MALE
        assert inf.n_x_het == n_spurious

    def test_intermediate_het_rate_is_ambiguous(self):
        """Neither model fits 4/25 hets well: no confident call, no declaration."""
        inf = infer_sex(_autosomes() + _x(25, 4), min_dp=20, min_gq=20)
        assert inf.sex is Sex.AMBIGUOUS
        assert inf.log10_lr is not None
        assert abs(inf.log10_lr) < SEX_LR_THRESHOLD
        assert inf.source == "inferred"
        assert inf.effective is Sex.AMBIGUOUS

    def test_too_few_sites_unavailable(self):
        inf = infer_sex(_autosomes() + _x(MIN_X_SITES - 1, 2), min_dp=20, min_gq=20)
        assert inf.sex is Sex.UNAVAILABLE
        assert inf.n_x_sites == MIN_X_SITES - 1
        assert inf.log10_lr is None
        assert inf.source == "unavailable"

    def test_no_chrx_unavailable(self):
        inf = infer_sex(_autosomes(), min_dp=20, min_gq=20)
        assert inf.sex is Sex.UNAVAILABLE
        assert inf.n_x_sites == 0
        assert inf.x_het_rate is None
        assert inf.effective is Sex.UNAVAILABLE

    def test_par_sites_do_not_count(self):
        """PAR1 chrX sites are diploid in both sexes and are ignored."""
        inf = infer_sex(_autosomes() + _x(25, 12, pos0=X_PAR1_POS), min_dp=20, min_gq=20)
        assert inf.sex is Sex.UNAVAILABLE
        assert inf.n_x_sites == 0

    def test_depth_filter(self):
        inf = infer_sex(_autosomes() + _x(25, 12, dp=10), min_dp=20, min_gq=20)
        assert inf.sex is Sex.UNAVAILABLE
        assert inf.n_x_sites == 0

    def test_gq_filter(self):
        inf = infer_sex(_autosomes() + _x(25, 12, gq=5), min_dp=20, min_gq=20)
        assert inf.sex is Sex.UNAVAILABLE
        # Absent GQ (e.g. bcftools output without the field) passes the filter.
        inf = infer_sex(_autosomes() + _x(25, 12, gq=None), min_dp=20, min_gq=20)
        assert inf.sex is Sex.FEMALE

    def test_uncalled_or_multiallelic_sites_skipped(self):
        bad = [_m("chrX", X_NONPAR_START + 1, (-1, -1)), _m("chrX", X_NONPAR_START + 2, (0, 2))]
        inf = infer_sex(_autosomes() + _x(25, 12) + bad, min_dp=20, min_gq=20)
        assert inf.n_x_sites == 25

    def test_y_sites_counted_only(self):
        y = [_m("chrY", 10_000_000 + i, (0, 0)) for i in range(3)]
        inf = infer_sex(_autosomes() + _x(25, 0) + y, min_dp=20, min_gq=20)
        assert inf.n_y_sites == 3
        assert inf.sex is Sex.MALE

    def test_no_autosomes_uses_default_female_rate(self):
        """Without autosomal sites the female model falls back to p=0.5."""
        inf = infer_sex(_x(25, 12), min_dp=20, min_gq=20)
        assert inf.sex is Sex.FEMALE

    def test_declared_resolves_ambiguous(self):
        inf = infer_sex(_autosomes() + _x(25, 4), min_dp=20, min_gq=20, declared=Sex.FEMALE)
        assert inf.sex is Sex.AMBIGUOUS
        assert inf.declared is Sex.FEMALE
        assert inf.source == "declared"
        assert inf.effective is Sex.FEMALE

    def test_declared_resolves_unavailable(self):
        inf = infer_sex(_autosomes(), min_dp=20, min_gq=20, declared=Sex.MALE)
        assert inf.sex is Sex.UNAVAILABLE
        assert inf.source == "declared"
        assert inf.effective is Sex.MALE

    def test_declared_agrees(self):
        inf = infer_sex(_autosomes() + _x(25, 12), min_dp=20, min_gq=20, declared=Sex.FEMALE)
        assert inf.source == "inferred+declared"
        assert inf.effective is Sex.FEMALE

    def test_conflict_keeps_inferred(self):
        """A confident inference against the label: conflict, data win for routing."""
        inf = infer_sex(_autosomes() + _x(25, 12), min_dp=20, min_gq=20, declared=Sex.MALE)
        assert inf.sex is Sex.FEMALE
        assert inf.source == "conflict"
        assert inf.effective is Sex.FEMALE


def _depths(rel: float, auto_depth: int = 1000, n_auto: int = 10) -> list[tuple[str, int, int]]:
    """Forced-pileup depths whose chrY relative depth is ``rel``.

    Three non-PAR chrY sites (SRY, ZFY, AMELY positions) with the max at
    ``rel * auto_depth``, plus a chrY PAR1 site at high depth that must be ignored.
    """
    rows = [("chr1", 1_000_000 + i, auto_depth) for i in range(n_auto)]
    top = round(rel * auto_depth)
    rows += [("chrY", 2_787_394, top), ("chrY", 2_979_984, top // 2), ("chrY", 6_869_956, 0)]
    rows.append(("chrY", 1_000_000, 10 * auto_depth))
    return rows


class TestChrYSecondary:
    @pytest.mark.parametrize(
        "rel, expected",
        [
            (None, None),
            (0.0, "F"),
            (CHRY_REL_DEPTH_FEMALE_WEAK, "F"),
            (0.1, "F*"),
            (CHRY_REL_DEPTH_MALE_WEAK, "F*"),
            (0.3, "M*"),
            (CHRY_REL_DEPTH_MALE, "M*"),
            (0.51, "M"),
            (1.2, "M"),
        ],
    )
    def test_grade_matches_lab_thresholds(self, rel, expected):
        assert grade_chry_depth(rel) == expected

    def test_no_depth_input(self):
        inf = infer_sex(_autosomes() + _x(25, 4), min_dp=20, min_gq=20)
        assert inf.chry_rel_depth is None
        assert inf.chry_call is None
        assert not inf.chry_resolved

    def test_depth_ratio_ignores_par(self):
        inf = infer_sex(_autosomes() + _x(25, 12), min_dp=20, min_gq=20, region_depth=_depths(0.8))
        assert inf.chry_rel_depth == pytest.approx(0.8)
        assert inf.chry_call == "M"

    def test_male_resolves_ambiguous(self):
        """The plan's known case: a male with spurious chrX hets is ambiguous on chrX alone."""
        markers = _autosomes() + _x(25, 4)
        assert infer_sex(markers, min_dp=20, min_gq=20).sex is Sex.AMBIGUOUS
        inf = infer_sex(markers, min_dp=20, min_gq=20, region_depth=_depths(0.8))
        assert inf.sex is Sex.MALE
        assert inf.chry_resolved
        assert inf.source == "inferred"
        assert inf.effective is Sex.MALE

    def test_male_resolves_unavailable(self):
        inf = infer_sex(_autosomes(), min_dp=20, min_gq=20, region_depth=_depths(0.8))
        assert inf.sex is Sex.MALE
        assert inf.chry_resolved
        assert inf.n_x_sites == 0

    @pytest.mark.parametrize("rel", [0.3, 0.1, 0.0])
    def test_weak_or_absent_chry_does_not_resolve(self, rel):
        """``M*`` is too weak; ``F*`` / ``F`` (absent chrY) are not evidence for a female."""
        inf = infer_sex(_autosomes() + _x(25, 4), min_dp=20, min_gq=20, region_depth=_depths(rel))
        assert inf.sex is Sex.AMBIGUOUS
        assert not inf.chry_resolved

    def test_never_overrules_confident_chrx(self):
        inf = infer_sex(_autosomes() + _x(25, 12), min_dp=20, min_gq=20, region_depth=_depths(0.8))
        assert inf.sex is Sex.FEMALE
        assert inf.chry_call == "M"
        assert not inf.chry_resolved

    def test_chry_male_against_declared_female_is_conflict(self):
        inf = infer_sex(
            _autosomes() + _x(25, 4),
            min_dp=20,
            min_gq=20,
            declared=Sex.FEMALE,
            region_depth=_depths(0.8),
        )
        assert inf.sex is Sex.MALE
        assert inf.source == "conflict"

    def test_assess_sex_aligns_donor_depths(self):
        host = _autosomes() + _x(25, 4)
        donors = [_autosomes() + _x(25, 4), _autosomes() + _x(25, 4)]
        res = assess_sex(
            host,
            donors,
            min_dp=20,
            min_gq=20,
            host_region_depth=_depths(0.8),
            donor_region_depths=[None, _depths(0.9)],
        )
        assert res.host.sex is Sex.MALE
        assert res.donors[0].sex is Sex.AMBIGUOUS
        assert res.donors[1].sex is Sex.MALE
        assert res.pair is PairStatus.UNKNOWN
        res = assess_sex(
            host,
            donors[1:],
            min_dp=20,
            min_gq=20,
            host_region_depth=_depths(0.8),
            donor_region_depths=[_depths(0.9)],
        )
        assert res.pair is PairStatus.MATCHED_MALE


class TestPairStatus:
    @pytest.mark.parametrize(
        "host, donors, expected",
        [
            (Sex.FEMALE, [Sex.FEMALE], PairStatus.MATCHED_FEMALE),
            (Sex.MALE, [Sex.MALE], PairStatus.MATCHED_MALE),
            (Sex.FEMALE, [Sex.MALE], PairStatus.MISMATCHED),
            (Sex.MALE, [Sex.FEMALE], PairStatus.MISMATCHED),
            (Sex.FEMALE, [Sex.FEMALE, Sex.FEMALE], PairStatus.MATCHED_FEMALE),
            (Sex.FEMALE, [Sex.FEMALE, Sex.MALE], PairStatus.MISMATCHED),
            (Sex.AMBIGUOUS, [Sex.FEMALE], PairStatus.UNKNOWN),
            (Sex.FEMALE, [Sex.UNAVAILABLE], PairStatus.UNKNOWN),
            (Sex.MALE, [Sex.MALE, Sex.AMBIGUOUS], PairStatus.UNKNOWN),
            (Sex.UNAVAILABLE, [Sex.UNAVAILABLE], PairStatus.UNKNOWN),
        ],
    )
    def test_matrix(self, host, donors, expected):
        assert pair_status(host, donors) is expected


class TestParseDeclaredSex:
    @pytest.mark.parametrize(
        "text, expected",
        [
            ("F", Sex.FEMALE),
            ("f", Sex.FEMALE),
            (" female ", Sex.FEMALE),
            ("FEMALE", Sex.FEMALE),
            ("M", Sex.MALE),
            ("male", Sex.MALE),
            ("Male", Sex.MALE),
            ("unknown", None),
            ("X", None),
            ("", None),
            (None, None),
        ],
    )
    def test_parse(self, text, expected):
        assert parse_declared_sex(text) is expected


class TestAssessSex:
    def test_mismatched_pair(self):
        host = _autosomes() + _x(25, 12)
        donor = _autosomes() + _x(25, 0)
        res = assess_sex(host, [donor], min_dp=20, min_gq=20)
        assert res.host.sex is Sex.FEMALE
        assert res.donors[0].sex is Sex.MALE
        assert res.pair is PairStatus.MISMATCHED

    def test_declared_donors_aligned(self):
        host = _autosomes() + _x(25, 12)
        donors = [_autosomes(), _autosomes() + _x(25, 12)]
        res = assess_sex(
            host,
            donors,
            min_dp=20,
            min_gq=20,
            declared_host=Sex.FEMALE,
            declared_donors=[Sex.FEMALE, None],
        )
        assert res.host.source == "inferred+declared"
        assert res.donors[0].source == "declared"
        assert res.donors[0].effective is Sex.FEMALE
        assert res.donors[1].source == "inferred"
        assert res.pair is PairStatus.MATCHED_FEMALE

    def test_autosome_only_panel_unknown(self):
        res = assess_sex(_autosomes(), [_autosomes()], min_dp=20, min_gq=20)
        assert res.host.sex is Sex.UNAVAILABLE
        assert res.donors[0].sex is Sex.UNAVAILABLE
        assert res.pair is PairStatus.UNKNOWN


class TestReexports:
    def test_enums_are_the_leaf_module_objects(self):
        """``qc.sex`` re-exports the ``sex_types`` enums (no duplicate definitions)."""
        assert Sex is SexLeaf
        assert PairStatus is PairStatusLeaf
