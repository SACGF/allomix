"""Tests for allomix.contigs: contig classification and the PAR mask."""

import pytest
from bioutils.par import get_par_map, get_par_meta

from allomix.contigs import (
    PAR_REGIONS,
    ContigClass,
    classify_contig,
    is_sex_chrom,
    normalize_chrom,
)
from allomix.genotype import is_sex_chrom as is_sex_chrom_reexport

# Per-build values as bioutils.par gives them (interbase), recorded in the module docstring.
GRCH37 = {
    "X": {"PAR1": [60000, 2699520], "PAR2": [154931043, 155260560]},
    "Y": {"PAR1": [10000, 2649520], "PAR2": [59034049, 59363566]},
}
GRCH38 = {
    "X": {"PAR1": [10000, 2781479], "PAR2": [155701382, 156030895]},
    "Y": {"PAR1": [10000, 2781479], "PAR2": [56887902, 57217415]},
}


# ---------------------------------------------------------------------------
# PAR boundaries
# ---------------------------------------------------------------------------


class TestParConstants:
    def test_shape(self):
        """Both builds' intervals per chromosome, sorted."""
        assert PAR_REGIONS == {
            "X": (
                (10000, 2781479),
                (60000, 2699520),
                (154931043, 155260560),
                (155701382, 156030895),
            ),
            "Y": (
                (10000, 2649520),
                (10000, 2781479),
                (56887902, 57217415),
                (59034049, 59363566),
            ),
        }

    def test_par1_overlaps_but_par2_does_not(self):
        """Document the build relationship behind the mask.

        PAR1 overlaps between builds. PAR2 does not overlap on either chromosome
        (GRCh38 moved the chrX PAR2 by 770 kb and the chrY PAR2 by 2.1 Mb), so
        the gap between the two PAR2 intervals is non-PAR.
        """

        def overlaps(a, b):
            return max(a[0], b[0]) < min(a[1], b[1])

        for chrom in ("X", "Y"):
            assert overlaps(GRCH37[chrom]["PAR1"], GRCH38[chrom]["PAR1"]), chrom
            assert not overlaps(GRCH37[chrom]["PAR2"], GRCH38[chrom]["PAR2"]), chrom

    @pytest.mark.parametrize(
        ("chrom", "pos", "expected"),
        [
            # The gap between the two builds' chrX PAR2 intervals is non-PAR in
            # both builds (non-PAR Xq28 in GRCh38; past the chromosome end in
            # GRCh37), so it must not be masked.
            ("chrX", 155_260_561, ContigClass.X_NONPAR),
            ("chrX", 155_500_000, ContigClass.X_NONPAR),
            ("chrX", 155_701_382, ContigClass.X_NONPAR),
            # Likewise the chrY PAR2 gap (non-PAR Yq in GRCh37; past the end in
            # GRCh38).
            ("chrY", 57_217_416, ContigClass.Y_NONPAR),
            ("chrY", 58_000_000, ContigClass.Y_NONPAR),
            ("chrY", 59_034_049, ContigClass.Y_NONPAR),
        ],
    )
    def test_par2_gap_is_nonpar(self, chrom, pos, expected):
        assert classify_contig(chrom, pos) is expected


class TestXParEdges:
    @pytest.mark.parametrize(
        ("pos", "expected"),
        [
            # PAR1 start (both builds' starts: GRCh38 10001 is the union start,
            # GRCh37 60001 falls inside the union).
            (10000, ContigClass.X_NONPAR),
            (10001, ContigClass.X_PAR),
            (60000, ContigClass.X_PAR),
            (60001, ContigClass.X_PAR),
            # PAR1 end: GRCh37 ends at 2699520, GRCh38 at 2781479. The band
            # between is PAR in the union.
            (2699520, ContigClass.X_PAR),
            (2699521, ContigClass.X_PAR),
            (2781479, ContigClass.X_PAR),
            (2781480, ContigClass.X_NONPAR),
            # Bulk of the chromosome.
            (100_000_000, ContigClass.X_NONPAR),
            # PAR2, GRCh37 interval 154931044-155260560 (both edges).
            (154931043, ContigClass.X_NONPAR),
            (154931044, ContigClass.X_PAR),
            (155260560, ContigClass.X_PAR),
            (155260561, ContigClass.X_NONPAR),
            # PAR2, GRCh38 interval 155701383-156030895 (both edges).
            (155701382, ContigClass.X_NONPAR),
            (155701383, ContigClass.X_PAR),
            (156030895, ContigClass.X_PAR),
            (156030896, ContigClass.X_NONPAR),
        ],
    )
    def test_boundaries(self, pos, expected):
        assert classify_contig("chrX", pos) is expected
        assert classify_contig("X", pos) is expected


class TestYParEdges:
    @pytest.mark.parametrize(
        ("pos", "expected"),
        [
            # PAR1: both builds start at 10001; GRCh37 ends 2649520, GRCh38
            # 2781479.
            (10000, ContigClass.Y_NONPAR),
            (10001, ContigClass.Y_PAR),
            (2649520, ContigClass.Y_PAR),
            (2649521, ContigClass.Y_PAR),
            (2781479, ContigClass.Y_PAR),
            (2781480, ContigClass.Y_NONPAR),
            # SRY sits 5.9 kb past PAR1 on GRCh38 and must stay non-PAR.
            (2787394, ContigClass.Y_NONPAR),
            (20_000_000, ContigClass.Y_NONPAR),
            # PAR2, GRCh38 interval 56887903-57217415 (both edges).
            (56887902, ContigClass.Y_NONPAR),
            (56887903, ContigClass.Y_PAR),
            (57217415, ContigClass.Y_PAR),
            (57217416, ContigClass.Y_NONPAR),
            # PAR2, GRCh37 interval 59034050-59363566 (both edges).
            (59034049, ContigClass.Y_NONPAR),
            (59034050, ContigClass.Y_PAR),
            (59363566, ContigClass.Y_PAR),
            (59363567, ContigClass.Y_NONPAR),
        ],
    )
    def test_boundaries(self, pos, expected):
        assert classify_contig("chrY", pos) is expected
        assert classify_contig("Y", pos) is expected


# ---------------------------------------------------------------------------
# Contig names
# ---------------------------------------------------------------------------


class TestNames:
    @pytest.mark.parametrize("n", range(1, 23))
    def test_autosomes(self, n):
        assert classify_contig(f"chr{n}", 1_000_000) is ContigClass.AUTOSOME
        assert classify_contig(str(n), 1_000_000) is ContigClass.AUTOSOME
        assert not is_sex_chrom(f"chr{n}")
        assert not is_sex_chrom(str(n))

    @pytest.mark.parametrize("name", ["chrX", "X", "chrx", "x", "CHRX", "Chrx"])
    def test_x_prefix_and_case(self, name):
        assert classify_contig(name, 50_000_000) is ContigClass.X_NONPAR
        assert classify_contig(name, 1_000_000) is ContigClass.X_PAR
        assert is_sex_chrom(name)

    @pytest.mark.parametrize("name", ["chrY", "Y", "chry", "y"])
    def test_y_prefix_and_case(self, name):
        assert classify_contig(name, 20_000_000) is ContigClass.Y_NONPAR
        assert classify_contig(name, 1_000_000) is ContigClass.Y_PAR
        assert is_sex_chrom(name)

    @pytest.mark.parametrize("name", ["chrM", "M", "chrMT", "MT", "chrm", "mt"])
    def test_mito(self, name):
        assert classify_contig(name, 1) is ContigClass.MT
        assert classify_contig(name, 16_569) is ContigClass.MT
        assert is_sex_chrom(name)

    @pytest.mark.parametrize(
        "name",
        [
            "chr1_KI270706v1_random",
            "chrUn_KN707606v1_decoy",
            "chrX_KI270880v1_alt",
            "chrY_KI270740v1_random",
            "HLA-A*01:01:01:01",
            "chrEBV",
            "GL000220.1",
            "chr23",
            "chr0",
            "",
        ],
    )
    def test_other(self, name):
        assert classify_contig(name, 1000) is ContigClass.OTHER
        assert is_sex_chrom(name)

    def test_position_irrelevant_off_xy(self):
        """Only X and Y are split by position."""
        for pos in (1, 10001, 1_000_000, 300_000_000):
            assert classify_contig("chr1", pos) is ContigClass.AUTOSOME
            assert classify_contig("chrM", pos) is ContigClass.MT
            assert classify_contig("chrUn_decoy", pos) is ContigClass.OTHER

    def test_normalize_chrom(self):
        assert normalize_chrom("chrX") == "X"
        assert normalize_chrom("chrx") == "X"
        assert normalize_chrom("X") == "X"
        assert normalize_chrom("chr1_KI270706v1_random") == "1_KI270706V1_RANDOM"
        # Only the leading "chr" is stripped.
        assert normalize_chrom("Xchr") == "XCHR"


class TestReexport:
    def test_genotype_reexports_is_sex_chrom(self):
        """Existing ``from allomix.genotype import is_sex_chrom`` keeps working."""
        assert is_sex_chrom_reexport is is_sex_chrom


# ---------------------------------------------------------------------------
# bioutils.par source data
# ---------------------------------------------------------------------------


def test_bioutils_par_matches_recorded_builds():
    """bioutils.par still carries the NCBI values recorded in the module docstring.

    ``PAR_REGIONS`` is built from ``bioutils.par`` at import, so this pins the
    source data: a change in a bioutils release fails here rather than moving
    the mask silently.
    """
    # classify_contig tests a 1-based POS as start_i < pos <= end_i, which
    # relies on bioutils keeping interbase coordinates.
    assert get_par_meta("GRCh38")["coordinates"].startswith("interbase")
    for build, expected in (("GRCh37", GRCH37), ("GRCh38", GRCH38)):
        assert get_par_map(build) == expected, build
