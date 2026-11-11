"""Tests for allomix.contigs: contig classification and the PAR mask."""

import pytest

from allomix.contigs import (
    ContigClass,
    classify_contig,
    is_sex_chrom,
    normalize_chrom,
)
from allomix.genotype import is_sex_chrom as is_sex_chrom_reexport

# ---------------------------------------------------------------------------
# PAR boundaries
# ---------------------------------------------------------------------------


class TestPar2Gap:
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
