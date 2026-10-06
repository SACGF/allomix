"""Tests for allomix.contigs: contig classification and the PAR mask."""

import pytest

from allomix.contigs import (
    PAR_REGIONS,
    ContigClass,
    classify_contig,
    is_sex_chrom,
    normalize_chrom,
)
from allomix.genotype import is_sex_chrom as is_sex_chrom_reexport

# Per-build source values, 1-based inclusive, as recorded in the module comment.
GRCH37 = {
    "X": {"PAR1": (60001, 2699520), "PAR2": (154931044, 155260560)},
    "Y": {"PAR1": (10001, 2649520), "PAR2": (59034050, 59363566)},
}
GRCH38 = {
    "X": {"PAR1": (10001, 2781479), "PAR2": (155701383, 156030895)},
    "Y": {"PAR1": (10001, 2781479), "PAR2": (56887903, 57217415)},
}


def _union(*par_maps):
    """Exact union of PAR intervals across builds: pool per chromosome, merge overlaps."""
    out = {}
    for chrom in sorted({c for m in par_maps for c in m}):
        ivs = sorted(tuple(iv) for m in par_maps for iv in m.get(chrom, {}).values())
        merged: list[tuple[int, int]] = []
        for s_, e_ in ivs:
            if merged and s_ <= merged[-1][1] + 1:
                merged[-1] = (merged[-1][0], max(merged[-1][1], e_))
            else:
                merged.append((s_, e_))
        out[chrom] = tuple(merged)
    return out


# ---------------------------------------------------------------------------
# PAR boundaries
# ---------------------------------------------------------------------------


class TestParConstants:
    def test_union_of_recorded_builds(self):
        """The module constants are exactly the union of the two builds' values."""
        assert PAR_REGIONS == _union(GRCH37, GRCH38)

    def test_shape(self):
        """PAR1 merges into one interval; each PAR2 stays as two."""
        assert PAR_REGIONS == {
            "X": ((10001, 2781479), (154931044, 155260560), (155701383, 156030895)),
            "Y": ((10001, 2781479), (56887903, 57217415), (59034050, 59363566)),
        }
        for ivs in PAR_REGIONS.values():
            assert list(ivs) == sorted(ivs)
            assert all(s_ <= e_ for s_, e_ in ivs)
            assert all(ivs[i][1] < ivs[i + 1][0] for i in range(len(ivs) - 1))

    def test_par1_overlaps_but_par2_does_not(self):
        """Document the build relationship behind the constants' shape.

        PAR1 overlaps between builds and merges to one interval. PAR2 does not
        overlap on either chromosome (GRCh38 moved the chrX PAR2 by 770 kb and
        the chrY PAR2 by 2.1 Mb), so each PAR2 is kept as two intervals and the
        gap between them is non-PAR.
        """

        def overlaps(a, b):
            return max(a[0], b[0]) <= min(a[1], b[1])

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
# Cross-check against bioutils.par (fork davmlaw/bioutils@add-par-regions)
# ---------------------------------------------------------------------------


def test_constants_match_bioutils_union():
    """Constants equal the union of the two NCBI builds as packaged in bioutils.par.

    Skipped unless the fork is installed:
    ``pip install "git+https://github.com/davmlaw/bioutils@add-par-regions"``.
    """
    bioutils_par = pytest.importorskip("bioutils.par")
    grch37 = bioutils_par.get_par_map("GRCh37")
    grch38 = bioutils_par.get_par_map("GRCh38")
    # bioutils returns {chrom: {par_name: [start_i, end_i]}} in interbase
    # (0-based, half-open) coordinates like bioutils.cytobands; allomix keeps
    # the NCBI 1-based inclusive values, so convert as (start_i + 1, end_i).
    assert bioutils_par.get_par_meta("GRCh38")["coordinates"].startswith("interbase")
    as_tuples37 = {c: {n: (iv[0] + 1, iv[1]) for n, iv in r.items()} for c, r in grch37.items()}
    as_tuples38 = {c: {n: (iv[0] + 1, iv[1]) for n, iv in r.items()} for c, r in grch38.items()}
    assert as_tuples37 == GRCH37
    assert as_tuples38 == GRCH38
    assert PAR_REGIONS == _union(as_tuples37, as_tuples38)
