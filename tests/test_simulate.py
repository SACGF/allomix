"""Tests for allomix.simulate — synthetic chimeric VCF generation."""

from __future__ import annotations

import random
import textwrap
from pathlib import Path

import pytest

from allomix.simulate import (
    alt_dose,
    blend_vcfs,
    expected_vaf,
    extract_depth,
    extract_gt,
    gt_from_counts,
    is_informative,
    parse_vcf,
    sample_allele_counts,
    write_vcf,
)

# ---------------------------------------------------------------------------
# Helpers for creating minimal VCF files
# ---------------------------------------------------------------------------

MINIMAL_HEADER = textwrap.dedent("""\
    ##fileformat=VCFv4.2
    ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
    ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
    ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
    ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
    ##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled likelihoods">
    ##FORMAT=<ID=AF,Number=A,Type=Float,Description="Variant allele frequency">
    ##contig=<ID=chr1,length=248956422>
""")


def _write_test_vcf(
    path: Path,
    sample_name: str,
    records: list[tuple[str, int, str, str, str, str]],
) -> None:
    """Write a minimal test VCF.

    Each record is (chrom, pos, ref, alt, format_str, sample_str).
    """
    with open(path, "w") as fh:
        fh.write(MINIMAL_HEADER)
        fh.write(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}\n")
        for chrom, pos, ref, alt, fmt, samp in records:
            fh.write(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t100\tPASS\t.\t{fmt}\t{samp}\n")


# ---------------------------------------------------------------------------
# Tests: alt_dose
# ---------------------------------------------------------------------------


class TestAltDose:
    def test_hom_ref(self) -> None:
        assert alt_dose((0, 0)) == 0

    def test_het(self) -> None:
        assert alt_dose((0, 1)) == 1
        assert alt_dose((1, 0)) == 1

    def test_hom_alt(self) -> None:
        assert alt_dose((1, 1)) == 2


# ---------------------------------------------------------------------------
# Tests: expected_vaf
# ---------------------------------------------------------------------------


class TestExpectedVaf:
    """Test expected VAF for all 9 genotype combinations."""

    @pytest.mark.parametrize(
        "host_gt, donor_gt, frac, exp_vaf",
        [
            # Both hom-ref -> always 0
            ((0, 0), (0, 0), 0.0, 0.0),
            ((0, 0), (0, 0), 0.5, 0.0),
            ((0, 0), (0, 0), 1.0, 0.0),
            # host 0/0, donor 0/1 -> f * 1 / 2
            ((0, 0), (0, 1), 0.0, 0.0),
            ((0, 0), (0, 1), 0.5, 0.25),
            ((0, 0), (0, 1), 1.0, 0.5),
            # host 0/0, donor 1/1 -> f * 2 / 2 = f
            ((0, 0), (1, 1), 0.0, 0.0),
            ((0, 0), (1, 1), 0.5, 0.5),
            ((0, 0), (1, 1), 1.0, 1.0),
            # host 0/1, donor 0/0 -> (1-f) * 1 / 2
            ((0, 1), (0, 0), 0.0, 0.5),
            ((0, 1), (0, 0), 0.5, 0.25),
            ((0, 1), (0, 0), 1.0, 0.0),
            # host 0/1, donor 0/1 -> always 0.5
            ((0, 1), (0, 1), 0.0, 0.5),
            ((0, 1), (0, 1), 0.5, 0.5),
            ((0, 1), (0, 1), 1.0, 0.5),
            # host 0/1, donor 1/1 -> ((1-f) + 2f) / 2
            ((0, 1), (1, 1), 0.0, 0.5),
            ((0, 1), (1, 1), 0.5, 0.75),
            ((0, 1), (1, 1), 1.0, 1.0),
            # host 1/1, donor 0/0 -> (1-f) * 2 / 2 = 1-f
            ((1, 1), (0, 0), 0.0, 1.0),
            ((1, 1), (0, 0), 0.5, 0.5),
            ((1, 1), (0, 0), 1.0, 0.0),
            # host 1/1, donor 0/1 -> ((1-f)*2 + f*1) / 2
            ((1, 1), (0, 1), 0.0, 1.0),
            ((1, 1), (0, 1), 0.5, 0.75),
            ((1, 1), (0, 1), 1.0, 0.5),
            # host 1/1, donor 1/1 -> always 1.0
            ((1, 1), (1, 1), 0.0, 1.0),
            ((1, 1), (1, 1), 0.5, 1.0),
            ((1, 1), (1, 1), 1.0, 1.0),
        ],
    )
    def test_expected_vaf(
        self,
        host_gt: tuple[int, int],
        donor_gt: tuple[int, int],
        frac: float,
        exp_vaf: float,
    ) -> None:
        result = expected_vaf(host_gt, donor_gt, frac)
        assert result == pytest.approx(exp_vaf, abs=1e-10)


# ---------------------------------------------------------------------------
# Tests: is_informative
# ---------------------------------------------------------------------------


class TestIsInformative:
    """Informative = different alt dose between host and donor."""

    @pytest.mark.parametrize(
        "host_gt, donor_gt, expected",
        [
            ((0, 0), (0, 0), False),  # same dose 0
            ((0, 0), (0, 1), True),  # 0 vs 1
            ((0, 0), (1, 1), True),  # 0 vs 2
            ((0, 1), (0, 0), True),  # 1 vs 0
            ((0, 1), (0, 1), False),  # same dose 1
            ((0, 1), (1, 1), True),  # 1 vs 2
            ((1, 1), (0, 0), True),  # 2 vs 0
            ((1, 1), (0, 1), True),  # 2 vs 1
            ((1, 1), (1, 1), False),  # same dose 2
        ],
    )
    def test_informativeness(
        self,
        host_gt: tuple[int, int],
        donor_gt: tuple[int, int],
        expected: bool,
    ) -> None:
        assert is_informative(host_gt, donor_gt) is expected


# ---------------------------------------------------------------------------
# Tests: sample_allele_counts
# ---------------------------------------------------------------------------


class TestSampleAlleleCounts:
    def test_zero_depth(self) -> None:
        ref, alt = sample_allele_counts(0.5, 0)
        assert ref == 0
        assert alt == 0

    def test_counts_sum_to_depth(self) -> None:
        rng = random.Random(42)
        for _ in range(50):
            depth = rng.randint(100, 5000)
            vaf = rng.random()
            ref, alt = sample_allele_counts(vaf, depth, rng)
            assert ref + alt == depth

    def test_vaf_zero_gives_all_ref(self) -> None:
        ref, alt = sample_allele_counts(0.0, 1000, random.Random(1))
        assert alt == 0
        assert ref == 1000

    def test_vaf_one_gives_all_alt(self) -> None:
        ref, alt = sample_allele_counts(1.0, 1000, random.Random(1))
        assert ref == 0
        assert alt == 1000

    def test_binomial_in_expected_range(self) -> None:
        """With depth=10000 and vaf=0.5, ALT count should be near 5000."""
        rng = random.Random(123)
        alt_counts = [sample_allele_counts(0.5, 10000, rng)[1] for _ in range(20)]
        mean_alt = sum(alt_counts) / len(alt_counts)
        # Should be within ~2% of 5000
        assert 4800 < mean_alt < 5200


# ---------------------------------------------------------------------------
# Tests: gt_from_counts
# ---------------------------------------------------------------------------


class TestGtFromCounts:
    def test_hom_ref(self) -> None:
        assert gt_from_counts(1000, 0) == "0/0"

    def test_het(self) -> None:
        assert gt_from_counts(500, 500) == "0/1"

    def test_hom_alt(self) -> None:
        assert gt_from_counts(0, 1000) == "1/1"

    def test_low_alt(self) -> None:
        # 4% ALT -> hom ref
        assert gt_from_counts(960, 40) == "0/0"

    def test_high_alt(self) -> None:
        # 96% ALT -> hom alt
        assert gt_from_counts(40, 960) == "1/1"

    def test_zero_depth(self) -> None:
        assert gt_from_counts(0, 0) == "./."


# ---------------------------------------------------------------------------
# Tests: extract_gt and extract_depth
# ---------------------------------------------------------------------------


class TestExtractGt:
    def test_simple_het(self, tmp_path: Path) -> None:
        _write_test_vcf(
            tmp_path / "test.vcf",
            "SAMPLE",
            [
                ("chr1", 100, "A", "T", "GT:AD:DP:GQ:PL:AF", "0/1:500,500:1000:99:100,0,100:0.5"),
            ],
        )
        _, records = parse_vcf(tmp_path / "test.vcf")
        assert extract_gt(records[0]) == (0, 1)

    def test_hom_alt(self, tmp_path: Path) -> None:
        _write_test_vcf(
            tmp_path / "test.vcf",
            "SAMPLE",
            [
                ("chr1", 100, "A", "T", "GT:AD:DP:GQ:PL:AF", "1/1:0,1000:1000:99:100,100,0:1.0"),
            ],
        )
        _, records = parse_vcf(tmp_path / "test.vcf")
        assert extract_gt(records[0]) == (1, 1)

    def test_nocall(self, tmp_path: Path) -> None:
        _write_test_vcf(
            tmp_path / "test.vcf",
            "SAMPLE",
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "./.:.:."),
            ],
        )
        _, records = parse_vcf(tmp_path / "test.vcf")
        assert extract_gt(records[0]) is None


class TestExtractDepth:
    def test_from_dp(self, tmp_path: Path) -> None:
        _write_test_vcf(
            tmp_path / "test.vcf",
            "SAMPLE",
            [
                ("chr1", 100, "A", "T", "GT:AD:DP:GQ:PL:AF", "0/1:500,500:1000:99:100,0,100:0.5"),
            ],
        )
        _, records = parse_vcf(tmp_path / "test.vcf")
        assert extract_depth(records[0]) == 1000

    def test_from_ad_fallback(self, tmp_path: Path) -> None:
        _write_test_vcf(
            tmp_path / "test.vcf",
            "SAMPLE",
            [
                ("chr1", 100, "A", "T", "GT:AD", "0/1:600,400"),
            ],
        )
        _, records = parse_vcf(tmp_path / "test.vcf")
        assert extract_depth(records[0]) == 1000


# ---------------------------------------------------------------------------
# Tests: blend_vcfs end-to-end
# ---------------------------------------------------------------------------


def _make_pair_vcfs(
    tmp_path: Path,
    host_records: list[tuple[str, int, str, str, str, str]],
    donor_records: list[tuple[str, int, str, str, str, str]],
) -> tuple[Path, Path]:
    """Create host and donor VCFs in tmp_path and return their paths."""
    host_path = tmp_path / "host.vcf"
    donor_path = tmp_path / "donor.vcf"
    _write_test_vcf(host_path, "HOST", host_records)
    _write_test_vcf(donor_path, "DONOR", donor_records)
    return host_path, donor_path


class TestBlendVcfs:
    def test_fraction_zero_matches_host(self, tmp_path: Path) -> None:
        """At f=0, output should match host genotypes."""
        host_path, donor_path = _make_pair_vcfs(
            tmp_path,
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "0/0:1000,0:1000"),
                ("chr1", 200, "A", "T", "GT:AD:DP", "0/1:500,500:1000"),
                ("chr1", 300, "A", "T", "GT:AD:DP", "1/1:0,1000:1000"),
            ],
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "1/1:0,1000:1000"),
                ("chr1", 200, "A", "T", "GT:AD:DP", "0/0:1000,0:1000"),
                ("chr1", 300, "A", "T", "GT:AD:DP", "0/1:500,500:1000"),
            ],
        )

        result = blend_vcfs(host_path, donor_path, 0.0, target_depth=2000, seed=42)
        assert result.num_markers == 3

        # Parse the output records to check GTs
        gts = []
        for line in result.records:
            sample = line.split("\t")[9]
            gt = sample.split(":")[0]
            gts.append(gt)
        assert gts == ["0/0", "0/1", "1/1"]

    def test_fraction_one_matches_donor(self, tmp_path: Path) -> None:
        """At f=1, output should match donor genotypes."""
        host_path, donor_path = _make_pair_vcfs(
            tmp_path,
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "0/0:1000,0:1000"),
                ("chr1", 200, "A", "T", "GT:AD:DP", "0/1:500,500:1000"),
                ("chr1", 300, "A", "T", "GT:AD:DP", "1/1:0,1000:1000"),
            ],
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "1/1:0,1000:1000"),
                ("chr1", 200, "A", "T", "GT:AD:DP", "0/0:1000,0:1000"),
                ("chr1", 300, "A", "T", "GT:AD:DP", "0/1:500,500:1000"),
            ],
        )

        result = blend_vcfs(host_path, donor_path, 1.0, target_depth=2000, seed=42)
        gts = []
        for line in result.records:
            sample = line.split("\t")[9]
            gt = sample.split(":")[0]
            gts.append(gt)
        assert gts == ["1/1", "0/0", "0/1"]

    def test_only_shared_loci(self, tmp_path: Path) -> None:
        """Only loci present in both VCFs should appear in output."""
        host_path, donor_path = _make_pair_vcfs(
            tmp_path,
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "0/1:500,500:1000"),
                ("chr1", 200, "A", "T", "GT:AD:DP", "0/0:1000,0:1000"),
                ("chr1", 999, "G", "C", "GT:AD:DP", "1/1:0,1000:1000"),
            ],
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "1/1:0,1000:1000"),
                ("chr1", 300, "A", "T", "GT:AD:DP", "0/1:500,500:1000"),
            ],
        )

        result = blend_vcfs(host_path, donor_path, 0.5, target_depth=1000, seed=1)
        # Only chr1:100 is shared
        assert result.num_markers == 1

    def test_informative_count(self, tmp_path: Path) -> None:
        """Informative markers = those where host and donor differ in alt dose."""
        host_path, donor_path = _make_pair_vcfs(
            tmp_path,
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "0/0:1000,0:1000"),  # vs 1/1 -> informative
                (
                    "chr1",
                    200,
                    "A",
                    "T",
                    "GT:AD:DP",
                    "0/1:500,500:1000",
                ),  # vs 0/1 -> not informative
                ("chr1", 300, "A", "T", "GT:AD:DP", "1/1:0,1000:1000"),  # vs 0/0 -> informative
            ],
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "1/1:0,1000:1000"),
                ("chr1", 200, "A", "T", "GT:AD:DP", "0/1:500,500:1000"),
                ("chr1", 300, "A", "T", "GT:AD:DP", "0/0:1000,0:1000"),
            ],
        )

        result = blend_vcfs(host_path, donor_path, 0.5, target_depth=1000, seed=1)
        assert result.num_markers == 3
        assert result.num_informative == 2

    def test_write_and_reparse(self, tmp_path: Path) -> None:
        """Write a blended VCF and verify it can be re-parsed."""
        host_path, donor_path = _make_pair_vcfs(
            tmp_path,
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "0/0:1000,0:1000"),
                ("chr1", 200, "A", "T", "GT:AD:DP", "0/1:500,500:1000"),
            ],
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "1/1:0,1000:1000"),
                ("chr1", 200, "A", "T", "GT:AD:DP", "0/0:1000,0:1000"),
            ],
        )

        result = blend_vcfs(
            host_path,
            donor_path,
            0.5,
            target_depth=1000,
            sample_name="test_blend",
            seed=99,
        )
        out_path = tmp_path / "blended.vcf"
        write_vcf(result, out_path)

        # Re-parse the written VCF
        header, records = parse_vcf(out_path)
        assert any("#CHROM" in line for line in header)
        assert len(records) == 2
        # Sample name should appear in the header
        chrom_line = [line for line in header if line.startswith("#CHROM")][0]
        assert "test_blend" in chrom_line
        # Each record should have parseable GT
        for rec in records:
            gt = extract_gt(rec)
            assert gt is not None

    def test_invalid_fraction_raises(self, tmp_path: Path) -> None:
        host_path, donor_path = _make_pair_vcfs(
            tmp_path,
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "0/0:1000,0:1000"),
            ],
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "1/1:0,1000:1000"),
            ],
        )
        with pytest.raises(ValueError, match="donor_fraction"):
            blend_vcfs(host_path, donor_path, 1.5)

    def test_reproducible_with_seed(self, tmp_path: Path) -> None:
        """Same seed should produce identical output."""
        host_path, donor_path = _make_pair_vcfs(
            tmp_path,
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "0/0:1000,0:1000"),
                ("chr1", 200, "A", "T", "GT:AD:DP", "0/1:500,500:1000"),
            ],
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "1/1:0,1000:1000"),
                ("chr1", 200, "A", "T", "GT:AD:DP", "0/0:1000,0:1000"),
            ],
        )

        r1 = blend_vcfs(host_path, donor_path, 0.5, target_depth=1000, seed=42)
        r2 = blend_vcfs(host_path, donor_path, 0.5, target_depth=1000, seed=42)
        assert r1.records == r2.records

    def test_ref_only_host_with_variant_donor(self, tmp_path: Path) -> None:
        """When host has ALT='.', donor's ALT allele should be used."""
        host_path, donor_path = _make_pair_vcfs(
            tmp_path,
            [
                ("chr1", 100, "A", ".", "GT:AD:DP", "0/0:1000:1000"),
            ],
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "1/1:0,1000:1000"),
            ],
        )

        result = blend_vcfs(host_path, donor_path, 0.5, target_depth=1000, seed=1)
        assert result.num_markers == 1
        # The ALT column should be T, not .
        alt_col = result.records[0].split("\t")[4]
        assert alt_col == "T"


# ---------------------------------------------------------------------------
# Tests: parse_vcf with real example data
# ---------------------------------------------------------------------------


class TestParseRealVcf:
    """Test parsing against the example VCF shipped with the project."""

    EXAMPLE_VCF = Path(__file__).resolve().parent.parent / "data" / "idt_rhampseq_sid_example.vcf"

    @pytest.mark.skipif(
        not EXAMPLE_VCF.exists(),
        reason="Example VCF not available",
    )
    def test_parse_example_vcf(self) -> None:
        header, records = parse_vcf(self.EXAMPLE_VCF)
        assert len(header) > 0
        assert any(line.startswith("##fileformat") for line in header)
        assert len(records) > 0

        # Every record should have a parseable GT
        for rec in records:
            gt = extract_gt(rec)
            assert gt is not None, f"Failed to parse GT at {rec.locus}"

        # Every record with a variant should have depth
        for rec in records:
            depth = extract_depth(rec)
            assert depth is not None and depth > 0, f"No depth at {rec.locus}"
