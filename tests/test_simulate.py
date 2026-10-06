"""Tests for allomix.simulate — synthetic chimeric VCF generation."""

import filecmp
import math
import os
import random
import statistics
import sys
import tempfile
import textwrap
from pathlib import Path

import numpy as np
import pytest

from allomix.genotype import InformativeMarker
from allomix.simulate import (
    HostAberration,
    alt_dose,
    assign_cnloh_aberrations,
    assign_cnv_aberrations,
    blend_from_genotype_dicts,
    blend_vcfs,
    build_joint_vcf,
    build_joint_vcf_from_genotype_dicts,
    chrx_copy_number,
    cn_weighted_vaf,
    expected_vaf,
    expected_vaf_multi,
    extract_depth,
    extract_gt,
    generate_paired_related_genotypes,
    generate_related_genotypes,
    generate_sex_chrom_genotypes,
    generate_sibling_trio_genotypes,
    gt_from_counts,
    is_chrx_par,
    is_informative,
    parse_text_vcf,
    sample_allele_counts,
    thin_informative_markers,
    write_genotype_vcf,
    write_joint_vcf,
    write_vcf,
)

# The sexchrom fixture generator lives in scripts/ (not installed).
SCRIPTS_DIR = Path(__file__).resolve().parent.parent / "scripts"
sys.path.insert(0, str(SCRIPTS_DIR))

import generate_sexchrom_test_data  # noqa: E402

SEXCHROM_FIXTURE_DIR = Path(__file__).resolve().parent / "test_data" / "sexchrom"

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
    with open(path, "w", encoding="utf-8") as fh:
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

    def test_overdispersion_inflates_variance(self) -> None:
        """Finite rho should widen the VAF spread well beyond binomial."""
        depth, n = 2000, 400
        rng_bin = random.Random(7)
        rng_bb = random.Random(7)
        binom = [sample_allele_counts(0.5, depth, rng_bin)[1] / depth for _ in range(n)]
        betab = [sample_allele_counts(0.5, depth, rng_bb, rho=50.0)[1] / depth for _ in range(n)]
        var_binom = statistics.pvariance(binom)
        var_betab = statistics.pvariance(betab)
        # Beta-binomial var(VAF) = p(1-p)/n * (n+rho)/(rho+1); at p=0.5, n=2000,
        # rho=50 that is ~40x the binomial variance. Mean stays at 0.5.
        assert var_betab > 10 * var_binom
        assert abs(statistics.mean(betab) - 0.5) < 0.02

    def test_infinite_rho_matches_binomial(self) -> None:
        """rho=inf (default) must reproduce the binomial draw exactly."""
        a = sample_allele_counts(0.3, 1500, random.Random(99))
        b = sample_allele_counts(0.3, 1500, random.Random(99), rho=float("inf"))
        assert a == b

    def test_rho_het_only_boundary_stays_binomial(self) -> None:
        """rho_marker_type='het_only' must leave VAF=0 and VAF=1 at the binomial
        error background; the donor-absent allele rate matches e/3 within
        sampling noise, with no extra-binomial inflation.
        """
        n_trials = 400
        depth = 2000
        e = 0.003
        rng0 = random.Random(11)
        rng1 = random.Random(22)
        alt_at_zero = [
            sample_allele_counts(
                0.0,
                depth,
                rng0,
                error_rate=e,
                rho=100.0,
                rho_marker_type="het_only",
            )[1]
            for _ in range(n_trials)
        ]
        ref_at_one = [
            sample_allele_counts(
                1.0,
                depth,
                rng1,
                error_rate=e,
                rho=100.0,
                rho_marker_type="het_only",
            )[0]
            for _ in range(n_trials)
        ]
        # Empirical mean = e/3 per direction; tolerance is ~3x binomial SE.
        expected = e / 3.0
        se = math.sqrt(expected * (1.0 - expected) / depth) / math.sqrt(n_trials)
        mean_zero = sum(alt_at_zero) / (n_trials * depth)
        mean_one = sum(ref_at_one) / (n_trials * depth)
        assert abs(mean_zero - expected) < 5 * se
        assert abs(mean_one - expected) < 5 * se
        # Variance should be at the binomial floor (no overdispersion). The
        # binomial var(VAF) = expected*(1-expected)/depth; allow a generous 2x
        # ceiling (still well below the ~12x beta-binomial inflation rho=100
        # would give if it were applied here).
        bin_var = expected * (1.0 - expected) / depth
        var_zero = statistics.pvariance([y / depth for y in alt_at_zero])
        assert var_zero < 2.5 * bin_var, (
            f"VAF=0 boundary variance {var_zero:.2e} > 2.5x binomial floor "
            f"{bin_var:.2e}; rho was applied at boundary"
        )

    def test_rho_het_only_intermediate_still_overdispersed(self) -> None:
        """rho_marker_type='het_only' must still inflate variance at VAF=0.5
        by the expected beta-binomial factor.
        """
        depth, n = 2000, 400
        rho = 100.0
        rng_bin = random.Random(33)
        rng_bb = random.Random(33)
        binom = [sample_allele_counts(0.5, depth, rng_bin)[1] / depth for _ in range(n)]
        betab = [
            sample_allele_counts(
                0.5,
                depth,
                rng_bb,
                rho=rho,
                rho_marker_type="het_only",
            )[1]
            / depth
            for _ in range(n)
        ]
        var_binom = statistics.pvariance(binom)
        var_betab = statistics.pvariance(betab)
        # Expected inflation factor 1 + (n-1)/(rho+1) at p=0.5, n=2000, rho=100
        # is ~20.8x. Empirical estimates of variance are noisy at n=400 reps
        # so allow a wide band around the expected ratio.
        expected_ratio = 1.0 + (depth - 1) / (rho + 1.0)
        ratio = var_betab / var_binom
        assert 0.4 * expected_ratio < ratio < 2.5 * expected_ratio, (
            f"var-inflation ratio {ratio:.1f} far from expected {expected_ratio:.1f}"
        )

    def test_rho_marker_type_invalid_raises(self) -> None:
        with pytest.raises(ValueError, match="rho_marker_type"):
            sample_allele_counts(
                0.5,
                1000,
                random.Random(1),
                rho=100.0,
                rho_marker_type="nope",
            )


class TestSampleAlleleCountsErrorModel:
    """Verify the 4-state (trinucleotide) error model in sample_allele_counts."""

    def test_error_model_pure_ref_floor(self) -> None:
        """With vaf=0.0 and error_rate=0.03, expected ALT rate is e/3 = 0.01."""
        rng = random.Random(42)
        n_trials = 200
        depth = 10000
        alt_counts = [
            sample_allele_counts(0.0, depth, rng, error_rate=0.03)[1] for _ in range(n_trials)
        ]
        mean_alt_rate = sum(alt_counts) / (n_trials * depth)
        # Expected: e/3 = 0.01.  Allow +/- 0.002 for sampling noise.
        assert abs(mean_alt_rate - 0.01) < 0.002, (
            f"Expected ALT rate ~0.01 (e/3), got {mean_alt_rate:.4f}"
        )

    def test_error_model_pure_alt_floor(self) -> None:
        """With vaf=1.0 and error_rate=0.03, expected REF rate is e/3 = 0.01."""
        rng = random.Random(42)
        n_trials = 200
        depth = 10000
        ref_counts = [
            sample_allele_counts(1.0, depth, rng, error_rate=0.03)[0] for _ in range(n_trials)
        ]
        mean_ref_rate = sum(ref_counts) / (n_trials * depth)
        assert abs(mean_ref_rate - 0.01) < 0.002, (
            f"Expected REF rate ~0.01 (e/3), got {mean_ref_rate:.4f}"
        )

    def test_error_model_matches_estimator(self) -> None:
        """The simulator's effective ALT probability should match the estimator.

        For vaf=0.3 and error_rate=0.01, the 4-state conditional model gives:
            p_alt = 0.3*(1-0.01) + 0.7*0.01/3 = 0.29933...
            p_binomial = p_alt / (1 - 2*0.01/3) = 0.29933 / 0.99333 = 0.30133...
        """
        rng = random.Random(123)
        n_trials = 500
        depth = 10000
        e = 0.01
        p_alt = 0.3 * (1 - e) + 0.7 * e / 3.0
        expected_p = p_alt / (1.0 - 2.0 * e / 3.0)
        alt_counts = [
            sample_allele_counts(0.3, depth, rng, error_rate=e)[1] for _ in range(n_trials)
        ]
        mean_alt_rate = sum(alt_counts) / (n_trials * depth)
        assert abs(mean_alt_rate - expected_p) < 0.001, (
            f"Expected {expected_p:.6f}, got {mean_alt_rate:.6f}"
        )


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
        _, records = parse_text_vcf(tmp_path / "test.vcf")
        assert extract_gt(records[0]) == (0, 1)

    def test_hom_alt(self, tmp_path: Path) -> None:
        _write_test_vcf(
            tmp_path / "test.vcf",
            "SAMPLE",
            [
                ("chr1", 100, "A", "T", "GT:AD:DP:GQ:PL:AF", "1/1:0,1000:1000:99:100,100,0:1.0"),
            ],
        )
        _, records = parse_text_vcf(tmp_path / "test.vcf")
        assert extract_gt(records[0]) == (1, 1)

    def test_nocall(self, tmp_path: Path) -> None:
        _write_test_vcf(
            tmp_path / "test.vcf",
            "SAMPLE",
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "./.:.:."),
            ],
        )
        _, records = parse_text_vcf(tmp_path / "test.vcf")
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
        _, records = parse_text_vcf(tmp_path / "test.vcf")
        assert extract_depth(records[0]) == 1000

    def test_from_ad_fallback(self, tmp_path: Path) -> None:
        _write_test_vcf(
            tmp_path / "test.vcf",
            "SAMPLE",
            [
                ("chr1", 100, "A", "T", "GT:AD", "0/1:600,400"),
            ],
        )
        _, records = parse_text_vcf(tmp_path / "test.vcf")
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
        header, records = parse_text_vcf(out_path)
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
# Tests: parse_text_vcf with the synthetic example VCF
# ---------------------------------------------------------------------------


class TestParseExampleVcf:
    """Test parsing the synthetic single-sample example VCF.

    Uses the committed synthetic fixture (made-up coordinates) rather than a
    real panel VCF, whose marker positions are proprietary.
    """

    EXAMPLE_VCF = Path(__file__).resolve().parent / "test_data" / "single_sample_example.vcf"

    def test_parse_example_vcf(self) -> None:
        header, records = parse_text_vcf(self.EXAMPLE_VCF)
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


class TestPairedRelatedGenotypes:
    """generate_paired_related_genotypes shares a host across relatedness levels."""

    LEVELS = ["unrelated", "cousin", "half-sibling", "sibling"]

    def test_host_and_palt_shared_across_levels(self):
        rng = random.Random(7)
        panels = generate_paired_related_genotypes(200, self.LEVELS, rng)
        ref = panels["unrelated"]
        for rel in self.LEVELS[1:]:
            assert len(panels[rel]) == len(ref)
            for a, b in zip(ref, panels[rel]):
                assert a["host_gt"] == b["host_gt"]
                assert a["p_alt"] == b["p_alt"]

    def test_informative_count_monotone_non_increasing(self):
        # Averaged over replicates, informative markers should fall as relatedness
        # rises. The paired design makes this hold per replicate too.
        n_markers = 300
        for rep in range(5):
            rng = random.Random(100 + rep)
            panels = generate_paired_related_genotypes(n_markers, self.LEVELS, rng)
            counts = [sum(m["informative"] for m in panels[rel]) for rel in self.LEVELS]
            for lo, hi in zip(counts, counts[1:]):
                assert hi <= lo, f"non-monotone informative counts {counts} (rep {rep})"

    def test_unknown_level_raises(self):
        rng = random.Random(1)
        with pytest.raises(ValueError):
            generate_paired_related_genotypes(10, ["unrelated", "bogus"], rng)


class TestBlendVcfLocusDropout:
    """num_markers should match len(records) when dropout occurs."""

    def test_num_markers_matches_records_with_dropout(self):
        rng = random.Random(42)
        geno = generate_related_genotypes(50, "unrelated", rng)

        with tempfile.TemporaryDirectory() as tmpdir:
            host_path = os.path.join(tmpdir, "host.vcf")
            donor_path = os.path.join(tmpdir, "donor.vcf")
            write_genotype_vcf(geno, host_path, "HOST", key="host_gt")
            write_genotype_vcf(geno, donor_path, "DONOR", key="donor_gt")

            result = blend_vcfs(
                host_path,
                donor_path,
                donor_fraction=0.20,
                target_depth=1000,
                seed=42,
                locus_dropout_rate=0.20,
            )
            assert result.num_markers == len(result.records), (
                f"num_markers={result.num_markers} but len(records)={len(result.records)}"
            )


# ---------------------------------------------------------------------------
# Tests: build_joint_vcf
# ---------------------------------------------------------------------------


class TestBuildJointVcf:
    """Test the multi-sample joint VCF builder."""

    def test_basic_structure(self, tmp_path: Path) -> None:
        """Joint VCF should have correct header and sample columns."""
        host_path, donor_path = _make_pair_vcfs(
            tmp_path,
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "0/0:1000,0:1000"),
                ("chr1", 200, "A", "T", "GT:AD:DP", "1/1:0,1000:1000"),
            ],
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "1/1:0,1000:1000"),
                ("chr1", 200, "A", "T", "GT:AD:DP", "0/0:1000,0:1000"),
            ],
        )
        result = build_joint_vcf(
            host_path=str(host_path),
            donor_paths=[str(donor_path)],
            admix_fractions=[0.0, 0.5],
            admix_sample_names=["TP1", "TP2"],
            target_depth=1000,
            seed=42,
        )
        assert result.num_markers == 2
        assert result.sample_names == ["HOST", "DONOR", "TP1", "TP2"]
        # Check header has all sample names
        chrom_line = [line for line in result.header if line.startswith("#CHROM")][0]
        assert "HOST" in chrom_line
        assert "DONOR" in chrom_line
        assert "TP1" in chrom_line
        assert "TP2" in chrom_line

    def test_write_and_parse(self, tmp_path: Path) -> None:
        """Write joint VCF and verify it can be parsed back."""
        host_path, donor_path = _make_pair_vcfs(
            tmp_path,
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "0/0:1000,0:1000"),
                ("chr1", 200, "A", "T", "GT:AD:DP", "1/1:0,1000:1000"),
            ],
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "1/1:0,1000:1000"),
                ("chr1", 200, "A", "T", "GT:AD:DP", "0/0:1000,0:1000"),
            ],
        )
        result = build_joint_vcf(
            host_path=str(host_path),
            donor_paths=[str(donor_path)],
            admix_fractions=[0.10],
            admix_sample_names=["ADMIX"],
            target_depth=2000,
            seed=42,
        )
        out = tmp_path / "joint.vcf"
        write_joint_vcf(result, out)

        # Re-parse with simulate.parse_text_vcf (text-based parser)
        header, records = parse_text_vcf(out)
        assert len(records) == 2
        chrom_line = [line for line in header if line.startswith("#CHROM")][0]
        assert "HOST" in chrom_line
        assert "DONOR" in chrom_line
        assert "ADMIX" in chrom_line

    def test_informative_count(self, tmp_path: Path) -> None:
        """Informative markers should be counted correctly."""
        host_path, donor_path = _make_pair_vcfs(
            tmp_path,
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "0/0:1000,0:1000"),  # informative
                ("chr1", 200, "A", "T", "GT:AD:DP", "0/1:500,500:1000"),  # not informative
                ("chr1", 300, "A", "T", "GT:AD:DP", "1/1:0,1000:1000"),  # informative
            ],
            [
                ("chr1", 100, "A", "T", "GT:AD:DP", "1/1:0,1000:1000"),
                ("chr1", 200, "A", "T", "GT:AD:DP", "0/1:500,500:1000"),
                ("chr1", 300, "A", "T", "GT:AD:DP", "0/0:1000,0:1000"),
            ],
        )
        result = build_joint_vcf(
            host_path=str(host_path),
            donor_paths=[str(donor_path)],
            admix_fractions=[0.10],
            admix_sample_names=["ADMIX"],
            target_depth=1000,
            seed=42,
        )
        assert result.num_markers == 3
        assert result.num_informative == 2

    def test_mismatched_lengths_raises(self, tmp_path: Path) -> None:
        """Mismatched fractions and names should raise ValueError."""
        host_path, donor_path = _make_pair_vcfs(
            tmp_path,
            [("chr1", 100, "A", "T", "GT:AD:DP", "0/0:1000,0:1000")],
            [("chr1", 100, "A", "T", "GT:AD:DP", "1/1:0,1000:1000")],
        )
        with pytest.raises(ValueError, match="admix_fractions length"):
            build_joint_vcf(
                host_path=str(host_path),
                donor_paths=[str(donor_path)],
                admix_fractions=[0.10, 0.50],
                admix_sample_names=["ONLY_ONE"],
                target_depth=1000,
                seed=42,
            )


# ---------------------------------------------------------------------------
# Tests: host CN-LoH aberrations
# ---------------------------------------------------------------------------


class TestCnWeightedVaf:
    def test_no_aberration_matches_diploid(self) -> None:
        """With no aberration the CN-weighted VAF equals the diploid model."""
        for host_gt in [(0, 0), (0, 1), (1, 1)]:
            for donor_gt in [(0, 0), (0, 1), (1, 1)]:
                for f in [0.0, 0.1, 0.5, 0.9]:
                    assert cn_weighted_vaf(host_gt, [donor_gt], [f], None) == pytest.approx(
                        expected_vaf(host_gt, donor_gt, f)
                    )

    def test_cnloh_pure_host_retains_alt(self) -> None:
        """Pure-host CN-LoH retaining ALT drives a germline het to VAF 1.0."""
        aberr = HostAberration(cn=2, alt_copies=2, clonal_fraction=1.0)
        assert cn_weighted_vaf((0, 1), [(0, 0)], [0.0], aberr) == pytest.approx(1.0)

    def test_cnloh_pure_host_retains_ref(self) -> None:
        """Pure-host CN-LoH retaining REF drives a germline het to VAF 0.0."""
        aberr = HostAberration(cn=2, alt_copies=0, clonal_fraction=1.0)
        assert cn_weighted_vaf((0, 1), [(0, 0)], [0.0], aberr) == pytest.approx(0.0)

    def test_cnloh_shifts_vaf_vs_diploid(self) -> None:
        """A CN-LoH het with a hom-ref donor shifts VAF above the diploid value."""
        host_gt, donor_gt, f = (0, 1), (0, 0), 0.2
        aberr = HostAberration(cn=2, alt_copies=2, clonal_fraction=1.0)
        diploid = expected_vaf(host_gt, donor_gt, f)  # 0.5 * (1 - f) = 0.4
        cnloh = cn_weighted_vaf(host_gt, [donor_gt], [f], aberr)  # (1 - f) = 0.8
        assert diploid == pytest.approx(0.4)
        assert cnloh == pytest.approx(0.8)

    def test_cnloh_partial_clone(self) -> None:
        """A 50% clone sits halfway between germline and full-clone VAF."""
        host_gt, donor_gt, f = (0, 1), (0, 0), 0.0
        aberr = HostAberration(cn=2, alt_copies=2, clonal_fraction=0.5)
        # half normal het (0.5) + half clone hom-alt (1.0) = 0.75
        assert cn_weighted_vaf(host_gt, [donor_gt], [f], aberr) == pytest.approx(0.75)

    def test_deletion_changes_denominator(self) -> None:
        """A host deletion (CN1) raises the apparent donor VAF at that locus."""
        # Host het loses the ALT homolog (CN1, 0 ALT copies), donor is hom-alt.
        host_gt, donor_gt, f = (0, 1), (1, 1), 0.1
        aberr = HostAberration(cn=1, alt_copies=0, clonal_fraction=1.0)
        # num = f*2 = 0.2 ; den = (1-f)*1 + f*2 = 0.9 + 0.2 = 1.1
        assert cn_weighted_vaf(host_gt, [donor_gt], [f], aberr) == pytest.approx(0.2 / 1.1)


class TestAssignCnlohAberrations:
    def test_only_hets_affected(self) -> None:
        """Homozygous markers never receive an aberration."""
        markers = [
            {"host_gt": (0, 0)},
            {"host_gt": (1, 1)},
            {"host_gt": (0, 1)},
        ]
        rng = random.Random(0)
        aberrs = assign_cnloh_aberrations(
            markers, fraction_affected=1.0, clonal_fraction=1.0, rng=rng
        )
        assert aberrs[0] is None
        assert aberrs[1] is None
        assert aberrs[2] is not None
        assert aberrs[2].cn == 2
        assert aberrs[2].alt_copies in (0, 2)

    def test_fraction_zero_assigns_none(self) -> None:
        markers = [{"host_gt": (0, 1)} for _ in range(20)]
        rng = random.Random(0)
        aberrs = assign_cnloh_aberrations(markers, 0.0, 1.0, rng)
        assert all(a is None for a in aberrs)

    def test_fraction_roughly_matches(self) -> None:
        markers = [{"host_gt": (0, 1)} for _ in range(2000)]
        rng = random.Random(1)
        aberrs = assign_cnloh_aberrations(markers, 0.3, 1.0, rng)
        frac = sum(1 for a in aberrs if a is not None) / len(aberrs)
        assert 0.25 < frac < 0.35


class TestBlendVcfsWithAberrations:
    def test_aberration_length_validated(self, tmp_path: Path) -> None:
        host_path, donor_path = _make_pair_vcfs(
            tmp_path,
            [("chr1", 200, "A", "T", "GT:AD:DP", "0/1:500,500:1000")],
            [("chr1", 200, "A", "T", "GT:AD:DP", "0/0:1000,0:1000")],
        )
        with pytest.raises(ValueError, match="host_aberrations length"):
            blend_vcfs(
                host_path,
                donor_path,
                0.1,
                target_depth=2000,
                seed=1,
                host_aberrations=[None, None],
            )

    def test_cnloh_raises_observed_vaf(self, tmp_path: Path) -> None:
        """CN-LoH at a host het (hom-ref donor) lifts the observed ALT VAF."""
        # Host het, donor hom-ref, small donor fraction. Diploid expects ~0.45;
        # CN-LoH retaining ALT in a pure clone expects ~0.9.
        host_path, donor_path = _make_pair_vcfs(
            tmp_path,
            [("chr1", 200, "A", "T", "GT:AD:DP", "0/1:500,500:1000")],
            [("chr1", 200, "A", "T", "GT:AD:DP", "0/0:1000,0:1000")],
        )
        aberr = [HostAberration(cn=2, alt_copies=2, clonal_fraction=1.0)]
        result = blend_vcfs(
            host_path,
            donor_path,
            0.1,
            target_depth=20000,
            seed=7,
            error_rate=0.0,
            host_aberrations=aberr,
        )
        sample = result.records[0].split("\t")[9]
        ad_ref, ad_alt = (int(x) for x in sample.split(":")[1].split(","))
        observed = ad_alt / (ad_ref + ad_alt)
        assert observed > 0.85


class TestAssignCnvAberrations:
    def test_invalid_kind(self) -> None:
        rng = random.Random(0)
        with pytest.raises(ValueError, match="kind must be"):
            assign_cnv_aberrations([{"host_gt": (0, 1)}], 1.0, 1.0, rng, kind="bogus")

    def test_deletion_cn1_all_genotypes(self) -> None:
        """Deletions affect hom markers too (DNA-mass effect), unlike CN-LoH."""
        markers = [{"host_gt": (0, 0)}, {"host_gt": (1, 1)}, {"host_gt": (0, 1)}]
        rng = random.Random(0)
        aberrs = assign_cnv_aberrations(markers, 1.0, 1.0, rng, kind="deletion")
        assert all(a is not None for a in aberrs)
        assert all(a.cn == 1 for a in aberrs)
        assert aberrs[0].alt_copies == 0  # hom-ref -> retains REF
        assert aberrs[1].alt_copies == 1  # hom-alt -> retains ALT
        assert aberrs[2].alt_copies in (0, 1)  # het -> retains one homolog

    def test_gain_cn3_all_genotypes(self) -> None:
        markers = [{"host_gt": (0, 0)}, {"host_gt": (1, 1)}, {"host_gt": (0, 1)}]
        rng = random.Random(0)
        aberrs = assign_cnv_aberrations(markers, 1.0, 1.0, rng, kind="gain")
        assert all(a is not None and a.cn == 3 for a in aberrs)
        assert aberrs[0].alt_copies == 0  # 0,0 + dup 0
        assert aberrs[1].alt_copies == 3  # 1,1 + dup 1
        assert aberrs[2].alt_copies in (1, 2)  # 0,1 + dup of one homolog

    def test_cnloh_only_hets(self) -> None:
        markers = [{"host_gt": (0, 0)}, {"host_gt": (0, 1)}]
        rng = random.Random(0)
        aberrs = assign_cnv_aberrations(markers, 1.0, 1.0, rng, kind="cnloh")
        assert aberrs[0] is None
        assert aberrs[1] is not None and aberrs[1].cn == 2

    def test_cnloh_wrapper_matches(self) -> None:
        markers = [{"host_gt": (0, 1)} for _ in range(50)]
        a = assign_cnv_aberrations(markers, 0.5, 0.8, random.Random(3), kind="cnloh")
        b = assign_cnloh_aberrations(markers, 0.5, 0.8, random.Random(3))
        assert [None if x is None else (x.cn, x.alt_copies) for x in a] == [
            None if y is None else (y.cn, y.alt_copies) for y in b
        ]

    def test_deletion_raises_donor_vaf(self) -> None:
        """A host deletion losing the host allele raises apparent donor VAF."""
        # Host hom-ref, donor hom-alt, small donor fraction. Diploid VAF = f.
        # Host deletion (CN1) halves host DNA, so donor VAF rises above f.
        host_gt, donor_gt, f = (0, 0), (1, 1), 0.1
        aberr = HostAberration(cn=1, alt_copies=0, clonal_fraction=1.0)
        diploid = expected_vaf(host_gt, donor_gt, f)
        deleted = cn_weighted_vaf(host_gt, [donor_gt], [f], aberr)
        assert diploid == pytest.approx(0.1)
        # num = f*2 = 0.2 ; den = (1-f)*1 + f*2 = 1.1 -> 0.1818
        assert deleted == pytest.approx(0.2 / 1.1)
        assert deleted > diploid

    def test_gain_lowers_donor_vaf(self) -> None:
        """A host gain adds host DNA, diluting apparent donor VAF."""
        host_gt, donor_gt, f = (0, 0), (1, 1), 0.1
        aberr = HostAberration(cn=3, alt_copies=0, clonal_fraction=1.0)
        diploid = expected_vaf(host_gt, donor_gt, f)
        gained = cn_weighted_vaf(host_gt, [donor_gt], [f], aberr)
        # num = f*2 = 0.2 ; den = (1-f)*3 + f*2 = 2.9 -> 0.069
        assert gained == pytest.approx(0.2 / 2.9)
        assert gained < diploid


# ---------------------------------------------------------------------------
# Depth thinning of informative markers (thin_informative_markers)
# ---------------------------------------------------------------------------


def _make_markers(
    n: int, mean_depth: int, depth_cv: float, vaf: float, seed: int
) -> list[InformativeMarker]:
    """Build a synthetic informative-marker list with a real-ish depth spread.

    Per-marker depths are drawn log-normal (mean ``mean_depth``, CV ``depth_cv``)
    and ALT counts binomial at ``vaf``, so the list looks like a parsed admix
    sample with locus-to-locus depth variation to thin.
    """
    rng = np.random.default_rng(seed)
    sigma2 = math.log(1 + depth_cv**2)
    mu = math.log(mean_depth) - sigma2 / 2
    sigma = math.sqrt(sigma2)
    markers = []
    for i in range(n):
        dp = max(1, int(round(math.exp(rng.normal(mu, sigma)))))
        alt = int(rng.binomial(dp, vaf))
        markers.append(
            InformativeMarker(
                chrom="chr1",
                pos=1000 + i,
                ref="A",
                alt="G",
                host_gt=(0, 0),
                donor_gts=[(0, 1)],
                marker_type=1,
                admix_ad_ref=dp - alt,
                admix_ad_alt=alt,
                admix_dp=dp,
            )
        )
    return markers


class TestThinInformativeMarkers:
    """Tests for global-rate binomial thinning of admix counts."""

    def test_rate_one_is_passthrough(self) -> None:
        """rate=1.0 returns markers with unchanged counts."""
        markers = _make_markers(50, 2000, 0.43, 0.25, seed=1)
        rng = np.random.default_rng(0)
        out = thin_informative_markers(markers, 1.0, rng)
        assert len(out) == len(markers)
        for a, b in zip(markers, out):
            assert (a.admix_ad_ref, a.admix_ad_alt, a.admix_dp) == (
                b.admix_ad_ref,
                b.admix_ad_alt,
                b.admix_dp,
            )

    def test_invalid_rate_raises(self) -> None:
        markers = _make_markers(5, 2000, 0.0, 0.25, seed=1)
        rng = np.random.default_rng(0)
        for bad in (0.0, -0.1, 1.5):
            with pytest.raises(ValueError):
                thin_informative_markers(markers, bad, rng)

    def test_dp_is_ref_plus_alt(self) -> None:
        """Thinned admix_dp equals the thinned ref + alt counts at every marker."""
        markers = _make_markers(200, 1500, 0.43, 0.3, seed=2)
        rng = np.random.default_rng(7)
        out = thin_informative_markers(markers, 0.3, rng)
        for m in out:
            assert m.admix_dp == m.admix_ad_ref + m.admix_ad_alt

    def test_mean_depth_scales_with_rate(self) -> None:
        """Mean thinned depth ~ rate * original mean over many markers."""
        markers = _make_markers(2000, 2000, 0.43, 0.25, seed=3)
        orig_mean = float(np.mean([m.admix_dp for m in markers]))
        rng = np.random.default_rng(11)
        rate = 0.25
        out = thin_informative_markers(markers, rate, rng)
        thin_mean = float(np.mean([m.admix_dp for m in out]))
        assert thin_mean == pytest.approx(rate * orig_mean, rel=0.03)

    def test_allele_ratio_preserved_in_expectation(self) -> None:
        """Pooled ALT fraction is preserved over many seeds (thinning is unbiased)."""
        markers = _make_markers(100, 2000, 0.43, 0.3, seed=4)
        orig_alt = sum(m.admix_ad_alt for m in markers)
        orig_tot = sum(m.admix_dp for m in markers)
        orig_frac = orig_alt / orig_tot
        rate = 0.2
        fracs = []
        for s in range(8):
            rng = np.random.default_rng(100 + s)
            out = thin_informative_markers(markers, rate, rng)
            alt = sum(m.admix_ad_alt for m in out)
            tot = sum(m.admix_dp for m in out)
            fracs.append(alt / tot)
        assert float(np.mean(fracs)) == pytest.approx(orig_frac, abs=0.01)

    def test_depth_cv_preserved(self) -> None:
        """Thinning at a global rate preserves the locus-to-locus depth CV.

        A per-marker normalisation would flatten this; the global rate must not.
        """
        markers = _make_markers(3000, 2000, 0.43, 0.25, seed=5)
        depths = np.array([m.admix_dp for m in markers], dtype=float)
        orig_cv = depths.std() / depths.mean()
        rng = np.random.default_rng(13)
        out = thin_informative_markers(markers, 0.25, rng)
        td = np.array([m.admix_dp for m in out], dtype=float)
        thin_cv = td.std() / td.mean()
        # Binomial thinning adds a little sampling variance at low counts but the
        # real depth spread dominates, so the CV is preserved to within ~15%.
        assert thin_cv == pytest.approx(orig_cv, rel=0.15)

    def test_deterministic_with_seeded_rng(self) -> None:
        """Same seeded generator gives identical thinned counts."""
        markers = _make_markers(100, 2000, 0.43, 0.25, seed=6)
        out_a = thin_informative_markers(markers, 0.4, np.random.default_rng(42))
        out_b = thin_informative_markers(markers, 0.4, np.random.default_rng(42))
        for a, b in zip(out_a, out_b):
            assert (a.admix_ad_ref, a.admix_ad_alt) == (b.admix_ad_ref, b.admix_ad_alt)

    def test_input_not_mutated(self) -> None:
        """The input markers are left unchanged (fresh copies returned)."""
        markers = _make_markers(50, 2000, 0.43, 0.25, seed=8)
        snapshot = [(m.admix_ad_ref, m.admix_ad_alt, m.admix_dp) for m in markers]
        thin_informative_markers(markers, 0.3, np.random.default_rng(0))
        after = [(m.admix_ad_ref, m.admix_ad_alt, m.admix_dp) for m in markers]
        assert snapshot == after


# ---------------------------------------------------------------------------
# Tests: sex-chromosome genotypes and copy-number-weighted blending
# ---------------------------------------------------------------------------

FRACTIONS = [0.0, 0.01, 0.05, 0.2, 0.5, 0.8, 1.0]
GTS = [(0, 0), (0, 1), (1, 1)]


class TestChrxCopyNumber:
    def test_autosome_and_unknown_sex_are_diploid(self) -> None:
        assert chrx_copy_number("chr1", 5, "M") == 2
        assert chrx_copy_number("chrX", 50_000_000, None) == 2

    def test_male_nonpar_is_hemizygous(self) -> None:
        assert chrx_copy_number("chrX", 50_000_000, "M") == 1
        assert chrx_copy_number("X", 50_000_000, "M") == 1
        assert chrx_copy_number("chrX", 50_000_000, "F") == 2

    def test_par_union_boundaries(self) -> None:
        # Fixture PAR positions and the union mask edges (GRCh38 PAR1 end,
        # GRCh37 PAR2 start).
        for pos in (1_000_000, 155_800_000, 10_001, 2_781_479, 154_931_044, 156_030_895):
            assert is_chrx_par(pos), pos
            assert chrx_copy_number("chrX", pos, "M") == 2
        for pos in (10_000, 2_781_480, 154_931_043, 10_000_000, 150_000_000):
            assert not is_chrx_par(pos), pos
            assert chrx_copy_number("chrX", pos, "M") == 1

    def test_alt_contig_is_diploid(self) -> None:
        assert chrx_copy_number("chrX_KI270880v1_alt", 100_000, "M") == 2

    def test_invalid_sex_raises(self) -> None:
        with pytest.raises(ValueError):
            chrx_copy_number("chrX", 50_000_000, "male")


class TestGenerateSexChromGenotypes:
    def test_male_never_het_without_spurious_rate(self) -> None:
        for host_sex, donor_sex in (("M", "M"), ("M", "F"), ("F", "M")):
            rng = random.Random(3)
            markers = generate_sex_chrom_genotypes(
                500, host_sex, donor_sex, rng, male_chrx_spurious_het_rate=0.0
            )
            for who, sex in (("host", host_sex), ("donor", donor_sex)):
                if sex != "M":
                    continue
                for m in markers:
                    assert m[f"{who}_gt"][0] == m[f"{who}_gt"][1], m
                    assert m[f"{who}_gt"] == m[f"{who}_true_gt"]
                    assert m[f"{who}_gt"] in ((0, 0), (1, 1))
                    assert m[f"{who}_cn"] == 1

    def test_female_is_ordinary_diploid(self) -> None:
        rng = random.Random(4)
        markers = generate_sex_chrom_genotypes(300, "F", "F", rng)
        assert all(m["host_cn"] == 2 and m["donor_cn"] == 2 for m in markers)
        assert all(m["host_gt"] == m["host_true_gt"] for m in markers)
        assert any(alt_dose(m["host_gt"]) == 1 for m in markers)
        assert all(m["chrom"] == "chrX" and not m["par"] for m in markers)

    def test_spurious_rate_one_corrupts_every_male_call(self) -> None:
        rng = random.Random(5)
        markers = generate_sex_chrom_genotypes(50, "M", "F", rng, male_chrx_spurious_het_rate=1.0)
        for m in markers:
            assert m["host_gt"] == (0, 1)
            assert m["host_true_gt"] in ((0, 0), (1, 1))
            assert m["donor_gt"] == m["donor_true_gt"]  # female untouched

    def test_spurious_rate_roughly_matches(self) -> None:
        rng = random.Random(6)
        markers = generate_sex_chrom_genotypes(3000, "M", "M", rng)  # default 0.04
        n_het = sum(1 for m in markers if m["host_gt"] == (0, 1))
        assert 0.025 < n_het / len(markers) < 0.055

    def test_true_genotypes_independent_of_spurious_rate(self) -> None:
        a = generate_sex_chrom_genotypes(
            200, "M", "M", random.Random(7), male_chrx_spurious_het_rate=0.0
        )
        b = generate_sex_chrom_genotypes(
            200, "M", "M", random.Random(7), male_chrx_spurious_het_rate=0.5
        )
        assert [m["host_true_gt"] for m in a] == [m["host_true_gt"] for m in b]
        assert [m["donor_true_gt"] for m in a] == [m["donor_true_gt"] for m in b]

    def test_par_markers_diploid_and_in_mask(self) -> None:
        rng = random.Random(8)
        markers = generate_sex_chrom_genotypes(5, "M", "M", rng, n_par=2)
        nonpar, par = markers[:5], markers[5:]
        assert all(not m["par"] and not is_chrx_par(m["pos"]) for m in nonpar)
        assert len(par) == 2
        assert all(m["par"] and is_chrx_par(m["pos"]) for m in par)
        assert all(m["host_cn"] == 2 and m["donor_cn"] == 2 for m in par)
        assert len({m["pos"] for m in markers}) == len(markers)

    def test_informative_uses_true_allele_fraction(self) -> None:
        rng = random.Random(9)
        markers = generate_sex_chrom_genotypes(300, "M", "F", rng)
        for m in markers:
            expected = alt_dose(m["host_true_gt"]) != alt_dose(m["donor_true_gt"])
            assert m["informative"] == expected

    def test_invalid_sex_raises(self) -> None:
        with pytest.raises(ValueError):
            generate_sex_chrom_genotypes(10, "male", "F", random.Random(0))


class TestCnWeightedVafSex:
    def test_male_minor_at_female_hom_site(self) -> None:
        """Male donor hom-alt at a female host hom-ref site: f / (2 - f)."""
        for f in FRACTIONS:
            vaf = cn_weighted_vaf((0, 0), [(1, 1)], [f], None, host_cn=2, donor_cns=[1])
            assert vaf == pytest.approx(f / (2.0 - f))

    def test_female_hom_alt_minor_at_male_hom_ref(self) -> None:
        """Female donor hom-alt at a male host hom-ref site: 2f / (1 + f)."""
        for f in FRACTIONS:
            vaf = cn_weighted_vaf((0, 0), [(1, 1)], [f], None, host_cn=1, donor_cns=[2])
            assert vaf == pytest.approx(2.0 * f / (1.0 + f))

    def test_ff_reduces_to_autosomal(self) -> None:
        for h in GTS:
            for d in GTS:
                for f in FRACTIONS:
                    assert cn_weighted_vaf(h, [d], [f], None, 2, [2]) == expected_vaf(h, d, f)

    def test_mm_hom_hom_ploidy_cancels(self) -> None:
        """Male/male hom contrasts: copy number cancels, the diploid math is exact."""
        for h in ((0, 0), (1, 1)):
            for d in ((0, 0), (1, 1)):
                for f in FRACTIONS:
                    assert cn_weighted_vaf(h, [d], [f], None, 1, [1]) == pytest.approx(
                        expected_vaf(h, d, f)
                    )

    def test_autosomal_two_donor_unchanged(self) -> None:
        """Copy number 2 everywhere returns expected_vaf_multi exactly (not approx)."""
        for h in GTS:
            for d1 in GTS:
                for d2 in GTS:
                    for f1, f2 in ((0.0, 0.0), (0.1, 0.05), (0.3, 0.3), (0.5, 0.5)):
                        assert cn_weighted_vaf(h, [d1, d2], [f1, f2], None, 2, [2, 2]) == (
                            expected_vaf_multi(h, [d1, d2], [f1, f2])
                        )

    def test_aberration_with_default_cn_matches_previous_formula(self) -> None:
        host_gt, donor_gt, f = (0, 1), (1, 1), 0.1
        aberr = HostAberration(cn=1, alt_copies=0, clonal_fraction=0.6)
        c = 0.6
        num = (1 - f) * (1 - c) * 1 + (1 - f) * c * 0 + f * 2
        den = (1 - f) * (1 - c) * 2 + (1 - f) * c * 1 + f * 2
        assert cn_weighted_vaf(host_gt, [donor_gt], [f], aberr) == pytest.approx(num / den)


def _write_pair_vcfs(markers: list[dict], tmp_path: Path) -> tuple[Path, Path]:
    host_path = tmp_path / "host.vcf"
    donor_path = tmp_path / "donor.vcf"
    write_genotype_vcf(markers, host_path, "HOST", key="host_gt")
    write_genotype_vcf(markers, donor_path, "DONOR", key="donor_gt")
    return host_path, donor_path


def _alt_frac(record_line: str) -> float:
    fields = record_line.split("\t")
    ad = dict(zip(fields[8].split(":"), fields[9].split(":")))["AD"]
    ref, alt = (int(x) for x in ad.split(","))
    return alt / (ref + alt)


class TestBlendingWithSex:
    def test_blend_from_genotype_dicts_autosomal_identical_with_explicit_cn(self) -> None:
        """Trio dicts with explicit cn=2 / true-gt keys blend byte-identically."""
        markers = generate_sibling_trio_genotypes(80, random.Random(11))
        tagged = [
            {
                **m,
                "host_cn": 2,
                "donor1_cn": 2,
                "donor2_cn": 2,
                "host_true_gt": m["host_gt"],
                "donor1_true_gt": m["donor1_gt"],
                "donor2_true_gt": m["donor2_gt"],
            }
            for m in markers
        ]
        plain = blend_from_genotype_dicts(markers, [0.1, 0.2], seed=5, depth_cv=0.43)
        with_cn = blend_from_genotype_dicts(tagged, [0.1, 0.2], seed=5, depth_cv=0.43)
        assert plain.records == with_cn.records
        assert plain.header == with_cn.header

    def test_blend_from_genotype_dicts_accepts_donor_gt_key(self) -> None:
        markers = generate_related_genotypes(40, "unrelated", random.Random(12))
        result = blend_from_genotype_dicts(markers, [0.3], seed=1)
        assert result.num_markers == 40
        assert result.num_informative == sum(m["informative"] for m in markers)

    def test_blend_vcfs_autosomal_identical_with_sexes(self, tmp_path: Path) -> None:
        """On autosomes, declaring sexes changes nothing (cn 2 everywhere)."""
        markers = generate_related_genotypes(60, "unrelated", random.Random(13))
        host_path, donor_path = _write_pair_vcfs(markers, tmp_path)
        kwargs = dict(donor_fraction=0.2, target_depth=1000, seed=3, depth_cv=0.43)
        base = blend_vcfs(host_path, donor_path, **kwargs)
        for sexes in (("F", "F"), ("M", "M"), ("M", "F")):
            same = blend_vcfs(
                host_path, donor_path, host_sex=sexes[0], donor_sex=sexes[1], **kwargs
            )
            assert same.records == base.records

    def test_blend_vcfs_male_host_female_donor_chrx(self, tmp_path: Path) -> None:
        """A female hom-alt donor at a male hom-ref chrX site reads 2f/(1+f), PAR stays diploid."""
        markers = generate_sex_chrom_genotypes(
            60, "M", "F", random.Random(14), n_par=2, male_chrx_spurious_het_rate=0.0
        )
        host_path, donor_path = _write_pair_vcfs(markers, tmp_path)
        f = 0.2
        result = blend_vcfs(
            host_path,
            donor_path,
            donor_fraction=f,
            target_depth=200_000,
            seed=2,
            error_rate=0.0,
            host_sex="M",
            donor_sex="F",
        )
        by_pos = {int(line.split("\t")[1]): line for line in result.records}
        checked_nonpar = checked_par = 0
        for m in markers:
            if m["host_true_gt"] != (0, 0) or m["donor_true_gt"] != (1, 1):
                continue
            observed = _alt_frac(by_pos[m["pos"]])
            if m["par"]:
                assert observed == pytest.approx(f, abs=0.01)
                checked_par += 1
            else:
                assert observed == pytest.approx(2 * f / (1 + f), abs=0.01)
                checked_nonpar += 1
        assert checked_nonpar > 0

    def test_blend_vcfs_male_het_at_nonpar_raises(self, tmp_path: Path) -> None:
        markers = generate_sex_chrom_genotypes(
            5, "M", "M", random.Random(15), male_chrx_spurious_het_rate=0.0
        )
        markers[2]["host_gt"] = (0, 1)
        host_path, donor_path = _write_pair_vcfs(markers, tmp_path)
        with pytest.raises(ValueError, match="hemizygous"):
            blend_vcfs(host_path, donor_path, donor_fraction=0.2, host_sex="M", donor_sex="M")
        # Without sexes the same VCFs blend as diploid (legacy behaviour).
        blend_vcfs(host_path, donor_path, donor_fraction=0.2)

    def test_build_joint_vcf_sex_aware_matches_formula(self, tmp_path: Path) -> None:
        markers = generate_sex_chrom_genotypes(
            40, "F", "M", random.Random(16), male_chrx_spurious_het_rate=0.0
        )
        host_path, donor_path = _write_pair_vcfs(markers, tmp_path)
        f = 0.1
        result = build_joint_vcf(
            host_path,
            [donor_path],
            [f],
            ["ADMIX"],
            target_depth=200_000,
            seed=4,
            error_rate=0.0,
            host_sex="F",
            donor_sexes=["M"],
        )
        checked = 0
        for m, line in zip(markers, result.records):
            if m["host_true_gt"] == (0, 0) and m["donor_true_gt"] == (1, 1):
                fields = line.split("\t")
                ref, alt = (int(x) for x in fields[-1].split(":")[1].split(","))
                assert alt / (ref + alt) == pytest.approx(f / (2 - f), abs=0.01)
                checked += 1
        assert checked > 0

    def test_build_joint_vcf_donor_sexes_length_validated(self, tmp_path: Path) -> None:
        markers = generate_related_genotypes(5, "unrelated", random.Random(17))
        host_path, donor_path = _write_pair_vcfs(markers, tmp_path)
        with pytest.raises(ValueError):
            build_joint_vcf(host_path, [donor_path], [0.1], ["A"], donor_sexes=["M", "F"])

    def test_joint_from_dicts_spurious_het_written_but_blended_true(self) -> None:
        """Host column shows the erroneous 0/1; admix reads come from the true allele."""
        markers = generate_sex_chrom_genotypes(
            6, "M", "M", random.Random(18), male_chrx_spurious_het_rate=0.0
        )
        target = markers[1]
        target["host_true_gt"] = (1, 1)
        target["host_gt"] = (0, 1)
        target["donor_true_gt"] = target["donor_gt"] = (0, 0)
        result = build_joint_vcf_from_genotype_dicts(
            markers, [0.0], ["PURE_HOST"], target_depth=100_000, seed=1, error_rate=0.0
        )
        fields = result.records[1].split("\t")
        assert result.sample_names == ["HOST", "DONOR", "PURE_HOST"]
        assert fields[9].startswith("0/1:")
        host_ref, host_alt = (int(x) for x in fields[9].split(":")[1].split(","))
        assert host_alt / (host_ref + host_alt) == pytest.approx(0.5, abs=0.02)
        admix_ref, admix_alt = (int(x) for x in fields[11].split(":")[1].split(","))
        assert admix_alt == 100_000 and admix_ref == 0

    def test_joint_from_dicts_sorts_written_gt(self) -> None:
        markers = generate_sex_chrom_genotypes(300, "F", "F", random.Random(19))
        result = build_joint_vcf_from_genotype_dicts(markers, [0.5], ["A"], target_depth=50, seed=0)
        assert any(m["host_gt"] == (1, 0) for m in markers)
        assert not any("\t1/0:" in line for line in result.records)


class TestSexchromFixtureGenerator:
    def test_deterministic_and_matches_committed(self, tmp_path: Path) -> None:
        """Two runs at the default seed agree with each other and with tests/test_data/sexchrom."""
        run_a = tmp_path / "a"
        run_b = tmp_path / "b"
        assert generate_sexchrom_test_data.main(["--outdir", str(run_a)]) == 0
        assert generate_sexchrom_test_data.main(["--outdir", str(run_b)]) == 0
        names = sorted(p.name for p in run_a.iterdir())
        assert names == ["joint_FF.vcf", "joint_MF.vcf", "joint_MM.vcf", "truth_table.tsv"]
        for name in names:
            assert filecmp.cmp(run_a / name, run_b / name, shallow=False), name
            assert filecmp.cmp(run_a / name, SEXCHROM_FIXTURE_DIR / name, shallow=False), (
                f"{name} differs from the committed fixture; regenerate with "
                "scripts/generate_sexchrom_test_data.py if the change is intended"
            )

    def test_mm_host_has_spurious_chrx_het(self) -> None:
        """The MM fixture gives the het-drop rule a male 0/1 at non-PAR chrX."""
        _, records = parse_text_vcf(SEXCHROM_FIXTURE_DIR / "joint_MM.vcf")
        hets = [
            r
            for r in records
            if r.chrom == "chrX" and not is_chrx_par(r.pos) and r.sample.startswith("0/1:")
        ]
        assert hets, "no spurious male het in the MM host column"
        # Every other male non-PAR chrX call is hom-encoded, never haploid.
        for r in records:
            if r.chrom == "chrX" and not is_chrx_par(r.pos):
                assert r.sample.split(":")[0] in ("0/0", "0/1", "1/1")

    def test_fixture_contents(self) -> None:
        for label in ("FF", "MM", "MF"):
            header, records = parse_text_vcf(SEXCHROM_FIXTURE_DIR / f"joint_{label}.vcf")
            assert header[-1].split("\t")[9:] == ["HOST", "DONOR", "ADMIX_F0.95", "ADMIX_F0.80"]
            chroms = [r.chrom for r in records]
            assert chroms.count("chrX_KI270880v1_alt") == 1
            x_par = [r for r in records if r.chrom == "chrX" and is_chrx_par(r.pos)]
            x_nonpar = [r for r in records if r.chrom == "chrX" and not is_chrx_par(r.pos)]
            assert len(x_par) == 2 and len(x_nonpar) == 20
            assert sum(1 for c in chroms if c not in ("chrX", "chrX_KI270880v1_alt")) == 60
