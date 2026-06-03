"""Tests for allomix.simulate — synthetic chimeric VCF generation."""

import math
import os
import random
import statistics
import tempfile
import textwrap
from pathlib import Path

import pytest

from allomix.simulate import (
    HostAberration,
    alt_dose,
    assign_cnloh_aberrations,
    assign_cnv_aberrations,
    blend_vcfs,
    build_joint_vcf,
    cn_weighted_vaf,
    expected_vaf,
    extract_depth,
    extract_gt,
    generate_related_genotypes,
    gt_from_counts,
    is_informative,
    parse_text_vcf,
    sample_allele_counts,
    write_genotype_vcf,
    write_joint_vcf,
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

    def test_binomial_in_expected_range(self) -> None:
        """With depth=10000 and vaf=0.5, ALT count should be near 5000."""
        rng = random.Random(123)
        alt_counts = [sample_allele_counts(0.5, 10000, rng)[1] for _ in range(20)]
        mean_alt = sum(alt_counts) / len(alt_counts)
        # Should be within ~2% of 5000
        assert 4800 < mean_alt < 5200

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
                0.0, depth, rng0, error_rate=e, rho=100.0, rho_marker_type="het_only",
            )[1]
            for _ in range(n_trials)
        ]
        ref_at_one = [
            sample_allele_counts(
                1.0, depth, rng1, error_rate=e, rho=100.0, rho_marker_type="het_only",
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
        binom = [
            sample_allele_counts(0.5, depth, rng_bin)[1] / depth for _ in range(n)
        ]
        betab = [
            sample_allele_counts(
                0.5, depth, rng_bb, rho=rho, rho_marker_type="het_only",
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
                0.5, 1000, random.Random(1), rho=100.0, rho_marker_type="nope",
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

    def test_not_symmetric_model(self) -> None:
        """Verify we are NOT using the old symmetric model.

        Under the old symmetric model, vaf=0.0 with error_rate=0.03 would give
        p_obs = 0.03.  Under the 4-state model it should be 0.01.
        """
        rng = random.Random(42)
        n_trials = 200
        depth = 10000
        alt_counts = [
            sample_allele_counts(0.0, depth, rng, error_rate=0.03)[1] for _ in range(n_trials)
        ]
        mean_alt_rate = sum(alt_counts) / (n_trials * depth)
        # Under old symmetric model this would be ~0.03.
        # Under 4-state model it should be ~0.01.
        assert mean_alt_rate < 0.02, (
            f"ALT rate {mean_alt_rate:.4f} is too high; looks like the old symmetric error model"
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

TEST_DATA_DIR = Path(__file__).resolve().parent / "test_data"


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

    def test_from_existing_test_data(self) -> None:
        """Build from the existing host/donor test data VCFs."""
        host_path = TEST_DATA_DIR / "host.vcf"
        donor_path = TEST_DATA_DIR / "donor.vcf"
        if not host_path.exists():
            pytest.skip("Test data not available")

        result = build_joint_vcf(
            host_path=str(host_path),
            donor_paths=[str(donor_path)],
            admix_fractions=[0.10],
            admix_sample_names=["TP1"],
            target_depth=2000,
            seed=42,
        )
        assert result.num_markers > 0
        assert result.num_informative > 0


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
        aberrs = assign_cnloh_aberrations(markers, fraction_affected=1.0, clonal_fraction=1.0, rng=rng)
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
                host_path, donor_path, 0.1, target_depth=2000, seed=1,
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
            host_path, donor_path, 0.1, target_depth=20000, seed=7,
            error_rate=0.0, host_aberrations=aberr,
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
