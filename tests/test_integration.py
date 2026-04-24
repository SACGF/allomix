"""Integration tests -- full pipeline from joint VCFs through to results.

Uses the joint-called test VCFs (produced by build_joint_vcf) and runs
genotype -> chimerism -> qc -> report, verifying the output.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from allomix.chimerism import estimate_single_donor_bb
from allomix.cli import main
from allomix.genotype import classify_markers, parse_vcf
from allomix.qc import assess_quality
from allomix.report import timeline_json, to_json, to_tsv
from allomix.simulate import build_joint_vcf, write_joint_vcf

TEST_DATA_DIR = Path(__file__).resolve().parent / "test_data"
JOINT_VCF = TEST_DATA_DIR / "joint_single_donor.vcf"
JOINT_MULTI_VCF = TEST_DATA_DIR / "joint_multi_donor.vcf"


def _run_pipeline(
    vcf_path,
    host_sample="HOST",
    donor_sample="DONOR",
    admix_sample="ADMIX_F0.10",
    min_dp=0,
    min_gq=0,
):
    """Run the full genotype -> chimerism -> qc pipeline from a joint VCF."""
    host = parse_vcf(vcf_path, sample=host_sample, min_dp=0, min_gq=0)
    donor = parse_vcf(vcf_path, sample=donor_sample, min_dp=0, min_gq=0)
    admix = parse_vcf(vcf_path, sample=admix_sample, min_dp=0, min_gq=0)

    genotypes = classify_markers(
        host,
        [donor],
        admix,
        min_dp=min_dp,
        min_gq=min_gq,
        pass_only=False,
    )
    genotypes.sample_name = admix_sample

    result = estimate_single_donor_bb(genotypes.informative, error_rate=0.01)
    qc = assess_quality(result, genotypes)
    return result, qc, genotypes


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestFullPipeline:
    """End-to-end: joint VCF -> genotype -> chimerism -> QC."""

    def test_pure_host(self):
        """f=0.0: estimate should be ~0%."""
        result, _, _ = _run_pipeline(JOINT_VCF, admix_sample="ADMIX_F0.00")
        assert result.donor_fraction < 0.02

    def test_pure_donor(self):
        """f=1.0: estimate should be ~100%."""
        result, _, _ = _run_pipeline(JOINT_VCF, admix_sample="ADMIX_F1.00")
        assert result.donor_fraction > 0.98

    def test_fifty_fifty(self):
        """f=0.5: estimate should be near 50%."""
        result, _, _ = _run_pipeline(JOINT_VCF, admix_sample="ADMIX_F0.50")
        assert 0.40 < result.donor_fraction < 0.60

    def test_ten_percent(self):
        """f=0.10: estimate should be near 10%."""
        result, _, _ = _run_pipeline(JOINT_VCF, admix_sample="ADMIX_F0.10")
        assert 0.05 < result.donor_fraction < 0.20

    def test_one_percent(self):
        """f=0.01: estimate should be near 1% (testing low-fraction sensitivity)."""
        result, _, _ = _run_pipeline(JOINT_VCF, admix_sample="ADMIX_F0.01")
        assert result.donor_fraction < 0.05

    def test_ci_contains_truth(self):
        """CI should contain the true fraction for well-behaved data.

        With 100 markers and depth 2000, the CI is narrow so we test f=0.50
        where sampling variability has the least relative effect.
        """
        result, _, _ = _run_pipeline(JOINT_VCF, admix_sample="ADMIX_F0.50")
        lo, hi = result.donor_fraction_ci
        assert lo <= 0.50 <= hi, (
            f"f_true=0.50: CI [{lo:.4f}, {hi:.4f}] "
            f"does not contain truth, estimate={result.donor_fraction:.4f}"
        )

    def test_qc_passes(self):
        """QC should pass for well-behaved synthetic data."""
        _, qc, _ = _run_pipeline(JOINT_VCF, admix_sample="ADMIX_F0.10")
        assert qc.pass_
        assert qc.n_informative > 0

    def test_per_marker_results(self):
        """Per-marker results should be populated."""
        result, _, _ = _run_pipeline(JOINT_VCF, admix_sample="ADMIX_F0.10")
        assert len(result.per_marker) > 0
        for mr in result.per_marker:
            assert mr.dp > 0
            assert 0.0 <= mr.observed_vaf <= 1.0


class TestReportIntegration:
    """Integration: pipeline results -> report output."""

    @pytest.fixture
    def pipeline_result(self):
        return _run_pipeline(JOINT_VCF, admix_sample="ADMIX_F0.10")

    def test_tsv_output(self, pipeline_result, tmp_path):
        result, qc, _ = pipeline_result
        out = tmp_path / "results.tsv"
        with open(out, "w", encoding="utf-8") as f:
            to_tsv(result, qc, f)
        content = out.read_text()
        lines = content.strip().split("\n")
        assert len(lines) == 2  # header + data
        assert "donor_pct" in lines[0]

    def test_tsv_verbose(self, pipeline_result, tmp_path):
        result, qc, _ = pipeline_result
        out = tmp_path / "results_verbose.tsv"
        with open(out, "w", encoding="utf-8") as f:
            to_tsv(result, qc, f, verbose=True)
        content = out.read_text()
        lines = content.strip().split("\n")
        # header + data + blank + marker_header + marker_lines
        assert len(lines) > 4

    def test_json_output(self, pipeline_result):
        result, qc, genotypes = pipeline_result
        data = to_json(result, qc, sample_name=genotypes.sample_name)
        assert "donor_pct" in data
        assert isinstance(data["donor_pct"], float)
        json.dumps(data)

    def test_timeline(self):
        results = []
        for sample_name in ["ADMIX_F0.00", "ADMIX_F0.10", "ADMIX_F0.50"]:
            result, qc, genotypes = _run_pipeline(JOINT_VCF, admix_sample=sample_name)
            results.append((genotypes.sample_name, result, qc))

        data = timeline_json(results)
        assert "timepoints" in data
        assert len(data["timepoints"]) == 3
        json.dumps(data)


class TestCLIIntegration:
    """Test CLI wiring runs without error."""

    def test_monitor_tsv(self, tmp_path):
        out = tmp_path / "cli_out.tsv"
        rc = main(
            [
                "monitor",
                "--vcf",
                str(JOINT_VCF),
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR",
                "--sample",
                "ADMIX_F0.10",
                "--output",
                str(out),
                "--min-dp",
                "0",
                "--min-gq",
                "0",
            ]
        )
        assert rc == 0
        content = out.read_text()
        assert "donor_pct" in content

    def test_monitor_json(self, tmp_path):
        out = tmp_path / "cli_out.json"
        rc = main(
            [
                "monitor",
                "--vcf",
                str(JOINT_VCF),
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR",
                "--sample",
                "ADMIX_F0.10",
                "--output",
                str(out),
                "--format",
                "json",
                "--min-dp",
                "0",
                "--min-gq",
                "0",
            ]
        )
        assert rc == 0
        data = json.loads(out.read_text())
        assert "donor_pct" in data

    def test_timeline(self, tmp_path):
        out = tmp_path / "timeline.json"
        rc = main(
            [
                "timeline",
                "--vcf",
                str(JOINT_VCF),
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR",
                "--sample",
                "ADMIX_F0.00",
                "--sample",
                "ADMIX_F0.10",
                "--output",
                str(out),
                "--min-dp",
                "0",
                "--min-gq",
                "0",
            ]
        )
        assert rc == 0
        data = json.loads(out.read_text())
        assert len(data["timepoints"]) == 2

    def test_monitor_multiple_samples(self, tmp_path):
        """Monitor with multiple admixture samples."""
        out = tmp_path / "multi.tsv"
        rc = main(
            [
                "monitor",
                "--vcf",
                str(JOINT_VCF),
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR",
                "--sample",
                "ADMIX_F0.00",
                "--sample",
                "ADMIX_F0.10",
                "--sample",
                "ADMIX_F0.50",
                "--output",
                str(out),
                "--min-dp",
                "0",
                "--min-gq",
                "0",
            ]
        )
        assert rc == 0

    def test_invalid_sample_name(self):
        """CLI should fail with clear error for bad sample name."""
        with pytest.raises(SystemExit):
            main(
                [
                    "monitor",
                    "--vcf",
                    str(JOINT_VCF),
                    "--host-sample",
                    "NONEXISTENT",
                    "--donor-sample",
                    "DONOR",
                    "--sample",
                    "ADMIX_F0.10",
                ]
            )


class TestDynamicJointVcf:
    """Test building a joint VCF on the fly and running pipeline through it."""

    def test_round_trip(self, tmp_path):
        """Build joint VCF, write it, run pipeline, verify estimate."""
        host_vcf = TEST_DATA_DIR / "host.vcf"
        donor_vcf = TEST_DATA_DIR / "donor.vcf"

        result = build_joint_vcf(
            host_path=str(host_vcf),
            donor_paths=[str(donor_vcf)],
            admix_fractions=[0.20],
            admix_sample_names=["TP1"],
            host_sample_name="H",
            donor_sample_names=["D"],
            target_depth=2000,
            seed=99,
        )
        joint_path = tmp_path / "joint.vcf"
        write_joint_vcf(result, joint_path)

        # Run through the CLI
        out = tmp_path / "out.tsv"
        rc = main(
            [
                "monitor",
                "--vcf",
                str(joint_path),
                "--host-sample",
                "H",
                "--donor-sample",
                "D",
                "--sample",
                "TP1",
                "--output",
                str(out),
                "--min-dp",
                "0",
                "--min-gq",
                "0",
            ]
        )
        assert rc == 0
        content = out.read_text()
        assert "donor_pct" in content
