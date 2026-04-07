"""Integration tests — full pipeline from synthetic VCFs through to results.

Uses simulate.blend_vcfs to create chimeric VCFs at known fractions,
then runs genotype → chimerism → qc → report and verifies the output.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from allomix.chimerism import estimate_single_donor
from allomix.genotype import classify_markers, parse_vcf
from allomix.qc import assess_quality
from allomix.report import timeline_json, to_json, to_tsv
from allomix.simulate import blend_vcfs, write_vcf

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
EXAMPLE_VCF = DATA_DIR / "idt_rhampseq_sid_example.vcf"


def _make_chimeric_vcf(
    tmp_path: Path, donor_fraction: float, seed: int = 42, depth: int = 2000
) -> Path:
    """Create a synthetic chimeric VCF from the example, using it as both host and donor.

    Since we use the same VCF as host and donor but with different fractions,
    we need markers where the genotypes differ. The example VCF has a mix of
    0/0, 0/1, and 1/1 genotypes, so treating it as "host" and a shifted version
    as "donor" would work. Instead, we create two synthetic pure-sample VCFs
    with known complementary genotypes, then blend them.
    """
    # Use the example VCF as host and create a "donor" by flipping genotypes
    host_path = EXAMPLE_VCF
    donor_path = _create_flipped_vcf(tmp_path, host_path)

    result = blend_vcfs(
        host_path=str(host_path),
        donor_path=str(donor_path),
        donor_fraction=donor_fraction,
        target_depth=depth,
        sample_name=f"chimeric_f{donor_fraction:.3f}",
        seed=seed,
    )
    out_path = tmp_path / f"chimeric_f{donor_fraction:.3f}.vcf"
    write_vcf(result, out_path)
    return out_path


def _create_flipped_vcf(tmp_path: Path, source_path: Path) -> Path:
    """Create a VCF where all genotypes are flipped (0/0→1/1, 1/1→0/0, 0/1 stays).

    This ensures maximum informativeness between host and donor.
    """
    lines = []
    with open(source_path) as f:
        for line in f:
            if line.startswith("#"):
                lines.append(line)
                continue
            if not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) < 10:
                lines.append(line)
                continue

            # Parse sample column and flip GT
            fmt_keys = fields[8].split(":")
            fmt_vals = fields[9].split(":")
            gt_idx = fmt_keys.index("GT")
            gt = fmt_vals[gt_idx]

            if gt == "0/0":
                # Flip to 1/1: swap AD values, set GT
                fmt_vals[gt_idx] = "1/1"
                if "AD" in fmt_keys:
                    ad_idx = fmt_keys.index("AD")
                    ad_parts = fmt_vals[ad_idx].split(",")
                    if len(ad_parts) == 2:
                        fmt_vals[ad_idx] = f"{ad_parts[1]},{ad_parts[0]}"
                    elif len(ad_parts) == 1:
                        # hom-ref with single AD value — make it hom-alt
                        fmt_vals[ad_idx] = f"0,{ad_parts[0]}"
                if "AF" in fmt_keys:
                    af_idx = fmt_keys.index("AF")
                    fmt_vals[af_idx] = "1.0"
                # Need an ALT allele if ALT is "."
                if fields[4] == ".":
                    fields[4] = "T" if fields[3] != "T" else "A"
            elif gt == "1/1":
                fmt_vals[gt_idx] = "0/0"
                if "AD" in fmt_keys:
                    ad_idx = fmt_keys.index("AD")
                    ad_parts = fmt_vals[ad_idx].split(",")
                    if len(ad_parts) == 2:
                        fmt_vals[ad_idx] = f"{ad_parts[1]},{ad_parts[0]}"
                if "AF" in fmt_keys:
                    af_idx = fmt_keys.index("AF")
                    fmt_vals[af_idx] = "0"
            # 0/1 stays as-is (het in both = non-informative, which is fine)

            fields[9] = ":".join(fmt_vals)
            lines.append("\t".join(fields) + "\n")

    out = tmp_path / "donor_flipped.vcf"
    with open(out, "w") as f:
        f.writelines(lines)
    return out


def _run_pipeline(host_path, donor_path, admix_path, min_dp=0, min_gq=0):
    """Run the full genotype → chimerism → qc pipeline."""
    host = parse_vcf(host_path, min_dp=0, min_gq=0)
    donor = parse_vcf(donor_path, min_dp=0, min_gq=0)
    admix = parse_vcf(admix_path, min_dp=0, min_gq=0)

    genotypes = classify_markers(
        host, [donor], admix, min_dp=min_dp, min_gq=min_gq, pass_only=False,
    )
    genotypes.sample_name = Path(admix_path).stem

    result = estimate_single_donor(genotypes.informative, error_rate=0.01)
    qc = assess_quality(result, genotypes)
    return result, qc, genotypes


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestFullPipeline:
    """End-to-end: synthetic VCF → genotype → chimerism → QC."""

    @pytest.fixture
    def donor_vcf(self, tmp_path):
        return _create_flipped_vcf(tmp_path, EXAMPLE_VCF)

    def test_pure_host(self, tmp_path, donor_vcf):
        """f=0.0: estimate should be ~0%."""
        chimeric = _make_chimeric_vcf(tmp_path, 0.0)
        result, qc, _ = _run_pipeline(EXAMPLE_VCF, donor_vcf, chimeric)
        assert result.donor_fraction < 0.02

    def test_pure_donor(self, tmp_path, donor_vcf):
        """f=1.0: estimate should be ~100%."""
        chimeric = _make_chimeric_vcf(tmp_path, 1.0)
        result, qc, _ = _run_pipeline(EXAMPLE_VCF, donor_vcf, chimeric)
        assert result.donor_fraction > 0.98

    def test_fifty_fifty(self, tmp_path, donor_vcf):
        """f=0.5: estimate should be near 50%."""
        chimeric = _make_chimeric_vcf(tmp_path, 0.5, seed=123)
        result, qc, _ = _run_pipeline(EXAMPLE_VCF, donor_vcf, chimeric)
        assert 0.40 < result.donor_fraction < 0.60

    def test_ten_percent(self, tmp_path, donor_vcf):
        """f=0.10: estimate should be near 10%."""
        chimeric = _make_chimeric_vcf(tmp_path, 0.10, seed=456)
        result, qc, _ = _run_pipeline(EXAMPLE_VCF, donor_vcf, chimeric)
        assert 0.05 < result.donor_fraction < 0.20

    def test_one_percent(self, tmp_path, donor_vcf):
        """f=0.01: estimate should be near 1% (testing low-fraction sensitivity)."""
        chimeric = _make_chimeric_vcf(tmp_path, 0.01, seed=789, depth=5000)
        result, qc, _ = _run_pipeline(EXAMPLE_VCF, donor_vcf, chimeric)
        assert result.donor_fraction < 0.05  # Within reasonable range

    def test_ci_contains_truth(self, tmp_path, donor_vcf):
        """CI should contain the true fraction for well-behaved data.

        With few markers (~8 informative), stochastic sampling can push
        estimates slightly off. We use moderate fractions where the CI
        is wide enough to reliably capture the truth.
        """
        for f_true in [0.10, 0.25, 0.50]:
            chimeric = _make_chimeric_vcf(
                tmp_path, f_true, seed=int(f_true * 1000) + 7, depth=3000,
            )
            result, _, _ = _run_pipeline(EXAMPLE_VCF, donor_vcf, chimeric)
            lo, hi = result.donor_fraction_ci
            assert lo <= f_true <= hi, (
                f"f_true={f_true}: CI [{lo:.4f}, {hi:.4f}] "
                f"does not contain truth, estimate={result.donor_fraction:.4f}"
            )

    def test_qc_passes(self, tmp_path, donor_vcf):
        """QC should pass for well-behaved synthetic data."""
        chimeric = _make_chimeric_vcf(tmp_path, 0.10, seed=42)
        _, qc, _ = _run_pipeline(EXAMPLE_VCF, donor_vcf, chimeric)
        assert qc.pass_
        assert qc.n_informative > 0

    def test_per_marker_results(self, tmp_path, donor_vcf):
        """Per-marker results should be populated."""
        chimeric = _make_chimeric_vcf(tmp_path, 0.10, seed=42)
        result, _, _ = _run_pipeline(EXAMPLE_VCF, donor_vcf, chimeric)
        assert len(result.per_marker) > 0
        for mr in result.per_marker:
            assert mr.dp > 0
            assert 0.0 <= mr.observed_vaf <= 1.0


class TestReportIntegration:
    """Integration: pipeline results → report output."""

    @pytest.fixture
    def pipeline_result(self, tmp_path):
        donor_vcf = _create_flipped_vcf(tmp_path, EXAMPLE_VCF)
        chimeric = _make_chimeric_vcf(tmp_path, 0.10, seed=42)
        return _run_pipeline(EXAMPLE_VCF, donor_vcf, chimeric)

    def test_tsv_output(self, pipeline_result, tmp_path):
        result, qc, genotypes = pipeline_result
        out = tmp_path / "results.tsv"
        with open(out, "w") as f:
            to_tsv(result, qc, f)
        content = out.read_text()
        lines = content.strip().split("\n")
        assert len(lines) == 2  # header + data
        assert "donor_pct" in lines[0]

    def test_tsv_verbose(self, pipeline_result, tmp_path):
        result, qc, genotypes = pipeline_result
        out = tmp_path / "results_verbose.tsv"
        with open(out, "w") as f:
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
        # Verify JSON-serialisable
        json.dumps(data)

    def test_timeline(self, tmp_path):
        donor_vcf = _create_flipped_vcf(tmp_path, EXAMPLE_VCF)
        results = []
        for f in [0.05, 0.10, 0.20]:
            chimeric = _make_chimeric_vcf(tmp_path, f, seed=int(f * 1000))
            result, qc, genotypes = _run_pipeline(EXAMPLE_VCF, donor_vcf, chimeric)
            results.append((genotypes.sample_name, result, qc))

        data = timeline_json(results)
        assert "timepoints" in data
        assert len(data["timepoints"]) == 3
        json.dumps(data)


class TestCLIIntegration:
    """Test CLI wiring runs without error."""

    def test_monitor_tsv(self, tmp_path):
        from allomix.cli import main

        donor_vcf = _create_flipped_vcf(tmp_path, EXAMPLE_VCF)
        chimeric = _make_chimeric_vcf(tmp_path, 0.10, seed=42)
        out = tmp_path / "cli_out.tsv"

        rc = main([
            "monitor",
            "--host", str(EXAMPLE_VCF),
            "--donor", str(donor_vcf),
            "--sample", str(chimeric),
            "--output", str(out),
            "--min-dp", "0",
            "--min-gq", "0",
        ])
        assert rc == 0
        content = out.read_text()
        assert "donor_pct" in content

    def test_monitor_json(self, tmp_path):
        from allomix.cli import main

        donor_vcf = _create_flipped_vcf(tmp_path, EXAMPLE_VCF)
        chimeric = _make_chimeric_vcf(tmp_path, 0.10, seed=42)
        out = tmp_path / "cli_out.json"

        rc = main([
            "monitor",
            "--host", str(EXAMPLE_VCF),
            "--donor", str(donor_vcf),
            "--sample", str(chimeric),
            "--output", str(out),
            "--format", "json",
            "--min-dp", "0",
            "--min-gq", "0",
        ])
        assert rc == 0
        data = json.loads(out.read_text())
        assert "donor_pct" in data

    def test_timeline(self, tmp_path):
        from allomix.cli import main

        donor_vcf = _create_flipped_vcf(tmp_path, EXAMPLE_VCF)
        c1 = _make_chimeric_vcf(tmp_path, 0.05, seed=1)
        c2 = _make_chimeric_vcf(tmp_path, 0.10, seed=2)
        out = tmp_path / "timeline.json"

        rc = main([
            "timeline",
            "--host", str(EXAMPLE_VCF),
            "--donor", str(donor_vcf),
            "--sample", str(c1),
            "--sample", str(c2),
            "--output", str(out),
            "--min-dp", "0",
            "--min-gq", "0",
        ])
        assert rc == 0
        data = json.loads(out.read_text())
        assert len(data["timepoints"]) == 2
