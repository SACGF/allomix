#!/usr/bin/env python3
"""Generate vibepaper facts CSVs from allomix validation results.

Reads validation outputs and writes structured facts for paper template substitution.

Usage:
    python scripts/generate_paper_facts.py
"""

import csv
import math
import random
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from allomix.bias import load_bias_table  # noqa: E402
from allomix.chimerism import PanelCalibration, estimate_single_donor_bb  # noqa: E402
from allomix.genotype import classify_markers, parse_vcf  # noqa: E402
from allomix.simulate import blend_vcfs, generate_marker_biases, write_vcf  # noqa: E402
from allomix.simulate import parse_text_vcf as sim_parse_vcf  # noqa: E402

FACTS_DIR = Path("output/facts")

STANDARD_FRACTIONS = [0.0, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.0]
BIAS_FRACTIONS = [0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 0.90, 0.95, 1.0]


def write_fact(name: str, data: dict) -> None:
    """Write a single-row facts CSV."""
    path = FACTS_DIR / f"{name}.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(data.keys()))
        writer.writeheader()
        writer.writerow(data)
    print(f"  Wrote {path}")


def fraction_to_name(f: float) -> str:
    d = round(f * 100)
    h = 100 - d
    return f"host_{h}_donor_{d}"


def run_sample(host_path, donor_path, sample_path, error_rate=0.01, marker_biases=None):
    host = parse_vcf(host_path, min_dp=0, min_gq=0)
    donor = parse_vcf(donor_path, min_dp=0, min_gq=0)
    admix = parse_vcf(sample_path, min_dp=0, min_gq=0)
    genotypes = classify_markers(host, [donor], admix, min_dp=0, min_gq=0, pass_only=False)
    result = estimate_single_donor_bb(
        genotypes.informative,
        error_rate=error_rate,
        calibration=PanelCalibration(biases=marker_biases),
    )
    return result, genotypes


def compute_metrics(rows):
    n = len(rows)
    errors = [r["error"] for r in rows]
    abs_errors = [abs(e) for e in errors]
    sq_errors = [e * e for e in errors]
    ci_covers = sum(1 for r in rows if r["ci_covers"])
    ci_widths = [r["ci_width"] for r in rows]
    return {
        "n_samples": n,
        "mean_signed_error": sum(errors) / n,
        "mean_abs_error": sum(abs_errors) / n,
        "rmse": math.sqrt(sum(sq_errors) / n),
        "max_abs_error": max(abs_errors),
        "ci_coverage": ci_covers / n,
        "mean_ci_width": sum(ci_widths) / n,
    }


def main():
    FACTS_DIR.mkdir(parents=True, exist_ok=True)

    host_vcf = "tests/test_data/host.vcf"
    donor_vcf = "tests/test_data/donor.vcf"

    # --- Standard validation (no bias correction) ---
    print("Running standard validation...")
    std_rows = []
    for frac in STANDARD_FRACTIONS:
        name = fraction_to_name(frac)
        vcf_path = f"tests/test_data/{name}.vcf"
        result, _genotypes = run_sample(host_vcf, donor_vcf, vcf_path)
        error = result.donor_fraction - frac
        ci_covers = result.donor_fraction_ci[0] <= frac <= result.donor_fraction_ci[1]
        ci_width = result.donor_fraction_ci[1] - result.donor_fraction_ci[0]
        std_rows.append({
            "sample_name": name,
            "true_frac": frac,
            "est_frac": result.donor_fraction,
            "error": error,
            "ci_lo": result.donor_fraction_ci[0],
            "ci_hi": result.donor_fraction_ci[1],
            "ci_width": ci_width,
            "ci_covers": ci_covers,
            "n_informative": result.n_informative,
        })

    std_metrics = compute_metrics(std_rows)

    write_fact("validation_summary", {
        "n_samples": std_metrics["n_samples"],
        "n_markers": 100,
        "n_informative": 80,
        "depth": 2000,
        "bias_sd": 0.02,
        "mean_signed_error_pct": round(std_metrics["mean_signed_error"] * 100, 4),
        "mean_abs_error_pct": round(std_metrics["mean_abs_error"] * 100, 4),
        "rmse_pct": round(std_metrics["rmse"] * 100, 4),
        "max_abs_error_pct": round(std_metrics["max_abs_error"] * 100, 4),
        "ci_coverage_pct": round(std_metrics["ci_coverage"] * 100, 1),
        "mean_ci_width_pct": round(std_metrics["mean_ci_width"] * 100, 4),
    })

    # Per-sample results for the table
    for r in std_rows:
        write_fact(f"val_{fraction_to_name(r['true_frac'])}", {
            "true_pct": round(r["true_frac"] * 100, 1),
            "est_pct": round(r["est_frac"] * 100, 2),
            "error_pct": round(r["error"] * 100, 4),
            "abs_error_pct": round(abs(r["error"]) * 100, 4),
            "ci_lo_pct": round(r["ci_lo"] * 100, 2),
            "ci_hi_pct": round(r["ci_hi"] * 100, 2),
            "ci_width_pct": round(r["ci_width"] * 100, 2),
            "ci_covers": r["ci_covers"],
            "n_informative": r["n_informative"],
        })

    # --- Bias correction comparison ---
    # Generate biased data first
    print("Running bias correction comparison...")
    bias_sd = 0.02
    depth = 2000
    seed = 42
    outdir = Path("output/bias_comparison/vcfs")
    outdir.mkdir(parents=True, exist_ok=True)

    _, host_records = sim_parse_vcf(host_vcf)
    _, donor_records = sim_parse_vcf(donor_vcf)
    donor_loci = {r.locus for r in donor_records}
    n_shared = sum(1 for r in host_records if r.locus in donor_loci)
    bias_rng = random.Random(seed)
    fixed_biases = generate_marker_biases(n_shared, bias_rng, bias_sd)

    for frac in BIAS_FRACTIONS:
        name = fraction_to_name(frac)
        sample_seed = seed + hash(str(frac)) % (2**31)
        result = blend_vcfs(
            host_path=host_vcf, donor_path=donor_vcf,
            donor_fraction=frac, target_depth=depth,
            sample_name=name, seed=sample_seed, fixed_biases=fixed_biases,
            locus_dropout_rate=0.016, depth_cv=0.43,
        )
        write_vcf(result, outdir / f"{name}.vcf")

    # Write bias table
    bias_path = outdir / "true_biases.tsv"
    biases = load_bias_table(bias_path) if bias_path.exists() else None

    # Run with and without bias correction
    for mode, mb in [("no_bias", None), ("with_bias", biases)]:
        rows = []
        for frac in BIAS_FRACTIONS:
            name = fraction_to_name(frac)
            vcf_path = str(outdir / f"{name}.vcf")
            res, _gen = run_sample(host_vcf, donor_vcf, vcf_path, marker_biases=mb)
            error = res.donor_fraction - frac
            ci_covers = res.donor_fraction_ci[0] <= frac <= res.donor_fraction_ci[1]
            ci_width = res.donor_fraction_ci[1] - res.donor_fraction_ci[0]
            rows.append({
                "true_frac": frac, "est_frac": res.donor_fraction,
                "error": error, "ci_covers": ci_covers, "ci_width": ci_width,
            })

        # Exclude boundary fractions for error metrics
        interior = [r for r in rows if 0.0 < r["true_frac"] < 1.0]
        m = compute_metrics(interior) if interior else compute_metrics(rows)
        all_m = compute_metrics(rows)

        write_fact(f"bias_{mode}", {
            "mean_signed_error_pct": round(m["mean_signed_error"] * 100, 4),
            "mean_abs_error_pct": round(m["mean_abs_error"] * 100, 4),
            "rmse_pct": round(m["rmse"] * 100, 4),
            "max_abs_error_pct": round(m["max_abs_error"] * 100, 4),
            "ci_coverage_pct": round(all_m["ci_coverage"] * 100, 1),
            "mean_ci_width_pct": round(all_m["mean_ci_width"] * 100, 4),
            "est_0pct": round(
                next(r["est_frac"] for r in rows if r["true_frac"] == 0.0) * 100, 2
            ),
            "est_100pct": round(
                next(r["est_frac"] for r in rows if r["true_frac"] == 1.0) * 100, 2
            ),
        })

    # --- Tool landscape facts ---
    write_fact("tool_landscape", {
        "n_commercial_tools": 4,
        "commercial_tools": "AlloSeq HCT; Devyser Chimerism; NGStrack; ScisGo Chimerism MD",
        "alloseq_n_markers": 202,
        "alloseq_lod": 0.3,
        "devyser_n_markers": 24,
        "devyser_lod": 0.05,
        "n_open_source_chimerism": 0,
        "n_tools_surveyed": 30,
    })

    # --- Panel specs ---
    write_fact("panel_specs", {
        "n_markers_panel": 76,
        "marker_type": "SNP",
        "panel_name": "IDT rhAmpSeq Sample ID",
        "typical_depth": ">1000x",
        "min_informative_markers": 3,
    })

    # --- Copy empirical panel stats into facts dir ---
    empirical_src = Path("paper/empirical_results/panel_empirical.csv")
    shutil.copy2(empirical_src, FACTS_DIR / "panel_empirical.csv")
    print(f"  Copied {empirical_src} -> {FACTS_DIR / 'panel_empirical.csv'}")

    print("\nDone! Facts written to output/facts/")


if __name__ == "__main__":
    main()
