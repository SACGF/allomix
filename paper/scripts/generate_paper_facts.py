#!/usr/bin/env python3
"""Generate vibepaper facts CSVs from allomix validation results.

Reads validation outputs and writes structured facts for paper template substitution.

Usage:
    python scripts/generate_paper_facts.py
"""

import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from allomix.calibration.bias import load_bias_table  # noqa: E402
from allomix.estimate.chimerism import estimate_single_donor_bb  # noqa: E402
from allomix.genotype import classify_markers, parse_vcf  # noqa: E402
from allomix.estimate.likelihood import PanelCalibration  # noqa: E402
from allomix.simulate import blend_vcfs, generate_marker_biases, write_vcf  # noqa: E402
from allomix.simulate import parse_text_vcf as sim_parse_vcf  # noqa: E402
from facts_io import read_fact_csv, write_fact_csv  # noqa: E402
from paper_seed import derive_seed  # noqa: E402

FACTS_DIR = Path("output/facts")

BIAS_FRACTIONS = [0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 0.90, 0.95, 1.0]


def write_fact(name: str, data: dict) -> None:
    """Write a typed facts CSV (see facts_io)."""
    path = write_fact_csv(FACTS_DIR / f"{name}.csv", data)
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

    # Bias correction comparison
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
        sample_seed = derive_seed(seed, frac)
        result = blend_vcfs(
            host_path=host_vcf, donor_path=donor_vcf,
            donor_fraction=frac, target_depth=depth,
            sample_name=name, seed=sample_seed, fixed_biases=fixed_biases,
            locus_dropout_rate=0.016, depth_cv=0.43,
        )
        write_vcf(result, outdir / f"{name}.vcf")

    bias_path = outdir / "true_biases.tsv"
    biases = load_bias_table(bias_path) if bias_path.exists() else None

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

    # Simulation calibration inputs shared across the validation scripts
    # (run_depth_validation.py, run_relatedness_validation.py, etc.), rounded
    # from empirical panel_empirical
    # values (depth CV 0.429, no-call 1.6%). Exposing them as a fact keeps the
    # paper prose in sync with what the simulator actually ran.
    write_fact("sim_calibration", {
        "depth_cv": 0.43,
        "locus_dropout_pct": 1.6,
        "seq_error_pct": 1,
    })

    # The commercial-tool figures (tool_landscape.*) are literature values, not
    # measured here: they live in paper/frozen_facts.csv with their citations.

    # Panel specs
    write_fact("panel_specs", {
        "n_markers_panel": 76,
        "marker_type": "SNP",
        "panel_name": "IDT rhAmpSeq Sample ID",
        "typical_depth": ">1000x",
        "min_informative_markers": 3,
    })

    # Empirical panel stats (measured on the internal data by
    # scripts/measure_panel_bias.py, committed as a one-row CSV) rewritten typed.
    empirical_src = Path("paper/empirical_results/panel_empirical.csv")
    write_fact_csv(FACTS_DIR / "panel_empirical.csv", read_fact_csv(empirical_src))
    print(f"  Wrote {empirical_src} -> {FACTS_DIR / 'panel_empirical.csv'} (typed)")

    print("\nDone! Facts written to output/facts/")


if __name__ == "__main__":
    main()
