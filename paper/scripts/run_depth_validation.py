#!/usr/bin/env python3
"""Run allomix validation across multiple sequencing depths.

Generates synthetic chimeric VCFs at each depth, runs allomix, and produces:
1. Per-depth validation results and summary metrics
2. Multi-depth comparison figures
3. Facts CSVs for vibepaper

Usage:
    python scripts/run_depth_validation.py
    python scripts/run_depth_validation.py --depths 50 100 200 500 1000 2000
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from allomix.chimerism import estimate_single_donor
from allomix.genotype import classify_markers, parse_vcf
from allomix.simulate import (
    blend_vcfs,
    generate_marker_biases_realistic,
    parse_vcf as sim_parse_vcf,
    write_vcf,
)

FRACTIONS = [0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 0.90, 0.95, 1.0]
DEFAULT_DEPTHS = [50, 100, 200, 500, 1000]
FACTS_DIR = Path("output/facts")


def fraction_to_name(f: float) -> str:
    d = round(f * 100)
    h = 100 - d
    return f"host_{h}_donor_{d}"


def generate_and_run(
    host_vcf: str,
    donor_vcf: str,
    depth: int,
    bias_sd: float,
    seed: int,
    outdir: Path,
    depth_cv: float = 0.43,
    locus_dropout_rate: float = 0.016,
) -> list[dict]:
    """Generate synthetic data at a given depth and run allomix on it."""
    vcf_dir = outdir / f"depth_{depth}"
    vcf_dir.mkdir(parents=True, exist_ok=True)

    # Generate consistent per-marker biases (same across depths, heavy-tailed)
    _, host_records = sim_parse_vcf(host_vcf)
    _, donor_records = sim_parse_vcf(donor_vcf)
    donor_loci = {r.locus for r in donor_records}
    n_shared = sum(1 for r in host_records if r.locus in donor_loci)
    bias_rng = random.Random(seed)
    fixed_biases = generate_marker_biases_realistic(n_shared, bias_rng)

    # Generate blended VCFs
    for frac in FRACTIONS:
        name = fraction_to_name(frac)
        sample_seed = seed + hash(str(frac)) % (2**31)
        result = blend_vcfs(
            host_path=host_vcf, donor_path=donor_vcf,
            donor_fraction=frac, target_depth=depth,
            sample_name=name, seed=sample_seed, fixed_biases=fixed_biases,
            locus_dropout_rate=locus_dropout_rate,
            depth_cv=depth_cv,
        )
        write_vcf(result, vcf_dir / f"{name}.vcf")

    # Run allomix on each
    rows = []
    for frac in FRACTIONS:
        name = fraction_to_name(frac)
        vcf_path = str(vcf_dir / f"{name}.vcf")

        host = parse_vcf(host_vcf, min_dp=0, min_gq=0)
        donor = parse_vcf(donor_vcf, min_dp=0, min_gq=0)
        admix = parse_vcf(vcf_path, min_dp=0, min_gq=0)

        genotypes = classify_markers(host, [donor], admix, min_dp=0, min_gq=0, pass_only=False)
        result = estimate_single_donor(genotypes.informative, error_rate=0.01)

        error = result.donor_fraction - frac
        ci_covers = result.donor_fraction_ci[0] <= frac <= result.donor_fraction_ci[1]
        ci_width = result.donor_fraction_ci[1] - result.donor_fraction_ci[0]

        rows.append({
            "depth": depth,
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

    return rows


def compute_metrics(rows: list[dict]) -> dict:
    """Compute aggregate metrics, excluding boundary fractions."""
    interior = [r for r in rows if 0.0 < r["true_frac"] < 1.0]
    n = len(interior)
    if n == 0:
        return {}

    errors = [r["error"] for r in interior]
    abs_errors = [abs(e) for e in errors]
    sq_errors = [e * e for e in errors]

    all_rows = rows
    ci_covers = sum(1 for r in all_rows if r["ci_covers"])
    ci_widths = [r["ci_width"] for r in all_rows]

    return {
        "n_samples": len(all_rows),
        "n_interior": n,
        "mean_signed_error": sum(errors) / n,
        "mean_abs_error": sum(abs_errors) / n,
        "rmse": math.sqrt(sum(sq_errors) / n),
        "max_abs_error": max(abs_errors),
        "ci_coverage": ci_covers / len(all_rows),
        "mean_ci_width": sum(ci_widths) / len(all_rows),
    }


def write_fact(name: str, data: dict) -> None:
    path = FACTS_DIR / f"{name}.csv"
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(data.keys()))
        writer.writeheader()
        writer.writerow(data)


def plot_results(all_results: dict[int, list[dict]], outdir: Path) -> None:
    """Generate multi-depth comparison figures."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
    except ImportError:
        print("matplotlib not available, skipping plots", file=sys.stderr)
        return

    depths = sorted(all_results.keys())
    colors = plt.cm.viridis_r([i / (len(depths) - 1) for i in range(len(depths))])

    # --- Figure 1: Multi-panel scatter (truth vs estimated) ---
    n_depths = len(depths)
    fig, axes = plt.subplots(1, n_depths, figsize=(4 * n_depths, 4), sharey=True)
    if n_depths == 1:
        axes = [axes]

    for ax, depth, color in zip(axes, depths, colors):
        rows = all_results[depth]
        truths = [r["true_frac"] * 100 for r in rows]
        ests = [r["est_frac"] * 100 for r in rows]

        ax.scatter(truths, ests, c=[color], s=40, edgecolors="white", linewidth=0.5, zorder=3)
        ax.plot([0, 100], [0, 100], "k--", alpha=0.4, linewidth=1)
        ax.set_xlabel("True donor %")
        ax.set_title(f"{depth}x", fontsize=13, fontweight="bold")
        ax.set_xlim(-2, 102)
        ax.set_ylim(-2, 102)
        ax.set_aspect("equal")
        ax.grid(True, alpha=0.2)

    axes[0].set_ylabel("Estimated donor %")
    fig.suptitle("In silico validation across sequencing depths", fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(outdir / "fig1_depth_scatter.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Figure 2: Residuals across depths ---
    fig, axes = plt.subplots(1, n_depths, figsize=(4 * n_depths, 3.5), sharey=True)
    if n_depths == 1:
        axes = [axes]

    for ax, depth, color in zip(axes, depths, colors):
        rows = all_results[depth]
        truths = [r["true_frac"] * 100 for r in rows]
        errors = [r["error"] * 100 for r in rows]

        ax.scatter(truths, errors, c=[color], s=40, edgecolors="white", linewidth=0.5, zorder=3)
        ax.axhline(0, color="k", linestyle="--", alpha=0.4)
        ax.set_xlabel("True donor %")
        ax.set_title(f"{depth}x", fontsize=13, fontweight="bold")
        ax.grid(True, alpha=0.2)

    axes[0].set_ylabel("Error (est - true) %")
    fig.suptitle("Estimation error across sequencing depths", fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(outdir / "fig2_depth_residuals.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Figure 3: Summary metrics vs depth ---
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    metrics_by_depth = {}
    for depth in depths:
        metrics_by_depth[depth] = compute_metrics(all_results[depth])

    mae_vals = [metrics_by_depth[d]["mean_abs_error"] * 100 for d in depths]
    rmse_vals = [metrics_by_depth[d]["rmse"] * 100 for d in depths]
    max_vals = [metrics_by_depth[d]["max_abs_error"] * 100 for d in depths]
    ci_cov_vals = [metrics_by_depth[d]["ci_coverage"] * 100 for d in depths]
    ci_width_vals = [metrics_by_depth[d]["mean_ci_width"] * 100 for d in depths]

    # Accuracy vs depth
    ax = axes[0]
    ax.plot(depths, mae_vals, "o-", color="steelblue", linewidth=2, markersize=8, label="MAE")
    ax.plot(depths, rmse_vals, "s-", color="darkorange", linewidth=2, markersize=8, label="RMSE")
    ax.plot(depths, max_vals, "^-", color="firebrick", linewidth=2, markersize=8, label="Max error")
    ax.set_xlabel("Sequencing depth (x)", fontsize=11)
    ax.set_ylabel("Error (%)", fontsize=11)
    ax.set_title("Accuracy vs depth", fontsize=12)
    ax.set_xscale("log")
    ax.set_xticks(depths)
    ax.set_xticklabels([str(d) for d in depths])
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # CI coverage vs depth
    ax = axes[1]
    ax.plot(depths, ci_cov_vals, "o-", color="steelblue", linewidth=2, markersize=8)
    ax.axhline(95, color="k", linestyle="--", alpha=0.4, label="Nominal 95%")
    ax.set_xlabel("Sequencing depth (x)", fontsize=11)
    ax.set_ylabel("CI coverage (%)", fontsize=11)
    ax.set_title("CI coverage vs depth", fontsize=12)
    ax.set_xscale("log")
    ax.set_xticks(depths)
    ax.set_xticklabels([str(d) for d in depths])
    ax.set_ylim(0, 105)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # CI width vs depth
    ax = axes[2]
    ax.plot(depths, ci_width_vals, "o-", color="steelblue", linewidth=2, markersize=8)
    ax.set_xlabel("Sequencing depth (x)", fontsize=11)
    ax.set_ylabel("Mean CI width (%)", fontsize=11)
    ax.set_title("CI width vs depth", fontsize=12)
    ax.set_xscale("log")
    ax.set_xticks(depths)
    ax.set_xticklabels([str(d) for d in depths])
    ax.grid(True, alpha=0.3)

    fig.suptitle("allomix performance vs sequencing depth", fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(outdir / "fig3_depth_summary.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    print(f"Figures saved to {outdir}/", file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run allomix validation across multiple sequencing depths.",
    )
    parser.add_argument("--host", default="tests/test_data/host.vcf")
    parser.add_argument("--donor", default="tests/test_data/donor.vcf")
    parser.add_argument(
        "--depths", type=int, nargs="+", default=DEFAULT_DEPTHS,
        help=f"Depths to test (default: {DEFAULT_DEPTHS})",
    )
    parser.add_argument("--bias-sd", type=float, default=0.02)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--outdir", default="output/depth_validation")
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    FACTS_DIR.mkdir(parents=True, exist_ok=True)

    all_results: dict[int, list[dict]] = {}

    for depth in sorted(args.depths):
        print(f"\n{'='*50}", file=sys.stderr)
        print(f"Depth: {depth}x", file=sys.stderr)
        print(f"{'='*50}", file=sys.stderr)

        rows = generate_and_run(
            args.host, args.donor, depth,
            bias_sd=args.bias_sd, seed=args.seed, outdir=outdir,
        )
        all_results[depth] = rows

        metrics = compute_metrics(rows)
        print(f"  MAE:          {metrics['mean_abs_error']*100:.4f}%", file=sys.stderr)
        print(f"  RMSE:         {metrics['rmse']*100:.4f}%", file=sys.stderr)
        print(f"  Max error:    {metrics['max_abs_error']*100:.4f}%", file=sys.stderr)
        print(f"  CI coverage:  {metrics['ci_coverage']:.1%}", file=sys.stderr)
        print(f"  CI width:     {metrics['mean_ci_width']*100:.4f}%", file=sys.stderr)

        # Write per-depth facts
        write_fact(f"depth_{depth}", {
            "depth": depth,
            "n_samples": metrics["n_samples"],
            "mean_abs_error_pct": round(metrics["mean_abs_error"] * 100, 4),
            "rmse_pct": round(metrics["rmse"] * 100, 4),
            "max_abs_error_pct": round(metrics["max_abs_error"] * 100, 4),
            "ci_coverage_pct": round(metrics["ci_coverage"] * 100, 1),
            "mean_ci_width_pct": round(metrics["mean_ci_width"] * 100, 4),
        })

    # Write combined summary table as facts
    depths_sorted = sorted(args.depths)
    summary_rows = []
    for depth in depths_sorted:
        m = compute_metrics(all_results[depth])
        summary_rows.append({
            "depth": depth,
            "mae": round(m["mean_abs_error"] * 100, 2),
            "rmse": round(m["rmse"] * 100, 2),
            "max_err": round(m["max_abs_error"] * 100, 2),
            "ci_cov": round(m["ci_coverage"] * 100, 1),
            "ci_width": round(m["mean_ci_width"] * 100, 2),
        })

    # Write summary table CSV (multi-row, not a vibepaper fact)
    summary_path = outdir / "depth_summary.tsv"
    with open(summary_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["depth", "mae", "rmse", "max_err", "ci_cov", "ci_width"],
                           delimiter="\t")
        w.writeheader()
        w.writerows(summary_rows)
    print(f"\nSummary table: {summary_path}", file=sys.stderr)

    # Write per-sample results for each depth
    for depth in depths_sorted:
        results_path = outdir / f"results_{depth}x.tsv"
        with open(results_path, "w", newline="") as f:
            fields = ["sample_name", "true_pct", "est_pct", "error_pct",
                       "ci_lo_pct", "ci_hi_pct", "ci_covers"]
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            w.writeheader()
            for r in all_results[depth]:
                w.writerow({
                    "sample_name": r["sample_name"],
                    "true_pct": f"{r['true_frac']*100:.1f}",
                    "est_pct": f"{r['est_frac']*100:.2f}",
                    "error_pct": f"{r['error']*100:.4f}",
                    "ci_lo_pct": f"{r['ci_lo']*100:.2f}",
                    "ci_hi_pct": f"{r['ci_hi']*100:.2f}",
                    "ci_covers": r["ci_covers"],
                })

    # Generate plots
    plot_results(all_results, outdir)

    # Copy figures to facts dir
    import shutil
    for fig in ["fig1_depth_scatter.png", "fig2_depth_residuals.png", "fig3_depth_summary.png"]:
        src = outdir / fig
        if src.exists():
            shutil.copy2(src, FACTS_DIR / fig)
            print(f"  Copied to {FACTS_DIR / fig}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
