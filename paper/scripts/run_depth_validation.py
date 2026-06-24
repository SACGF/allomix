#!/usr/bin/env python3
"""Run allomix validation across multiple sequencing depths with replicates.

Generates synthetic chimeric VCFs at each depth with N independent replicates
(different random seeds for per-marker bias and sampling noise), runs allomix,
and produces:
1. Per-depth validation results and summary metrics (mean ± SD across replicates)
2. Multi-depth comparison figures including boxplots
3. Facts CSVs for vibepaper

Usage:
    python scripts/run_depth_validation.py
    python scripts/run_depth_validation.py --depths 50 100 200 500 1000 --n-replicates 5
"""

import argparse
import csv
import math
import random
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from allomix.chimerism import estimate_single_donor_bb  # noqa: E402
from allomix.genotype import classify_markers, parse_vcf  # noqa: E402
from allomix.likelihood import PanelCalibration  # noqa: E402
from allomix.simulate import (  # noqa: E402
    blend_vcfs,
    generate_marker_biases_realistic,
    write_vcf,
)
from allomix.simulate import (  # noqa: E402
    parse_text_vcf as sim_parse_vcf,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paper_quick import qval  # noqa: E402  (also patches savefig for the watermark)

FRACTIONS = [0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 0.90, 0.95, 1.0]
DEFAULT_DEPTHS = [50, 100, 200, 500, 1000]
DEFAULT_N_REPLICATES = qval(5, 2)
FACTS_DIR = Path("output/facts")


def fraction_to_name(f: float) -> str:
    d = round(f * 100)
    h = 100 - d
    return f"host_{h}_donor_{d}"


def generate_and_run(
    host_vcf: str,
    donor_vcf: str,
    depth: int,
    seed: int,
    outdir: Path,
    depth_cv: float = 0.43,
    locus_dropout_rate: float = 0.016,
    bias_correction: bool = False,
) -> list[dict]:
    """Generate synthetic data at a given depth and run allomix on it."""
    vcf_dir = outdir / f"depth_{depth}" / f"seed_{seed}"
    vcf_dir.mkdir(parents=True, exist_ok=True)

    # Generate consistent per-marker biases (same across depths, heavy-tailed)
    _, host_records = sim_parse_vcf(host_vcf)
    _, donor_records = sim_parse_vcf(donor_vcf)
    donor_loci = {r.locus for r in donor_records}
    n_shared = sum(1 for r in host_records if r.locus in donor_loci)
    bias_rng = random.Random(seed)
    fixed_biases = generate_marker_biases_realistic(n_shared, bias_rng)

    # Generate blended VCFs and capture bias mapping
    bias_dict = None
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
        # Capture bias dict from first blend (same biases for all fractions)
        if bias_dict is None and result.marker_biases is not None:
            bias_dict = {(c, p, r, a): b for c, p, r, a, b in result.marker_biases}

    # Run allomix on each
    marker_biases = bias_dict if bias_correction else None
    rows = []
    for frac in FRACTIONS:
        name = fraction_to_name(frac)
        vcf_path = str(vcf_dir / f"{name}.vcf")

        host = parse_vcf(host_vcf, min_dp=0, min_gq=0)
        donor = parse_vcf(donor_vcf, min_dp=0, min_gq=0)
        admix = parse_vcf(vcf_path, min_dp=0, min_gq=0)

        genotypes = classify_markers(host, [donor], admix, min_dp=0, min_gq=0, pass_only=False)
        result = estimate_single_donor_bb(
            genotypes.informative,
            error_rate=0.01,
            calibration=PanelCalibration(biases=marker_biases),
        )

        error = result.donor_fraction - frac
        ci_covers = result.donor_fraction_ci[0] <= frac <= result.donor_fraction_ci[1]
        ci_width = result.donor_fraction_ci[1] - result.donor_fraction_ci[0]

        rows.append({
            "depth": depth,
            "seed": seed,
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
    """Compute aggregate metrics, excluding boundary fractions for error metrics."""
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


def aggregate_replicate_metrics(replicate_metrics: list[dict]) -> dict:
    """Compute mean and SD of metrics across replicates."""
    n = len(replicate_metrics)
    keys = ["mean_abs_error", "rmse", "max_abs_error", "ci_coverage", "mean_ci_width"]
    agg = {}
    for key in keys:
        vals = [m[key] for m in replicate_metrics]
        mean = sum(vals) / n
        sd = math.sqrt(sum((v - mean) ** 2 for v in vals) / (n - 1)) if n > 1 else 0.0
        agg[f"{key}_mean"] = mean
        agg[f"{key}_sd"] = sd
    agg["n_replicates"] = n
    return agg


def write_fact(name: str, data: dict) -> None:
    path = FACTS_DIR / f"{name}.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(data.keys()))
        writer.writeheader()
        writer.writerow(data)


def plot_results(
    all_results: dict[int, list[list[dict]]],
    all_metrics: dict[int, list[dict]],
    outdir: Path,
) -> None:
    """Generate multi-depth comparison figures with replicate support."""
    depths = sorted(all_results.keys())
    n_depths = len(depths)

    # --- Figure 1: Multi-panel scatter (truth vs estimated), all replicates overlaid ---
    fig, axes = plt.subplots(1, n_depths, figsize=(4 * n_depths, 4), sharey=True)
    if n_depths == 1:
        axes = [axes]

    colors = plt.get_cmap("viridis_r")([i / (len(depths) - 1) for i in range(len(depths))])

    for ax, depth, color in zip(axes, depths, colors):
        for rep_rows in all_results[depth]:
            truths = [r["true_frac"] * 100 for r in rep_rows]
            ests = [r["est_frac"] * 100 for r in rep_rows]
            ax.scatter(
                truths, ests, c=[color], s=25, alpha=0.5,
                edgecolors="white", linewidth=0.3, zorder=3,
            )
        ax.plot([0, 100], [0, 100], "k--", alpha=0.4, linewidth=1)
        ax.set_xlabel("True donor %")
        ax.set_title(f"{depth}x", fontsize=13, fontweight="bold")
        ax.set_xlim(-2, 102)
        ax.set_ylim(-2, 102)
        ax.set_aspect("equal")
        ax.grid(True, alpha=0.2)

    axes[0].set_ylabel("Estimated donor %")
    n_reps = len(next(iter(all_results.values())))
    fig.suptitle(
        f"In silico validation across sequencing depths (N={n_reps} replicates)",
        fontsize=14, y=1.02,
    )
    fig.tight_layout()
    fig.savefig(outdir / "fig1_depth_scatter.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Figure 2: Boxplots of absolute error by depth ---
    fig, ax = plt.subplots(figsize=(8, 5))

    box_data = []
    box_labels = []
    for depth in depths:
        # Collect absolute errors from interior fractions across all replicates
        abs_errors = []
        for rep_rows in all_results[depth]:
            for r in rep_rows:
                if 0.0 < r["true_frac"] < 1.0:
                    abs_errors.append(abs(r["error"]) * 100)
        box_data.append(abs_errors)
        box_labels.append(f"{depth}x")

    ax.boxplot(
        box_data, tick_labels=box_labels, patch_artist=True,
        boxprops=dict(facecolor="steelblue", alpha=0.6),
        medianprops=dict(color="firebrick", linewidth=2),
        whiskerprops=dict(color="grey"),
        capprops=dict(color="grey"),
        flierprops=dict(marker="o", markerfacecolor="steelblue", alpha=0.4, markersize=4),
    )
    ax.set_xlabel("Sequencing depth", fontsize=12)
    ax.set_ylabel("Absolute error (%)", fontsize=12)
    ax.set_title(
        f"Estimation error distribution by depth (N={n_reps} replicates)",
        fontsize=13,
    )
    ax.grid(True, alpha=0.3, axis="y")
    fig.tight_layout()
    fig.savefig(outdir / "fig2_depth_boxplots.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Figure 3: Summary metrics vs depth with error bars ---
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    mae_means = []
    mae_sds = []
    rmse_means = []
    rmse_sds = []
    max_means = []
    max_sds = []
    ci_cov_means = []
    ci_cov_sds = []
    ci_width_means = []
    ci_width_sds = []

    for depth in depths:
        agg = aggregate_replicate_metrics(all_metrics[depth])
        mae_means.append(agg["mean_abs_error_mean"] * 100)
        mae_sds.append(agg["mean_abs_error_sd"] * 100)
        rmse_means.append(agg["rmse_mean"] * 100)
        rmse_sds.append(agg["rmse_sd"] * 100)
        max_means.append(agg["max_abs_error_mean"] * 100)
        max_sds.append(agg["max_abs_error_sd"] * 100)
        ci_cov_means.append(agg["ci_coverage_mean"] * 100)
        ci_cov_sds.append(agg["ci_coverage_sd"] * 100)
        ci_width_means.append(agg["mean_ci_width_mean"] * 100)
        ci_width_sds.append(agg["mean_ci_width_sd"] * 100)

    # Accuracy vs depth
    ax = axes[0]
    ax.errorbar(
        depths, mae_means, yerr=mae_sds, fmt="o-",
        color="steelblue", linewidth=2, markersize=8, capsize=4, label="MAE",
    )
    ax.errorbar(
        depths, rmse_means, yerr=rmse_sds, fmt="s-",
        color="darkorange", linewidth=2, markersize=8, capsize=4, label="RMSE",
    )
    ax.errorbar(
        depths, max_means, yerr=max_sds, fmt="^-",
        color="firebrick", linewidth=2, markersize=8, capsize=4, label="Max error",
    )
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
    ax.errorbar(
        depths, ci_cov_means, yerr=ci_cov_sds, fmt="o-",
        color="steelblue", linewidth=2, markersize=8, capsize=4,
    )
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
    ax.errorbar(
        depths, ci_width_means, yerr=ci_width_sds, fmt="o-",
        color="steelblue", linewidth=2, markersize=8, capsize=4,
    )
    ax.set_xlabel("Sequencing depth (x)", fontsize=11)
    ax.set_ylabel("Mean CI width (%)", fontsize=11)
    ax.set_title("CI width vs depth", fontsize=12)
    ax.set_xscale("log")
    ax.set_xticks(depths)
    ax.set_xticklabels([str(d) for d in depths])
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        f"allomix performance vs sequencing depth (N={n_reps}, mean ± SD)",
        fontsize=14, y=1.02,
    )
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
    parser.add_argument(
        "--n-replicates", "-n", type=int, default=DEFAULT_N_REPLICATES,
        help=f"Number of replicates per depth (default: {DEFAULT_N_REPLICATES})",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--outdir", default="output/depth_validation")
    parser.add_argument(
        "--bias-correction", action="store_true",
        help="Pass known marker biases to the estimator for bias correction",
    )
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    FACTS_DIR.mkdir(parents=True, exist_ok=True)

    # all_results[depth] = list of replicate row-lists
    # all_metrics[depth] = list of per-replicate metric dicts
    all_results: dict[int, list[list[dict]]] = {}
    all_metrics: dict[int, list[dict]] = {}

    for depth in sorted(args.depths):
        print(f"\n{'='*50}", file=sys.stderr)
        print(f"Depth: {depth}x  ({args.n_replicates} replicates)", file=sys.stderr)
        print(f"{'='*50}", file=sys.stderr)

        all_results[depth] = []
        all_metrics[depth] = []

        for rep in range(args.n_replicates):
            rep_seed = args.seed + rep * 1000
            rows = generate_and_run(
                args.host, args.donor, depth, seed=rep_seed, outdir=outdir,
                bias_correction=args.bias_correction,
            )
            all_results[depth].append(rows)

            metrics = compute_metrics(rows)
            all_metrics[depth].append(metrics)
            print(
                f"  Rep {rep}: MAE={metrics['mean_abs_error']*100:.4f}%  "
                f"RMSE={metrics['rmse']*100:.4f}%  "
                f"Max={metrics['max_abs_error']*100:.4f}%  "
                f"CI={metrics['ci_coverage']:.1%}",
                file=sys.stderr,
            )

        agg = aggregate_replicate_metrics(all_metrics[depth])
        print(
            f"  Mean MAE: {agg['mean_abs_error_mean']*100:.4f}% "
            f"± {agg['mean_abs_error_sd']*100:.4f}%",
            file=sys.stderr,
        )
        print(
            f"  Mean RMSE: {agg['rmse_mean']*100:.4f}% "
            f"± {agg['rmse_sd']*100:.4f}%",
            file=sys.stderr,
        )

        # Write per-depth facts (mean ± SD across replicates)
        write_fact(f"depth_{depth}", {
            "depth": depth,
            "n_replicates": args.n_replicates,
            "n_samples_per_rep": agg.get("n_replicates", args.n_replicates),
            "mean_abs_error_pct": round(agg["mean_abs_error_mean"] * 100, 2),
            "mean_abs_error_sd_pct": round(agg["mean_abs_error_sd"] * 100, 2),
            "rmse_pct": round(agg["rmse_mean"] * 100, 2),
            "rmse_sd_pct": round(agg["rmse_sd"] * 100, 2),
            "max_abs_error_pct": round(agg["max_abs_error_mean"] * 100, 2),
            "max_abs_error_sd_pct": round(agg["max_abs_error_sd"] * 100, 2),
            "ci_coverage_pct": round(agg["ci_coverage_mean"] * 100, 1),
            "ci_coverage_sd_pct": round(agg["ci_coverage_sd"] * 100, 1),
            "mean_ci_width_pct": round(agg["mean_ci_width_mean"] * 100, 2),
            "mean_ci_width_sd_pct": round(agg["mean_ci_width_sd"] * 100, 2),
        })

    # Write combined summary table as TSV
    depths_sorted = sorted(args.depths)
    summary_rows = []
    for depth in depths_sorted:
        agg = aggregate_replicate_metrics(all_metrics[depth])
        summary_rows.append({
            "depth": depth,
            "n_replicates": args.n_replicates,
            "mae_mean": round(agg["mean_abs_error_mean"] * 100, 2),
            "mae_sd": round(agg["mean_abs_error_sd"] * 100, 2),
            "rmse_mean": round(agg["rmse_mean"] * 100, 2),
            "rmse_sd": round(agg["rmse_sd"] * 100, 2),
            "max_err_mean": round(agg["max_abs_error_mean"] * 100, 2),
            "max_err_sd": round(agg["max_abs_error_sd"] * 100, 2),
            "ci_cov_mean": round(agg["ci_coverage_mean"] * 100, 1),
            "ci_width_mean": round(agg["mean_ci_width_mean"] * 100, 2),
        })

    summary_path = outdir / "depth_summary.tsv"
    with open(summary_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "depth", "n_replicates", "mae_mean", "mae_sd", "rmse_mean", "rmse_sd",
                "max_err_mean", "max_err_sd", "ci_cov_mean", "ci_width_mean",
            ],
            delimiter="\t",
        )
        w.writeheader()
        w.writerows(summary_rows)
    print(f"\nSummary table: {summary_path}", file=sys.stderr)

    # Write per-sample results for each depth and replicate
    for depth in depths_sorted:
        for rep_idx, rep_rows in enumerate(all_results[depth]):
            results_path = outdir / f"results_{depth}x_rep{rep_idx}.tsv"
            with open(results_path, "w", newline="", encoding="utf-8") as f:
                fields = [
                    "sample_name", "true_pct", "est_pct", "error_pct",
                    "ci_lo_pct", "ci_hi_pct", "ci_covers",
                ]
                w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
                w.writeheader()
                for r in rep_rows:
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
    plot_results(all_results, all_metrics, outdir)

    # Copy figures to facts dir
    for fig in ["fig1_depth_scatter.png", "fig2_depth_boxplots.png", "fig3_depth_summary.png"]:
        src = outdir / fig
        if src.exists():
            shutil.copy2(src, FACTS_DIR / fig)
            print(f"  Copied to {FACTS_DIR / fig}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
