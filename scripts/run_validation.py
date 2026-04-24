#!/usr/bin/env python3
"""Run allomix on synthetic test data and produce a validation report.

Reads a truth table (TSV with sample_name and true_donor_fraction), runs
allomix on each sample, and produces:

1. validation_results.tsv — per-sample truth vs estimate with metrics
2. validation_summary.tsv — aggregate accuracy metrics
3. validation_scatter.png — truth vs estimated scatter plot
4. validation_residuals.png — Bland-Altman style residual plot
5. validation_ci.png — CI coverage visualisation

Usage:
    python scripts/run_validation.py \
        --host tests/test_data/host.vcf \
        --donor tests/test_data/donor.vcf \
        --truth tests/test_data/truth_table.tsv \
        --vcf-dir tests/test_data \
        --outdir validation_output

    # Timeline data:
    python scripts/run_validation.py \
        --host tests/test_data/host.vcf \
        --donor tests/test_data/donor.vcf \
        --truth tests/test_data/timeline/truth_table.tsv \
        --vcf-dir tests/test_data/timeline \
        --outdir validation_output_timeline
"""

from __future__ import annotations

import argparse
import csv
import logging
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from allomix.chimerism import estimate_single_donor_bb  # noqa: E402
from allomix.genotype import classify_markers, parse_vcf  # noqa: E402
from allomix.qc import assess_quality  # noqa: E402

log = logging.getLogger(__name__)


def run_sample(
    host_path: str,
    donor_path: str,
    sample_path: str,
    min_dp: int = 0,
    min_gq: int = 0,
    error_rate: float = 0.01,
) -> dict:
    """Run allomix pipeline on one sample. Returns a results dict."""
    host = parse_vcf(host_path, min_dp=0, min_gq=0)
    donor = parse_vcf(donor_path, min_dp=0, min_gq=0)
    admix = parse_vcf(sample_path, min_dp=0, min_gq=0)

    genotypes = classify_markers(
        host, [donor], admix, min_dp=min_dp, min_gq=min_gq, pass_only=False,
    )
    genotypes.sample_name = Path(sample_path).stem

    result = estimate_single_donor_bb(genotypes.informative, error_rate=error_rate)
    qc = assess_quality(result, genotypes)

    return {
        "estimated_donor_fraction": result.donor_fraction,
        "ci_lo": result.donor_fraction_ci[0],
        "ci_hi": result.donor_fraction_ci[1],
        "n_informative": result.n_informative,
        "n_used": result.n_markers_used,
        "log_likelihood": result.log_likelihood,
        "qc_pass": qc.pass_,
    }


def compute_metrics(rows: list[dict]) -> dict:
    """Compute aggregate validation metrics from per-sample results."""
    n = len(rows)
    if n == 0:
        return {}

    errors = [r["error"] for r in rows]
    abs_errors = [abs(e) for e in errors]
    sq_errors = [e * e for e in errors]

    ci_covers = sum(1 for r in rows if r["ci_covers_truth"])
    ci_widths = [r["ci_width"] for r in rows]

    mean_error = sum(errors) / n
    mean_abs_error = sum(abs_errors) / n
    rmse = math.sqrt(sum(sq_errors) / n)
    max_abs_error = max(abs_errors)
    ci_coverage_rate = ci_covers / n
    mean_ci_width = sum(ci_widths) / n

    return {
        "n_samples": n,
        "mean_signed_error_pct": mean_error * 100,
        "mean_abs_error_pct": mean_abs_error * 100,
        "rmse_pct": rmse * 100,
        "max_abs_error_pct": max_abs_error * 100,
        "ci_coverage_rate": ci_coverage_rate,
        "mean_ci_width_pct": mean_ci_width * 100,
    }


def write_results_tsv(rows: list[dict], path: Path) -> None:
    """Write per-sample validation results."""
    fields = [
        "sample_name", "true_donor_pct", "estimated_donor_pct",
        "error_pct", "abs_error_pct", "ci_lo_pct", "ci_hi_pct",
        "ci_width_pct", "ci_covers_truth", "n_informative", "n_used", "qc_pass",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for r in rows:
            writer.writerow({
                "sample_name": r["sample_name"],
                "true_donor_pct": f"{r['true_donor_fraction'] * 100:.2f}",
                "estimated_donor_pct": f"{r['estimated_donor_fraction'] * 100:.2f}",
                "error_pct": f"{r['error'] * 100:.4f}",
                "abs_error_pct": f"{abs(r['error']) * 100:.4f}",
                "ci_lo_pct": f"{r['ci_lo'] * 100:.2f}",
                "ci_hi_pct": f"{r['ci_hi'] * 100:.2f}",
                "ci_width_pct": f"{r['ci_width'] * 100:.2f}",
                "ci_covers_truth": str(r["ci_covers_truth"]),
                "n_informative": r["n_informative"],
                "n_used": r["n_used"],
                "qc_pass": str(r["qc_pass"]),
            })


def write_summary_tsv(metrics: dict, path: Path) -> None:
    """Write aggregate validation summary."""
    with open(path, "w", encoding="utf-8") as f:
        f.write("metric\tvalue\n")
        f.write(f"n_samples\t{metrics['n_samples']}\n")
        f.write(f"mean_signed_error_pct\t{metrics['mean_signed_error_pct']:.4f}\n")
        f.write(f"mean_abs_error_pct\t{metrics['mean_abs_error_pct']:.4f}\n")
        f.write(f"rmse_pct\t{metrics['rmse_pct']:.4f}\n")
        f.write(f"max_abs_error_pct\t{metrics['max_abs_error_pct']:.4f}\n")
        f.write(f"ci_coverage_rate\t{metrics['ci_coverage_rate']:.4f}\n")
        f.write(f"mean_ci_width_pct\t{metrics['mean_ci_width_pct']:.4f}\n")


def try_plot(rows: list[dict], outdir: Path) -> None:
    """Generate validation plots."""
    truths = [r["true_donor_fraction"] * 100 for r in rows]
    estimates = [r["estimated_donor_fraction"] * 100 for r in rows]
    errors = [r["error"] * 100 for r in rows]
    ci_los = [r["ci_lo"] * 100 for r in rows]
    ci_his = [r["ci_hi"] * 100 for r in rows]
    covers = [r["ci_covers_truth"] for r in rows]

    # --- Scatter: truth vs estimated ---
    fig, ax = plt.subplots(figsize=(7, 7))
    ax.scatter(truths, estimates, c="steelblue", s=50, edgecolors="white", linewidth=0.5, zorder=3)
    ax.plot([0, 100], [0, 100], "k--", alpha=0.5, label="Identity")
    ax.set_xlabel("True donor %", fontsize=12)
    ax.set_ylabel("Estimated donor %", fontsize=12)
    ax.set_title("allomix Validation: Truth vs Estimated", fontsize=13)
    ax.set_xlim(-2, 102)
    ax.set_ylim(-2, 102)
    ax.set_aspect("equal")
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / "validation_scatter.png", dpi=150)
    plt.close(fig)

    # --- Residuals (Bland-Altman style) ---
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.scatter(truths, errors, c="steelblue", s=50, edgecolors="white", linewidth=0.5, zorder=3)
    ax.axhline(0, color="k", linestyle="--", alpha=0.5)
    ax.set_xlabel("True donor %", fontsize=12)
    ax.set_ylabel("Error (estimated − true) %", fontsize=12)
    ax.set_title("allomix Validation: Residuals", fontsize=13)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / "validation_residuals.png", dpi=150)
    plt.close(fig)

    # --- CI coverage ---
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, r in enumerate(rows):
        color = "steelblue" if covers[i] else "firebrick"
        ax.plot([truths[i], truths[i]], [ci_los[i], ci_his[i]], color=color, linewidth=2, alpha=0.7)
    ax.scatter(truths, estimates, c="black", s=30, zorder=5, label="Estimate")
    ax.plot([0, 100], [0, 100], "k--", alpha=0.3)
    ax.set_xlabel("True donor %", fontsize=12)
    ax.set_ylabel("Estimated donor % (with 95% CI)", fontsize=12)
    ax.set_title("allomix Validation: Confidence Interval Coverage", fontsize=13)
    n_cover = sum(covers)
    ax.legend(
        [f"CI covers truth ({n_cover}/{len(rows)})", "CI misses truth"],
        fontsize=9,
    )
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(outdir / "validation_ci.png", dpi=150)
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Run allomix validation against known truth.",
    )
    parser.add_argument("--host", required=True, help="Host genotype VCF")
    parser.add_argument("--donor", required=True, help="Donor genotype VCF")
    parser.add_argument("--truth", required=True, help="Truth table TSV")
    parser.add_argument("--vcf-dir", required=True, help="Directory containing chimeric VCFs")
    parser.add_argument(
        "--outdir", default="output/validation",
        help="Output directory for results (default: output/validation)",
    )
    parser.add_argument("--min-dp", type=int, default=0, help="Minimum depth filter")
    parser.add_argument("--min-gq", type=int, default=0, help="Minimum GQ filter")
    parser.add_argument(
        "--error-rate", type=float, default=0.01, help="Sequencing error rate",
    )
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    vcf_dir = Path(args.vcf_dir)

    # Read truth table
    with open(args.truth, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        truth_rows = list(reader)

    rows: list[dict] = []

    for tr in truth_rows:
        sample_name = tr["sample_name"]
        true_frac = float(tr["true_donor_fraction"])
        vcf_path = vcf_dir / f"{sample_name}.vcf"

        if not vcf_path.exists():
            log.warning("%s not found, skipping", vcf_path)
            continue

        log.info("Running %s (truth=%.1f%%) ...", sample_name, true_frac * 100)

        result = run_sample(
            args.host, args.donor, str(vcf_path),
            min_dp=args.min_dp, min_gq=args.min_gq,
            error_rate=args.error_rate,
        )

        error = result["estimated_donor_fraction"] - true_frac
        ci_covers = result["ci_lo"] <= true_frac <= result["ci_hi"]
        ci_width = result["ci_hi"] - result["ci_lo"]

        rows.append({
            "sample_name": sample_name,
            "true_donor_fraction": true_frac,
            "estimated_donor_fraction": result["estimated_donor_fraction"],
            "error": error,
            "ci_lo": result["ci_lo"],
            "ci_hi": result["ci_hi"],
            "ci_width": ci_width,
            "ci_covers_truth": ci_covers,
            "n_informative": result["n_informative"],
            "n_used": result["n_used"],
            "qc_pass": result["qc_pass"],
        })

    if not rows:
        log.error("No samples processed")
        return 1

    # Write per-sample results
    results_path = outdir / "validation_results.tsv"
    write_results_tsv(rows, results_path)
    log.info("Per-sample results: %s", results_path)

    # Compute and write summary
    metrics = compute_metrics(rows)
    summary_path = outdir / "validation_summary.tsv"
    write_summary_tsv(metrics, summary_path)
    log.info("Summary metrics: %s", summary_path)

    # Log summary
    log.info("")
    log.info("=" * 50)
    log.info("VALIDATION SUMMARY")
    log.info("=" * 50)
    log.info("  Samples:            %d", metrics["n_samples"])
    log.info("  Mean signed error:  %+.4f%%", metrics["mean_signed_error_pct"])
    log.info("  Mean abs error:     %.4f%%", metrics["mean_abs_error_pct"])
    log.info("  RMSE:               %.4f%%", metrics["rmse_pct"])
    log.info("  Max abs error:      %.4f%%", metrics["max_abs_error_pct"])
    log.info("  CI coverage rate:   %.1f%%", metrics["ci_coverage_rate"] * 100)
    log.info("  Mean CI width:      %.4f%%", metrics["mean_ci_width_pct"])
    log.info("=" * 50)

    try_plot(rows, outdir)
    log.info("Plots saved to %s/", outdir)

    return 0


if __name__ == "__main__":
    sys.exit(main())
