#!/usr/bin/env python3
"""Compare allomix accuracy with and without bias correction on simulated biased data.

Generates synthetic chimeric VCFs with per-marker capture bias, then runs allomix
both with and without bias correction. Prints a side-by-side comparison of accuracy
metrics to show the effect of bias correction.

Usage:
    python scripts/compare_bias_correction.py
    python scripts/compare_bias_correction.py --bias-sd 0.03 --depth 1000
"""

import argparse
import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

import random  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from allomix.bias import load_bias_table  # noqa: E402
from allomix.chimerism import PanelCalibration, estimate_single_donor_bb  # noqa: E402
from allomix.genotype import classify_markers, parse_vcf  # noqa: E402
from allomix.simulate import (  # noqa: E402
    blend_vcfs,
    generate_marker_biases,
    write_vcf,
)
from allomix.simulate import (  # noqa: E402
    parse_vcf as sim_parse_vcf,
)

FRACTIONS = [0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 0.90, 0.95, 1.0]


def _fraction_to_name(f: float) -> str:
    d = round(f * 100)
    h = 100 - d
    return f"host_{h}_donor_{d}"


def generate_biased_data(
    host_vcf: str,
    donor_vcf: str,
    outdir: Path,
    bias_sd: float,
    depth: int,
    seed: int,
) -> tuple[Path, Path]:
    """Generate biased test data with consistent per-marker biases.

    All fractions share the same per-marker biases (as in real data where
    biases are a property of the capture panel, not the sample).

    Returns (truth_table_path, bias_table_path).
    """
    outdir.mkdir(parents=True, exist_ok=True)

    # Pre-generate panel biases: parse VCFs to count shared markers,
    # then generate biases once with a fixed seed.
    _, host_records = sim_parse_vcf(host_vcf)
    _, donor_records = sim_parse_vcf(donor_vcf)
    donor_loci = {r.locus for r in donor_records}
    n_shared = sum(1 for r in host_records if r.locus in donor_loci)

    bias_rng = random.Random(seed)
    fixed_biases = generate_marker_biases(n_shared, bias_rng, bias_sd)

    truth_rows = []

    # Generate one blend to capture marker identities for the bias table
    first_result = None

    for frac in FRACTIONS:
        name = _fraction_to_name(frac)
        # Each fraction gets a different seed for allele count sampling,
        # but uses the SAME fixed_biases for panel capture bias.
        sample_seed = seed + hash(str(frac)) % (2**31)

        result = blend_vcfs(
            host_path=host_vcf,
            donor_path=donor_vcf,
            donor_fraction=frac,
            target_depth=depth,
            sample_name=name,
            seed=sample_seed,
            fixed_biases=fixed_biases,
            locus_dropout_rate=0.016,
            depth_cv=0.43,
        )
        write_vcf(result, outdir / f"{name}.vcf")

        if first_result is None:
            first_result = result

        truth_rows.append({
            "sample_name": name,
            "true_donor_fraction": f"{frac:.6f}",
        })

    # Write truth table
    truth_path = outdir / "truth_table.tsv"
    with open(truth_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, ["sample_name", "true_donor_fraction"], delimiter="\t")
        w.writeheader()
        w.writerows(truth_rows)

    # Write bias table from the first result's marker_biases
    bias_path = outdir / "true_biases.tsv"
    with open(bias_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["chrom", "pos", "ref", "alt", "bias", "n_het"])
        if first_result and first_result.marker_biases:
            for chrom, pos, ref, alt, bias in first_result.marker_biases:
                w.writerow([chrom, pos, ref, alt, f"{bias:.6f}", 1])

    return truth_path, bias_path


def run_validation(
    host_vcf: str,
    donor_vcf: str,
    truth_path: Path,
    vcf_dir: Path,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
) -> list[dict]:
    """Run allomix on all samples and compute error metrics."""
    with open(truth_path, encoding="utf-8") as f:
        truth_rows = list(csv.DictReader(f, delimiter="\t"))

    rows = []
    for tr in truth_rows:
        name = tr["sample_name"]
        true_frac = float(tr["true_donor_fraction"])
        vcf_path = vcf_dir / f"{name}.vcf"

        host = parse_vcf(host_vcf, min_dp=0, min_gq=0)
        donor = parse_vcf(donor_vcf, min_dp=0, min_gq=0)
        admix = parse_vcf(str(vcf_path), min_dp=0, min_gq=0)

        genotypes = classify_markers(host, [donor], admix, min_dp=0, min_gq=0, pass_only=False)
        result = estimate_single_donor_bb(
            genotypes.informative,
            error_rate=0.01,
            calibration=PanelCalibration(biases=marker_biases or {}),
        )

        error = result.donor_fraction - true_frac
        ci_covers = result.donor_fraction_ci[0] <= true_frac <= result.donor_fraction_ci[1]
        ci_width = result.donor_fraction_ci[1] - result.donor_fraction_ci[0]

        rows.append({
            "sample_name": name,
            "true_frac": true_frac,
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
    """Compute aggregate metrics."""
    n = len(rows)
    if n == 0:
        return {}

    # Exclude boundary fractions (0% and 100%) from error metrics
    # since they have no room for bias to manifest
    interior = [r for r in rows if 0.0 < r["true_frac"] < 1.0]
    ni = len(interior)

    errors = [r["error"] for r in interior]
    abs_errors = [abs(e) for e in errors]
    sq_errors = [e * e for e in errors]

    ci_covers = sum(1 for r in rows if r["ci_covers"])
    ci_widths = [r["ci_width"] for r in rows]

    return {
        "n_samples": n,
        "n_interior": ni,
        "mean_signed_error": sum(errors) / ni if ni else 0,
        "mean_abs_error": sum(abs_errors) / ni if ni else 0,
        "rmse": math.sqrt(sum(sq_errors) / ni) if ni else 0,
        "max_abs_error": max(abs_errors) if abs_errors else 0,
        "ci_coverage": ci_covers / n if n else 0,
        "mean_ci_width": sum(ci_widths) / n if n else 0,
    }


def print_comparison(m_no: dict, m_yes: dict, bias_sd: float) -> None:
    """Print side-by-side comparison table."""
    print(f"\n{'=' * 70}")
    print(f"BIAS CORRECTION COMPARISON  (bias_sd = {bias_sd})")
    print(f"{'=' * 70}")
    print(f"{'Metric':<30} {'No correction':>18} {'With correction':>18}")
    print(f"{'-' * 30} {'-' * 18} {'-' * 18}")

    def row(label, key, fmt=".4f", scale=100, unit="%"):
        v1 = m_no[key] * scale
        v2 = m_yes[key] * scale
        delta = v2 - v1
        sign = "+" if delta >= 0 else ""
        print(f"{label:<30} {v1:>17{fmt}}{unit} {v2:>17{fmt}}{unit}  ({sign}{delta:{fmt}}{unit})")

    row("Mean signed error", "mean_signed_error")
    row("Mean absolute error", "mean_abs_error")
    row("RMSE", "rmse")
    row("Max absolute error", "max_abs_error")

    # CI coverage as percentage, not *100
    v1 = m_no["ci_coverage"]
    v2 = m_yes["ci_coverage"]
    delta = v2 - v1
    sign = "+" if delta >= 0 else ""
    print(f"{'CI coverage rate':<30} {v1:>17.1%}  {v2:>17.1%}   ({sign}{delta:.1%})")

    row("Mean CI width", "mean_ci_width")
    print(f"{'=' * 70}")


def print_per_sample(rows_no: list[dict], rows_yes: list[dict]) -> None:
    """Print per-sample comparison."""
    print(f"\n{'Sample':<20} {'Truth':>6} {'No bias':>8} {'Err':>7} "
          f"{'Corrected':>10} {'Err':>7} {'CI(no)':>12} {'CI(yes)':>12}")
    print("-" * 100)
    for rn, ry in zip(rows_no, rows_yes):
        t = rn["true_frac"] * 100
        en = rn["est_frac"] * 100
        ey = ry["est_frac"] * 100
        err_n = rn["error"] * 100
        err_y = ry["error"] * 100
        ci_n = f"[{rn['ci_lo']*100:.2f},{rn['ci_hi']*100:.2f}]"
        ci_y = f"[{ry['ci_lo']*100:.2f},{ry['ci_hi']*100:.2f}]"
        cover_n = "*" if rn["ci_covers"] else " "
        cover_y = "*" if ry["ci_covers"] else " "
        print(f"{rn['sample_name']:<20} {t:>5.1f}% {en:>7.2f}% {err_n:>+6.2f}% "
              f"{ey:>9.2f}% {err_y:>+6.2f}% {ci_n:>12}{cover_n} {ci_y:>12}{cover_y}")
    print("\n  * = CI covers truth")


def try_plot(rows_no: list[dict], rows_yes: list[dict], outdir: Path, bias_sd: float) -> None:
    """Generate comparison plot."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    truths = [r["true_frac"] * 100 for r in rows_no]
    est_no = [r["est_frac"] * 100 for r in rows_no]
    est_yes = [r["est_frac"] * 100 for r in rows_yes]
    err_no = [r["error"] * 100 for r in rows_no]
    err_yes = [r["error"] * 100 for r in rows_yes]

    # Scatter: truth vs estimated
    ax = axes[0]
    ax.scatter(truths, est_no, c="firebrick", s=50, label="No correction", alpha=0.7, zorder=3)
    ax.scatter(truths, est_yes, c="steelblue", s=50, label="Bias-corrected", alpha=0.7, zorder=3)
    ax.plot([0, 100], [0, 100], "k--", alpha=0.4)
    ax.set_xlabel("True donor %")
    ax.set_ylabel("Estimated donor %")
    ax.set_title("Truth vs Estimated")
    ax.legend(fontsize=9)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.3)

    # Residuals
    ax = axes[1]
    ax.scatter(truths, err_no, c="firebrick", s=50, label="No correction", alpha=0.7, zorder=3)
    ax.scatter(truths, err_yes, c="steelblue", s=50, label="Bias-corrected", alpha=0.7, zorder=3)
    ax.axhline(0, color="k", linestyle="--", alpha=0.4)
    ax.set_xlabel("True donor %")
    ax.set_ylabel("Error (est - true) %")
    ax.set_title("Residuals")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # CI coverage
    ax = axes[2]
    for i, r in enumerate(rows_no):
        color = "firebrick"
        alpha = 0.5
        ax.plot([truths[i] - 0.3, truths[i] - 0.3],
                [r["ci_lo"] * 100, r["ci_hi"] * 100], color=color, linewidth=2, alpha=alpha)
    for i, r in enumerate(rows_yes):
        color = "steelblue"
        alpha = 0.5
        ax.plot([truths[i] + 0.3, truths[i] + 0.3],
                [r["ci_lo"] * 100, r["ci_hi"] * 100], color=color, linewidth=2, alpha=alpha)
    ax.plot([0, 100], [0, 100], "k--", alpha=0.3)
    ax.set_xlabel("True donor %")
    ax.set_ylabel("Estimated % with CI")
    ax.set_title("Confidence Intervals")
    ax.grid(True, alpha=0.3)

    fig.suptitle(f"Bias Correction Comparison (bias_sd = {bias_sd})", fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(outdir / "bias_correction_comparison.png", dpi=150, bbox_inches="tight")
    plt.close(fig)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Compare allomix with and without bias correction on simulated data.",
    )
    parser.add_argument(
        "--host", default="tests/test_data/host.vcf",
        help="Host genotype VCF (default: tests/test_data/host.vcf)",
    )
    parser.add_argument(
        "--donor", default="tests/test_data/donor.vcf",
        help="Donor genotype VCF (default: tests/test_data/donor.vcf)",
    )
    parser.add_argument(
        "--bias-sd", type=float, default=0.02,
        help="Per-marker capture bias SD (default: 0.02)",
    )
    parser.add_argument("--depth", type=int, default=2000, help="Target depth (default: 2000)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed (default: 42)")
    parser.add_argument(
        "--outdir", default="output/bias_comparison",
        help="Output directory (default: output/bias_comparison)",
    )
    args = parser.parse_args(argv)

    outdir = Path(args.outdir)

    print(f"Generating biased test data (bias_sd={args.bias_sd}, "
          f"depth={args.depth}) ...", file=sys.stderr)

    truth_path, bias_path = generate_biased_data(
        args.host, args.donor, outdir / "vcfs",
        bias_sd=args.bias_sd, depth=args.depth, seed=args.seed,
    )

    # Load the true bias table
    biases = load_bias_table(bias_path)
    print(f"  Bias table: {len(biases)} markers", file=sys.stderr)

    # Run without bias correction
    print("Running allomix WITHOUT bias correction ...", file=sys.stderr)
    rows_no = run_validation(
        args.host, args.donor, truth_path, outdir / "vcfs",
        marker_biases=None,
    )

    # Run with bias correction
    print("Running allomix WITH bias correction ...", file=sys.stderr)
    rows_yes = run_validation(
        args.host, args.donor, truth_path, outdir / "vcfs",
        marker_biases=biases,
    )

    # Compute and display metrics
    m_no = compute_metrics(rows_no)
    m_yes = compute_metrics(rows_yes)

    print_per_sample(rows_no, rows_yes)
    print_comparison(m_no, m_yes, args.bias_sd)

    outdir.mkdir(parents=True, exist_ok=True)
    try_plot(rows_no, rows_yes, outdir, args.bias_sd)
    print(f"\nPlot saved to {outdir}/bias_correction_comparison.png", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
