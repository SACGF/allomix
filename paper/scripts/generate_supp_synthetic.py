#!/usr/bin/env python3
"""Generate supplementary figures for synthetic data and bias model validation.

Produces Figures S1-S6 showing:
  - Empirical vs simulated bias distributions
  - Depth distribution comparison
  - Het VAF comparison
  - Noise component ablation study
  - CI calibration across replicates
  - Per-marker residual analysis

Usage:
    python paper/scripts/generate_supp_synthetic.py
    python paper/scripts/generate_supp_synthetic.py --calibration-batch 0 10
"""

import argparse
import csv
import math
import random
import sys
import tempfile
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from scipy.stats import gaussian_kde, norm  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from allomix.chimerism import estimate_single_donor_bb  # noqa: E402
from allomix.genotype import classify_markers, parse_vcf  # noqa: E402
from allomix.simulate import (  # noqa: E402
    blend_vcfs,
    generate_marker_biases,
    generate_marker_biases_realistic,
    sample_marker_depths,
    write_vcf,
)

# --- Config ---
HOST_VCF = "tests/test_data/host.vcf"
DONOR_VCF = "tests/test_data/donor.vcf"
EMPIRICAL_PER_MARKER = "paper/empirical_results/panel_per_marker.tsv"
FACTS_DIR = Path("output/facts")
CALIBRATION_DIR = Path("output/calibration")
FRACTIONS = [0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 0.90, 0.95, 1.0]
N_ABLATION_REPS = 10
N_CALIBRATION_REPS = 100
N_CALIBRATION_BATCHES = 10

# Ablation conditions: each defines blend_vcfs kwargs and whether to pass
# the true bias table to the estimator (simulating bias correction).
CONDITIONS = {
    "Ideal": {
        "blend": {
            "marker_bias_sd": 0.0,
            "depth_cv": 0.0,
            "error_rate": 0.0,
            "locus_dropout_rate": 0.0,
            "realistic_biases": False,
        },
        "correct_bias": False,
    },
    "Bias only": {
        "blend": {
            "marker_bias_sd": 0.0,
            "depth_cv": 0.0,
            "error_rate": 0.0,
            "locus_dropout_rate": 0.0,
            "realistic_biases": True,
        },
        "correct_bias": False,
    },
    "Bias corrected": {
        "blend": {
            "marker_bias_sd": 0.0,
            "depth_cv": 0.0,
            "error_rate": 0.0,
            "locus_dropout_rate": 0.0,
            "realistic_biases": True,
        },
        "correct_bias": True,
    },
    "Depth only": {
        "blend": {
            "marker_bias_sd": 0.0,
            "depth_cv": 0.43,
            "error_rate": 0.0,
            "locus_dropout_rate": 0.0,
            "realistic_biases": False,
        },
        "correct_bias": False,
    },
    "Error only": {
        "blend": {
            "marker_bias_sd": 0.0,
            "depth_cv": 0.0,
            "error_rate": 0.01,
            "locus_dropout_rate": 0.0,
            "realistic_biases": False,
        },
        "correct_bias": False,
    },
    "Full": {
        "blend": {
            "marker_bias_sd": 0.0,
            "depth_cv": 0.43,
            "error_rate": 0.01,
            "locus_dropout_rate": 0.016,
            "realistic_biases": True,
        },
        "correct_bias": True,
    },
}


def load_empirical_per_marker() -> list[dict]:
    """Load per-marker stats from panel_per_marker.tsv."""
    rows = []
    with open(EMPIRICAL_PER_MARKER, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(
                {
                    "median_bias": float(row["median_bias"]),
                    "mean_bias": float(row["mean_bias"]),
                    "sd_within": float(row["sd_within"]),
                    "mean_depth": float(row["mean_depth"]),
                    "depth_cv": float(row["depth_cv"]),
                    "call_rate": float(row["call_rate"]),
                    "het_ratio_vs_hwe": float(row["het_ratio_vs_hwe"]),
                }
            )
    return rows


def _build_bias_table(
    blend_result,
) -> dict[tuple[str, int, str, str], float] | None:
    """Convert BlendResult.marker_biases into the dict format expected by
    estimate_single_donor_bb(marker_biases=...).
    """
    if not blend_result.marker_biases:
        return None
    return {
        (chrom, pos, ref, alt): bias for chrom, pos, ref, alt, bias in blend_result.marker_biases
    }


# ---------------------------------------------------------------------------
# Figure S1: Empirical vs simulated bias distributions
# ---------------------------------------------------------------------------


def plot_bias_distributions(empirical_biases, ax_hist, ax_cdf, rng, n_draws=10_000):
    """Figure S1: empirical vs simulated bias distributions."""
    sim_simple = generate_marker_biases(n_draws, rng, bias_sd=0.0175)
    sim_mixture = generate_marker_biases_realistic(n_draws, rng)

    # Panel A: histogram + KDE
    bins = np.linspace(-0.12, 0.12, 40)
    ax_hist.hist(
        empirical_biases,
        bins=bins,
        density=True,
        alpha=0.5,
        color="0.4",
        label="Empirical (71 markers)",
        edgecolor="white",
    )
    xs = np.linspace(-0.12, 0.12, 300)
    kde_simple = gaussian_kde(sim_simple)
    kde_mixture = gaussian_kde(sim_mixture)
    ax_hist.plot(xs, kde_simple(xs), color="firebrick", linewidth=2, label="Gaussian (SD=0.0175)")
    ax_hist.plot(xs, kde_mixture(xs), color="steelblue", linewidth=2, label="Mixture model")
    # Vertical dashed lines at empirical P5, P95
    p5, p95 = np.percentile(empirical_biases, [5, 95])
    ax_hist.axvline(p5, color="0.6", linestyle=":", alpha=0.5)
    ax_hist.axvline(p95, color="0.6", linestyle=":", alpha=0.5)
    ax_hist.set_xlabel("Per-marker bias (VAF shift)")
    ax_hist.set_ylabel("Density")
    ax_hist.set_title("A", fontsize=12, fontweight="bold", loc="left")
    ax_hist.legend(fontsize=9)

    # Panel B: CDF of |bias|
    emp_abs = sorted(abs(b) for b in empirical_biases)
    sim_simple_abs = sorted(abs(b) for b in sim_simple)
    sim_mixture_abs = sorted(abs(b) for b in sim_mixture)
    ax_cdf.plot(
        emp_abs,
        np.linspace(0, 1, len(emp_abs)),
        color="0.3",
        linewidth=2.5,
        label="Empirical",
    )
    ax_cdf.plot(
        sim_simple_abs,
        np.linspace(0, 1, len(sim_simple_abs)),
        color="firebrick",
        linewidth=1.5,
        alpha=0.7,
        label="Gaussian",
    )
    ax_cdf.plot(
        sim_mixture_abs,
        np.linspace(0, 1, len(sim_mixture_abs)),
        color="steelblue",
        linewidth=1.5,
        alpha=0.7,
        label="Mixture",
    )
    ax_cdf.axvline(0.03, color="0.6", linestyle=":", alpha=0.5)
    ax_cdf.set_xlabel("|Bias|")
    ax_cdf.set_ylabel("Cumulative proportion")
    ax_cdf.set_title("B", fontsize=12, fontweight="bold", loc="left")
    ax_cdf.legend(fontsize=9)


# ---------------------------------------------------------------------------
# Figure S2: Depth distribution
# ---------------------------------------------------------------------------


def plot_depth_distributions(marker_depths, marker_depth_cvs, ax_hist, ax_cv):
    """Figure S2: empirical depth distribution vs log-normal model."""
    rng = random.Random(42)
    sim_depths = sample_marker_depths(10_000, mean_depth=1732, depth_cv=0.429, rng=rng)

    # Panel A: depth histogram
    bins = np.linspace(0, 4000, 30)
    ax_hist.hist(
        marker_depths,
        bins=bins,
        density=True,
        alpha=0.5,
        color="0.4",
        label="Empirical (71 markers)",
        edgecolor="white",
    )
    kde = gaussian_kde([float(d) for d in sim_depths])
    xs = np.linspace(0, 4000, 300)
    ax_hist.plot(xs, kde(xs), color="steelblue", linewidth=2, label="Log-normal model")
    ax_hist.set_xlabel("Mean depth per marker")
    ax_hist.set_ylabel("Density")
    ax_hist.set_title("A", fontsize=12, fontweight="bold", loc="left")
    ax_hist.legend(fontsize=9)

    # Panel B: per-marker depth CV (sorted)
    sorted_cvs = sorted(marker_depth_cvs)
    ax_cv.bar(range(len(sorted_cvs)), sorted_cvs, color="steelblue", alpha=0.6)
    ax_cv.axhline(
        np.mean(marker_depth_cvs),
        color="firebrick",
        linestyle="--",
        linewidth=1.5,
        label=f"Mean CV = {np.mean(marker_depth_cvs):.3f}",
    )
    ax_cv.set_xlabel("Marker (sorted by CV)")
    ax_cv.set_ylabel("Within-marker depth CV")
    ax_cv.set_title("B", fontsize=12, fontweight="bold", loc="left")
    ax_cv.legend(fontsize=9)


# ---------------------------------------------------------------------------
# Figure S3: Het VAF violin (supplementary)
# ---------------------------------------------------------------------------


def plot_het_vaf_violin(empirical_median_bias, ax_violin):
    """Figure S3: het VAF distribution, empirical vs simulated."""
    rng = random.Random(42)
    n = len(empirical_median_bias)
    sim_biases = generate_marker_biases_realistic(n, rng)

    emp_het_vafs = [0.5 + b for b in empirical_median_bias]
    sim_het_vafs = [0.5 + b for b in sim_biases]

    parts = ax_violin.violinplot(
        [emp_het_vafs, sim_het_vafs], positions=[1, 2], showmeans=True, showmedians=True
    )
    # Style the violin bodies
    for pc in parts["bodies"]:
        pc.set_facecolor("steelblue")
        pc.set_alpha(0.5)
    ax_violin.set_xticks([1, 2])
    ax_violin.set_xticklabels(["Empirical", "Simulated"])
    ax_violin.set_ylabel("Median het VAF per marker")
    ax_violin.axhline(0.5, color="0.6", linestyle="--", alpha=0.5)


# ---------------------------------------------------------------------------
# Main paper figure: Bias stability
# ---------------------------------------------------------------------------


def plot_bias_stability(empirical_median_bias, empirical_sd_within, ax):
    """Main paper figure: per-marker bias magnitude vs within-marker VAF SD.

    Validates the fixed-bias-per-marker assumption. If within-marker SD is
    roughly constant regardless of bias magnitude, the fixed-bias model is
    appropriate.
    """
    abs_bias = np.abs(empirical_median_bias)
    ax.scatter(
        abs_bias,
        empirical_sd_within,
        s=30,
        alpha=0.7,
        color="steelblue",
        edgecolors="white",
        linewidth=0.5,
    )
    ax.set_xlabel("|Median bias| per marker")
    ax.set_ylabel("Within-marker VAF SD")
    r = np.corrcoef(abs_bias, empirical_sd_within)[0, 1]
    ax.set_title(f"r = {r:.2f}", fontsize=10)
    ax.grid(True, alpha=0.2)


# ---------------------------------------------------------------------------
# Figure S4: Noise component ablation study
# ---------------------------------------------------------------------------


def run_ablation(host_vcf: str, donor_vcf: str, tmpdir: Path) -> dict[str, list[dict]]:
    """Run all conditions x fractions x replicates, return results dict."""
    results: dict[str, list[dict]] = {}

    for cond_name, cond in CONDITIONS.items():
        cond_results = []
        print(f"  Ablation: {cond_name}", file=sys.stderr)
        for rep in range(N_ABLATION_REPS):
            seed = 1000 * rep + 1
            for frac in FRACTIONS:
                frac_seed = seed + hash(str(frac)) % (2**31)
                blend = blend_vcfs(
                    host_vcf,
                    donor_vcf,
                    frac,
                    target_depth=500,
                    seed=frac_seed,
                    **cond["blend"],
                )
                vcf_path = tmpdir / f"{cond_name}_{rep}_{frac}.vcf"
                write_vcf(blend, vcf_path)

                bias_table = _build_bias_table(blend) if cond["correct_bias"] else None

                host = parse_vcf(host_vcf, min_dp=0, min_gq=0)
                donor = parse_vcf(donor_vcf, min_dp=0, min_gq=0)
                admix = parse_vcf(str(vcf_path), min_dp=0, min_gq=0)
                markers = classify_markers(
                    host, [donor], admix, min_dp=0, min_gq=0, pass_only=False
                )
                est = estimate_single_donor_bb(
                    markers.informative, error_rate=0.01, marker_biases=bias_table
                )

                cond_results.append(
                    {
                        "true_frac": frac,
                        "est_frac": est.donor_fraction,
                        "error": est.donor_fraction - frac,
                        "replicate": rep,
                    }
                )
        results[cond_name] = cond_results

    return results


def plot_ablation(results, ax_rmse, ax_per_frac):
    """Figure S4: ablation study."""
    colors = {
        "Ideal": "#4CAF50",
        "Bias only": "#F44336",
        "Bias corrected": "#E91E63",
        "Depth only": "#FF9800",
        "Error only": "#9C27B0",
        "Full": "#2196F3",
    }
    linestyles = {c: "--" if "corrected" in c.lower() else "-" for c in CONDITIONS}

    # Panel A: RMSE per condition
    rmses = {}
    for cond, rows in results.items():
        interior = [r for r in rows if 0.0 < r["true_frac"] < 1.0]
        rmse = math.sqrt(sum(r["error"] ** 2 for r in interior) / len(interior))
        rmses[cond] = rmse

    x_pos = range(len(rmses))
    ax_rmse.bar(
        x_pos,
        [rmses[c] * 100 for c in CONDITIONS],
        color=[colors[c] for c in CONDITIONS],
        alpha=0.7,
    )
    ax_rmse.set_xticks(list(x_pos))
    ax_rmse.set_xticklabels(CONDITIONS.keys(), rotation=35, ha="right", fontsize=8)
    ax_rmse.set_ylabel("RMSE (%)")
    ax_rmse.set_title("A", fontsize=12, fontweight="bold", loc="left")
    ax_rmse.grid(True, alpha=0.2, axis="y")

    # Panel B: per-fraction mean |error|
    for cond in CONDITIONS:
        rows = results[cond]
        fracs = sorted(set(r["true_frac"] for r in rows))
        mean_abs = []
        for f in fracs:
            f_rows = [r for r in rows if r["true_frac"] == f]
            mean_abs.append(np.mean([abs(r["error"]) for r in f_rows]) * 100)
        ax_per_frac.plot(
            [f * 100 for f in fracs],
            mean_abs,
            "o-",
            color=colors[cond],
            linestyle=linestyles[cond],
            label=cond,
            markersize=5,
            linewidth=1.5,
        )

    ax_per_frac.set_xlabel("True donor fraction (%)")
    ax_per_frac.set_ylabel("Mean |error| (%)")
    ax_per_frac.set_title("B", fontsize=12, fontweight="bold", loc="left")
    ax_per_frac.legend(fontsize=8, ncol=2)
    ax_per_frac.grid(True, alpha=0.2)


# ---------------------------------------------------------------------------
# Figure S5: CI calibration
# ---------------------------------------------------------------------------


def run_calibration_batch(
    host_vcf: str,
    donor_vcf: str,
    batch_idx: int,
    reps_per_batch: int,
    outpath: str | Path,
) -> None:
    """Run one batch of calibration replicates, write results to CSV.

    Called via: python generate_supp_synthetic.py --calibration-batch 0 10
    """
    params = CONDITIONS["Full"]["blend"]
    rows = []
    start_rep = batch_idx * reps_per_batch

    for rep in range(start_rep, start_rep + reps_per_batch):
        seed = 5000 + rep
        print(f"    Calibration rep {rep}", file=sys.stderr)
        with tempfile.TemporaryDirectory() as tmpdir:
            for frac in FRACTIONS:
                frac_seed = seed + hash(str(frac)) % (2**31)
                blend = blend_vcfs(
                    host_vcf,
                    donor_vcf,
                    frac,
                    target_depth=500,
                    seed=frac_seed,
                    **params,
                )
                vcf_path = Path(tmpdir) / f"cal_{rep}_{frac}.vcf"
                write_vcf(blend, vcf_path)

                bias_table = _build_bias_table(blend)

                host = parse_vcf(host_vcf, min_dp=0, min_gq=0)
                donor = parse_vcf(donor_vcf, min_dp=0, min_gq=0)
                admix = parse_vcf(str(vcf_path), min_dp=0, min_gq=0)
                markers = classify_markers(
                    host, [donor], admix, min_dp=0, min_gq=0, pass_only=False
                )
                est = estimate_single_donor_bb(
                    markers.informative, error_rate=0.01, marker_biases=bias_table
                )

                ci_lo, ci_hi = est.donor_fraction_ci
                rows.append(
                    {
                        "true_frac": frac,
                        "est_frac": est.donor_fraction,
                        "ci_lo": ci_lo,
                        "ci_hi": ci_hi,
                        "ci_width": ci_hi - ci_lo,
                        "ci_covers": ci_lo <= frac <= ci_hi,
                        "replicate": rep,
                    }
                )

    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with open(outpath, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def load_calibration_batches() -> list[dict]:
    """Merge all batch CSVs from output/calibration/ into one list."""
    results = []
    for p in sorted(CALIBRATION_DIR.glob("batch_*.csv")):
        with open(p, encoding="utf-8") as f:
            for row in csv.DictReader(f):
                results.append(
                    {
                        "true_frac": float(row["true_frac"]),
                        "est_frac": float(row["est_frac"]),
                        "ci_lo": float(row["ci_lo"]),
                        "ci_hi": float(row["ci_hi"]),
                        "ci_width": float(row["ci_width"]),
                        "ci_covers": row["ci_covers"] == "True",
                        "replicate": int(row["replicate"]),
                    }
                )
    return results


def plot_calibration(results, ax_coverage, ax_width):
    """Figure S5: CI calibration."""
    fracs = sorted(set(r["true_frac"] for r in results))

    coverages = []
    widths_mean = []
    widths_sd = []

    for f in fracs:
        f_rows = [r for r in results if r["true_frac"] == f]
        cov = sum(r["ci_covers"] for r in f_rows) / len(f_rows) * 100
        coverages.append(cov)
        ws = [r["ci_width"] * 100 for r in f_rows]
        widths_mean.append(np.mean(ws))
        widths_sd.append(np.std(ws))

    frac_pcts = [f * 100 for f in fracs]

    # Panel A: coverage
    ax_coverage.plot(frac_pcts, coverages, "o-", color="steelblue", markersize=6, linewidth=1.5)
    ax_coverage.axhline(95, color="firebrick", linestyle="--", linewidth=1, alpha=0.7)
    ax_coverage.fill_between(frac_pcts, 90, 100, color="firebrick", alpha=0.05)
    ax_coverage.set_xlabel("True donor fraction (%)")
    ax_coverage.set_ylabel("CI coverage rate (%)")
    ax_coverage.set_ylim(70, 105)
    ax_coverage.set_title("A", fontsize=12, fontweight="bold", loc="left")
    ax_coverage.grid(True, alpha=0.2)

    # Panel B: CI width
    ax_width.errorbar(
        frac_pcts,
        widths_mean,
        yerr=widths_sd,
        fmt="o-",
        color="steelblue",
        capsize=3,
        markersize=5,
        linewidth=1.5,
    )
    ax_width.set_xlabel("True donor fraction (%)")
    ax_width.set_ylabel("95% CI width (%)")
    ax_width.set_title("B", fontsize=12, fontweight="bold", loc="left")
    ax_width.grid(True, alpha=0.2)


# ---------------------------------------------------------------------------
# Figure S6: Per-marker residuals
# ---------------------------------------------------------------------------


def plot_residuals(host_vcf, donor_vcf, ax_hist, ax_scatter, tmpdir):
    """Figure S6: per-marker residuals under full simulation."""
    frac = 0.30
    result = blend_vcfs(
        host_vcf,
        donor_vcf,
        frac,
        target_depth=500,
        seed=42,
        realistic_biases=True,
        depth_cv=0.43,
        error_rate=0.01,
        locus_dropout_rate=0.016,
    )
    vcf_path = tmpdir / "residual_check.vcf"
    write_vcf(result, vcf_path)

    host = parse_vcf(host_vcf, min_dp=0, min_gq=0)
    donor = parse_vcf(donor_vcf, min_dp=0, min_gq=0)
    admix = parse_vcf(str(vcf_path), min_dp=0, min_gq=0)
    markers = classify_markers(host, [donor], admix, min_dp=0, min_gq=0, pass_only=False)

    expected_vafs = []
    observed_vafs = []
    for m in markers.informative:
        host_alt_dose = m.host_gt[0] + m.host_gt[1]
        donor_alt_dose = m.donor_gts[0][0] + m.donor_gts[0][1]
        exp_vaf = ((1 - frac) * host_alt_dose + frac * donor_alt_dose) / 2.0
        total = m.admix_ad_ref + m.admix_ad_alt
        obs_vaf = m.admix_ad_alt / total if total > 0 else 0
        expected_vafs.append(exp_vaf)
        observed_vafs.append(obs_vaf)

    residuals = [o - e for o, e in zip(observed_vafs, expected_vafs)]

    # Panel A: histogram
    ax_hist.hist(
        residuals,
        bins=30,
        density=True,
        alpha=0.5,
        color="steelblue",
        edgecolor="white",
        label="Observed residuals",
    )
    xs = np.linspace(-0.15, 0.15, 200)
    sd_residual = np.std(residuals)
    ax_hist.plot(
        xs,
        norm.pdf(xs, 0, sd_residual),
        color="firebrick",
        linewidth=2,
        label=f"Normal fit (SD={sd_residual:.3f})",
    )
    ax_hist.set_xlabel("Residual (observed - expected VAF)")
    ax_hist.set_ylabel("Density")
    ax_hist.set_title("A", fontsize=12, fontweight="bold", loc="left")
    ax_hist.legend(fontsize=9)

    # Panel B: residuals vs expected VAF
    ax_scatter.scatter(
        expected_vafs,
        residuals,
        s=25,
        alpha=0.6,
        color="steelblue",
        edgecolors="white",
        linewidth=0.5,
    )
    ax_scatter.axhline(0, color="0.4", linestyle="--", linewidth=1)
    ax_scatter.set_xlabel("Expected VAF")
    ax_scatter.set_ylabel("Residual")
    ax_scatter.set_title("B", fontsize=12, fontweight="bold", loc="left")


# ---------------------------------------------------------------------------
# Facts CSV
# ---------------------------------------------------------------------------


def write_supp_facts(ablation_results, cal_results, empirical_biases, empirical_sd_within):
    """Write summary facts CSV for supplementary text template variables."""
    facts = {}

    # Ablation RMSE per condition
    for cond, rows in ablation_results.items():
        interior = [r for r in rows if 0.0 < r["true_frac"] < 1.0]
        rmse = math.sqrt(sum(r["error"] ** 2 for r in interior) / len(interior))
        key = cond.lower().replace(" ", "_")
        facts[f"ablation_rmse_{key}_pct"] = f"{rmse * 100:.3f}"

    # Calibration overall coverage
    if cal_results:
        n_covered = sum(1 for r in cal_results if r["ci_covers"])
        facts["cal_coverage_pct"] = f"{n_covered / len(cal_results) * 100:.1f}"

    # Empirical bias summary
    abs_biases = [abs(b) for b in empirical_biases]
    facts["n_empirical_markers"] = str(len(empirical_biases))
    facts["empirical_p95_abs_bias"] = f"{np.percentile(abs_biases, 95):.4f}"

    # Bias stability correlation
    r = np.corrcoef(np.abs(empirical_biases), empirical_sd_within)[0, 1]
    facts["bias_stability_r"] = f"{r:.2f}"

    path = FACTS_DIR / "supp_synthetic.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(facts.keys()))
        writer.writeheader()
        writer.writerow(facts)
    print(f"  Wrote {path}", file=sys.stderr)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    """Entry point. Supports two modes:

    1. Default (no args): generate all figures (S1-S4, S6, main-paper bias
       stability) and merge pre-computed calibration batches for S5.
    2. --calibration-batch BATCH_IDX REPS_PER_BATCH: run one calibration
       batch and write its CSV. Called by Snakemake in parallel.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--calibration-batch",
        nargs=2,
        type=int,
        metavar=("BATCH", "REPS"),
        help="Run calibration batch BATCH with REPS replicates, then exit.",
    )
    args = parser.parse_args()

    if args.calibration_batch is not None:
        batch_idx, reps = args.calibration_batch
        outpath = CALIBRATION_DIR / f"batch_{batch_idx}.csv"
        print(f"Running calibration batch {batch_idx} ({reps} reps)...", file=sys.stderr)
        run_calibration_batch(HOST_VCF, DONOR_VCF, batch_idx, reps, outpath)
        print(f"  Wrote {outpath}", file=sys.stderr)
        return

    FACTS_DIR.mkdir(parents=True, exist_ok=True)
    markers = load_empirical_per_marker()
    rng = random.Random(42)
    empirical_biases = [m["median_bias"] for m in markers]

    # --- Figure S1: Bias distributions ---
    print("Generating Figure S1: Bias distributions...", file=sys.stderr)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    plot_bias_distributions(empirical_biases, ax1, ax2, rng)
    fig.tight_layout()
    fig.savefig(FACTS_DIR / "figS1_bias_distributions.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Figure S2: Depth distributions ---
    print("Generating Figure S2: Depth distributions...", file=sys.stderr)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    marker_depths = [m["mean_depth"] for m in markers]
    marker_depth_cvs = [m["depth_cv"] for m in markers]
    plot_depth_distributions(marker_depths, marker_depth_cvs, ax1, ax2)
    fig.tight_layout()
    fig.savefig(FACTS_DIR / "figS2_depth_distributions.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Figure S3: Het VAF violin (supplementary) ---
    print("Generating Figure S3: Het VAF violin...", file=sys.stderr)
    fig, ax = plt.subplots(1, 1, figsize=(5, 4.5))
    plot_het_vaf_violin([m["median_bias"] for m in markers], ax)
    fig.tight_layout()
    fig.savefig(FACTS_DIR / "figS3_het_vaf.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Main paper figure: Bias stability ---
    print("Generating bias stability figure...", file=sys.stderr)
    fig, ax = plt.subplots(1, 1, figsize=(5, 4.5))
    plot_bias_stability(
        [m["median_bias"] for m in markers],
        [m["sd_within"] for m in markers],
        ax,
    )
    fig.tight_layout()
    fig.savefig(FACTS_DIR / "fig_bias_stability.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Figure S4: Ablation study ---
    print("Generating Figure S4: Ablation study...", file=sys.stderr)
    with tempfile.TemporaryDirectory() as tmpdir:
        ablation_results = run_ablation(HOST_VCF, DONOR_VCF, Path(tmpdir))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    plot_ablation(ablation_results, ax1, ax2)
    fig.tight_layout()
    fig.savefig(FACTS_DIR / "figS4_ablation.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Figure S5: CI calibration (merge pre-computed batches) ---
    print("Generating Figure S5: CI calibration...", file=sys.stderr)
    cal_results = load_calibration_batches()
    if not cal_results:
        print(
            "WARNING: No calibration batches found in output/calibration/. "
            "Run calibration batches first.",
            file=sys.stderr,
        )
    else:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
        plot_calibration(cal_results, ax1, ax2)
        fig.tight_layout()
        fig.savefig(FACTS_DIR / "figS5_ci_calibration.png", dpi=150, bbox_inches="tight")
        plt.close(fig)

    # --- Figure S6: Residuals ---
    print("Generating Figure S6: Residuals...", file=sys.stderr)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
    with tempfile.TemporaryDirectory() as tmpdir:
        plot_residuals(HOST_VCF, DONOR_VCF, ax1, ax2, Path(tmpdir))
    fig.tight_layout()
    fig.savefig(FACTS_DIR / "figS6_residuals.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # --- Write facts CSV ---
    write_supp_facts(
        ablation_results,
        cal_results,
        empirical_biases,
        [m["sd_within"] for m in markers],
    )

    print("Done. Figures saved to output/facts/figS*.png", file=sys.stderr)


if __name__ == "__main__":
    main()
