#!/usr/bin/env python3
"""Generate multi-donor validation figure for the allomix paper.

Panel A: Per-donor accuracy scatter (true vs estimated) with identity line and CIs.
Panel B: 2D log-likelihood contour for one representative mixture.

Usage:
    python paper/scripts/generate_multidonor_figure.py
"""

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from scipy.stats import chi2  # noqa: E402

from allomix.chimerism import estimate_multi_donor, total_log_likelihood_multi_bb  # noqa: E402
from allomix.genotype import classify_markers, parse_vcf  # noqa: E402

DATA_DIR = Path("tests/test_data/multidonor")
FACTS_DIR = Path("output/facts")

# Representative interior mixture for contour plot
CONTOUR_SAMPLE = "host_60_d1_30_d2_10"  # 60% host, 30% d1, 10% d2
CONTOUR_TRUE_F1 = 0.30
CONTOUR_TRUE_F2 = 0.10


def run_all_samples():
    """Run multi-donor estimation on all chimeric VCFs."""
    host_vcf = str(DATA_DIR / "host.vcf")
    donor1_vcf = str(DATA_DIR / "donor1.vcf")
    donor2_vcf = str(DATA_DIR / "donor2.vcf")

    host = parse_vcf(host_vcf, min_dp=0, min_gq=0)
    donor1 = parse_vcf(donor1_vcf, min_dp=0, min_gq=0)
    donor2 = parse_vcf(donor2_vcf, min_dp=0, min_gq=0)

    truth_rows = []
    with open(DATA_DIR / "truth_table.tsv", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            truth_rows.append(row)

    results = []
    for t in truth_rows:
        sample_path = str(DATA_DIR / f"{t['sample_name']}.vcf")
        admix = parse_vcf(sample_path, min_dp=0, min_gq=0)
        genotypes = classify_markers(
            host, [donor1, donor2], admix, min_dp=0, min_gq=0, pass_only=False
        )
        result = estimate_multi_donor(genotypes.informative, n_donors=2, error_rate=0.01)

        results.append(
            {
                "sample_name": t["sample_name"],
                "true_f1": float(t["true_donor1_fraction"]),
                "true_f2": float(t["true_donor2_fraction"]),
                "est_f1": result.donor_fractions[0],
                "est_f2": result.donor_fractions[1],
                "ci1_lo": result.donor_fraction_cis[0][0],
                "ci1_hi": result.donor_fraction_cis[0][1],
                "ci2_lo": result.donor_fraction_cis[1][0],
                "ci2_hi": result.donor_fraction_cis[1][1],
                "markers": genotypes.informative,
            }
        )

    return results


def plot_figure(results: list[dict]) -> None:
    """Generate the 2-panel figure."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))

    # --- Panel A: Per-donor accuracy scatter ---
    ax = axes[0]

    for r in results:
        # Donor 1 (blue circles)
        ax.errorbar(
            r["true_f1"] * 100,
            r["est_f1"] * 100,
            yerr=[
                [max(0, (r["est_f1"] - r["ci1_lo"]) * 100)],
                [max(0, (r["ci1_hi"] - r["est_f1"]) * 100)],
            ],
            fmt="o",
            color="steelblue",
            markersize=5,
            capsize=2,
            elinewidth=0.8,
            alpha=0.8,
            zorder=3,
        )
        # Donor 2 (orange triangles)
        ax.errorbar(
            r["true_f2"] * 100,
            r["est_f2"] * 100,
            yerr=[
                [max(0, (r["est_f2"] - r["ci2_lo"]) * 100)],
                [max(0, (r["ci2_hi"] - r["est_f2"]) * 100)],
            ],
            fmt="^",
            color="darkorange",
            markersize=5,
            capsize=2,
            elinewidth=0.8,
            alpha=0.8,
            zorder=3,
        )

    ax.plot([0, 100], [0, 100], "k--", alpha=0.4, linewidth=1, label="Identity")
    ax.set_xlabel("True donor fraction (%)", fontsize=11)
    ax.set_ylabel("Estimated donor fraction (%)", fontsize=11)
    ax.set_title("A. Per-donor estimation accuracy", fontsize=12, fontweight="bold", loc="left")
    ax.set_xlim(-3, 103)
    ax.set_ylim(-3, 103)
    ax.set_aspect("equal")
    ax.grid(True, alpha=0.2)

    # Custom legend
    legend_elements = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="steelblue",
            markersize=8,
            label="Donor 1",
        ),
        Line2D(
            [0],
            [0],
            marker="^",
            color="w",
            markerfacecolor="darkorange",
            markersize=8,
            label="Donor 2",
        ),
        Line2D([0], [0], color="k", linestyle="--", alpha=0.4, label="Identity"),
    ]
    ax.legend(handles=legend_elements, fontsize=9, loc="upper left")

    # --- Panel B: 2D log-likelihood contour ---
    ax = axes[1]

    # Find the representative sample's markers
    contour_result = None
    for r in results:
        if r["sample_name"] == CONTOUR_SAMPLE:
            contour_result = r
            break

    if contour_result is None:
        print(f"Warning: {CONTOUR_SAMPLE} not found, skipping panel B", file=sys.stderr)
    else:
        markers = contour_result["markers"]

        # Evaluate LL on a grid
        grid_n = 200
        f1_range = np.linspace(0, 1, grid_n)
        f2_range = np.linspace(0, 1, grid_n)
        ll_grid = np.full((grid_n, grid_n), np.nan)

        for i, f1 in enumerate(f1_range):
            for j, f2 in enumerate(f2_range):
                if f1 + f2 <= 1.0:
                    ll_grid[j, i] = total_log_likelihood_multi_bb(markers, [f1, f2], error_rate=0.01)

        # Find maximum LL
        ll_max = np.nanmax(ll_grid)
        delta_ll = ll_max - ll_grid

        # 95% joint CI contour: chi2(df=2) / 2
        threshold_95 = chi2.ppf(0.95, df=2) / 2  # ~2.996

        # Grey out infeasible region
        infeasible = np.zeros_like(ll_grid, dtype=bool)
        for i in range(grid_n):
            for j in range(grid_n):
                if f1_range[i] + f2_range[j] > 1.0:
                    infeasible[j, i] = True

        # Mask infeasible region
        delta_ll_masked = np.where(infeasible, np.nan, delta_ll)

        # Filled contour of delta-LL
        levels = [0, 0.5, 1, threshold_95, 5, 10, 20, 50]
        cf = ax.contourf(
            f1_range * 100,
            f2_range * 100,
            delta_ll_masked,
            levels=levels,
            cmap="Blues_r",
            alpha=0.8,
            extend="max",
        )

        # 95% CI contour line
        ax.contour(
            f1_range * 100,
            f2_range * 100,
            delta_ll_masked,
            levels=[threshold_95],
            colors="firebrick",
            linewidths=1.5,
            linestyles="--",
        )

        # Grey out infeasible region
        f1_fill = np.array([0, 100, 100])
        f2_fill = np.array([100, 100, 0])
        ax.fill(f1_fill, f2_fill, color="0.85", alpha=0.6, zorder=2)
        ax.plot([0, 100], [100, 0], color="0.5", linewidth=1, linestyle="-", zorder=2)

        # Mark true point and MLE
        ax.plot(
            CONTOUR_TRUE_F1 * 100,
            CONTOUR_TRUE_F2 * 100,
            "*",
            color="gold",
            markersize=14,
            markeredgecolor="k",
            markeredgewidth=0.8,
            zorder=5,
            label="True",
        )
        ax.plot(
            contour_result["est_f1"] * 100,
            contour_result["est_f2"] * 100,
            "o",
            color="firebrick",
            markersize=8,
            markeredgecolor="k",
            markeredgewidth=0.8,
            zorder=5,
            label="MLE",
        )

        cbar = fig.colorbar(cf, ax=ax, shrink=0.8)
        cbar.set_label("Δ log-likelihood from maximum", fontsize=9)

        ax.set_xlabel("Donor 1 fraction (%)", fontsize=11)
        ax.set_ylabel("Donor 2 fraction (%)", fontsize=11)
        ax.set_title("B. Log-likelihood surface", fontsize=12, fontweight="bold", loc="left")
        ax.set_xlim(0, 80)
        ax.set_ylim(0, 80)
        ax.set_aspect("equal")
        ax.grid(True, alpha=0.2)

        # Legend for panel B
        legend_b = [
            Line2D(
                [0],
                [0],
                marker="*",
                color="w",
                markerfacecolor="gold",
                markeredgecolor="k",
                markersize=12,
                label="True",
            ),
            Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                markerfacecolor="firebrick",
                markeredgecolor="k",
                markersize=8,
                label="MLE",
            ),
            Line2D([0], [0], color="firebrick", linestyle="--", linewidth=1.5, label="95% CI"),
        ]
        ax.legend(handles=legend_b, fontsize=9, loc="upper right")

    fig.tight_layout()

    # Save
    FACTS_DIR.mkdir(parents=True, exist_ok=True)
    fig_path = FACTS_DIR / "fig_multidonor.png"
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {fig_path}", file=sys.stderr)


def main() -> int:
    print("Generating multi-donor figure...", file=sys.stderr)
    results = run_all_samples()
    plot_figure(results)
    return 0


if __name__ == "__main__":
    sys.exit(main())
