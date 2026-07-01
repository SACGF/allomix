"""Contamination-floor dose-response figure for SRP434573 (issue #19).

Reads the long-format per-site export written by probe_contam_median_srp434573.py
(output/facts/srp_contam_persite.csv: one row per pooled consensus-homozygous site
with n_carriers, n_alleles, minor_frac) and renders the contamination-floor boxplot:

  output/facts/fig_srp_contam.png   Supplementary Figure S13. x = n_carriers
                                       (0..5): each other co-pooled individual
                                       carrying the minor allele counts 1.

The boxplot shows per-site minor fraction per dose bin on a log y-axis (exact-zero
sites drawn at a small floor), the per-bin median overlaid as a connected line, and
n_sites annotated under each box. The monotonic rise from the no-carrier floor is
the evidence that the floor is real co-pooled material (index hopping), not flat
sequencing error.

An alternative dose-weighted x-axis (n_alleles, het=1/hom=2, 0..8+) is kept behind
the --allele-dose flag for review; it renders to output/srp_contam_alleles.png and
is not part of the paper build (the two x-axes tell the same story; the carrier
count was chosen for cleaner, denser bins).

House style follows plot_srp434573.py.
"""

import argparse
import csv
from pathlib import Path
from statistics import median

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import paper_quick  # noqa: E402, F401  -- quick-build watermark (import for side effect)
from matplotlib.ticker import FuncFormatter  # noqa: E402

OUT = Path("output")
FACTS_DIR = OUT / "facts"
PERSITE = FACTS_DIR / "srp_contam_persite.csv"
YFLOOR_PCT = 0.001  # percent; where exact-zero sites are drawn on the log axis
MEDIAN_COLOR = "#d62728"
BOX_COLOR = "#1f77b4"


def _fmt_pct(v: float, _pos: int) -> str:
    if v >= 1:
        return f"{v:g}%"
    return f"{v:g}%".rstrip("0").rstrip(".")


def _read() -> list[tuple[int, int, float]]:
    rows = []
    with open(PERSITE, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rows.append((int(r["n_carriers"]), int(r["n_alleles"]), float(r["minor_frac"])))
    return rows


def _floored_pct(frac: float) -> float:
    """Per-site minor fraction as a percent, with exact zeros at the log-axis floor."""
    return max(frac * 100.0, YFLOOR_PCT)


def plot_dose(
    rows: list[tuple[int, int, float]],
    col: int,
    bins: list[int],
    labels: list[str],
    merge_last: bool,
    xlabel: str,
    title: str,
    out_path: Path,
) -> None:
    """Boxplot of per-site minor fraction by co-pooled dose bin.

    Args:
        col: index into each row of the binning variable (0 carriers, 1 alleles).
        merge_last: if True, the final bin collects all values >= bins[-1] ("k+").
    """
    by_bin: dict[int, list[float]] = {b: [] for b in bins}
    for row in rows:
        d = row[col]
        if merge_last and d >= bins[-1]:
            by_bin[bins[-1]].append(_floored_pct(row[2]))
        elif d in by_bin:
            by_bin[d].append(_floored_pct(row[2]))

    positions = list(range(len(bins)))
    data = [by_bin[b] for b in bins]
    counts = [len(d) for d in data]
    medians = [median(d) if d else YFLOOR_PCT for d in data]

    fig, ax = plt.subplots(figsize=(8.0, 5.6))
    bp = ax.boxplot(
        data, positions=positions, widths=0.6, showfliers=True,
        patch_artist=True, medianprops=dict(color="white", linewidth=1.2),
        flierprops=dict(marker="o", markersize=2.5, markerfacecolor="0.5",
                        markeredgecolor="none", alpha=0.4),
        whiskerprops=dict(color="0.4"), capprops=dict(color="0.4"),
    )
    for box in bp["boxes"]:
        box.set(facecolor=BOX_COLOR, alpha=0.55, edgecolor="0.3")

    ax.plot(positions, medians, marker="D", markersize=7, color=MEDIAN_COLOR,
            linewidth=1.8, markeredgecolor="white", markeredgewidth=0.7, zorder=5,
            label="per-bin median")

    ax.set_yscale("log")
    ax.axhline(YFLOOR_PCT, color="0.6", linewidth=0.8, linestyle=":", zorder=0)
    ax.set_ylim(YFLOOR_PCT * 0.7, max(max(d) for d in data if d) * 1.4)
    ax.yaxis.set_major_formatter(FuncFormatter(_fmt_pct))
    ax.grid(True, axis="y", which="both", alpha=0.2)

    ax.set_xticks(positions)
    ax.set_xticklabels(labels)
    ax.set_xlim(-0.7, len(bins) - 0.3)
    ax.set_xlabel(xlabel, fontsize=12, labelpad=26)
    ax.set_ylabel("per-site minor fraction (log scale)", fontsize=12)
    ax.set_title(title, fontsize=12.5, fontweight="bold")

    for pos, n in zip(positions, counts):
        ax.annotate(f"n={n}", xy=(pos, 0.0), xytext=(0, -20), textcoords="offset points",
                    xycoords=("data", "axes fraction"), ha="center", va="top",
                    fontsize=8.5, color="0.3")

    ax.legend(fontsize=9.5, loc="upper left", framealpha=0.92)
    fig.tight_layout()
    fig.subplots_adjust(bottom=0.20)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--allele-dose", action="store_true",
        help="also render the dose-weighted (n_alleles) x-axis view to "
             "output/srp_contam_alleles.png (review-only, not in the paper build)",
    )
    args = parser.parse_args()

    FACTS_DIR.mkdir(parents=True, exist_ok=True)
    rows = _read()

    # Supplementary Figure S13: carrier-count x-axis (the paper figure).
    plot_dose(
        rows, col=0,
        bins=[0, 1, 2, 3, 4, 5],
        labels=["0", "1", "2", "3", "4", "5"],
        merge_last=False,
        xlabel="co-pooled individuals carrying the minor allele (het or hom = 1)",
        title="SRP434573 contamination floor by co-pooled carrier count",
        out_path=FACTS_DIR / "fig_srp_contam.png",
    )

    if args.allele_dose:
        plot_dose(
            rows, col=1,
            bins=[0, 1, 2, 3, 4, 5, 6, 7, 8],
            labels=["0", "1", "2", "3", "4", "5", "6", "7", "8+"],
            merge_last=True,
            xlabel="co-pooled minor alleles (het = 1, hom = 2)",
            title="SRP434573 contamination floor by co-pooled allele dose",
            out_path=OUT / "srp_contam_alleles.png",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
