#!/usr/bin/env python3
"""Plot the combined 2x2 LoD grid (Figure 1): test x relatedness.

Columns are the two readouts allomix runs on the same sample (MLE magnitude
estimate, host-presence detection test); rows are donor/host relatedness
(unrelated, full sibling). Reading left to right within a row shows the
quantify-vs-detect gap; reading top to bottom within a column shows the
identity-by-descent penalty (siblings share more genotype, so fewer markers are
informative). All four panels share log-log axes so every cell is directly
comparable.

Both readouts come from the same simulated sweep (matched depths, panel sizes,
donor/host pairs, and sequencing-error model), so the side-by-side comparison is
fair. The MLE panels also carry the limit-of-blank (LoB) as a faint dashed line;
the presence test has no blank (its null is the sequencing-error background).

Inputs:
  output/facts/lod_summary.csv                 # MLE LoD + LoB
  output/facts/presence_lod_curve_summary.csv  # presence-test LoD

Output:
  output/facts/fig_lod_curves.png
"""

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import paper_quick  # noqa: E402, F401  -- quick-build watermark (import for side effect)
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator  # noqa: E402

FACTS_DIR = Path("output/facts")

RELATEDNESS_ORDER = ["unrelated", "sibling"]
ROW_TITLES = {"unrelated": "Unrelated", "sibling": "Full sibling"}
THRESHOLDS = [0.5, 1.0]  # percent (reference action-zone lines)
# 25-marker panels are not realistic clinically; drop them from the plot (they
# stay in the summary CSVs so the headline facts are unaffected).
MIN_PLOT_MARKERS = 50
Y_LIMITS = (0.02, 5.0)
Y_TICKS = [0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 5]


def _format_pct(v: float, _pos: int) -> str:
    """Format a percent tick with the minimum decimals needed."""
    if v >= 1:
        return f"{v:g}%"
    return f"{v:g}%".rstrip("0").rstrip(".")


def _read_mle(path: Path) -> dict[str, dict[int, list[tuple]]]:
    """Read lod_summary.csv into by_rel[rel][depth] -> [(nm, lod, lo, hi, lob)]."""
    by_rel: dict[str, dict[int, list[tuple]]] = defaultdict(lambda: defaultdict(list))
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                depth = int(r["depth"])
                nm = int(r["n_markers"])

                def _f(key: str) -> float:
                    return float(r[key]) if r[key] not in ("", "nan") else float("nan")

                row = (nm, _f("lod_pct"), _f("lod_pct_ci_lo"), _f("lod_pct_ci_hi"), _f("lob_pct"))
            except (ValueError, KeyError):
                continue
            if nm < MIN_PLOT_MARKERS:
                continue
            by_rel[r["relatedness"]][depth].append(row)
    return by_rel


def _read_presence(path: Path) -> dict[str, dict[int, list[tuple]]]:
    """Read presence_lod_curve_summary.csv into by_rel[rel][depth] -> [(nm, lod, lo, hi)]."""
    by_rel: dict[str, dict[int, list[tuple]]] = defaultdict(lambda: defaultdict(list))
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                depth = int(r["depth"])
                nm = int(r["n_markers"])

                def _f(key: str) -> float:
                    return float(r[key]) if r[key] not in ("", "nan") else float("nan")

                row = (nm, _f("lod_pct"), _f("lod_pct_ci_lo"), _f("lod_pct_ci_hi"))
            except (ValueError, KeyError):
                continue
            if nm < MIN_PLOT_MARKERS:
                continue
            by_rel[r["relatedness"]][depth].append(row)
    return by_rel


def _finite_band(cell: list[tuple]) -> tuple[list[int], list[float], list[float], list[float]]:
    """Sort a [(nm, lod, lo, hi, ...)] cell by panel size, drop non-positive LoD."""
    out_x, out_y, out_lo, out_hi = [], [], [], []
    for row in sorted(cell, key=lambda t: t[0]):
        x, y, lo, hi = row[0], row[1], row[2], row[3]
        if y is None or not math.isfinite(y) or y <= 0:
            continue
        out_x.append(x)
        out_y.append(y)
        out_lo.append(lo if (lo is not None and math.isfinite(lo) and lo > 0) else y)
        out_hi.append(hi if (hi is not None and math.isfinite(hi) and hi > 0) else y)
    return out_x, out_y, out_lo, out_hi


def _depth_colors(depths: list[int]) -> dict[int, tuple]:
    """viridis_r cropped to [0.35, 1.0]: deepest depth -> dark purple, shallowest -> teal."""
    cmap = plt.get_cmap("viridis_r")
    if len(depths) > 1:
        return {d: cmap(0.35 + 0.65 * i / (len(depths) - 1)) for i, d in enumerate(depths)}
    return {depths[0]: cmap(1.0)}


def _style_axis(ax, nmarkers: list[int]) -> None:
    """Apply the shared log-log panel styling."""
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(nmarkers)
    ax.set_xticklabels([str(n) for n in nmarkers])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.grid(True, which="both", alpha=0.2)
    ax.yaxis.set_major_locator(FixedLocator(Y_TICKS))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_formatter(FuncFormatter(_format_pct))
    ax.set_ylim(*Y_LIMITS)
    for thr in THRESHOLDS:
        ax.axhline(thr, color="black", linestyle="--", linewidth=0.7, alpha=0.35)


def plot(mle_path: Path, presence_path: Path, out_path: Path) -> None:
    mle = _read_mle(mle_path)
    presence = _read_presence(presence_path)
    if not mle:
        raise SystemExit(f"No rows in {mle_path}")
    if not presence:
        raise SystemExit(f"No rows in {presence_path}")

    depths = sorted(
        {d for rel in mle.values() for d in rel} | {d for rel in presence.values() for d in rel}
    )
    nmarkers = sorted(
        {row[0] for rel in mle.values() for cell in rel.values() for row in cell}
        | {row[0] for rel in presence.values() for cell in rel.values() for row in cell}
    )
    colors = _depth_colors(depths)

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.4), sharex=True, sharey=True)

    for i, rel in enumerate(RELATEDNESS_ORDER):
        # MLE column (left): LoD solid + band, LoB faint dashed
        ax = axes[i][0]
        for depth in depths:
            cell = mle.get(rel, {}).get(depth, [])
            x_f, y_f, lo_f, hi_f = _finite_band(cell)
            if x_f:
                ax.plot(x_f, y_f, "o-", color=colors[depth], linewidth=1.8, markersize=6)
                ax.fill_between(x_f, lo_f, hi_f, color=colors[depth], alpha=0.15, linewidth=0)
            # LoB dashed (column 5 of each MLE row).
            lob = [(r[0], r[4]) for r in sorted(cell, key=lambda t: t[0])]
            lob = [(nm, v) for nm, v in lob if v is not None and math.isfinite(v) and v > 0]
            if lob:
                ax.plot(
                    [p[0] for p in lob],
                    [p[1] for p in lob],
                    "--",
                    color=colors[depth],
                    linewidth=1.0,
                    alpha=0.55,
                )
        _style_axis(ax, nmarkers)

        # Presence column (right): LoD solid square + band
        ax = axes[i][1]
        for depth in depths:
            cell = presence.get(rel, {}).get(depth, [])
            x_f, y_f, lo_f, hi_f = _finite_band(cell)
            if not x_f:
                continue
            ax.plot(
                x_f, y_f, "s-", color=colors[depth], linewidth=1.8, markersize=6, label=f"{depth}x"
            )
            ax.fill_between(x_f, lo_f, hi_f, color=colors[depth], alpha=0.15, linewidth=0)
        _style_axis(ax, nmarkers)

    axes[0][0].set_title("MLE (quantify)", fontsize=12, fontweight="bold")
    axes[0][1].set_title("Presence test (detect)", fontsize=12, fontweight="bold")
    for i, rel in enumerate(RELATEDNESS_ORDER):
        axes[i][0].set_ylabel("Limit of detection (%)", fontsize=11)
        axes[i][0].annotate(
            ROW_TITLES[rel],
            xy=(-0.27, 0.5),
            xycoords="axes fraction",
            rotation=90,
            ha="center",
            va="center",
            fontsize=13,
            fontweight="bold",
        )
    for j in range(2):
        axes[1][j].set_xlabel("Panel size (markers)", fontsize=11)

    axes[0][1].legend(title="Depth", fontsize=9, loc="upper right", framealpha=0.9)
    style_handles = [
        plt.Line2D([0], [0], color="grey", linewidth=1.8, marker="o", label="LoD (CI excludes 0)"),
        plt.Line2D(
            [0], [0], color="grey", linewidth=1.0, linestyle="--", alpha=0.7, label="LoB (dashed)"
        ),
    ]
    axes[0][0].legend(handles=style_handles, fontsize=8.5, loc="lower left", framealpha=0.9)

    fig.suptitle(
        "Limit of detection as a function of panel size and sequencing depth",
        fontsize=13,
        y=1.0,
    )
    fig.text(
        0.5,
        0.965,
        "simulated (in silico) data; MLE quantification (left) vs recipient-presence detection (right), "
        "matched sweep",
        ha="center",
        va="top",
        fontsize=9,
        color="0.35",
    )
    fig.tight_layout(rect=(0.03, 0, 1, 0.96))
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}", file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--mle-summary", default=str(FACTS_DIR / "lod_summary.csv"))
    parser.add_argument(
        "--presence-summary", default=str(FACTS_DIR / "presence_lod_curve_summary.csv")
    )
    parser.add_argument("--out", default=str(FACTS_DIR / "fig_lod_curves.png"))
    args = parser.parse_args(argv)
    plot(Path(args.mle_summary), Path(args.presence_summary), Path(args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
