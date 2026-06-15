#!/usr/bin/env python3
"""Plot the real-data subsample LoD curves from subsample_lod_summary.csv.

Two standalone figures (real SRP434573 data only, no simulated overlay), each
styled like ``plot_lod_curves.py`` so they read as Figure-1-style panels: LoD (%)
vs panel size on log-log axes, one curve per mean depth, the 10th-90th percentile
band across mixtures shaded, dashed reference lines at 0.5% and 1%.

  output/facts/fig_subsample_lod_mle.png       magnitude (MLE) LoD, with the
      analytical lod_fraction overlaid as a faint dashed line (both MLE-side, so
      they share one figure).
  output/facts/fig_subsample_lod_presence.png  presence-test LoD.

The curves are pseudo-replicates from one real read draw per mixture (sub-sampled
reads, not independent low-depth libraries); see
``claude/public_data_subsample_plan.md`` section 7 for the framing the caption
must carry. The low-fraction plateau is the co-pooled contamination floor of this
dataset, not allomix's intrinsic limit.
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

THRESHOLDS = [0.5, 1.0]  # percent (reference action-zone lines)
# 25-marker panels are not realistic clinically; drop them from the plot (they
# stay in the summary CSV so the headline facts are unaffected).
MIN_PLOT_MARKERS = 50


def _read_summary(path: Path) -> list[dict]:
    """Read subsample_lod_summary.csv, coercing numeric columns."""
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                r["depth"] = int(r["depth"])
                r["n_markers"] = int(r["n_markers"])
                for k in ("lod_pct", "lod_pct_ci_lo", "lod_pct_ci_hi"):
                    r[k] = float(r[k]) if r[k] not in ("", "nan") else float("nan")
            except (ValueError, KeyError):
                continue
            rows.append(r)
    return rows


def _by_depth(rows: list[dict], test: str) -> dict[int, list[tuple]]:
    """Group one test's rows into depth -> [(n_markers, lod, lo, hi), ...]."""
    out: dict[int, list[tuple]] = defaultdict(list)
    for r in rows:
        if r["test"] != test or r["n_markers"] < MIN_PLOT_MARKERS:
            continue
        out[r["depth"]].append(
            (r["n_markers"], r["lod_pct"], r["lod_pct_ci_lo"], r["lod_pct_ci_hi"])
        )
    return out


def _finite_xy(cell: list[tuple]) -> tuple[list[int], list[float], list[float], list[float]]:
    """Sort by panel size and drop non-finite / non-positive LoD points."""
    out_x, out_y, out_lo, out_hi = [], [], [], []
    for x, y, lo, hi in sorted(cell, key=lambda t: t[0]):
        if y is None or not math.isfinite(y) or y <= 0:
            continue
        out_x.append(x)
        out_y.append(y)
        out_lo.append(lo if (lo is not None and math.isfinite(lo) and lo > 0) else y)
        out_hi.append(hi if (hi is not None and math.isfinite(hi) and hi > 0) else y)
    return out_x, out_y, out_lo, out_hi


def _depth_colors(depths: list[int]) -> dict[int, tuple]:
    """viridis_r palette: deepest depth -> dark purple, shallowest -> teal."""
    cmap = plt.get_cmap("viridis_r")
    if len(depths) > 1:
        return {d: cmap(0.35 + 0.65 * i / (len(depths) - 1)) for i, d in enumerate(depths)}
    return {depths[0]: cmap(1.0)}


def _style_axis(ax, nmarkers: list[int]) -> None:
    """Apply the shared log-log panel styling used by plot_lod_curves.py."""
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(nmarkers)
    ax.set_xticklabels([str(n) for n in nmarkers])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xlabel("Panel size (markers)", fontsize=11)
    ax.grid(True, which="both", alpha=0.2)
    ax.yaxis.set_major_locator(FixedLocator([0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 5, 10]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_formatter(FuncFormatter(_format_pct))
    ax.set_ylim(0.05, 10.0)
    for thr in THRESHOLDS:
        ax.axhline(thr, color="black", linestyle="--", linewidth=0.7, alpha=0.35)
    ax.set_ylabel("Limit of detection (% minor)", fontsize=11)


def _format_pct(v: float, _pos: int) -> str:
    """Format a percent tick with the minimum decimals needed."""
    if v >= 1:
        return f"{v:g}%"
    return f"{v:g}%".rstrip("0").rstrip(".")


def plot_mle(rows: list[dict], out_path: Path) -> None:
    """MLE magnitude LoD vs panel size, with the analytical LoD overlaid."""
    empirical = _by_depth(rows, "mle")
    analytical = _by_depth(rows, "mle_analytical")
    if not empirical:
        raise SystemExit("No mle rows to plot")

    depths = sorted(empirical)
    nmarkers = sorted({nm for cell in empirical.values() for nm, *_ in cell})
    colors = _depth_colors(depths)

    fig, ax = plt.subplots(figsize=(7, 5.4))
    for depth in depths:
        x, y, lo, hi = _finite_xy(empirical[depth])
        if not x:
            continue
        ax.plot(x, y, "o-", color=colors[depth], linewidth=1.8, markersize=6, label=f"{depth}x")
        ax.fill_between(x, lo, hi, color=colors[depth], alpha=0.15, linewidth=0)
        # Faint dashed analytical (Fisher-info) LoD on the same axes.
        ax_x, ax_y, _, _ = _finite_xy(analytical.get(depth, []))
        if ax_x:
            ax.plot(ax_x, ax_y, "--", color=colors[depth], linewidth=1.0, alpha=0.55)

    _style_axis(ax, nmarkers)
    ax.set_title("Real-data LoD: magnitude (MLE) estimate", fontsize=12, fontweight="bold")

    depth_legend = ax.legend(title="Mean depth", fontsize=9, loc="upper right", framealpha=0.9)
    ax.add_artist(depth_legend)
    style_handles = [
        plt.Line2D(
            [0],
            [0],
            color="grey",
            linewidth=1.8,
            linestyle="-",
            label="empirical MLE LoD (CI excludes 0)",
        ),
        plt.Line2D(
            [0],
            [0],
            color="grey",
            linewidth=1.0,
            linestyle="--",
            alpha=0.7,
            label="analytical LoD (Fisher info)",
        ),
    ]
    ax.legend(handles=style_handles, fontsize=8.5, loc="lower left", framealpha=0.9)

    fig.text(
        0.5,
        0.005,
        "SRP434573 sub-sampled reads (pseudo-replicates); low-fraction plateau is "
        "this dataset's contamination floor",
        ha="center",
        va="bottom",
        fontsize=8,
        color="0.35",
    )
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}", file=sys.stderr)


def plot_presence(rows: list[dict], out_path: Path) -> None:
    """Presence-test LoD vs panel size."""
    presence = _by_depth(rows, "presence")
    if not presence:
        raise SystemExit("No presence rows to plot")

    depths = sorted(presence)
    nmarkers = sorted({nm for cell in presence.values() for nm, *_ in cell})
    colors = _depth_colors(depths)

    fig, ax = plt.subplots(figsize=(7, 5.4))
    for depth in depths:
        x, y, lo, hi = _finite_xy(presence[depth])
        if not x:
            continue
        ax.plot(x, y, "s-", color=colors[depth], linewidth=1.8, markersize=6, label=f"{depth}x")
        ax.fill_between(x, lo, hi, color=colors[depth], alpha=0.15, linewidth=0)

    _style_axis(ax, nmarkers)
    ax.set_title("Real-data LoD: presence test (LRT)", fontsize=12, fontweight="bold")
    ax.legend(title="Mean depth", fontsize=9, loc="upper right", framealpha=0.9)

    fig.text(
        0.5,
        0.005,
        "SRP434573 sub-sampled reads (pseudo-replicates); detection = host-presence LRT p < 0.05",
        ha="center",
        va="bottom",
        fontsize=8,
        color="0.35",
    )
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}", file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--summary", default=str(FACTS_DIR / "subsample_lod_summary.csv"))
    parser.add_argument("--out-mle", default=str(FACTS_DIR / "fig_subsample_lod_mle.png"))
    parser.add_argument("--out-presence", default=str(FACTS_DIR / "fig_subsample_lod_presence.png"))
    args = parser.parse_args(argv)

    rows = _read_summary(Path(args.summary))
    if not rows:
        raise SystemExit(f"No rows in {args.summary}")
    plot_mle(rows, Path(args.out_mle))
    plot_presence(rows, Path(args.out_presence))
    return 0


if __name__ == "__main__":
    sys.exit(main())
