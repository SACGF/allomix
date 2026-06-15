#!/usr/bin/env python3
"""Plot the simulated presence-test LoD curves (issue #25, companion to Figure 1).

Reads ``presence_lod_curve_summary.csv`` from
``run_presence_lod_validation.py`` and draws a two-facet figure (unrelated,
sibling) parallel to the MLE Figure 1 (``plot_lod_curves.py``): LoD (%) vs panel
size on log-log axes, one curve per sequencing depth, the 10th-90th percentile
band across donor/host pairs shaded, dashed reference lines at 0.5% and 1%.

Detection here is the host-presence LRT at p < 0.05 (the test's own null is the
sequencing-error background, so no blank/LoB is needed). Per-marker bias is off in
the underlying sweep to keep that null calibrated.

Output: output/facts/fig_presence_lod_curves.png
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
FACET_TITLES = {"unrelated": "Unrelated", "sibling": "Full sibling"}
THRESHOLDS = [0.5, 1.0]  # percent (reference action-zone lines)
# 25-marker panels are not realistic clinically; drop them from the plot (they
# stay in the summary CSV so the headline facts are unaffected).
MIN_PLOT_MARKERS = 50


def _format_pct(v: float, _pos: int) -> str:
    """Format a percent tick with the minimum decimals needed."""
    if v >= 1:
        return f"{v:g}%"
    return f"{v:g}%".rstrip("0").rstrip(".")


def _read_summary(path: Path) -> list[dict]:
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


def _finite_xy(
    cell: list[tuple],
) -> tuple[list[int], list[float], list[float], list[float]]:
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


def plot(summary_path: Path, out_path: Path) -> None:
    rows = _read_summary(summary_path)
    if not rows:
        raise SystemExit(f"No rows in {summary_path}")

    by_rel: dict[str, dict[int, list[tuple]]] = defaultdict(lambda: defaultdict(list))
    depths_seen, nmarkers_seen = set(), set()
    for r in rows:
        if r["n_markers"] < MIN_PLOT_MARKERS:
            continue
        by_rel[r["relatedness"]][r["depth"]].append(
            (r["n_markers"], r["lod_pct"], r["lod_pct_ci_lo"], r["lod_pct_ci_hi"])
        )
        depths_seen.add(r["depth"])
        nmarkers_seen.add(r["n_markers"])

    depths = sorted(depths_seen)
    nmarkers = sorted(nmarkers_seen)
    # viridis_r cropped to [0.35, 1.0]: deepest depth -> dark purple, shallowest
    # -> green-teal (matches plot_lod_curves.py).
    cmap = plt.get_cmap("viridis_r")
    if len(depths) > 1:
        colors = {d: cmap(0.35 + 0.65 * i / (len(depths) - 1)) for i, d in enumerate(depths)}
    else:
        colors = {depths[0]: cmap(1.0)}

    rels = [r for r in RELATEDNESS_ORDER if r in by_rel]
    fig, axes = plt.subplots(1, len(rels), figsize=(6 * len(rels), 5.2), sharex=True, sharey=True)
    if len(rels) == 1:
        axes = [axes]

    for ax, rel in zip(axes, rels):
        for depth in depths:
            cell = by_rel.get(rel, {}).get(depth, [])
            x_f, y_f, lo_f, hi_f = _finite_xy(cell)
            if not x_f:
                continue
            ax.plot(
                x_f, y_f, "s-", color=colors[depth], linewidth=1.8, markersize=6, label=f"{depth}x"
            )
            ax.fill_between(x_f, lo_f, hi_f, color=colors[depth], alpha=0.15, linewidth=0)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks(nmarkers)
        ax.set_xticklabels([str(n) for n in nmarkers])
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_xlabel("Panel size (markers)", fontsize=11)
        ax.set_title(FACET_TITLES.get(rel, rel), fontsize=12, fontweight="bold")
        ax.grid(True, which="both", alpha=0.2)
        ax.yaxis.set_major_locator(FixedLocator([0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 5]))
        ax.yaxis.set_minor_locator(NullLocator())
        ax.yaxis.set_major_formatter(FuncFormatter(_format_pct))
        ax.set_ylim(0.05, 5.0)
        for thr in THRESHOLDS:
            ax.axhline(thr, color="black", linestyle="--", linewidth=0.7, alpha=0.35)

    axes[0].set_ylabel("Limit of detection (% minor)", fontsize=11)
    axes[-1].legend(title="Depth", fontsize=9, loc="upper right", framealpha=0.9)

    fig.suptitle(
        "Presence-test limit of detection by panel size and sequencing depth",
        fontsize=13,
        y=1.02,
    )
    fig.text(
        0.5,
        0.985,
        "simulated; detection = host-presence LRT p < 0.05 (per-marker bias off "
        "to keep the null calibrated)",
        ha="center",
        va="top",
        fontsize=9,
        color="0.35",
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}", file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--summary", default=str(FACTS_DIR / "presence_lod_curve_summary.csv"))
    parser.add_argument("--out", default=str(FACTS_DIR / "fig_presence_lod_curves.png"))
    args = parser.parse_args(argv)
    plot(Path(args.summary), Path(args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
