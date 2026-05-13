#!/usr/bin/env python3
"""Plot LoD curves from output/facts/lod_summary.csv.

Two facets (unrelated, sibling). Each facet plots LoD (%) vs panel size on
log-log axes, with one coloured curve per sequencing depth and shaded
bootstrap CIs. A vertical reference line marks the 76-marker rhAmpSeq panel.

Output: output/facts/fig5_lod_curves.png
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

FACTS_DIR = Path("output/facts")

RELATEDNESS_ORDER = ["unrelated", "sibling"]
FACET_TITLES = {"unrelated": "Unrelated", "sibling": "Full sibling"}
REFERENCE_PANEL = 76  # our lab's IDT rhAmpSeq panel
THRESHOLDS = [0.1, 0.5, 1.0]  # percent


def _read_summary(path: Path) -> list[dict]:
    rows = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                r["depth"] = int(r["depth"])
                r["n_markers"] = int(r["n_markers"])
                r["lob_pct"] = float(r["lob_pct"]) if r["lob_pct"] else float("nan")
                r["lod_pct"] = float(r["lod_pct"]) if r["lod_pct"] else float("nan")
                r["lod_pct_ci_lo"] = (float(r["lod_pct_ci_lo"])
                                      if r["lod_pct_ci_lo"] else float("nan"))
                r["lod_pct_ci_hi"] = (float(r["lod_pct_ci_hi"])
                                      if r["lod_pct_ci_hi"] else float("nan"))
            except ValueError:
                continue
            rows.append(r)
    return rows


def _finite_xy(xs: list[int], ys: list[float],
               lo: list[float], hi: list[float]
               ) -> tuple[list[int], list[float], list[float], list[float]]:
    out_x, out_y, out_lo, out_hi = [], [], [], []
    for x, y, l, h in zip(xs, ys, lo, hi):
        if y is None or not math.isfinite(y) or y <= 0:
            continue
        out_x.append(x)
        out_y.append(y)
        out_lo.append(l if math.isfinite(l) and l > 0 else y)
        out_hi.append(h if math.isfinite(h) and h > 0 else y)
    return out_x, out_y, out_lo, out_hi


def plot(summary_path: Path, out_path: Path) -> None:
    rows = _read_summary(summary_path)
    if not rows:
        raise SystemExit(f"No rows in {summary_path}")

    # Group: by_rel[rel][depth] -> list of (n_markers, lod_pct, ci_lo, ci_hi) sorted
    by_rel: dict[str, dict[int, list[tuple]]] = defaultdict(lambda: defaultdict(list))
    depths_seen = set()
    nmarkers_seen = set()
    for r in rows:
        by_rel[r["relatedness"]][r["depth"]].append(
            (r["n_markers"], r["lod_pct"], r["lod_pct_ci_lo"], r["lod_pct_ci_hi"])
        )
        depths_seen.add(r["depth"])
        nmarkers_seen.add(r["n_markers"])

    depths = sorted(depths_seen)
    nmarkers = sorted(nmarkers_seen)
    cmap = plt.get_cmap("viridis_r")
    colors = {d: cmap(i / max(len(depths) - 1, 1)) for i, d in enumerate(depths)}

    fig, axes = plt.subplots(1, len(RELATEDNESS_ORDER), figsize=(12, 5.2),
                             sharex=True, sharey=True)
    if len(RELATEDNESS_ORDER) == 1:
        axes = [axes]

    for ax, rel in zip(axes, RELATEDNESS_ORDER):
        for depth in depths:
            cell = sorted(by_rel.get(rel, {}).get(depth, []), key=lambda t: t[0])
            xs = [c[0] for c in cell]
            ys = [c[1] for c in cell]
            lo = [c[2] for c in cell]
            hi = [c[3] for c in cell]
            x_f, y_f, lo_f, hi_f = _finite_xy(xs, ys, lo, hi)
            if not x_f:
                continue
            ax.plot(x_f, y_f, "o-", color=colors[depth], linewidth=1.8,
                    markersize=6, label=f"{depth}x")
            ax.fill_between(x_f, lo_f, hi_f, color=colors[depth], alpha=0.15,
                            linewidth=0)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks(nmarkers)
        ax.set_xticklabels([str(n) for n in nmarkers])
        ax.set_xlabel("Panel size (markers)", fontsize=11)
        ax.set_title(FACET_TITLES.get(rel, rel), fontsize=12, fontweight="bold")
        ax.grid(True, which="both", alpha=0.2)

        ax.axvline(REFERENCE_PANEL, color="grey", linestyle=":", linewidth=1.2,
                   alpha=0.8)
        ymin, ymax = ax.get_ylim()
        # Threshold lines on the right-most facet only get labels.
        for thr in THRESHOLDS:
            ax.axhline(thr, color="black", linestyle="--", linewidth=0.7,
                       alpha=0.35)

    axes[0].set_ylabel("Limit of detection (% donor)", fontsize=11)
    # Annotate reference panel on the left facet.
    ymin, ymax = axes[0].get_ylim()
    axes[0].text(REFERENCE_PANEL, ymax * 0.7, f" {REFERENCE_PANEL}-marker panel",
                 fontsize=8.5, color="grey", rotation=90,
                 verticalalignment="top")

    # Threshold labels on the right facet.
    for thr in THRESHOLDS:
        axes[-1].text(nmarkers[-1] * 1.05, thr, f" {thr}%", fontsize=8.5,
                      color="black", verticalalignment="center")

    axes[-1].legend(title="Depth", fontsize=9, loc="upper right",
                    framealpha=0.9)

    fig.suptitle(
        "Limit of detection as a function of panel size and sequencing depth",
        fontsize=13, y=1.02,
    )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}", file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--summary", default=str(FACTS_DIR / "lod_summary.csv"))
    parser.add_argument("--out", default=str(FACTS_DIR / "fig5_lod_curves.png"))
    args = parser.parse_args(argv)
    plot(Path(args.summary), Path(args.out))
    return 0


if __name__ == "__main__":
    sys.exit(main())
