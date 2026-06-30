#!/usr/bin/env python3
"""Plot how recipient CNV/LoH affects the limit of detection, both directions (issue #13).

Reads output/facts/cnv_loh_summary.csv and draws a 2-row grid (one row per
detection direction) x 3 columns (aberration kind). The y axis is the LoD of the
minor component as a %, log-scaled with percent ticks, the same style as the
depth x markers LoD curves (plot_lod_grid.py). Standard (solid) vs robust
refit (dashed), one colour per relatedness, baseline (no aberration) at burden 0.

The two rows are the two low-fraction detection tasks:
  - "host": detect a recipient relapse (the minor component, carrying the
    aberration) against a clean donor background. Early-warning use; little
    affected by recipient CN-LoH.
  - "donor": detect the donor (minor) against a recipient CN-LoH background
    (mixed-chimerism / substantial recipient). Donor LoD badly inflated; robust
    refit recovers it.

A LoD above the 20% probed ceiling (component undetectable) is drawn above the
dotted line.
"""

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import paper_quick  # noqa: E402, F401  -- quick-build watermark (import for side effect)
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator  # noqa: E402

FACTS_DIR = Path("output/facts")
KINDS = ["cnloh", "deletion", "gain"]
KIND_LABELS = {"cnloh": "CN-LoH (copy-neutral)", "deletion": "Deletion (CN1)", "gain": "Gain (CN3)"}
# Relapse detection first (primary use), then donor detection.
MODES = ["host", "donor"]
MODE_YLABEL = {
    "host": "Relapse (recipient) LoD (%)",
    "donor": "Donor LoD (%)",
}
COLORS = {"unrelated": "#1f77b4", "sibling": "#d62728"}
MAX_PROBED_PCT = 20.0
CEILING_PCT = 32.0
Y_TICKS = [0.3, 0.5, 1, 2, 5, 10, 20]


def _format_pct(v: float, _pos: int) -> str:
    if v >= 1:
        return f"{v:g}%"
    return f"{v:g}%".rstrip("0").rstrip(".")


def _pct(x: str) -> float:
    """Summary LoD cell (a fraction, or 'inf'/'') -> % for plotting."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return float("nan")
    if v != v:
        return float("nan")
    if v == float("inf"):
        return CEILING_PCT
    return v * 100.0


def load_summary(path: Path) -> list[dict]:
    rows = []
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            rows.append(
                {
                    "mode": r["mode"],
                    "relatedness": r["relatedness"],
                    "kind": r["kind"],
                    "burden": float(r["burden"]),
                    # depth / n_markers present once the depth x markers grid is swept.
                    "depth": int(r["depth"]) if r.get("depth") not in (None, "") else None,
                    "n_markers": int(r["n_markers"]) if r.get("n_markers") not in (None, "") else None,
                    "lod_std": _pct(r["lod_std"]),
                    "lod_robust": _pct(r["lod_robust"]),
                }
            )
    return rows


def series(rows: list[dict], mode: str, kind: str, field: str) -> dict[str, dict[float, float]]:
    """Per-relatedness {burden: LoD %} for one mode/kind. Burden 0 = baseline."""
    acc: dict[str, dict[float, list[float]]] = defaultdict(lambda: defaultdict(list))
    for r in rows:
        if r["mode"] != mode:
            continue
        v = r[field]
        if v != v:  # NaN: fit failed
            continue
        if r["kind"] == "baseline":
            acc[r["relatedness"]][0.0].append(v)
        elif r["kind"] == kind:
            acc[r["relatedness"]][r["burden"]].append(v)
    return {rel: {b: statistics.fmean(v) for b, v in bur.items()} for rel, bur in acc.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=FACTS_DIR / "cnv_loh_summary.csv")
    parser.add_argument("--out", type=Path, default=FACTS_DIR / "fig_cnv_loh.png")
    parser.add_argument("--depth", type=int, default=None,
                        help="Depth slice for the burden figure (default: largest present)")
    parser.add_argument("--markers", type=int, default=None,
                        help="Panel-size slice for the burden figure (default: largest present)")
    args = parser.parse_args()

    rows = load_summary(args.summary)
    # The burden figure is a slice at one operating point; pick the largest swept
    # depth / panel unless overridden (single-point sweeps have just one).
    depths = [r["depth"] for r in rows if r["depth"] is not None]
    markers = [r["n_markers"] for r in rows if r["n_markers"] is not None]
    sel_depth = args.depth if args.depth is not None else (max(depths) if depths else None)
    sel_markers = args.markers if args.markers is not None else (max(markers) if markers else None)
    if sel_depth is not None:
        rows = [r for r in rows if r["depth"] in (None, sel_depth)
                and r["n_markers"] in (None, sel_markers)]
    modes = [m for m in MODES if any(r["mode"] == m for r in rows)]

    fig, axes = plt.subplots(
        len(modes), len(KINDS), figsize=(4.6 * len(KINDS), 4.4 * len(modes)),
        sharey=True, squeeze=False,
    )
    for i, mode in enumerate(modes):
        for j, kind in enumerate(KINDS):
            std = series(rows, mode, kind, "lod_std")
            rob = series(rows, mode, kind, "lod_robust")
            ax = axes[i][j]
            for rel in sorted(std):
                c = COLORS.get(rel)
                bs = sorted(std[rel])
                ax.plot(bs, [std[rel][b] for b in bs], marker="o", color=c, label=f"{rel} (std)")
                if rel in rob:
                    br = sorted(rob[rel])
                    ax.plot(br, [rob[rel][b] for b in br], marker="s", ls="--", color=c,
                            label=f"{rel} (robust)")
            ax.set_yscale("log")
            ax.set_ylim(0.2, CEILING_PCT * 1.25)
            ax.yaxis.set_major_locator(FixedLocator(Y_TICKS))
            ax.yaxis.set_minor_locator(NullLocator())
            ax.yaxis.set_major_formatter(FuncFormatter(_format_pct))
            ax.axhline(MAX_PROBED_PCT, color="0.7", ls=":", lw=1)
            ax.grid(True, which="both", alpha=0.2)
            if i == 0:
                ax.set_title(KIND_LABELS[kind])
            if i == len(modes) - 1:
                ax.set_xlabel("CN-aberration burden (fraction of eligible markers)")
            if j == 0:
                ax.set_ylabel(MODE_YLABEL.get(mode, mode))
                ax.legend(fontsize=8, loc="lower right")

    fig.suptitle(
        "Recipient CNV/LoH and the detection limit: relapse detection (top) is "
        "robust; donor detection in mixed chimerism (bottom) is degraded but "
        "recovered by robust refit"
    )
    panel_desc = (
        f"{sel_markers or '?'} markers, {sel_depth or '?'}x depth, pure clone (clonal fraction 1.0)"
    )
    fig.text(
        0.5, 0.005,
        f"{panel_desc}. Depth/markers set the baseline (burden 0) LoD per the depth x markers "
        "curves; the CN-LoH/deletion LoD floor is a systematic-bias limit (largely "
        "depth-independent). More markers also improve the robust recovery.",
        ha="center", va="bottom", fontsize=7, color="0.35", wrap=True,
    )
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
