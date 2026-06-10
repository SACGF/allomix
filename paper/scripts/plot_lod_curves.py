#!/usr/bin/env python3
"""Plot LoD curves from output/facts/lod_summary.csv.

Two facets (unrelated, sibling). Each facet plots LoD (%) vs panel size on
log-log axes, with one coloured curve per sequencing depth. The curve is the
median LoD across donor/host pairs and the shaded band is the 10th-90th
percentile across pairs (the IBD-driven spread, widest for siblings at small
panels).

Output: output/facts/fig5_lod_curves.png
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
from matplotlib.ticker import FuncFormatter, NullLocator  # noqa: E402


def _format_pct(v: float, _pos: int) -> str:
    """Format a percent value with the minimum number of decimals needed."""
    if v >= 1:
        return f"{v:g}%"
    return f"{v:g}%".rstrip("0").rstrip(".")

FACTS_DIR = Path("output/facts")

RELATEDNESS_ORDER = ["unrelated", "sibling"]
FACET_TITLES = {"unrelated": "Unrelated", "sibling": "Full sibling"}
THRESHOLDS = [0.5, 1.0]  # percent (reference action-zone lines)
# 25-marker panels aren't realistic in clinical practice; drop them from the
# plot to give the rest of the data more horizontal real estate. The data
# stays in lod_summary.csv (60 rows including 25-marker cells) so the
# headline facts are unaffected.
MIN_PLOT_MARKERS = 50


def _read_presence(path: Path, error_rate: float) -> dict[str, dict[int, list[tuple]]]:
    """Read presence_lod summary into by_rel[rel][depth] -> [(n_markers, lod_pct)].

    Filters to a single error rate (the presence sweep has an error-rate axis the
    chimerism sweep does not). ``presence_lod`` is a host fraction, converted to
    percent to share the chimerism LoD axis. Non-finite / unresolved cells are
    skipped.
    """
    by_rel: dict[str, dict[int, list[tuple]]] = defaultdict(lambda: defaultdict(list))
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                if abs(float(r["error_rate"]) - error_rate) > 1e-12:
                    continue
                nm = int(r["n_markers"])
                lod = float(r["presence_lod"])
            except (ValueError, KeyError):
                continue
            if not math.isfinite(lod) or lod <= 0:
                continue
            by_rel[r["relatedness"]][int(r["depth"])].append((nm, lod * 100.0))
    return by_rel


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


def plot(
    summary_path: Path,
    out_path: Path,
    presence_path: Path | None = None,
    presence_error_rate: float = 0.001,
) -> None:
    rows = _read_summary(summary_path)
    if not rows:
        raise SystemExit(f"No rows in {summary_path}")

    presence = (
        _read_presence(presence_path, presence_error_rate)
        if presence_path is not None and presence_path.exists()
        else None
    )

    # Group: by_rel[rel][depth] -> list of (n_markers, lod_pct, ci_lo, ci_hi, lob_pct) sorted
    by_rel: dict[str, dict[int, list[tuple]]] = defaultdict(lambda: defaultdict(list))
    depths_seen = set()
    nmarkers_seen = set()
    for r in rows:
        if r["n_markers"] < MIN_PLOT_MARKERS:
            continue
        by_rel[r["relatedness"]][r["depth"]].append(
            (r["n_markers"], r["lod_pct"], r["lod_pct_ci_lo"],
             r["lod_pct_ci_hi"], r["lob_pct"])
        )
        depths_seen.add(r["depth"])
        nmarkers_seen.add(r["n_markers"])

    depths = sorted(depths_seen)
    nmarkers = sorted(nmarkers_seen)
    # viridis_r cropped to [0.35, 1.0] avoids the bright-yellow end of the
    # default palette (low contrast on white). The deepest sequencing depth
    # (2000x, the most clinically interesting) maps to viridis_r(1.0) = dark
    # purple, the shallowest (100x) to viridis_r(0.35) = green-teal.
    cmap = plt.get_cmap("viridis_r")
    if len(depths) > 1:
        colors = {
            d: cmap(0.35 + 0.65 * i / (len(depths) - 1))
            for i, d in enumerate(depths)
        }
    else:
        colors = {depths[0]: cmap(1.0)}

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
            lob = [c[4] for c in cell]
            x_f, y_f, lo_f, hi_f = _finite_xy(xs, ys, lo, hi)
            if not x_f:
                continue
            ax.plot(x_f, y_f, "o-", color=colors[depth], linewidth=1.8,
                    markersize=6, label=f"{depth}x")
            ax.fill_between(x_f, lo_f, hi_f, color=colors[depth], alpha=0.15,
                            linewidth=0)
            # Plot LoB as a faint dashed line beneath LoD for the same depth.
            # The LoB curve (across-pair median of the per-pair blank 95th
            # percentile) shows the noise floor and is monotone in panel size
            # and depth; plotting it beside LoD makes the underlying
            # monotonicity of the estimator visible.
            lob_xy = [(x, y) for x, y in zip(xs, lob)
                      if y is not None and math.isfinite(y) and y > 0]
            if lob_xy:
                lob_x = [p[0] for p in lob_xy]
                lob_y = [p[1] for p in lob_xy]
                ax.plot(lob_x, lob_y, "--", color=colors[depth], linewidth=1.0,
                        alpha=0.55)

        # Overlay presence-test LoD (dotted, square markers) at matched depths,
        # same colour map. The presence sweep carries an error-rate axis; this
        # is the single rate chosen in main(). Markers shared with the chimerism
        # curve let you read the presence-vs-MLE LoD gap at each panel size.
        if presence is not None:
            for depth in depths:
                cell = sorted(presence.get(rel, {}).get(depth, []),
                              key=lambda t: t[0])
                pts = [(nm, lod) for nm, lod in cell
                       if nm >= MIN_PLOT_MARKERS and math.isfinite(lod) and lod > 0]
                if not pts:
                    continue
                ax.plot([p[0] for p in pts], [p[1] for p in pts], "s:",
                        color=colors.get(depth, "grey"), linewidth=1.6,
                        markersize=5, alpha=0.9)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xticks(nmarkers)
        ax.set_xticklabels([str(n) for n in nmarkers])
        # Suppress matplotlib's default log-axis minor ticks ("6×10¹", "3×10²")
        # that otherwise bleed through alongside our custom panel-size labels.
        ax.xaxis.set_minor_locator(NullLocator())
        ax.set_xlabel("Panel size (markers)", fontsize=11)
        ax.set_title(FACET_TITLES.get(rel, rel), fontsize=12, fontweight="bold")
        ax.grid(True, which="both", alpha=0.2)

        # Show y-axis tick labels as percentages instead of matplotlib's default
        # 10^-1 / 10^0 scientific format. Extra ticks at 0.3/0.4 because that
        # band is the clinically interesting action zone for deeper / larger
        # panels — the standard 1-2-5 convention would skip them.
        ax.yaxis.set_major_locator(
            plt.FixedLocator([0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 5])
        )
        ax.yaxis.set_minor_locator(NullLocator())
        ax.yaxis.set_major_formatter(FuncFormatter(_format_pct))
        # Trim y-range to focus on the action zone. The 25-marker cells (already
        # filtered above) and the worst-case sibling-100x-50-marker cell sit
        # above 2%. The lower bound reaches down to the deepest LoB band
        # (~0.023% at 2000x/400-marker) so the bottom LoD and LoB curves stay
        # inside the axes instead of running off the bottom.
        ax.set_ylim(0.02, 5.0)

        for thr in THRESHOLDS:
            ax.axhline(thr, color="black", linestyle="--", linewidth=0.7,
                       alpha=0.35)

    axes[0].set_ylabel("Limit of detection (% donor)", fontsize=11)

    axes[-1].legend(title="Depth", fontsize=9, loc="upper right",
                    framealpha=0.9)
    # Style key on the left facet: solid = LoD, dashed = LoB.
    style_handles = [
        plt.Line2D([0], [0], color="grey", linewidth=1.8, linestyle="-",
                   label="chimerism MLE LoD"),
        plt.Line2D([0], [0], color="grey", linewidth=1.0, linestyle="--",
                   alpha=0.7, label="LoB (dashed)"),
    ]
    if presence is not None:
        style_handles.append(
            plt.Line2D([0], [0], color="grey", linewidth=1.6, linestyle=":",
                       marker="s", markersize=5, label="presence LoD")
        )
    axes[0].legend(handles=style_handles, fontsize=8.5, loc="lower left",
                   framealpha=0.9)

    fig.suptitle(
        "Limit of detection as a function of panel size and sequencing depth",
        fontsize=13, y=1.02,
    )
    if presence is not None:
        # Make the error assumption explicit: the presence test is far more
        # error-sensitive than the MLE, so the two curves are only comparable
        # when the chimerism sweep was run at this same error rate.
        fig.text(
            0.5, 0.985,
            f"presence LoD overlaid at error rate {presence_error_rate:g}; "
            "valid only if the chimerism sweep used the same rate",
            ha="center", va="top", fontsize=9, color="0.35",
        )
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}", file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--summary", default=str(FACTS_DIR / "lod_summary.csv"))
    parser.add_argument("--out", default=str(FACTS_DIR / "fig5_lod_curves.png"))
    parser.add_argument(
        "--presence-summary", default=None,
        help="Optional presence_lod summary to overlay (e.g. "
             "output/facts/presence_lod_bypanel_summary.csv).",
    )
    parser.add_argument(
        "--presence-error-rate", type=float, default=0.001,
        help="Error-rate slice of the presence sweep to overlay (default 0.001).",
    )
    args = parser.parse_args(argv)
    plot(
        Path(args.summary),
        Path(args.out),
        Path(args.presence_summary) if args.presence_summary else None,
        args.presence_error_rate,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
