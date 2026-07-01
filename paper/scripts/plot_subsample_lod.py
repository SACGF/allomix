#!/usr/bin/env python3
"""Plot the real-data subsample LoD grid from subsample_lod_summary.csv.

One 2x2 figure (real SRP434573 data only), laid out like the simulated Figure 1
(``plot_lod_grid.py``):

  columns: MLE magnitude estimate (left) vs host-presence detection test (right)
  rows:    the mixtures that stop at 1% (top) vs the mixtures titrated to 0.5%
           (bottom), two disjoint sets

Each panel plots LoD (%) vs panel size on log-log axes, one curve per mean depth,
the 10th-90th percentile band across mixtures shaded, reference lines at 0.5% and
1%. The two rows are disjoint: only 3 of the 10 mixtures were titrated below 1%,
so the top row holds the 7 that stop at 1% and the bottom row the 3 that resolve
down to 0.5%. Keeping them apart stops the 0.5%-reaching mixtures from being
buried in a top-row median that the 1%-floored mixtures would otherwise pin at 1%.

A cell is censored when its median mixture detected every titration point it was
given: the true LoD is then at or below that mixture set's lowest dilution (1%
top row, 0.5% bottom row), drawn as an X (LoD <= marker), not a resolved value. The
flat regions are the bottom of the dilution grid, not allomix's intrinsic limit.

The per-mixture LoD is monotone-constrained across panel size upstream (the
panels are nested, so more markers cannot raise the LoD), which removes the
threshold-read-off noise that otherwise wobbles the median; see
run_subsample_lod.py. Points are jittered horizontally per depth so curves that
collapse onto the same value stay readable.

  output/facts/fig_subsample_lod_grid.png

The analytical (Fisher-info) LoD overlay is off by default (it crowds the panels);
pass --show-analytical to bring it back on the MLE column.
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

COLUMNS = [("mle", "MLE (quantify)", "o"), ("presence", "Presence test (detect)", "s")]
# Max fractional horizontal offset applied per depth so curves that collapse onto
# the same censored value (e.g. all depths at 1% top row, 0.5% bottom row) fan out
# and stay readable. Multiplicative because the panel-size axis is log-scaled.
JITTER = 0.05


def _format_pct(v: float, _pos: int) -> str:
    """Format a percent tick with the minimum decimals needed."""
    if v >= 1:
        return f"{v:g}%"
    return f"{v:g}%".rstrip("0").rstrip(".")


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
                r["censored"] = str(r.get("censored", "")).strip().lower() in ("true", "1")
                r["mixture_set"] = r.get("mixture_set") or "all"
            except (ValueError, KeyError):
                continue
            rows.append(r)
    return rows


def _by_depth(rows: list[dict], test: str, mixture_set: str = "all") -> dict[int, list[tuple]]:
    """Group one test/mixture_set's rows into depth -> [(n_markers, lod, lo, hi, censored), ...]."""
    out: dict[int, list[tuple]] = defaultdict(list)
    for r in rows:
        if (
            r["test"] != test
            or r.get("mixture_set", "all") != mixture_set
            or r["n_markers"] < MIN_PLOT_MARKERS
        ):
            continue
        out[r["depth"]].append(
            (
                r["n_markers"],
                r["lod_pct"],
                r["lod_pct_ci_lo"],
                r["lod_pct_ci_hi"],
                bool(r.get("censored", False)),
            )
        )
    return out


def _finite_xy(
    cell: list[tuple],
) -> tuple[list[int], list[float], list[float], list[float], list[bool]]:
    """Sort by panel size and drop non-finite / non-positive LoD points."""
    out_x, out_y, out_lo, out_hi, out_c = [], [], [], [], []
    for x, y, lo, hi, c in sorted(cell, key=lambda t: t[0]):
        if y is None or not math.isfinite(y) or y <= 0:
            continue
        out_x.append(x)
        out_y.append(y)
        out_lo.append(lo if (lo is not None and math.isfinite(lo) and lo > 0) else y)
        out_hi.append(hi if (hi is not None and math.isfinite(hi) and hi > 0) else y)
        out_c.append(bool(c))
    return out_x, out_y, out_lo, out_hi, out_c


def _draw_curve(ax, x, y, lo, hi, censored, color, marker, x_jitter: float = 1.0) -> None:
    """Draw one depth's LoD curve, distinguishing resolved from censored points.

    Resolved points get a filled marker and the 10-90% band shaded across them.
    Censored points (true LoD at or below the value) get an X, which stays legible
    where several depths pile onto the same censored value. ``x_jitter`` is a
    per-depth multiplicative offset that fans overlapping curves apart.
    """
    x = [xi * x_jitter for xi in x]
    ax.plot(x, y, "-", color=color, linewidth=1.8, zorder=3)

    res = [i for i, c in enumerate(censored) if not c]
    cen = [i for i, c in enumerate(censored) if c]

    if res:
        ax.plot(
            [x[i] for i in res],
            [y[i] for i in res],
            marker,
            color=color,
            markersize=6,
            linestyle="none",
            zorder=4,
        )
        bx = [x[i] for i in res]
        if len(bx) >= 2:
            ax.fill_between(
                bx,
                [lo[i] for i in res],
                [hi[i] for i in res],
                color=color,
                alpha=0.15,
                linewidth=0,
            )
    if cen:
        ax.plot(
            [x[i] for i in cen],
            [y[i] for i in cen],
            "X",
            color=color,
            markersize=8,
            markeredgecolor="white",
            markeredgewidth=0.6,
            linestyle="none",
            zorder=5,
        )


def _depth_colors(depths: list[int]) -> dict[int, tuple]:
    """viridis_r palette: deepest depth -> dark purple, shallowest -> teal."""
    cmap = plt.get_cmap("viridis_r")
    if len(depths) > 1:
        return {d: cmap(0.35 + 0.65 * i / (len(depths) - 1)) for i, d in enumerate(depths)}
    return {depths[0]: cmap(1.0)}


def _style_axis(ax, nmarkers: list[int]) -> None:
    """Apply the shared log-log panel styling (axis labels set by the caller)."""
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xticks(nmarkers)
    ax.set_xticklabels([str(n) for n in nmarkers])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.grid(True, which="both", alpha=0.2)
    ax.yaxis.set_major_locator(FixedLocator([0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 5, 10]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_formatter(FuncFormatter(_format_pct))
    ax.set_ylim(0.1, 10.0)
    for thr in THRESHOLDS:
        ax.axhline(thr, color="black", linestyle="--", linewidth=0.7, alpha=0.35)


def _set_pct(mixture_set: str) -> float:
    """Titration floor (percent) encoded in a set name like 'to_1pct'; inf for 'all'."""
    if mixture_set == "all":
        return float("inf")
    try:
        return float(mixture_set.replace("to_", "").replace("pct", ""))
    except ValueError:
        return float("inf")


def _set_nmix(rows: list[dict], mixture_set: str) -> int | None:
    """Mixture count recorded for a mixture set (for row labels)."""
    for r in rows:
        if r.get("mixture_set", "all") == mixture_set:
            try:
                return int(r.get("n_mixtures") or 0)
            except ValueError:
                return None
    return None


def _row_label(rows: list[dict], mixture_set: str) -> str:
    """Human row label, e.g. 'All mixtures (n=10)' or 'Titrated to 0.5% (n=3)'."""
    n = _set_nmix(rows, mixture_set)
    suffix = f" (n={n})" if n else ""
    if mixture_set == "all":
        return f"All mixtures{suffix}"
    frac = mixture_set.replace("to_", "").replace("pct", "")
    return f"Titrated to {frac}%{suffix}"


def _depth_jitter(depths: list[int]) -> dict[int, float]:
    """Per-depth multiplicative x-offset in [1-JITTER, 1+JITTER], centred on 1."""
    n = len(depths)
    if n < 2:
        return {depths[0]: 1.0} if depths else {}
    return {d: (1.0 + JITTER) ** (2 * i / (n - 1) - 1) for i, d in enumerate(depths)}


def _draw_panel(ax, rows, test, mixture_set, depths, colors, marker, nmarkers, show_analytical):
    """Draw all depth curves for one (test, mixture_set) panel."""
    jitter = _depth_jitter(depths)
    by_d = _by_depth(rows, test, mixture_set)
    for depth in depths:
        x, y, lo, hi, cens = _finite_xy(by_d.get(depth, []))
        if not x:
            continue
        _draw_curve(ax, x, y, lo, hi, cens, colors[depth], marker, x_jitter=jitter[depth])
    if show_analytical and test == "mle":
        analytical = _by_depth(rows, "mle_analytical", mixture_set)
        for depth in depths:
            ax_x, ax_y, *_ = _finite_xy(analytical.get(depth, []))
            if ax_x:
                ax.plot(ax_x, ax_y, "--", color=colors[depth], linewidth=1.0, alpha=0.55)
    _style_axis(ax, nmarkers)


def plot_grid(rows: list[dict], out_path: Path, show_analytical: bool = False) -> None:
    """Draw the 2x2 real-data LoD grid (test x mixture set)."""
    depths = sorted({r["depth"] for r in rows})
    nmarkers = sorted({r["n_markers"] for r in rows if r["n_markers"] >= MIN_PLOT_MARKERS})
    if not depths or not nmarkers:
        raise SystemExit("No plottable rows in summary")
    colors = _depth_colors(depths)

    # Two disjoint rows, ordered by titration floor descending: the mixtures that
    # stop at 1% (top) and those titrated to 0.5% (bottom). The "all" set stays in
    # the CSV for the headline facts but is not drawn, so the rows do not overlap.
    row_sets = sorted(
        {r.get("mixture_set", "all") for r in rows} - {"all"},
        key=_set_pct,
        reverse=True,
    )
    if not row_sets:
        row_sets = ["all"]

    fig, axes = plt.subplots(
        len(row_sets),
        2,
        figsize=(11.5, 4.7 * len(row_sets)),
        sharex=True,
        sharey=True,
        squeeze=False,
    )

    for i, mset in enumerate(row_sets):
        for j, (test, _title, marker) in enumerate(COLUMNS):
            _draw_panel(
                axes[i][j], rows, test, mset, depths, colors, marker, nmarkers, show_analytical
            )

    for j, (_test, title, _marker) in enumerate(COLUMNS):
        axes[0][j].set_title(title, fontsize=12, fontweight="bold")
    for i, mset in enumerate(row_sets):
        axes[i][0].set_ylabel("Limit of detection (% minor)", fontsize=11)
        axes[i][0].annotate(
            _row_label(rows, mset),
            xy=(-0.27, 0.5),
            xycoords="axes fraction",
            rotation=90,
            ha="center",
            va="center",
            fontsize=13,
            fontweight="bold",
        )
    for j in range(2):
        axes[-1][j].set_xlabel("Panel size (markers)", fontsize=11)

    depth_handles = [
        plt.Line2D([0], [0], color=colors[d], linewidth=1.8, label=f"{d}x") for d in depths
    ]
    axes[0][1].legend(
        handles=depth_handles, title="Mean depth", fontsize=9, loc="upper right", framealpha=0.9
    )
    style_handles = [
        plt.Line2D([0], [0], color="grey", marker="o", linewidth=1.8, label="resolved LoD"),
        plt.Line2D(
            [0], [0], color="grey", marker="X", linestyle="none", label="LoD ≤ marker (below grid)"
        ),
    ]
    if show_analytical:
        style_handles.append(
            plt.Line2D(
                [0], [0], color="grey", linestyle="--", linewidth=1.0, alpha=0.7,
                label="analytical LoD (Fisher info)",
            )
        )
    axes[0][0].legend(handles=style_handles, fontsize=8.5, loc="lower left", framealpha=0.9)

    fig.suptitle(
        "Real-data limit of detection (SRP434573 titrated mixtures)", fontsize=13, y=1.0
    )
    fig.text(
        0.5,
        0.965,
        "sub-sampled reads (pseudo-replicates); X marks LoD at or below the lowest titration "
        "(1% top row, 0.5% bottom row), not resolved lower; curves monotone-constrained "
        "(nested panels)",
        ha="center",
        va="top",
        fontsize=8.5,
        color="0.35",
    )
    fig.tight_layout(rect=(0.03, 0, 1, 0.96))
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}", file=sys.stderr)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--summary", default=str(FACTS_DIR / "subsample_lod_summary.csv"))
    parser.add_argument("--out", default=str(FACTS_DIR / "fig_subsample_lod_grid.png"))
    parser.add_argument(
        "--show-analytical",
        action="store_true",
        help="Overlay the analytical (Fisher-info) LoD on the MLE column (off by default).",
    )
    args = parser.parse_args(argv)

    rows = _read_summary(Path(args.summary))
    if not rows:
        raise SystemExit(f"No rows in {args.summary}")
    plot_grid(rows, Path(args.out), show_analytical=args.show_analytical)
    return 0


if __name__ == "__main__":
    sys.exit(main())
