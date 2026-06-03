#!/usr/bin/env python3
"""Plot simulated host-presence detection probability against spiked level.

Internal SA Path validation plot (not part of the allomix package).

Reads the in-silico presence-LOD facts written by the paper build
(``output/facts/presence_lod_summary.csv`` and its ``_rho_`` counterpart) and
draws detection probability as a function of the spiked host fraction. Each
summary row carries a semicolon-joined ``pos_fractions`` list and matching
``presence_det_rates`` (fraction of replicates the presence test calls
positive), so one row becomes one curve.

The point of the figure is the low-level falloff: how far down in host fraction
the test keeps detecting, and how much beta-binomial overdispersion (rho) moves
that limit. The plain summary is pure binomial sampling (optimistic); the
``_rho_`` summary applies the calibrated overdispersion, which is the dominant
control on the real limit. Both are drawn so the gap is visible.

Usage:
    python scripts/plot_presence_lod_curve.py \
        --summary output/facts/presence_lod_summary.csv \
        --rho-summary output/facts/presence_lod_rho_summary.csv \
        --depth 1000 --relatedness sibling \
        --output output/presence_lod_curve.png
"""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# Detection-probability target used to read an LOD off the curve.
DET_TARGET = 0.95

# One colour per error rate, assigned in ascending order.
ERR_COLORS = {
    0.0003: "#1b7837",
    0.001: "#2c7fb8",
    0.003: "#d95f0e",
    0.01: "#c0392b",
}


def _floats(field: str) -> list[float]:
    return [float(x) for x in field.split(";") if x]


def read_curves(path: Path, depth: int, relatedness: str) -> list[dict]:
    """Read summary rows matching depth and relatedness into curve dicts.

    Returns:
        List of dicts with error_rate, fractions (as %), det_rates and
        fp_rate (false-positive rate of the LRT at the zero-host blanks).
    """
    out: list[dict] = []
    with open(path) as f:
        for r in csv.DictReader(f):
            if int(r["depth"]) != depth or r["relatedness"] != relatedness:
                continue
            out.append(
                {
                    "error_rate": float(r["error_rate"]),
                    "frac_pct": [v * 100.0 for v in _floats(r["pos_fractions"])],
                    "det_rates": _floats(r["presence_det_rates"]),
                    "fp_rate": float(r["fp_rate_lrt"]),
                }
            )
    return out


def _draw_panel(ax, curves: list[dict], title: str) -> None:
    """Draw one model's detection curves (one line per error rate).

    The x position is the spiked host fraction on a log scale, but the axis is
    labelled in donor % (= 100 - host %, log-spaced by distance from 100%) to
    match the field convention and the other allomix plots. Full donor sits at
    the left; detection climbs as the donor fraction drops.
    """
    for c in sorted(curves, key=lambda d: d["error_rate"]):
        color = ERR_COLORS.get(c["error_rate"], "#666666")
        ax.plot(
            c["frac_pct"],
            c["det_rates"],
            "-",
            color=color,
            marker="o",
            ms=4,
            lw=1.7,
            label=f"err {c['error_rate'] * 100:g}%",
        )
    ax.axhline(DET_TARGET, color="#888888", lw=0.9, ls=":", zorder=1)
    ax.set_xscale("log")
    # Decade ticks keep the donor-% labels (99.999 / 99.99 / 99.9) from
    # crowding; the per-fraction markers still show the sampled levels.
    ax.set_xticks([0.001, 0.01, 0.1])
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{100 - v:g}"))
    ax.set_xlabel("Donor fraction (%), log-spaced by distance from 100%")
    ax.set_ylim(-0.03, 1.03)
    ax.set_title(title, fontsize=10)
    ax.grid(color="#eeeeee", lw=0.6)


def plot(
    binom: list[dict],
    rho: list[dict],
    depth: int,
    relatedness: str,
    output: Path,
) -> None:
    """Draw binomial vs overdispersed detection panels and save to ``output``."""
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 5.5), sharey=True)

    _draw_panel(axL, binom, "Binomial sampling (optimistic)")
    _draw_panel(axR, rho, "Beta-binomial overdispersion (rho)")
    axL.set_ylabel("Detection probability")
    axR.text(
        axR.get_xlim()[0],
        DET_TARGET,
        f" {DET_TARGET:.0%} detection",
        color="#666666",
        fontsize=8,
        va="bottom",
    )
    axL.legend(fontsize=9, loc="upper left", title="background error")

    fig.suptitle(
        f"Host-presence detection vs donor fraction — {relatedness} donor, {depth}x depth",
        fontsize=12,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(output, dpi=150)
    print(f"Wrote {output}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--summary", type=Path, default=Path("output/facts/presence_lod_summary.csv"))
    ap.add_argument(
        "--rho-summary", type=Path, default=Path("output/facts/presence_lod_rho_summary.csv")
    )
    ap.add_argument("--depth", type=int, default=1000, help="Depth slice to plot (default 1000)")
    ap.add_argument(
        "--relatedness",
        default="sibling",
        choices=["sibling", "unrelated"],
        help="Relatedness slice (default sibling, the HSCT-relevant harder case)",
    )
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    binom = read_curves(args.summary, args.depth, args.relatedness)
    rho = read_curves(args.rho_summary, args.depth, args.relatedness)
    if not binom:
        raise SystemExit(f"No rows for depth={args.depth} relatedness={args.relatedness}")
    plot(binom, rho, args.depth, args.relatedness, args.output)


if __name__ == "__main__":
    main()
