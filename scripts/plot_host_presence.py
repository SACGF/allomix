#!/usr/bin/env python3
"""Plot allomix host-presence detection: effect size and evidence together.

Internal SA Path validation plot (not part of the allomix package).

Takes an allomix ``batch.tsv`` (as produced by ``scripts/run_csv_batch.py``)
and draws a per-sample summary of the host-presence test. Each sample
contributes two linked quantities:

  - the estimated host fraction (effect size) with its confidence interval,
  - the presence-test p-value (evidence the fraction is above noise).

These have very different scales (donor % sitting just below 100 vs a
p-value spanning many orders of magnitude), so the effect size is drawn as a
forest and the p-value is listed as text rather than forced onto a second
axis. Donor % is on the Y axis, log-spaced by distance from 100% donor and
inverted so full donor is at the top (the same style as
plot_chimerism_comparison.py), with samples across X and the exact presence
p-value under each. Samples are coloured by whether the test calls host
present.

A p-value bar was rejected on purpose: its magnitude saturates (1e-44 and
1e-124 both just mean "certainly there"), so a bar wastes its range on a
distinction nobody acts on and crushes the only one that matters, whether p
crosses alpha. Listing the value keeps it exact; colour carries the call.

Convention: the effect axis is reported as % DONOR (= 100 - host fraction),
the standard of the field and consistent with every other allomix output.
The presence test itself estimates the host fraction; donor % is just
100 - that.

Usage:
    python scripts/plot_host_presence.py \
        output/validation_run9/batch.tsv \
        --output output/host_presence_run9.png
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# Significance threshold for calling host present. Samples at or below this
# p-value are drawn as "detected"; above it as "not detected".
ALPHA = 0.05

# p-values can underflow to exactly 0 (evidence beyond float range). Floor
# them here so -log10 is finite, and flag the bar so a capped value is not
# read as a precise one.
P_FLOOR = 1e-300

# Floor for the donor axis, expressed as host % (= 100 - donor %). Pure-donor
# samples have an estimated host fraction of 0, which a log axis cannot show;
# clamp them to this floor so they plot at the full-donor end. Matches the
# 0.02 default in plot_chimerism_comparison.py (donor 99.98%).
HOST_FLOOR_PCT = 0.02

COLOR_DETECTED = "#1b7837"  # green
COLOR_NOT = "#999999"  # grey


def short_label(name: str) -> str:
    """Patient code plus the leading numeric id, e.g. "REDACTED #5".

    The id keeps repeated codes (two REDACTED timepoints, two REDACTED panels)
    distinct. Falls back to the raw name if no upper-case code is found.
    """
    leading = name.split("_", 1)[0]
    codes = [t for t in name.split("_") if t.isalpha() and t.isupper() and len(t) >= 3]
    code = codes[-1] if codes else name
    return f"{code} #{leading}" if leading.isdigit() else code


def read_rows(path: Path) -> list[dict]:
    """Read batch.tsv rows that carry host-presence columns, parsed numeric.

    Returns:
        List of dicts with label, f_pct, ci_lo_pct, ci_hi_pct, p, markers,
        capped (whether p hit the floor) and detected. Rows without a
        host_f_est value are skipped.
    """
    rows: list[dict] = []
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            raw = r.get("host_f_est", "")
            if raw in (None, "", "NA"):
                continue
            p_raw = float(r["host_present_p"])
            capped = p_raw <= P_FLOOR
            p = max(p_raw, P_FLOOR)
            rows.append(
                {
                    "label": short_label(r["sample"]),
                    "f_pct": float(raw) * 100.0,
                    "ci_lo_pct": float(r["host_f_ci_lo"]) * 100.0,
                    "ci_hi_pct": float(r["host_f_ci_hi"]) * 100.0,
                    "p": p,
                    "capped": capped,
                    "markers": int(r["host_detect_markers"]),
                    "detected": p_raw <= ALPHA,
                }
            )
    return rows


def fmt_p(d: dict) -> str:
    """Format a presence p-value for the listing: exact, readable, honest.

    Underflowed values are shown as a bound, not a fake number; values near
    1 as a decimal; small values in scientific notation.
    """
    if d["capped"]:
        return "< 1e-300"
    p = d["p"]
    return f"{p:.2f}" if p >= 0.01 else f"{p:.1e}"


def plot(rows: list[dict], output: Path, title: str | None) -> None:
    """Draw the vertical host-presence forest and save it to ``output``.

    Donor % is on the Y axis (log-spaced by distance from 100%, inverted so
    full donor is at the top), samples across X, matching the house style of
    plot_chimerism_comparison.py. The exact presence p-value is listed along
    the bottom under each sample; colour carries the detected/not call.
    """
    # Sort so full donor (smallest host fraction) sits at the left.
    rows = sorted(rows, key=lambda d: d["f_pct"])
    n = len(rows)
    x = list(range(n))
    colors = [COLOR_DETECTED if d["detected"] else COLOR_NOT for d in rows]

    fig, ax = plt.subplots(figsize=(0.85 * n + 2.5, 7.5))

    # Plot internally as host % (= 100 - donor %) on a log axis so the action
    # near full donor is not compressed; the 100 - v formatter labels the axis
    # back in donor %. Pure-donor samples (host 0) clamp to the floor.
    def clamp(v: float) -> float:
        return max(v, HOST_FLOOR_PCT)

    for xi, d, c in zip(x, rows, colors):
        lo, hi = clamp(d["ci_lo_pct"]), clamp(d["ci_hi_pct"])
        ax.plot([xi, xi], [lo, hi], color=c, lw=1.6, solid_capstyle="round", zorder=2)
        ax.plot(xi, clamp(d["f_pct"]), "o", color=c, ms=7, zorder=3)

    ax.set_yscale("log")
    ax.set_ylim(105, HOST_FLOOR_PCT * 0.7)  # inverted: 100% donor at top
    ax.set_yticks([HOST_FLOOR_PCT, 0.1, 1, 10, 100])
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{100 - v:g}"))
    ax.set_ylabel("Donor fraction (%), log-spaced by distance from 100%")
    ax.set_xlim(-0.6, n - 0.4)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{d['label']}\n(m={d['markers']})" for d in rows], fontsize=8)
    ax.grid(axis="y", color="#eeeeee", lw=0.6, zorder=0)

    # Presence p-value listed along the bottom, one per sample column, coloured
    # by the call. Placed just below the axis in a blended transform (data x,
    # axes-fraction y) so it reads as a labelled row, not data.
    trans = ax.get_xaxis_transform()
    ax.text(
        -0.6, -0.16, "p =", fontsize=8, fontweight="bold", ha="right", va="center", transform=trans
    )
    for xi, d, c in zip(x, rows, colors):
        ax.text(xi, -0.16, fmt_p(d), fontsize=8, color=c, ha="center", va="center", transform=trans)

    handles = [
        plt.Line2D(
            [], [], marker="o", ls="", color=COLOR_DETECTED, label=f"host detected (p≤{ALPHA})"
        ),
        plt.Line2D([], [], marker="o", ls="", color=COLOR_NOT, label="not detected"),
    ]
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=9)

    if title:
        fig.suptitle(title, fontsize=12)
    fig.tight_layout(rect=(0, 0.08, 1, 0.97 if title else 1.0))
    fig.savefig(output, dpi=150)
    print(f"Wrote {output}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("batch_tsv", type=Path, help="allomix batch.tsv with host-presence columns")
    ap.add_argument("--output", type=Path, required=True, help="output PNG path")
    ap.add_argument("--title", default=None, help="optional figure title")
    args = ap.parse_args()

    rows = read_rows(args.batch_tsv)
    if not rows:
        raise SystemExit("No rows with host-presence columns found")
    plot(rows, args.output, args.title)


if __name__ == "__main__":
    main()
