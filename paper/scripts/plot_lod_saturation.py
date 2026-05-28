#!/usr/bin/env python3
"""LoD-vs-depth saturation: reconcile the simulated and the real per-sample LoD.

The per-sample LoD is set by the variance of the donor-fraction estimate, which
for a beta-binomial model is

    var(VAF) = p(1-p) * (n + rho) / (n * (rho + 1))   (n = depth, rho = overdispersion)

As depth grows this does not keep falling: it approaches p(1-p)/(rho+1), so the
LoD saturates and the effective depth is capped near rho+1. A pure-binomial
model (rho -> inf) keeps improving as 1/sqrt(n).

allomix's simulator draws read counts from a binomial (see simulate.py), so its
empirical LoD behaves like a high-rho (near-binomial) curve and reaches ~0.1-0.5%
at deep coverage. Real samples carry extra-binomial noise, so the estimator fits
a much smaller rho and the per-sample LoD lands near ~1% and barely improves with
depth. This script fits the LoD-vs-depth model

    LoD(n, M) = (A / sqrt(M)) * sqrt((n + rho) / (n * (rho + 1)))

(M = informative markers) to the in-silico sweep in lod_summary.csv, draws the
fitted curves and their saturation floor, and reports the fitted effective rho.
With --overlay-batch it also places real per-sample LoDs on the same axes and
fits their (much lower) rho, so the simulated and real limits can be compared
directly. That single overlay rho is a crude lumped estimate: it assumes the
real samples share the simulated scale A and one relatedness, so mixing sibling
and unrelated donors (or low-depth outliers) skews it. For the authoritative
per-sample overdispersion use scripts/diagnose_sample.py, which reads each
sample's fitted rho. The overlay here is for a quick visual reconciliation, not
a paper number, and real patient points should stay out of the in-silico paper
figure (run it without --overlay-batch for that).

Output:
    output/facts/fig_lod_saturation.png
    output/facts/lod_saturation.csv   (fitted A, rho, floor per relatedness/overlay)
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402
from scipy.optimize import curve_fit  # noqa: E402

FACTS_DIR = Path("output/facts")

# Colour per simulated panel size (informative markers vary, so key on n_markers).
_PANEL_CMAP = plt.get_cmap("viridis")


def _lod_model(X: np.ndarray, A: float, rho: float) -> np.ndarray:
    """LoD (donor %) as a function of depth n and marker count M.

    Args:
        X: 2-row array ``[n, M]`` (depth, informative markers).
        A: Overall scale (absorbs z, p(1-p), marker slopes).
        rho: Beta-binomial overdispersion (large = near binomial).

    Returns:
        Predicted LoD in percent.
    """
    n, m = X
    return (A / np.sqrt(m)) * np.sqrt((n + rho) / (n * (rho + 1.0)))


def _fit(n: np.ndarray, m: np.ndarray, lod: np.ndarray, fixed_a: float | None = None):
    """Fit (A, rho) of the LoD model. If ``fixed_a`` is set, fit rho only.

    Returns:
        ``(A, rho)``.
    """
    if fixed_a is not None:
        popt, _ = curve_fit(
            lambda X, rho: _lod_model(X, fixed_a, rho),
            np.vstack([n, m]),
            lod,
            p0=[1000.0],
            bounds=(1.0, 1e7),
            maxfev=20000,
        )
        return fixed_a, float(popt[0])
    a0 = float(np.median(lod * np.sqrt(m * n)))
    popt, _ = curve_fit(
        _lod_model,
        np.vstack([n, m]),
        lod,
        p0=[a0, 1000.0],
        bounds=([1e-6, 1.0], [1e6, 1e7]),
        maxfev=20000,
    )
    return float(popt[0]), float(popt[1])


def _read_lod_summary(
    path: Path, relatedness: str
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Read depth, marker count, informative count, and LoD% for one relatedness."""
    depth, n_markers, minf, lod = [], [], [], []
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r["relatedness"] != relatedness:
                continue
            lod_v = r["lod_pct"]
            if lod_v in ("", "NA"):
                continue
            depth.append(float(r["depth"]))
            n_markers.append(int(r["n_markers"]))
            minf.append(float(r["mean_n_informative"]))
            lod.append(float(lod_v))
    return np.array(depth), np.array(n_markers), np.array(minf), np.array(lod)


def _read_batch(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read mean_depth, informative markers, and LoD% from a monitor batch.tsv.

    Skips QC-FAIL rows and rows without a numeric lod_pct.
    """
    depth, minf, lod = [], [], []
    with path.open(encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            status = (r.get("qc_status") or r.get("qc_pass") or "PASS").upper()
            if status == "FAIL":
                continue
            lod_v = r.get("lod_pct")
            if lod_v in (None, "", "NA"):
                continue
            depth.append(float(r["mean_depth"]))
            minf.append(float(r["n_informative"]))
            lod.append(float(lod_v))
    return np.array(depth), np.array(minf), np.array(lod)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--lod-summary", type=Path, default=FACTS_DIR / "lod_summary.csv")
    parser.add_argument(
        "--relatedness",
        default="unrelated",
        choices=("unrelated", "sibling"),
        help="Which simulated relatedness to fit/plot (default unrelated).",
    )
    parser.add_argument(
        "--overlay-batch",
        type=Path,
        default=None,
        help="Optional real monitor batch.tsv to overlay (e.g. output/validation_run3/batch.tsv). "
        "Internal/real-data use; leave off for the in-silico paper figure.",
    )
    parser.add_argument(
        "--overlay-label", default="real (run3)", help="Legend label for the overlay."
    )
    parser.add_argument("--out", type=Path, default=FACTS_DIR / "fig_lod_saturation.png")
    parser.add_argument("--facts", type=Path, default=FACTS_DIR / "lod_saturation.csv")
    args = parser.parse_args()

    n, nm, minf, lod = _read_lod_summary(args.lod_summary, args.relatedness)
    if len(n) == 0:
        sys.exit(f"No '{args.relatedness}' rows with a numeric lod_pct in {args.lod_summary}")
    a_sim, rho_sim = _fit(n, minf, lod)

    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    depth_grid = np.logspace(math.log10(n.min() * 0.8), math.log10(n.max() * 1.25), 200)

    panels = sorted(set(zip(nm.tolist(), minf.tolist())))
    for n_mark, m_inf in panels:
        sel = nm == n_mark
        color = _PANEL_CMAP((panels.index((n_mark, m_inf)) + 0.5) / len(panels))
        ax.scatter(n[sel], lod[sel], color=color, s=36, zorder=3)
        ax.plot(
            depth_grid,
            _lod_model(np.vstack([depth_grid, np.full_like(depth_grid, m_inf)]), a_sim, rho_sim),
            color=color,
            lw=1.5,
            label=f"sim panel {n_mark} (M≈{m_inf:.0f})",
        )

    # Saturation floor + binomial reference for the largest panel.
    big_m = max(m for _, m in panels)
    floor_big = (a_sim / math.sqrt(big_m)) / math.sqrt(rho_sim + 1.0)
    ax.axhline(floor_big, color="0.5", ls=":", lw=1.0)
    ax.text(
        depth_grid[0],
        floor_big,
        f" sim floor {floor_big:.3f}% (M≈{big_m:.0f})",
        color="0.4",
        fontsize=8,
        va="bottom",
    )
    ax.plot(
        depth_grid,
        (a_sim / math.sqrt(big_m)) / np.sqrt(depth_grid),
        color="0.6",
        ls="--",
        lw=1.0,
        label="binomial (rho→∞)",
    )

    rho_real = None
    if args.overlay_batch:
        rd, rminf, rlod = _read_batch(args.overlay_batch)
        if len(rd):
            ax.scatter(
                rd,
                rlod,
                marker="D",
                s=70,
                color="#d62728",
                edgecolors="white",
                linewidths=0.8,
                alpha=0.9,
                zorder=6,
                label=args.overlay_label,
            )
            if len(rd) >= 4:
                _, rho_real = _fit(rd, rminf, rlod, fixed_a=a_sim)
                med_m = float(np.median(rminf))
                grid_m = np.full_like(depth_grid, med_m)
                ax.plot(
                    depth_grid,
                    _lod_model(np.vstack([depth_grid, grid_m]), a_sim, rho_real),
                    color="#d62728",
                    lw=1.5,
                    ls="-.",
                    label=f"{args.overlay_label} fit (M≈{med_m:.0f})",
                )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{v:g}"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{v:g}"))
    ax.set_xlabel("Mean depth (x)")
    ax.set_ylabel("Limit of detection (donor %)")
    title = f"LoD vs depth ({args.relatedness}): simulated rho≈{rho_sim:.0f}"
    if rho_real is not None:
        title += f", {args.overlay_label} rho≈{rho_real:.0f}"
    ax.set_title(title)
    ax.legend(fontsize=7.5, loc="upper right", framealpha=0.9)
    ax.grid(True, which="both", alpha=0.2)
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150)
    print(f"Wrote {args.out}")

    with args.facts.open("w", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(
            ["dataset", "relatedness", "A", "rho", "eff_depth_cap", "floor_pct_largest_panel"]
        )
        w.writerow(
            [
                "simulated",
                args.relatedness,
                f"{a_sim:.4f}",
                f"{rho_sim:.1f}",
                f"{rho_sim + 1:.0f}",
                f"{floor_big:.4f}",
            ]
        )
        if rho_real is not None:
            floor_real = (a_sim / math.sqrt(float(np.median(rminf)))) / math.sqrt(rho_real + 1.0)
            w.writerow(
                [
                    args.overlay_label,
                    "real",
                    f"{a_sim:.4f}",
                    f"{rho_real:.1f}",
                    f"{rho_real + 1:.0f}",
                    f"{floor_real:.4f}",
                ]
            )
    print(f"Wrote {args.facts}")


if __name__ == "__main__":
    main()
