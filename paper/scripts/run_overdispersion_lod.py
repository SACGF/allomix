#!/usr/bin/env python3
"""Quantify how much sequencing overdispersion (rho) moves the limit of detection.

allomix's simulator can now draw reads beta-binomial (simulate.sample_allele_counts
takes a ``rho``); pure-binomial sampling (rho -> inf) was hiding the extra-binomial
noise that real data carries, making the in-silico LoD optimistic. This script
sweeps rho at a fixed depth and panel size, simulates blank (pure-host) replicates
through the real estimator, and reports three LoD views per rho:

  - empirical LoB: the 95th percentile of the estimated donor fraction across
    blank replicates (the EP17 blank-based limit of blank);
  - tool LoD:      the estimator's own per-sample lod_fraction (median over reps),
    using the rho it fits from each blank;
  - analytic LoD:  detection_limit() on the panel at the true rho.

All three rise sharply as rho falls, so the figure shows directly how much the
choice of overdispersion matters (and the rho -> inf binomial reference shows how
optimistic the old simulator was).

Output:
    output/facts/fig_overdispersion_lod.png
    output/facts/overdispersion_lod.csv
"""

import argparse
import csv
import math
import random
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from paper_quick import qval  # noqa: E402  (also patches savefig for the watermark)

from allomix.chimerism import detection_limit, estimate_single_donor_bb  # noqa: E402
from allomix.genotype import InformativeMarker  # noqa: E402
from allomix.simulate import alt_dose, sample_allele_counts  # noqa: E402

FACTS_DIR = Path("output/facts")
ERROR_RATE = 0.01


def build_panel(n_markers: int, rng: random.Random) -> list[dict]:
    """Draw an unrelated host/donor panel of informative markers (host != donor).

    ``rng`` fixes the genotypes once; they are reused across rho/reps. Each dict
    carries host_gt, donor_gt, and host_vaf (the pure-host VAF at f=0).
    """
    panel: list[dict] = []
    pos = 0
    while len(panel) < n_markers:
        pos += 1
        p_alt = rng.uniform(0.2, 0.8)
        host_gt = _draw(p_alt, rng)
        donor_gt = _draw(p_alt, rng)
        if alt_dose(host_gt) == alt_dose(donor_gt):
            continue  # not informative
        panel.append(
            {
                "chrom": "chr1",
                "pos": pos,
                "host_gt": host_gt,
                "donor_gt": donor_gt,
                "host_vaf": alt_dose(host_gt) / 2.0,  # pure-host (blank) VAF
            }
        )
    return panel


def _draw(p_alt: float, rng: random.Random) -> tuple[int, int]:
    """Draw a diploid genotype from Hardy-Weinberg equilibrium."""
    return tuple(sorted(1 if rng.random() < p_alt else 0 for _ in range(2)))  # type: ignore[return-value]


def _markers_for_panel(panel: list[dict], depth: int) -> list[InformativeMarker]:
    """InformativeMarker objects at the given depth, for the analytic detection_limit.

    fraction_se reads the Fisher information from n = admix_ad_ref + admix_ad_alt,
    so the counts must sum to ``depth`` (their split is irrelevant at f=0); a
    zero-depth marker carries no information and forces an infinite LoD.
    """
    return [
        InformativeMarker(
            chrom=m["chrom"],
            pos=m["pos"],
            ref="A",
            alt="T",
            host_gt=m["host_gt"],
            donor_gts=[m["donor_gt"]],
            marker_type=0,
            admix_ad_ref=depth,
            admix_ad_alt=0,
            admix_dp=depth,
        )
        for m in panel
    ]


def simulate_blank(
    panel: list[dict], depth: int, rho: float, rng: random.Random
) -> list[InformativeMarker]:
    """Build one blank (pure-host) admixture replicate at the given overdispersion."""
    markers = []
    for m in panel:
        ref_c, alt_c = sample_allele_counts(m["host_vaf"], depth, rng, ERROR_RATE, rho)
        markers.append(
            InformativeMarker(
                chrom=m["chrom"],
                pos=m["pos"],
                ref="A",
                alt="T",
                host_gt=m["host_gt"],
                donor_gts=[m["donor_gt"]],
                marker_type=0,
                admix_ad_ref=ref_c,
                admix_ad_alt=alt_c,
                admix_dp=ref_c + alt_c,
            )
        )
    return markers


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--depth", type=int, default=1000, help="Mean depth (default 1000)")
    parser.add_argument(
        "--markers", type=int, default=150, help="Informative markers (default 150)"
    )
    parser.add_argument(
        "--n-reps", type=int, default=qval(60, 10),
        help="Blank replicates per rho (default 60; 10 in quick-build mode)",
    )
    parser.add_argument(
        "--rhos",
        type=float,
        nargs="+",
        default=[30, 100, 300, 1000, 3000, 10000],
        help="Finite rho grid to sweep (binomial rho->inf is added as a reference).",
    )
    parser.add_argument("--seed", type=int, default=20260528)
    parser.add_argument("--out", type=Path, default=FACTS_DIR / "fig_overdispersion_lod.png")
    parser.add_argument("--facts", type=Path, default=FACTS_DIR / "overdispersion_lod.csv")
    parser.add_argument(
        "--headline", type=Path, default=FACTS_DIR / "overdispersion_lod_headline.csv"
    )
    args = parser.parse_args()

    rng = random.Random(args.seed)
    panel = build_panel(args.markers, rng)
    info_markers = _markers_for_panel(panel, args.depth)

    rhos = [*sorted(args.rhos), float("inf")]
    rows = []
    for rho in rhos:
        est_fracs, tool_lods = [], []
        for rep in range(args.n_reps):
            rep_rng = random.Random(hash((args.seed, rho, rep)) & 0xFFFFFFFF)
            markers = simulate_blank(panel, args.depth, rho, rep_rng)
            result = estimate_single_donor_bb(markers, error_rate=ERROR_RATE)
            est_fracs.append(result.donor_fraction * 100.0)
            if math.isfinite(result.lod_fraction):
                tool_lods.append(result.lod_fraction * 100.0)
        emp_lob = float(np.percentile(est_fracs, 95))
        tool_lod = float(np.median(tool_lods)) if tool_lods else float("nan")
        _, an_lod = detection_limit(info_markers, ERROR_RATE, rho, None)
        rows.append(
            {
                "rho": rho,
                "eff_depth_cap": "inf" if math.isinf(rho) else f"{rho + 1:.0f}",
                "emp_lob_pct": emp_lob,
                "tool_lod_pct": tool_lod,
                "analytic_lod_pct": an_lod * 100.0,
            }
        )
        tag = "inf" if math.isinf(rho) else f"{rho:.0f}"
        print(
            f"rho={tag:>6}  emp_LoB={emp_lob:.3f}%  tool_LoD={tool_lod:.3f}%  "
            f"analytic_LoD={an_lod * 100:.3f}%"
        )

    _plot(rows, args)
    _write_facts(rows, args)
    _write_headline(rows, args)


def _plot(rows: list[dict], args: argparse.Namespace) -> None:
    finite = [r for r in rows if math.isfinite(r["rho"])]
    rho_x = [r["rho"] for r in finite]
    fig, ax = plt.subplots(figsize=(8.0, 6.0))
    ax.plot(
        rho_x, [r["analytic_lod_pct"] for r in finite], "-", color="#1f77b4", label="analytic LoD"
    )
    ax.plot(
        rho_x, [r["tool_lod_pct"] for r in finite], "o-", color="#d62728", label="tool LoD (sim)"
    )
    ax.plot(
        rho_x,
        [r["emp_lob_pct"] for r in finite],
        "s--",
        color="#2ca02c",
        label="empirical LoB (sim)",
    )

    binom = next(r for r in rows if math.isinf(r["rho"]))
    ax.axhline(binom["analytic_lod_pct"], color="0.5", ls=":", lw=1.0)
    ax.text(
        rho_x[0],
        binom["analytic_lod_pct"],
        f" binomial (rho→∞) LoD {binom['analytic_lod_pct']:.3f}%",
        color="0.4",
        fontsize=8,
        va="bottom",
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.invert_xaxis()  # more overdispersion (small rho) on the right
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{v:g}"))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _p: f"{v:g}"))
    ax.set_xlabel("Overdispersion concentration rho  (← more overdispersion)")
    ax.set_ylabel("Limit of detection (donor %)")
    ax.set_title(
        f"How much overdispersion matters: LoD vs rho "
        f"(depth={args.depth}x, M={args.markers}, unrelated)"
    )
    ax.legend(fontsize=8, loc="upper left", framealpha=0.9)
    ax.grid(True, which="both", alpha=0.2)
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150)
    print(f"Wrote {args.out}")


def _write_facts(rows: list[dict], args: argparse.Namespace) -> None:
    with args.facts.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "depth",
                "n_markers",
                "rho",
                "eff_depth_cap",
                "emp_lob_pct",
                "tool_lod_pct",
                "analytic_lod_pct",
            ]
        )
        for r in rows:
            w.writerow(
                [
                    args.depth,
                    args.markers,
                    "inf" if math.isinf(r["rho"]) else f"{r['rho']:.0f}",
                    r["eff_depth_cap"],
                    f"{r['emp_lob_pct']:.4f}",
                    f"{r['tool_lod_pct']:.4f}",
                    f"{r['analytic_lod_pct']:.4f}",
                ]
            )
    print(f"Wrote {args.facts}")


def _write_headline(rows: list[dict], args: argparse.Namespace) -> None:
    """Single-row headline facts for the paper (tool LoD at a few reference rhos)."""

    def lod_at(rho_target: float) -> float:
        return min(rows, key=lambda r: abs(r["rho"] - rho_target))["tool_lod_pct"]

    binom = next(r for r in rows if math.isinf(r["rho"]))["tool_lod_pct"]
    lod100 = lod_at(100)
    with args.headline.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "depth",
                "n_markers",
                "lod_binomial_pct",
                "lod_rho1000_pct",
                "lod_rho100_pct",
                "lod_rho30_pct",
                "fold_rho100_vs_binomial",
            ]
        )
        w.writerow(
            [
                args.depth,
                args.markers,
                f"{binom:.3f}",
                f"{lod_at(1000):.3f}",
                f"{lod100:.3f}",
                f"{lod_at(30):.3f}",
                f"{lod100 / binom:.1f}",
            ]
        )
    print(f"Wrote {args.headline}")


if __name__ == "__main__":
    main()
