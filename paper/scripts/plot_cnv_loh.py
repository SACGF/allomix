#!/usr/bin/env python3
"""Plot the effect of host copy-number aberrations on chimerism (issue #13).

Reads output/facts/cnv_loh_summary.csv and draws a grid: one column per
aberration kind (CN-LoH, deletion, gain), top row MAE vs burden, bottom row
95% CI coverage vs burden. Each kind's curve starts from the shared
no-aberration baseline at burden 0. Pure-clone cells (clonal_fraction = 1.0).

Output:
    output/facts/fig_cnv_loh.png
"""

from __future__ import annotations

import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

FACTS_DIR = Path("output/facts")
KINDS = ["cnloh", "deletion", "gain"]
KIND_LABELS = {"cnloh": "CN-LoH (copy-neutral)", "deletion": "Deletion (CN1)", "gain": "Gain (CN3)"}
COLORS = {"unrelated": "#1f77b4", "sibling": "#d62728"}


def load_summary(path: Path) -> list[dict]:
    rows = []
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            rows.append(
                {
                    "relatedness": r["relatedness"],
                    "kind": r["kind"],
                    "clonal_fraction": float(r["clonal_fraction"]),
                    "burden": float(r["burden"]),
                    "true_frac": float(r["true_frac"]),
                    "mae": float(r["mae"]),
                    "ci_coverage": float(r["ci_coverage"]),
                    "mae_robust": float(r["mae_robust"]),
                    "ci_coverage_robust": float(r["ci_coverage_robust"]),
                }
            )
    return rows


def series(rows: list[dict], kind: str, field: str, clonal: float) -> dict[str, dict[float, float]]:
    """Per-relatedness {burden: mean(field) over true_frac} for one kind.

    Burden 0 is taken from the shared baseline rows (kind == 'baseline').
    """
    acc: dict[str, dict[float, list[float]]] = defaultdict(lambda: defaultdict(list))
    for r in rows:
        if r["kind"] == "baseline":
            acc[r["relatedness"]][0.0].append(r[field])
        elif r["kind"] == kind and r["clonal_fraction"] == clonal:
            acc[r["relatedness"]][r["burden"]].append(r[field])
    return {rel: {b: statistics.fmean(v) for b, v in bur.items()} for rel, bur in acc.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=FACTS_DIR / "cnv_loh_summary.csv")
    parser.add_argument("--out", type=Path, default=FACTS_DIR / "fig_cnv_loh.png")
    parser.add_argument("--clonal", type=float, default=1.0)
    args = parser.parse_args()

    rows = load_summary(args.summary)

    fig, axes = plt.subplots(2, len(KINDS), figsize=(4.2 * len(KINDS), 8.0), sharex=True)
    mae_max = 0.0
    for j, kind in enumerate(KINDS):
        mae = series(rows, kind, "mae", args.clonal)
        cov = series(rows, kind, "ci_coverage", args.clonal)
        mae_r = series(rows, kind, "mae_robust", args.clonal)
        cov_r = series(rows, kind, "ci_coverage_robust", args.clonal)
        ax_mae, ax_cov = axes[0][j], axes[1][j]
        for rel in sorted(mae):
            bs = sorted(mae[rel])
            c = COLORS.get(rel)
            ax_mae.plot(bs, [mae[rel][b] for b in bs], marker="o", color=c, label=f"{rel} (std)")
            ax_mae.plot(bs, [mae_r[rel][b] for b in bs], marker="s", ls="--", color=c,
                        label=f"{rel} (robust)")
            ax_cov.plot(bs, [cov[rel][b] for b in bs], marker="o", color=c, label=f"{rel} (std)")
            ax_cov.plot(bs, [cov_r[rel][b] for b in bs], marker="s", ls="--", color=c,
                        label=f"{rel} (robust)")
            mae_max = max(mae_max, max(mae[rel].values()))
        ax_mae.set_title(KIND_LABELS[kind])
        ax_mae.grid(True, alpha=0.3)
        ax_cov.axhline(0.95, color="gray", linestyle="--", linewidth=1)
        ax_cov.set_ylim(0, 1.02)
        ax_cov.set_xlabel("burden (fraction of eligible markers)")
        ax_cov.grid(True, alpha=0.3)
        if j == 0:
            ax_mae.set_ylabel("MAE of donor fraction")
            ax_cov.set_ylabel("95% CI coverage")
            ax_mae.legend(fontsize=8)

    for j in range(len(KINDS)):
        axes[0][j].set_ylim(0, mae_max * 1.08)

    fig.suptitle(
        f"Host copy-number aberration impact, std vs robust refit "
        f"(pure clone, clonal={args.clonal:g})"
    )
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
