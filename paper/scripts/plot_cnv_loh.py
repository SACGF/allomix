#!/usr/bin/env python3
"""Plot the donor LoD inflation from host copy-number aberrations (issue #13).

Reads output/facts/cnv_loh_summary.csv and draws, per aberration kind, the donor
limit of detection (LoD, a donor fraction on a log axis, consistent with the
other LoD figures) versus the host CN-aberration burden. Standard estimator
(solid) vs robust refit (dashed), one colour per relatedness, with the
no-aberration baseline at burden 0.

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
# Donor fractions are probed up to 20%; a LoD above that is "undetectable here"
# and drawn at a ceiling marker above the plotted range.
MAX_PROBED_PCT = 20.0
CEILING_PCT = 35.0


def _f(x: str) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def load_summary(path: Path) -> list[dict]:
    rows = []
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            rows.append(
                {
                    "relatedness": r["relatedness"],
                    "kind": r["kind"],
                    "burden": float(r["burden"]),
                    "lod_std": _f(r["lod_std"]),
                    "lod_robust": _f(r["lod_robust"]),
                }
            )
    return rows


def series(rows: list[dict], kind: str, field: str) -> dict[str, dict[float, float]]:
    """Per-relatedness {burden: LoD in %}. Burden 0 comes from the shared baseline.

    Above-range LoD (inf) is mapped to a ceiling so it still plots; NaN (failed
    fit) is skipped.
    """
    acc: dict[str, dict[float, list[float]]] = defaultdict(lambda: defaultdict(list))
    for r in rows:
        v = r[field]
        if v != v:  # NaN: fit failed, skip
            continue
        pct = CEILING_PCT / 100.0 if v == float("inf") else v
        if r["kind"] == "baseline":
            acc[r["relatedness"]][0.0].append(pct)
        elif r["kind"] == kind:
            acc[r["relatedness"]][r["burden"]].append(pct)
    return {rel: {b: statistics.fmean(v) for b, v in bur.items()} for rel, bur in acc.items()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=FACTS_DIR / "cnv_loh_summary.csv")
    parser.add_argument("--out", type=Path, default=FACTS_DIR / "fig_cnv_loh.png")
    args = parser.parse_args()

    rows = load_summary(args.summary)

    fig, axes = plt.subplots(1, len(KINDS), figsize=(4.4 * len(KINDS), 4.6), sharey=True)
    for j, kind in enumerate(KINDS):
        std = series(rows, kind, "lod_std")
        rob = series(rows, kind, "lod_robust")
        ax = axes[j]
        for rel in sorted(std):
            c = COLORS.get(rel)
            bs = sorted(std[rel])
            ax.plot(bs, [std[rel][b] * 100 for b in bs], marker="o", color=c, label=f"{rel} (std)")
            if rel in rob:
                br = sorted(rob[rel])
                ax.plot(br, [rob[rel][b] * 100 for b in br], marker="s", ls="--", color=c,
                        label=f"{rel} (robust)")
        ax.set_yscale("log")
        ax.axhline(MAX_PROBED_PCT, color="gray", ls=":", lw=1)
        ax.set_ylim(top=CEILING_PCT * 1.4)
        ax.set_title(KIND_LABELS[kind])
        ax.set_xlabel("CN-aberration burden (fraction of eligible markers)")
        ax.grid(True, which="both", alpha=0.3)
        if j == 0:
            ax.set_ylabel("Donor limit of detection (%)")
            ax.legend(fontsize=8)
        if j == len(KINDS) - 1:
            ax.text(0.5, MAX_PROBED_PCT * 1.05, "above probed range (undetectable)",
                    fontsize=7, color="gray", ha="right", va="bottom")

    fig.suptitle("Host CNV/LoH inflates the donor LoD; robust refit recovers it (pure clone)")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150)
    print(f"Wrote {args.out}")


if __name__ == "__main__":
    main()
