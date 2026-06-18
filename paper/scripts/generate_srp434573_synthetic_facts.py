"""Summarise the semi-synthetic sub-0.5% SRP434573 mixtures for the paper (issue #5).

The public SRP434573 titration bottoms out at 0.5% host. ``scripts/mix_bams.sh``
blends each pair's two pure reference BAMs with ``samtools view --subsample`` to
make lower points (host fractions 0.1-0.5%), joint-called the same way. Those are
*semi-synthetic*: real reads, real noise, real GATK/bcftools path, artificial
mixing ratio. ``paper/scripts/run_srp434573_allomix.py`` runs allomix on them and
writes ``output/srp434573_synthetic.tsv``.

This script reads that TSV and writes:

  output/facts/srp434573_synthetic.csv     headline facts (template variables)
  output/facts/fig_srp434573_synthetic.png MLE + presence host% vs known fraction,
                                           points styled as semi-synthetic, with
                                           the real public 0.5% points overlaid.

If the synthetic TSV is missing or empty (a fresh checkout before the TAU-side
generation step has run), it writes an ``n_points=0`` stub CSV and a placeholder
figure so the Snakemake rule is still satisfied and the build stays green.

These points are always labelled "semi-synthetic" so a synthetic fraction is
never presented as a measured one.
"""

import csv
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import paper_quick  # noqa: E402, F401  -- quick-build watermark (import for side effect)
from matplotlib.ticker import FuncFormatter  # noqa: E402

OUT = Path("output")
FACTS_DIR = OUT / "facts"
SYN_TSV = OUT / "srp434573_synthetic.tsv"
TWO_TSV = OUT / "srp434573_two_person.tsv"

# A presence test counts as a detection when the host-present p-value clears 0.05.
DETECT_ALPHA = 0.05
# Log-axis stand-in for a "not detected" (0) estimate.
FLOOR = 0.02


def _read(path: Path) -> list[dict]:
    if not path.exists():
        return []
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _fmt_pct(v: float, _pos: int) -> str:
    if v >= 1:
        return f"{v:g}%"
    return f"{v:g}%".rstrip("0").rstrip(".")


def _median(xs: list[float]) -> float | None:
    return float(np.median(xs)) if xs else None


def real_05_medians(two: list[dict]) -> tuple[float | None, float | None]:
    """Median real-data MLE and presence host % at the 0.5% dilution."""
    mle = [_f(r["mle_pct"]) for r in two if _f(r.get("known_pct")) == 0.5]
    pres = [_f(r["presence_pct"]) for r in two if _f(r.get("known_pct")) == 0.5]
    return _median([m for m in mle if m is not None]), _median(
        [p for p in pres if p is not None]
    )


def compute_facts(syn: list[dict], two: list[dict]) -> dict:
    facts: dict[str, str] = {}
    fracs = sorted({_f(r["frac_pct"]) for r in syn if _f(r.get("frac_pct")) is not None})
    pairs = {r["mixture"] for r in syn}
    seeds = {r.get("seed") for r in syn}

    facts["n_points"] = str(len(syn))
    facts["n_pairs"] = str(len(pairs))
    facts["n_fractions"] = str(len(fracs))
    facts["n_seeds"] = str(len([s for s in seeds if s]))
    facts["frac_min_pct"] = f"{min(fracs):g}" if fracs else "0"
    facts["frac_max_pct"] = f"{max(fracs):g}" if fracs else "0"
    facts["frac_ladder"] = ", ".join(f"{x:g}%" for x in fracs)

    # Per-fraction median MLE / presence host % and presence detection rate.
    for fr in fracs:
        rows = [r for r in syn if _f(r["frac_pct"]) == fr]
        mle = [_f(r["mle_pct"]) for r in rows if _f(r["mle_pct"]) is not None]
        pres = [_f(r["presence_pct"]) for r in rows if _f(r["presence_pct"]) is not None]
        pvals = [_f(r["presence_p"]) for r in rows if _f(r["presence_p"]) is not None]
        detected = sum(1 for p in pvals if p is not None and p < DETECT_ALPHA)
        tag = f"{fr:g}".replace(".", "p")  # e.g. 0.1 -> 0p1 (CSV-safe key)
        mm = _median(mle)
        pm = _median(pres)
        facts[f"mle_med_{tag}"] = f"{mm:.3f}" if mm is not None else ""
        facts[f"presence_med_{tag}"] = f"{pm:.3f}" if pm is not None else ""
        facts[f"detect_rate_{tag}"] = (
            f"{detected / len(pvals):.2f}" if pvals else ""
        )
        facts[f"n_{tag}"] = str(len(rows))

    # Synthetic-vs-real cross-check at 0.5% (the anchored pairs have a real 0.5%).
    real_mle, real_pres = real_05_medians(two)
    syn_05 = [_f(r["mle_pct"]) for r in syn if _f(r["frac_pct"]) == 0.5]
    syn_05 = [m for m in syn_05 if m is not None]
    syn_05_med = _median(syn_05)
    facts["real_05_mle_med_pct"] = f"{real_mle:.3f}" if real_mle is not None else ""
    facts["syn_05_mle_med_pct"] = f"{syn_05_med:.3f}" if syn_05_med is not None else ""
    if real_mle is not None and syn_05_med is not None:
        facts["syn_minus_real_05_mle_pct"] = f"{syn_05_med - real_mle:.3f}"
    else:
        facts["syn_minus_real_05_mle_pct"] = ""
    return facts


def make_figure(syn: list[dict], two: list[dict], out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.6, 5.6))

    if not syn:
        ax.text(0.5, 0.5, "No semi-synthetic mixtures generated yet.\n"
                "Run scripts/make_semisynthetic_srp434573.py (TAU-side),\n"
                "then rebuild.", ha="center", va="center", fontsize=11,
                transform=ax.transAxes, color="0.4")
        ax.set_axis_off()
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Wrote {out_path} (placeholder, no synthetic data)", file=sys.stderr)
        return

    known = np.array([_f(r["frac_pct"]) for r in syn], dtype=float)
    mle = np.array([(_f(r["mle_pct"]) or 0.0) for r in syn], dtype=float)
    pres = np.array([(_f(r["presence_pct"]) or 0.0) for r in syn], dtype=float)
    # Small horizontal jitter so overlapping replicates at each fraction separate.
    rng = np.random.default_rng(0)
    jit = known * np.exp(rng.uniform(-0.05, 0.05, size=known.shape))

    lo = min(0.08, float(known.min()) * 0.8)
    lims = [lo, 1.2]
    ax.plot(lims, lims, color="0.4", ls="--", lw=1.1, zorder=1,
            label="perfect recovery (y = x)")
    ax.scatter(jit, np.where(mle > 0, mle, FLOOR), s=42, facecolor="none",
               edgecolor="#1f77b4", linewidth=1.2, marker="o", zorder=3,
               label="semi-synthetic MLE (100 − donor%)")
    ax.scatter(jit, np.where(pres > 0, pres, FLOOR), s=42, facecolor="none",
               edgecolor="#d62728", linewidth=1.2, marker="s", zorder=3,
               label="semi-synthetic presence-test")

    # Real public 0.5% points overlaid (filled) as the synthetic-vs-real anchor.
    real_mle, real_pres = real_05_medians(two)
    if real_mle is not None:
        ax.scatter([0.5], [real_mle], s=85, color="#1f77b4", edgecolor="black",
                   linewidth=0.8, marker="o", zorder=4, label="real 0.5% MLE (median)")
    if real_pres is not None and real_pres > 0:
        ax.scatter([0.5], [real_pres], s=85, color="#d62728", edgecolor="black",
                   linewidth=0.8, marker="s", zorder=4,
                   label="real 0.5% presence (median)")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*lims)
    ax.set_ylim(FLOOR * 0.7, 1.2)
    ticks = [0.1, 0.2, 0.3, 0.5, 1.0]
    ax.set_xticks(ticks)
    ax.set_yticks([FLOOR, 0.1, 0.2, 0.5, 1.0])
    ax.xaxis.set_major_formatter(FuncFormatter(_fmt_pct))
    ax.yaxis.set_major_formatter(
        FuncFormatter(lambda v, p: "0 (n.d.)" if abs(v - FLOOR) < 1e-9 else _fmt_pct(v, p))
    )
    ax.grid(True, which="both", alpha=0.2)
    ax.set_xlabel("Known host fraction (semi-synthetic mix)", fontsize=11)
    ax.set_ylabel("allomix estimated host %", fontsize=11)
    ax.set_title("Semi-synthetic sub-0.5% mixtures\n(subsampled real BAMs)",
                 fontsize=12, fontweight="bold", loc="left")
    ax.legend(fontsize=8, loc="upper left", framealpha=0.92)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path} ({len(syn)} synthetic points)", file=sys.stderr)


def main() -> int:
    FACTS_DIR.mkdir(parents=True, exist_ok=True)
    syn = _read(SYN_TSV)
    two = _read(TWO_TSV)

    facts = compute_facts(syn, two)
    path = FACTS_DIR / "srp434573_synthetic.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(facts.keys()))
        writer.writeheader()
        writer.writerow(facts)
    print(f"Wrote {path} (n_points={facts['n_points']})", file=sys.stderr)

    make_figure(syn, two, FACTS_DIR / "fig_srp434573_synthetic.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
