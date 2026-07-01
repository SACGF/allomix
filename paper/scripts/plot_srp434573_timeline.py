"""Real-data monitoring trajectory figure for the allomix paper.

Reads output/srp434573_two_person.tsv (written by run_srp434573_allomix.py) and
plots one mixture's titration ladder as a declining-chimerism monitoring
trajectory on real reads. The minor contributor is assigned to the host role
(see run_srp434573_allomix.py), so the ladder from 10% down to 0.5% host reads
as residual host declining over serial timepoints, the direction that matters
for relapse surveillance. These are independent titration samples, not serial
timepoints from one patient, so the figure is labelled a titration series
presented as a monitoring trajectory.

Outputs:
  output/facts/fig_srp434573_timeline.png   the trajectory figure (main-text Figure 2)
  output/facts/srp434573_timeline.csv       caption facts

Usage:
    python paper/scripts/plot_srp434573_timeline.py
"""

import csv
from pathlib import Path
from statistics import median

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import paper_quick  # noqa: E402, F401  -- quick-build watermark (import for side effect)
from matplotlib.ticker import FixedLocator, FuncFormatter, NullLocator  # noqa: E402

OUT = Path("output")
FACTS_DIR = OUT / "facts"

# The representative dilution ladder. F2 into M1 is the same mixture used for the
# committed example reports (docs/examples/srp434573_dilution_series.html): a full
# six-rung ladder from 10% down to 0.5% host, both endpoints present.
MIXTURE = "mix_F2_into_M1"


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _fmt_pct(v: float, _pos: int) -> str:
    if v >= 1:
        return f"{v:g}%"
    return f"{v:g}%".rstrip("0").rstrip(".")


def _read(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def main() -> int:
    rows = [r for r in _read(OUT / "srp434573_two_person.tsv") if r["mixture"] == MIXTURE]
    # The titrated rungs carry a known host fraction; the pure-donor / pure-host
    # endpoints do not. Order high host -> low host, the surveillance direction.
    ladder = sorted(
        (r for r in rows if _f(r["known_pct"]) is not None),
        key=lambda r: -_f(r["known_pct"]),
    )
    known = [_f(r["known_pct"]) for r in ladder]
    mle = [_f(r["mle_pct"]) or 0.0 for r in ladder]
    ci_lo = [_f(r["mle_ci_lo"]) or 0.0 for r in ladder]
    ci_hi = [_f(r["mle_ci_hi"]) or 0.0 for r in ladder]
    pres_p = [_f(r.get("presence_p")) for r in ladder]
    detected = [p is not None and p < 0.05 for p in pres_p]
    x = list(range(len(ladder)))

    # Mixture contamination floor: the in-data consensus-homozygous level, median
    # across this mixture's admixtures (the level below which a corrected estimate
    # is not separable from contamination).
    cfracs = [_f(r.get("contamination_frac")) for r in rows]
    cfracs = [c * 100 for c in cfracs if c is not None]
    floor_pct = median(cfracs) if cfracs else None

    fig, ax = plt.subplots(figsize=(7, 4.5))

    ax.plot(
        x, known, "s--", color="0.4", linewidth=1.5, markersize=6,
        label="Known host fraction", zorder=3,
    )
    yerr = [
        [max(m - lo, 0.0) for m, lo in zip(mle, ci_lo)],
        [max(hi - m, 0.0) for m, hi in zip(mle, ci_hi)],
    ]
    ax.errorbar(
        x, mle, yerr=yerr, fmt="o-", color="steelblue", linewidth=2.2,
        markersize=7, markeredgecolor="white", markeredgewidth=0.8,
        capsize=3, ecolor="steelblue", elinewidth=1.2,
        label="allomix estimate (95% CI)", zorder=5,
    )

    if floor_pct is not None:
        ax.axhline(floor_pct, color="firebrick", linestyle=":", linewidth=1.2, zorder=2)
        ax.annotate(
            f"contamination floor ~{floor_pct:.2g}%",
            xy=(x[-1], floor_pct), xytext=(0, 4), textcoords="offset points",
            fontsize=8, color="firebrick", ha="right", va="bottom",
        )

    # Mark every timepoint that returned a positive residual-host presence call.
    if all(detected):
        ax.text(
            0.02, 0.06,
            "residual-host presence test positive at all timepoints",
            transform=ax.transAxes, fontsize=8, color="steelblue",
            ha="left", va="bottom",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="white", edgecolor="none", alpha=0.8),
        )

    ax.set_yscale("log")
    ax.yaxis.set_major_locator(FixedLocator([0.1, 0.2, 0.5, 1, 2, 5, 10]))
    ax.yaxis.set_major_formatter(FuncFormatter(_fmt_pct))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_xticks(x)
    ax.set_xticklabels([f"{k:g}%" for k in known])
    ax.set_xlabel("Serial monitoring timepoint (titrated host fraction)", fontsize=11)
    ax.set_ylabel("Host fraction (%)", fontsize=11)
    ax.set_title(
        "Monitoring trajectory on real titrated mixtures (SRP434573)",
        fontsize=12, fontweight="bold", loc="left",
    )
    ax.grid(True, which="major", alpha=0.2)
    ax.legend(fontsize=9, loc="upper right")
    fig.tight_layout()

    FACTS_DIR.mkdir(parents=True, exist_ok=True)
    fig_path = FACTS_DIR / "fig_srp434573_timeline.png"
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    # Caption facts. MAE over the ladder rungs, plus the span and floor.
    abs_err = [abs(m - k) for m, k in zip(mle, known)]
    facts = {
        "mixture_label": MIXTURE.replace("mix_", "").replace("_into_", " into "),
        "n_timepoints": len(ladder),
        "host_start_pct": f"{known[0]:g}",
        "host_end_pct": f"{known[-1]:g}",
        "donor_start_pct": f"{100 - known[0]:g}",
        "donor_end_pct": f"{100 - known[-1]:g}",
        "mae_pct": round(sum(abs_err) / len(abs_err), 2),
        "presence_n_detected": sum(detected),
        "floor_pct": round(floor_pct, 2) if floor_pct is not None else "",
    }
    facts_path = FACTS_DIR / "srp434573_timeline.csv"
    with open(facts_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(facts))
        w.writeheader()
        w.writerow(facts)

    print(f"Saved {fig_path}")
    print(f"Wrote {facts_path}: {facts}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
