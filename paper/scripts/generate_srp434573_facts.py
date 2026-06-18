"""Summarise the SRP434573 real-data mixture results for the paper.

Reads the two TSVs produced by ``paper/scripts/run_srp434573_allomix.py`` (which runs
allomix on the joint-called SRP434573 VCFs) and writes:

  output/facts/srp434573.csv     headline facts (template variables)
  output/facts/fig_srp434573.png two-panel figure: dilution-series accuracy
                                 (known vs estimated host %) and the one
                                 three-person mixture

This script does no sequencing or joint calling: it only summarises the
already-computed allomix output. If the input TSVs are missing, regenerate them
first with ``paper/scripts/run_srp434573_allomix.py`` (which reads the committed
genotype snapshot in ``paper/public_data/SRP434573/genotypes``, or a freshly
joint-called ``output/genotypes/SRP434573`` if present).

The minor (titrated) contributor is mapped to HOST, so the reported quantity is
the host fraction and the dilution series (10% down to 0.5%) reads as a
declining-chimerism / relapse-monitoring series. See
``paper/public_data/SRP434573/README.md`` and ``doc/joint_calling.md``.
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
TWO_TSV = OUT / "srp434573_two_person.tsv"
THREE_TSV = OUT / "srp434573_three_person.tsv"
PANEL_BED = Path("paper/public_data/SRP434573/SRP434573.bed")

# Fractions at and above which the ~0.2% co-pooled contamination floor is
# negligible relative to the true host fraction (Methods / Discussion).
RELIABLE_MIN_PCT = 2.5


def _read(path: Path) -> list[dict]:
    if not path.exists():
        sys.exit(
            f"Missing {path}. Regenerate it with paper/scripts/run_srp434573_allomix.py "
            "(reads paper/public_data/SRP434573/genotypes)."
        )
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _f(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _count_intervals(bed: Path) -> dict[str, int]:
    """Count reconstructed panel intervals in the BED, split by chromosome class."""
    if not bed.exists():
        sys.exit(f"Missing {bed}; needed for panel interval counts.")
    total = autosomal = chrx = 0
    with open(bed, encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            chrom = line.split("\t", 1)[0]
            total += 1
            if chrom == "chrX":
                chrx += 1
            elif chrom.removeprefix("chr").isdigit():
                autosomal += 1
    return {"total": total, "autosomal": autosomal, "chrx": chrx}


def _fmt_pct(v: float, _pos: int) -> str:
    if v >= 1:
        return f"{v:g}%"
    return f"{v:g}%".rstrip("0").rstrip(".")


def compute_facts(two: list[dict], three: list[dict]) -> dict:
    facts: dict[str, str] = {}

    # Titration accuracy is computed on the dilution series only. The pure
    # host/donor endpoints (0% and 100% host) carry no titration known fraction
    # and serve as detection anchors in the figS12 plot, so drop them here.
    two = [r for r in two if _f(r["known_pct"]) is not None]

    known = np.array([_f(r["known_pct"]) for r in two])
    mle = np.array([(_f(r["mle_pct"]) or 0.0) for r in two])
    depth = np.array([_f(r["mean_depth"]) for r in two])
    qc = [r["qc"] for r in two]

    facts["n_mixtures"] = str(len({r["mixture"] for r in two}))
    facts["n_timepoints"] = str(len(two))
    facts["n_individuals"] = "7"
    facts["depth_median"] = f"{np.median(depth):.0f}"
    facts["depth_min"] = f"{depth.min():.0f}"
    facts["depth_max"] = f"{depth.max():.0f}"
    facts["n_review"] = str(sum(1 for q in qc if q == "REVIEW"))
    facts["n_pass"] = str(sum(1 for q in qc if q == "PASS"))

    # Panel and platform metadata (paper/public_data/SRP434573/README.md). The
    # SNP count is the thesis-stated 1,062; the intervals are reconstructed from
    # the aligned reads (1,052 = 1,025 autosomal + 27 chrX, derived from the BED).
    facts["panel_n_snps"] = "1062"
    facts["platform"] = "Illumina HiSeq 3000"
    facts["raw_depth_min"] = "1000"
    facts["raw_depth_max"] = "1900"
    iv = _count_intervals(PANEL_BED)
    facts["n_intervals"] = str(iv["total"])
    facts["n_intervals_autosomal"] = str(iv["autosomal"])
    facts["n_intervals_chrx"] = str(iv["chrx"])

    # Dilution ladder (minor-contributor %), derived from the known fractions
    # so the listed levels stay in step with the data.
    ladder = sorted({float(k) for k in known}, reverse=True)
    facts["n_dilution_levels"] = str(len(ladder))
    facts["dilution_max_pct"] = f"{ladder[0]:g}"
    facts["dilution_min_pct"] = f"{ladder[-1]:g}"
    facts["dilution_ladder"] = ", ".join(f"{v:g}%" for v in ladder)

    # Informative markers actually used in the MLE after GT-quality and depth
    # filtering (median and range across the two-person timepoints), distinct
    # from the panel size (1,062 SNPs) and the reconstructed intervals (1,052).
    n_used = np.array([_f(r["n_used"]) for r in two if _f(r["n_used"]) is not None])
    facts["markers_used_median"] = f"{np.median(n_used):.0f}"
    facts["markers_used_min"] = f"{n_used.min():.0f}"
    facts["markers_used_max"] = f"{n_used.max():.0f}"

    # Accuracy on the reliable range (known >= 2.5%), where the co-pooled
    # contamination floor (~0.2%) is small relative to the true fraction.
    rel = known >= RELIABLE_MIN_PCT
    rel_err = mle[rel] - known[rel]
    facts["n_reliable"] = str(int(rel.sum()))
    facts["mae_reliable_pct"] = f"{np.mean(np.abs(rel_err)):.2f}"
    facts["max_abs_err_reliable_pct"] = f"{np.max(np.abs(rel_err)):.2f}"
    # Concordance on the reliable range (Pearson r and r^2 of estimated vs known).
    r = np.corrcoef(known[rel], mle[rel])[0, 1]
    facts["r2_reliable"] = f"{r * r:.3f}"
    facts["r_reliable"] = f"{r:.3f}"
    # Log-scale concordance across the whole series (spans 0.5%-10%).
    pos = mle > 0
    rlog = np.corrcoef(np.log10(known[pos]), np.log10(mle[pos]))[0, 1]
    facts["r2_log_all"] = f"{rlog * rlog:.3f}"

    # The lowest titration level (0.5% host): the contamination floor competes
    # with the true fraction, so the estimate scatters around it.
    lo = np.isclose(known, 0.5)
    facts["n_lowest"] = str(int(lo.sum()))
    facts["lowest_known_pct"] = "0.5"
    facts["mle_lowest_min_pct"] = f"{mle[lo].min():.2f}"
    facts["mle_lowest_max_pct"] = f"{mle[lo].max():.2f}"
    facts["mle_lowest_mean_pct"] = f"{np.mean(mle[lo]):.2f}"

    # 1% host level (six observations across mixtures), the clinically relevant
    # residual-disease threshold.
    one = np.isclose(known, 1.0)
    facts["n_onepct"] = str(int(one.sum()))
    facts["mle_onepct_min_pct"] = f"{mle[one].min():.2f}"
    facts["mle_onepct_max_pct"] = f"{mle[one].max():.2f}"

    # Residual-host presence test (separate from the magnitude MLE): donor-
    # homozygous markers where the host carries the donor-absent allele. Report
    # how often it fires, how many markers it reads, and its host-fraction
    # estimate at the clinically relevant 1% level.
    pres_p = np.array([_f(r["presence_p"]) for r in two])
    pres_pct = np.array([(_f(r["presence_pct"]) or 0.0) for r in two])
    pres_mk = np.array([(_f(r["presence_markers"]) or 0.0) for r in two])
    detected = pres_p < 0.05
    facts["presence_n_detected"] = str(int(detected.sum()))
    facts["presence_n_total"] = str(len(two))
    facts["presence_markers_median"] = f"{np.median(pres_mk):.0f}"
    facts["presence_onepct_min_pct"] = f"{pres_pct[one].min():.2f}"
    facts["presence_onepct_max_pct"] = f"{pres_pct[one].max():.2f}"

    # Three-person mixture (1:3:5 of F2:M1:M2; host = F2, donors = M1, M2).
    by_comp = {r["component"]: r for r in three}
    for comp in ("F2", "M1", "M2"):
        r3 = by_comp.get(comp)
        if r3:
            facts[f"three_{comp.lower()}_known_pct"] = f"{_f(r3['known_pct']):.1f}"
            facts[f"three_{comp.lower()}_est_pct"] = f"{_f(r3['est_pct']):.1f}"

    return facts


def make_figure(two: list[dict], three: list[dict], out_path: Path) -> None:
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.5, 5.4))

    # Panel A: dilution-series accuracy (log-log scatter, both estimators).
    # Endpoints (0%/100% host) have no titration known fraction; exclude them.
    two = [r for r in two if _f(r["known_pct"]) is not None]
    known = [_f(r["known_pct"]) for r in two]
    mle = [(_f(r["mle_pct"]) or 0.0) for r in two]
    pres = [(_f(r["presence_pct"]) or 0.0) for r in two]
    floor = 0.04  # log-axis stand-in for "not detected" (0)
    mle_p = [m if m > 0 else floor for m in mle]
    pres_p = [p if p > 0 else floor for p in pres]

    # Axis floor sits just below the lowest estimate (~0.26% at the 0.5% dilution)
    # so the low-end scatter is kept without leaving a large empty lower-left.
    lims = [0.2, 13]
    axA.plot(lims, lims, color="0.4", ls="--", lw=1.2, zorder=1, label="perfect recovery (y = x)")
    axA.scatter(
        known,
        mle_p,
        s=55,
        color="#1f77b4",
        edgecolor="white",
        linewidth=0.6,
        marker="o",
        zorder=3,
        label="MLE (100 − donor%)",
    )
    axA.scatter(
        known,
        pres_p,
        s=55,
        facecolor="none",
        edgecolor="#d62728",
        linewidth=1.4,
        marker="s",
        zorder=3,
        label="presence-test",
    )
    axA.set_xscale("log")
    axA.set_yscale("log")
    axA.set_xlim(*lims)
    axA.set_ylim(*lims)
    axA.set_aspect("equal")
    ticks = [0.5, 1, 2.5, 5, 10]
    axA.set_xticks(ticks)
    axA.set_yticks(ticks)
    axA.xaxis.set_major_formatter(FuncFormatter(_fmt_pct))
    axA.yaxis.set_major_formatter(FuncFormatter(_fmt_pct))
    axA.grid(True, which="both", alpha=0.2)
    axA.set_xlabel("Known host fraction", fontsize=11)
    axA.set_ylabel("allomix estimated host %", fontsize=11)
    axA.set_title("A  Two-person dilution series", fontsize=12, fontweight="bold", loc="left")
    axA.legend(fontsize=8.5, loc="upper left", framealpha=0.92)

    # Panel B: three-person mixture.
    order = [("F2", "host"), ("M1", "donor"), ("M2", "donor")]
    by_comp = {r["component"]: r for r in three}
    labels, kk, ee, elo, ehi = [], [], [], [], []
    for comp, role in order:
        r3 = by_comp.get(comp)
        if not r3:
            continue
        labels.append(f"{comp}\n({role})")
        kk.append(_f(r3["known_pct"]))
        e = _f(r3["est_pct"])
        ee.append(e)
        lo, hi = _f(r3["ci_lo"]), _f(r3["ci_hi"])
        elo.append((e - lo) if (e is not None and lo is not None) else 0.0)
        ehi.append((hi - e) if (e is not None and hi is not None) else 0.0)
    x = range(len(labels))
    w = 0.38
    axB.bar([i - w / 2 for i in x], kk, width=w, color="0.7", label="known")
    axB.bar(
        [i + w / 2 for i in x],
        ee,
        width=w,
        color="#2c7fb8",
        yerr=[elo, ehi],
        capsize=4,
        label="allomix",
    )
    for i, (k, e) in enumerate(zip(kk, ee)):
        axB.text(i - w / 2, k + 1.2, f"{k:.1f}%", ha="center", fontsize=9, color="0.3")
        if e is not None:
            axB.text(
                i + w / 2,
                e + 1.2,
                f"{e:.1f}%",
                ha="center",
                fontsize=9,
                color="#2c7fb8",
                fontweight="bold",
            )
    axB.set_xticks(list(x))
    axB.set_xticklabels(labels, fontsize=10)
    axB.set_ylabel("fraction of sample (%)", fontsize=11)
    axB.set_ylim(0, max(kk + [e for e in ee if e is not None]) * 1.18)
    axB.set_title(
        "B  Three-person mixture (1:3:5 F2:M1:M2)", fontsize=12, fontweight="bold", loc="left"
    )
    axB.legend(fontsize=10)
    axB.grid(True, axis="y", alpha=0.2)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_path}", file=sys.stderr)


def main() -> int:
    FACTS_DIR.mkdir(parents=True, exist_ok=True)
    two = _read(TWO_TSV)
    three = _read(THREE_TSV)

    facts = compute_facts(two, three)
    path = FACTS_DIR / "srp434573.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(facts.keys()))
        writer.writeheader()
        writer.writerow(facts)
    print(f"Wrote {path}", file=sys.stderr)

    make_figure(two, three, FACTS_DIR / "fig_srp434573.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
