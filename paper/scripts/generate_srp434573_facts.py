"""Summarise the SRP434573 real-data mixture results for the paper.

Reads the two TSVs produced by ``paper/scripts/run_srp434573_allomix.py`` and writes
``output/facts/srp434573.csv`` (template variables) and
``output/facts/fig_srp434573.png`` (two-panel figure: dilution-series accuracy and the
three-person mixture). It only summarises already-computed allomix output. If the input
TSVs are missing, regenerate them with ``run_srp434573_allomix.py`` (reads the committed
snapshot in ``paper/public_data/SRP434573/genotypes``, or a freshly joint-called
``output/genotypes/SRP434573`` if present).

The minor (titrated) contributor is mapped to HOST, so the reported quantity is the host
fraction and the dilution series (10% down to 0.5%) reads as a declining-chimerism /
relapse-monitoring series. See ``paper/public_data/SRP434573/README.md`` and
``docs/joint_calling.md``.
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
# Baseline (no Step 30) two-person results, for the before/after contamination facts.
TWO_BASELINE_TSV = OUT / "srp434573_two_person_baseline.tsv"
THREE_TSV = OUT / "srp434573_three_person.tsv"
CONTAM_TABLE_DIR = OUT / "contam_tables"
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

    # Titration accuracy uses the dilution series only. The pure host/donor endpoints
    # (0%/100% host) carry no known fraction (serve as detection anchors in figS12).
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

    # Panel and platform metadata (paper/public_data/SRP434573/README.md). SNP count is
    # the thesis-stated 1,062; intervals are reconstructed from aligned reads (1,052 =
    # 1,025 autosomal + 27 chrX, derived from the BED).
    facts["panel_n_snps"] = "1062"
    facts["platform"] = "Illumina HiSeq 3000"
    facts["raw_depth_min"] = "1000"
    facts["raw_depth_max"] = "1900"
    iv = _count_intervals(PANEL_BED)
    facts["n_intervals"] = str(iv["total"])
    facts["n_intervals_autosomal"] = str(iv["autosomal"])
    facts["n_intervals_chrx"] = str(iv["chrx"])

    # Dilution ladder (minor-contributor %), derived from the known fractions so the
    # listed levels stay in step with the data.
    ladder = sorted({float(k) for k in known}, reverse=True)
    facts["n_dilution_levels"] = str(len(ladder))
    facts["dilution_max_pct"] = f"{ladder[0]:g}"
    facts["dilution_min_pct"] = f"{ladder[-1]:g}"
    facts["dilution_ladder"] = ", ".join(f"{v:g}%" for v in ladder)

    # Informative markers used in the MLE after GT-quality and depth filtering, distinct
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

    # 1% host level, the clinically relevant residual-disease threshold.
    one = np.isclose(known, 1.0)
    facts["n_onepct"] = str(int(one.sum()))
    facts["mle_onepct_min_pct"] = f"{mle[one].min():.2f}"
    facts["mle_onepct_max_pct"] = f"{mle[one].max():.2f}"

    # Residual-host presence test (separate from the magnitude MLE): donor-homozygous
    # markers where the host carries the donor-absent allele.
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


def _endpoint_floors(rows: list[dict]) -> list[float]:
    """MLE host % at the true-0%-host endpoints (pure-donor samples).

    These rows have ``known_pct`` blank; the < 50% cut excludes the pure-host endpoint
    (~100%). These are the floating floors Step 30 targets.
    """
    out = []
    for r in rows:
        if _f(r["known_pct"]) is not None:
            continue
        m = _f(r.get("mle_pct"))
        if m is not None and m < 50.0:
            out.append(m)
    return out


def contamination_facts(two: list[dict], two_base: list[dict] | None) -> dict:
    """Contamination-line and Step 30 before/after facts.

    ``two`` is the headline (Step 30) run, ``two_base`` the baseline (no
    correction) run; ``two_base`` may be None (then before/after keys are
    omitted). The contamination level is Step-30-independent (consensus-hom
    markers), so it is read off the headline run.
    """
    facts: dict[str, str] = {}

    # In-data contamination level (median per-site minor fraction over the error floor),
    # read off the headline run; identical under the correction.
    contam = [(_f(r.get("contamination_frac")) or 0.0) * 100 for r in two]
    contam = [c for c in contam if c > 0]
    if contam:
        facts["contam_floor_median_pct"] = f"{np.median(contam):.2f}"
        facts["contam_floor_max_pct"] = f"{np.max(contam):.2f}"

    # Per-mixture contamination level (median across that mixture's timepoints):
    # the height of each line in Figure S12.
    by_mix: dict[str, list[float]] = {}
    for r in two:
        c = _f(r.get("contamination_frac"))
        if c is not None and c > 0:
            by_mix.setdefault(r["mixture"], []).append(c * 100)
    if by_mix:
        per_mix = {m: float(np.median(v)) for m, v in by_mix.items()}
        facts["contam_line_min_pct"] = f"{min(per_mix.values()):.2f}"
        facts["contam_line_max_pct"] = f"{max(per_mix.values()):.2f}"

    # Step 30 gate outcome and slope range, from the saved per-mixture tables.
    if CONTAM_TABLE_DIR.exists():
        from allomix.calibration.contamination_table import load_contamination_table

        gated, slopes = 0, []
        tables = sorted(CONTAM_TABLE_DIR.glob("*.contam.tsv"))
        for t in tables:
            corr = load_contamination_table(t)
            if corr.gated:
                gated += 1
                slopes.append(corr.slope * 100)
        facts["correction_n_mixtures"] = str(len(tables))
        facts["correction_n_gated"] = str(gated)
        if slopes:
            facts["correction_slope_min_pct"] = f"{min(slopes):.3f}"
            facts["correction_slope_max_pct"] = f"{max(slopes):.3f}"

    # Zero-host endpoint floor, before vs after Step 30: the floating MLE at true 0%
    # host pulled toward 0.
    s30_floor = _endpoint_floors(two)
    if s30_floor:
        facts["endpoint_floor_max_corrected_pct"] = f"{max(s30_floor):.3f}"
        facts["endpoint_floor_median_corrected_pct"] = f"{np.median(s30_floor):.3f}"
    if two_base is not None:
        base_floor = _endpoint_floors(two_base)
        if base_floor:
            facts["endpoint_floor_max_baseline_pct"] = f"{max(base_floor):.3f}"
            facts["endpoint_floor_median_baseline_pct"] = f"{np.median(base_floor):.3f}"
        # 0.5% and 1% before/after means, on the titration rungs only.
        base = [r for r in two_base if _f(r["known_pct"]) is not None]
        bk = np.array([_f(r["known_pct"]) for r in base])
        bm = np.array([(_f(r["mle_pct"]) or 0.0) for r in base])
        if bk.size:
            lo = np.isclose(bk, 0.5)
            if lo.any():
                facts["mle_lowest_mean_baseline_pct"] = f"{np.mean(bm[lo]):.2f}"
            rel = bk >= RELIABLE_MIN_PCT
            if rel.any():
                facts["mae_reliable_baseline_pct"] = f"{np.mean(np.abs(bm[rel] - bk[rel])):.2f}"

    return facts


def zero_host_facts(two: list[dict]) -> dict:
    """True-0%-host detection specificity from the pure-donor endpoints.

    Each two-person pair's admix VCF carries the donor's own reads as an endpoint
    sample, piled through the same forced ``bcftools mpileup`` path as the titrated
    mixtures. Feeding the donor as the admix is a genuine 0%-host case (real reads,
    the real co-pooled contamination floor, no synthetic mixing), so these rows
    measure host-detection specificity at true zero: the MLE host floor, and whether
    the residual-host presence test correctly stays negative. Donor endpoints are the
    blank-``known_pct`` rows whose MLE host is < 50% (the pure-host endpoints sit near
    100%). Presence call is "absent" (correct) when its p-value is >= 0.05.
    """
    facts: dict[str, str] = {}
    donor_ep = [
        r for r in two
        if _f(r["known_pct"]) is None and (_f(r.get("mle_pct")) or 0.0) < 50.0
    ]
    if not donor_ep:
        return facts
    mle = [(_f(r.get("mle_pct")) or 0.0) for r in donor_ep]
    pres_pct = [(_f(r.get("presence_pct")) or 0.0) for r in donor_ep]
    p_vals = [_f(r.get("presence_p")) for r in donor_ep]
    absent_n = sum(1 for p in p_vals if p is not None and p >= 0.05)
    called_n = sum(1 for p in p_vals if p is not None and p < 0.05)

    facts["zero_host_n"] = str(len(donor_ep))
    facts["zero_host_mle_max_pct"] = f"{max(mle):.3f}"
    facts["zero_host_mle_median_pct"] = f"{float(np.median(mle)):.3f}"
    facts["zero_host_presence_absent_n"] = str(absent_n)
    facts["zero_host_presence_falsepos_n"] = str(called_n)
    facts["zero_host_presence_max_pct"] = f"{max(pres_pct):.3f}"
    return facts


def make_figure(two: list[dict], three: list[dict], out_path: Path) -> None:
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(12.5, 5.4))

    # Panel A: dilution-series accuracy (log-log scatter, both estimators). Endpoints
    # (0%/100% host) have no known fraction; exclude them.
    two = [r for r in two if _f(r["known_pct"]) is not None]
    known = [_f(r["known_pct"]) for r in two]
    mle = [(_f(r["mle_pct"]) or 0.0) for r in two]
    pres = [(_f(r["presence_pct"]) or 0.0) for r in two]
    floor = 0.04  # log-axis stand-in for "not detected" (0)
    mle_p = [m if m > 0 else floor for m in mle]
    pres_p = [p if p > 0 else floor for p in pres]

    # Axis floor sits just below the lowest estimate (~0.26% at the 0.5% dilution) so the
    # low-end scatter is kept without leaving a large empty lower-left.
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
    # The semi-synthetic sub-0.5% points (issue #5) are deliberately kept off this
    # real-data panel: they live on their own zoomed figure
    # (fig_srp434573_synthetic.png) so real and synthetic are never conflated.
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
    axA.set_xlabel("Known recipient fraction", fontsize=11)
    axA.set_ylabel("allomix estimated recipient %", fontsize=11)
    axA.set_title("A  Two-person dilution series", fontsize=12, fontweight="bold", loc="left")
    axA.legend(fontsize=8.5, loc="upper left", framealpha=0.92)

    # Panel B: three-person mixture.
    order = [("F2", "recipient"), ("M1", "donor"), ("M2", "donor")]
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

    two_base = _read(TWO_BASELINE_TSV) if TWO_BASELINE_TSV.exists() else None

    facts = compute_facts(two, three)
    facts.update(contamination_facts(two, two_base))
    facts.update(zero_host_facts(two))
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
