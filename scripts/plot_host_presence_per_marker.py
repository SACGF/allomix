#!/usr/bin/env python3
"""Per-marker host-presence structure for selected admixture samples.

Internal SA Path diagnostic (not part of the allomix package).

The cohort presence figure (plot_host_presence.py) gives one host fraction per
sample. This one opens the box: at each donor-homozygous marker (where every
donor is the same homozygote and the host carries the donor-absent allele),
it shows the host fraction that marker alone implies, so you can see whether
the host signal is spread evenly across markers or carried by a few.

For each marker the detector counts ``y`` donor-absent reads out of depth ``n``,
and the host carries that allele at dose ``h`` (1 if host het, 2 if host hom).
Under a single host fraction f the expected donor-absent VAF is f*(h/2), so the
marker's implied host fraction is (y/n)/(h/2). Plotting that per marker, sorted,
with each marker's binomial CI and the pooled MLE line, answers three things at
once: how many markers carry any host signal, whether the implied fractions
agree (even) or scatter beyond their CIs (overdispersed / a few drivers), and
how the per-marker picture relates to the single pooled estimate.

Points come from the local VCFs (no coordinates are surfaced). The pooled MLE
host fraction and presence p-value are read from the run's batch.tsv so the
reference line matches the reported result exactly.

Convention: the y axis is host fraction %, because this figure is about where
the (small) host signal sits; the cohort figure stays in donor %.

Markers the detector drops as alignment artifacts (strand / soft-clip /
read-position bias, e.g. the TP53 intron-3 indel site) are drawn as black X
marks and left out of the pooled MLE line and the marker counts, so the figure
shows both where the artifact markers sat and the cleaned-up pooled estimate.

Usage:
    python scripts/plot_host_presence_per_marker.py \
        --vcf-dir output/joint_called \
        --samples-csv-dir pipeline/sample_csvs \
        --batch output/validation_run10/batch.tsv \
        --run-label run10 \
        --panel-suffix .union_sid_haem_vendor_probes.vcf.gz \
        --admix 29_MO_HP_FULL_NDAD_2600807502K 18_MO_HP_FULL_BHOA_2532810149I \
                37_MO_HP_FULL_PCAH_2534210988I 2_MO_HP_FULL_QUDO_2526508341H \
        --output output/host_presence_per_marker_run10.png
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from allomix.analysis import analyse_sample
from allomix.genotype import parse_vcf

# Wilson-interval z for a 95% CI.
Z = 1.959963984540054

# Host dose -> colour. Dose 1 = host het (contributes the allele at half rate),
# dose 2 = host hom (full rate).
DOSE_COLOR = {1: "#7fcdbb", 2: "#225ea8"}


def short_label(name: str) -> str:
    """Patient code plus leading id, matching plot_host_presence.py."""
    leading = name.split("_", 1)[0]
    codes = [t for t in name.split("_") if t.isalpha() and t.isupper() and len(t) >= 3]
    code = codes[-1] if codes else name
    return f"{code} #{leading}" if leading.isdigit() else code


def _wilson(y: int, n: int) -> tuple[float, float, float]:
    """Wilson point and 95% interval for a binomial proportion.

    Returns (p_hat, lo, hi). Handles y == 0 (lower bound 0) without blowing up.
    """
    if n <= 0:
        return 0.0, 0.0, 0.0
    p = y / n
    denom = 1.0 + Z * Z / n
    center = (p + Z * Z / (2 * n)) / denom
    half = Z * math.sqrt(p * (1 - p) / n + Z * Z / (4 * n * n)) / denom
    return p, max(0.0, center - half), min(1.0, center + half)


def index_csvs(csv_dir: str) -> dict[str, dict]:
    """Map each ADMIX sample id to its patient stem, host and donor ids."""
    out: dict[str, dict] = {}
    for path in sorted(glob.glob(os.path.join(csv_dir, "*.csv"))):
        stem = Path(path).stem
        host = donor = None
        admixes: list[str] = []
        with open(path) as f:
            for r in csv.DictReader(f):
                stype = r["sample_type"].strip().upper()
                sid = r["sample_id"].strip()
                if stype == "HOST":
                    host = sid
                elif stype == "DONOR":
                    donor = sid
                elif stype == "ADMIX":
                    admixes.append(sid)
        for a in admixes:
            out[a] = {"stem": stem, "host": host, "donor": donor}
    return out


def read_batch_presence(path: Path) -> dict[str, dict]:
    """Map sample -> {f_pct, p, capped, n} from a batch.tsv with host columns."""
    out: dict[str, dict] = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if not r.get("host_f_est"):
                continue
            p_raw = float(r["host_present_p"])
            out[r["sample"]] = {
                "f_pct": float(r["host_f_est"]) * 100.0,
                "p": p_raw,
                "capped": p_raw <= 1e-300,
                "n": int(r["host_detect_markers"]),
            }
    return out


def per_marker_rows(
    vcf_dir: str,
    panel_suffix: str,
    meta: dict,
    admix_sample: str,
    min_gq: int,
    min_dp: int,
    error_rate: float,
    use_sex_chroms: bool,
) -> tuple[list[dict], dict]:
    """Compute per-marker implied host fraction and the pooled MLE for a sample.

    The per-marker implied fraction subtracts the same per-direction background
    the detector uses (``error_rate / 3`` under the global fallback), so a
    background-only marker lands near 0 and the points sit on the same footing
    as the pooled MLE line. Returns (per-marker dicts, pooled dict).
    """
    panel = os.path.join(vcf_dir, meta["stem"] + panel_suffix)
    admix_vcf = os.path.join(vcf_dir, meta["stem"] + ".admix.vcf.gz")
    host = parse_vcf(panel, sample=meta["host"], min_gq=min_gq, gt_ad_consistency=True)
    donor = parse_vcf(panel, sample=meta["donor"], min_gq=min_gq, gt_ad_consistency=True)
    admix = parse_vcf(admix_vcf, sample=admix_sample, min_dp=0)
    # One shared analysis path (allomix.analysis): same classify -> presence ->
    # donor-hom selection the `monitor` CLI runs, so the points and the pooled
    # line here match the reported batch exactly.
    analysis = analyse_sample(
        host,
        [donor],
        admix,
        min_dp=min_dp,
        min_gq=min_gq,
        error_rate=error_rate,
        use_sex_chroms=use_sex_chroms,
    )
    hp = analysis.result.host_presence
    e = error_rate / 3.0  # per-direction background under the global fallback

    def implied(p: float, coef: float) -> float:
        return 100.0 * max(0.0, p - e) / coef

    out: list[dict] = []
    # donor_hom_markers carries the same artifact flag host_presence_test uses,
    # so the artifact points are drawn but excluded from the pooled line/counts.
    for m in analysis.donor_hom_markers:
        p_hat, lo, hi = _wilson(m.y, m.n)
        out.append(
            {
                "f_pct": implied(p_hat, m.coef),
                "lo_pct": implied(lo, m.coef),
                "hi_pct": implied(hi, m.coef),
                "n": m.n,
                "h": m.h,
                "y": m.y,
                "artifact": m.artifact,
            }
        )
    pooled = {
        "f_pct": hp.f_host_mle * 100.0,
        "p": hp.lrt_pval,
        "capped": hp.lrt_pval <= 1e-300,
        "n": hp.n_markers,
    }
    return out, pooled


def _draw_panel(ax, label: str, markers: list[dict], pooled: dict | None) -> None:
    """Draw one sample's per-marker caterpillar with the pooled MLE line.

    Markers whose 95% CI excludes the pooled MLE are ringed: these are the
    ones that disagree with a single shared host fraction. If the host signal
    were evenly distributed and binomial, only about 5% would land outside by
    chance; many more means the spread is overdispersed or driven by a few
    markers, which is the same thing the chimerism goodness-of-fit flags.
    """
    markers = sorted(markers, key=lambda d: d["f_pct"])
    x = list(range(len(markers)))
    kept = [d for d in markers if not d["artifact"]]
    n_art = sum(1 for d in markers if d["artifact"])
    # "Has host presence" = signal above background: the background-subtracted
    # 95% CI lower bound clears 0. Counting y > 0 only measures coverage,
    # because at this depth nearly every marker catches a background read.
    # Counts use the kept (non-artifact) markers, matching the pooled line.
    n_lit = sum(1 for d in kept if d["lo_pct"] > 0)
    n_tot = len(kept)
    f_mle = pooled["f_pct"] if pooled is not None else None

    n_inconsistent = 0
    for xi, d in zip(x, markers):
        if d["artifact"]:
            # Filtered as an alignment artifact: black X, no CI bar, and left
            # out of the pooled line and all counts. This is the same set the
            # Manhattan plot rings with an X.
            ax.plot(xi, d["f_pct"], "x", color="#000000", ms=7.0, mew=1.7, zorder=5)
            continue
        c = DOSE_COLOR.get(d["h"], "#888888")
        excludes = f_mle is not None and (d["hi_pct"] < f_mle or d["lo_pct"] > f_mle)
        n_inconsistent += excludes
        ax.plot([xi, xi], [d["lo_pct"], d["hi_pct"]], color=c, lw=0.8, alpha=0.55, zorder=2)
        ax.plot(
            xi,
            d["f_pct"],
            "o",
            color=c,
            ms=4.0,
            zorder=3,
            markeredgecolor="#c0392b" if excludes else "none",
            markeredgewidth=1.0 if excludes else 0.0,
        )

    if pooled is not None:
        ax.axhline(pooled["f_pct"], color="#c0392b", lw=1.3, zorder=4)
        p = pooled["p"]
        p_str = "< 1e-300" if pooled["capped"] else (f"{p:.2f}" if p >= 0.01 else f"{p:.1e}")
        frac_inc = n_inconsistent / n_tot if n_tot else 0.0
        ax.text(
            0.97,
            0.96,
            f"pooled MLE = {pooled['f_pct']:.3f}%\n"
            f"p = {p_str}\n"
            f"{n_inconsistent}/{n_tot} markers off the line ({frac_inc:.0%})",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=8.5,
            color="#c0392b",
        )

    ax.axhline(0, color="#cccccc", lw=0.8, zorder=1)
    art_note = f", {n_art} filtered" if n_art else ""
    ax.set_title(f"{label}   ({n_lit}/{n_tot} markers above background{art_note})", fontsize=11)
    ax.set_xlabel("donor-homozygous marker (sorted by implied host fraction)")
    ax.set_ylabel("implied host fraction (%)")
    ax.grid(axis="y", color="#eeeeee", lw=0.6, zorder=0)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf-dir", default="output/joint_called")
    ap.add_argument("--samples-csv-dir", default="pipeline/sample_csvs")
    ap.add_argument("--batch", type=Path, default=Path("output/validation_run10/batch.tsv"))
    ap.add_argument(
        "--run-label",
        default="run10",
        help="Label shown in the figure title (default run10)",
    )
    ap.add_argument("--panel-suffix", default=".union_sid_haem_vendor_probes.vcf.gz")
    ap.add_argument("--admix", nargs="+", required=True, help="ADMIX sample ids to plot")
    ap.add_argument("--min-gq", type=int, default=20)
    ap.add_argument("--min-dp", type=int, default=20)
    ap.add_argument(
        "--use-sex-chroms",
        action="store_true",
        help="Keep sex/MT markers. Default off, matching the monitor run "
        "default; leave off for sex-mismatched pairs.",
    )
    ap.add_argument(
        "--error-rate",
        type=float,
        default=0.01,
        help="Global sequencing error rate, matching the run (default 0.01)",
    )
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()

    csv_index = index_csvs(args.samples_csv_dir)
    batch = read_batch_presence(args.batch)

    panels: list[tuple[str, list[dict], dict | None]] = []
    for admix_sample in args.admix:
        meta = csv_index.get(admix_sample)
        if meta is None:
            raise SystemExit(f"ADMIX sample {admix_sample!r} not found in {args.samples_csv_dir}")
        markers, pooled = per_marker_rows(
            args.vcf_dir,
            args.panel_suffix,
            meta,
            admix_sample,
            args.min_gq,
            args.min_dp,
            args.error_rate,
            args.use_sex_chroms,
        )
        # Cross-check the recomputed MLE against the reported run value.
        b = batch.get(admix_sample)
        if b is not None:
            print(
                f"{short_label(admix_sample):12} recomputed f={pooled['f_pct']:.3f}% "
                f"batch f={b['f_pct']:.3f}%  (n={pooled['n']} vs {b['n']})"
            )
        panels.append((short_label(admix_sample), markers, pooled))

    n = len(panels)
    ncol = 2 if n > 1 else 1
    nrow = math.ceil(n / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.8 * ncol, 5.2 * nrow), squeeze=False)
    flat = [ax for row in axes for ax in row]
    for ax, (label, markers, pooled) in zip(flat, panels):
        _draw_panel(ax, label, markers, pooled)
    for ax in flat[n:]:
        ax.axis("off")

    handles = [
        plt.Line2D([], [], marker="o", ls="", color=DOSE_COLOR[2], label="host hom (dose 2)"),
        plt.Line2D([], [], marker="o", ls="", color=DOSE_COLOR[1], label="host het (dose 1)"),
        plt.Line2D([], [], color="#c0392b", lw=1.3, label="pooled MLE host fraction"),
        plt.Line2D(
            [], [], marker="x", ls="", color="#000000", mew=1.7, label="filtered (artifact)"
        ),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=9)
    fig.suptitle(f"Per-marker host-presence structure — {args.run_label}", fontsize=12)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    fig.savefig(args.output, dpi=150)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
