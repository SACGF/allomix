#!/usr/bin/env python3
"""Manhattan plot of per-marker host-presence signal along the genome.

Diagnostic script (run from `scripts/`, not part of the installed allomix package).

Unlike the per-marker caterpillar (plot_host_presence_per_marker.py), which sorts
by magnitude, this keeps genomic position on the x axis, so a host CNV/LOH region
shows up as a contiguous run of markers lifted above the sample-wide pooled host
fraction while a lone artifact stays isolated.

y is the dose-normalised, background-subtracted implied host fraction at each
donor-homozygous marker, with the pooled MLE drawn as a horizontal line. Points
are coloured by chromosome; UPREG markers (significantly above the pooled
fraction, Bonferroni one-sided binomial) are ringed.

The x axis shows genomic coordinates, so the figure is written to a LOCAL file
only (see CLAUDE.md) and not surfaced to stdout, even though the marker sites are
public panel SNPs.

Usage:
    python scripts/host_presence_manhattan.py \
        --admix SAMPLE_AAAA SAMPLE_BBBB SAMPLE_CCCC \
        --out output/host_presence_manhattan.png
"""

import argparse
import csv
import glob
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import binom

from allomix.analysis import analyse_sample
from allomix.genotype import parse_vcf

ALPHA = 0.05

# hg38 chromosome lengths (bp), for a proportional genome x axis.
HG38 = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr3": 198295559,
    "chr4": 190214555,
    "chr5": 181538259,
    "chr6": 170805979,
    "chr7": 159345973,
    "chr8": 145138636,
    "chr9": 138394717,
    "chr10": 133797422,
    "chr11": 135086622,
    "chr12": 133275309,
    "chr13": 114364328,
    "chr14": 107043718,
    "chr15": 101991189,
    "chr16": 90338345,
    "chr17": 83257441,
    "chr18": 80373285,
    "chr19": 58617616,
    "chr20": 64444167,
    "chr21": 46709983,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
}
CHROM_ORDER = list(HG38.keys())
# Cumulative start offset of each chromosome along the concatenated axis.
_OFFSET: dict[str, int] = {}
_acc = 0
for _c in CHROM_ORDER:
    _OFFSET[_c] = _acc
    _acc += HG38[_c]
GENOME_LEN = _acc

ALT_COLORS = ["#3b528b", "#5ec962"]  # two-tone chromosome banding


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


def short_label(name: str) -> str:
    leading = name.split("_", 1)[0]
    codes = [t for t in name.split("_") if t.isalpha() and t.isupper() and len(t) >= 3]
    code = codes[-1] if codes else name
    return f"{code} #{leading}" if leading.isdigit() else code


def load_genes(path: Path) -> dict[str, list[tuple[int, int, str]]]:
    """Load a protein-coding gene BED into chrom -> [(start, end, name), ...]."""
    genes: dict[str, list[tuple[int, int, str]]] = {}
    if not path.exists():
        return genes
    with open(path) as f:
        for line in f:
            chrom, start, end, name = line.rstrip("\n").split("\t")
            genes.setdefault(chrom, []).append((int(start), int(end), name))
    for v in genes.values():
        v.sort()
    return genes


def nearest_gene(genes: dict, chrom: str, pos: int) -> tuple[str, int]:
    """Return (gene_name, distance_bp); distance 0 if the position is inside a gene."""
    best_name, best_dist = "", None
    for start, end, name in genes.get(chrom, []):
        dist = 0 if start <= pos <= end else min(abs(pos - start), abs(pos - end))
        if best_dist is None or dist < best_dist:
            best_name, best_dist = name, dist
    return best_name, (best_dist if best_dist is not None else -1)


def build_markers(
    meta: dict,
    admix_sample: str,
    vcf_dir: str,
    panel_suffix: str,
    min_gq: int,
    min_dp: int,
    error_rate: float,
    genes: dict,
) -> tuple[list[dict], dict]:
    """Per-marker implied host fraction with genome coordinate and UPREG flag.

    UPREG markers are annotated with their nearest protein-coding gene for the
    manual driver-gene check.
    """
    panel = os.path.join(vcf_dir, meta["stem"] + panel_suffix)
    admix_vcf = os.path.join(vcf_dir, meta["stem"] + ".admix.vcf.gz")
    host = parse_vcf(panel, sample=meta["host"], min_gq=min_gq, gt_ad_consistency=True)
    donor = parse_vcf(panel, sample=meta["donor"], min_gq=min_gq, gt_ad_consistency=True)
    admix = parse_vcf(admix_vcf, sample=admix_sample, min_dp=0)
    # Shared analysis path (allomix.analysis); sex chroms kept here so the
    # genomic view shows chrX/Y markers for investigation.
    analysis = analyse_sample(
        host,
        [donor],
        admix,
        min_dp=min_dp,
        min_gq=min_gq,
        error_rate=error_rate,
        use_sex_chroms=True,
    )
    hp = analysis.result.host_presence
    f = hp.f_host_mle
    e = error_rate / 3.0
    markers = analysis.donor_hom_markers
    bonf = ALPHA / len(markers) if markers else ALPHA

    out: list[dict] = []
    for m in markers:
        if m.chrom not in _OFFSET:
            continue
        vaf = m.y / m.n if m.n else 0.0
        p0 = min(max(e + f * m.coef, 1e-9), 0.5)
        p_up = float(binom.sf(m.y - 1, m.n, p0)) if m.y > 0 else 1.0
        # Artifact markers (strand/soft-clip/read-position bias) the detector
        # drops: marked, not counted as upregulated discoveries, and excluded
        # from the per-chromosome means.
        upreg = (p_up < bonf) and not m.artifact
        gene = ""
        if upreg or m.artifact:
            gname, gdist = nearest_gene(genes, m.chrom, m.pos)
            gene = gname if gdist == 0 else (f"{gname}~{gdist // 1000}kb" if gname else "")
        out.append(
            {
                "chrom": m.chrom,
                "x": _OFFSET[m.chrom] + m.pos,
                "implied_pct": 100.0 * max(0.0, vaf - e) / m.coef,
                "upreg": upreg,
                "artifact": m.artifact,
                "gene": gene,
            }
        )
    pooled = {"f_pct": f * 100.0, "p": hp.lrt_pval, "n": len(markers)}
    return out, pooled


def _draw_panel(ax, label: str, markers: list[dict], pooled: dict) -> None:
    """Draw one sample's panel: markers, pooled line, per-chromosome means."""
    by_chrom: dict[str, list[dict]] = {}
    for m in markers:
        if m.get("artifact"):
            # Artifact-filtered marker: grey x, not part of the per-chrom mean.
            ax.plot(m["x"], m["implied_pct"], "x", ms=7, color="#999999", mew=1.6, zorder=4)
            if m["gene"]:
                ax.annotate(
                    f"{m['gene']} (filtered)",
                    (m["x"], m["implied_pct"]),
                    textcoords="offset points",
                    xytext=(4, 4),
                    fontsize=6.5,
                    color="#777777",
                    fontstyle="italic",
                    zorder=5,
                )
            continue
        by_chrom.setdefault(m["chrom"], []).append(m)
        ci = CHROM_ORDER.index(m["chrom"]) % 2
        ax.plot(m["x"], m["implied_pct"], "o", ms=5, color=ALT_COLORS[ci], zorder=2)

    # Per-chromosome mean line: horizontal segment over each chromosome's span
    # at its markers' mean implied host fraction. A whole chromosome (or arm)
    # above the pooled line is the CNV/LOH signature, distinct from a lone spike.
    for chrom, ms in by_chrom.items():
        mean = sum(d["implied_pct"] for d in ms) / len(ms)
        x0 = _OFFSET[chrom]
        ax.plot([x0, x0 + HG38[chrom]], [mean, mean], color="black", lw=1.1, alpha=0.7, zorder=3)

    ax.axhline(pooled["f_pct"], color="#c0392b", lw=1.1, zorder=1)

    for m in markers:
        if m["upreg"]:
            ax.plot(
                m["x"], m["implied_pct"], "o", ms=9, mfc="none", mec="#c0392b", mew=1.5, zorder=4
            )
            if m["gene"]:
                ax.annotate(
                    m["gene"],
                    (m["x"], m["implied_pct"]),
                    textcoords="offset points",
                    xytext=(4, 4),
                    fontsize=7,
                    color="#c0392b",
                    fontstyle="italic",
                    zorder=5,
                )

    n_up = sum(1 for m in markers if m["upreg"])
    n_art = sum(1 for m in markers if m.get("artifact"))
    ax.set_title(
        f"{label}   pooled host {pooled['f_pct']:.3f}%   "
        f"{pooled['n']} markers, {n_up} upregulated (ringed), "
        f"{n_art} artifact-filtered (grey x)",
        fontsize=10,
    )
    ax.set_ylabel("implied host %")
    ax.set_xlim(0, GENOME_LEN)
    ax.grid(axis="y", color="#eeeeee", lw=0.6, zorder=0)

    # Chromosome labels on every panel (not just the bottom).
    ticks = [_OFFSET[c] + HG38[c] / 2 for c in CHROM_ORDER]
    ax.set_xticks(ticks)
    ax.set_xticklabels([c[3:] for c in CHROM_ORDER], fontsize=6.5)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf-dir", default="output/joint_called")
    ap.add_argument("--samples-csv-dir", default="pipeline/sample_csvs")
    ap.add_argument("--panel-suffix", default=".union_sid_haem_vendor_probes.vcf.gz")
    ap.add_argument("--admix", nargs="+", required=True)
    ap.add_argument("--min-gq", type=int, default=20)
    ap.add_argument("--min-dp", type=int, default=20)
    ap.add_argument("--error-rate", type=float, default=0.01)
    ap.add_argument("--genes-bed", type=Path, default=Path("output/refseq109_genes.bed"))
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument(
        "--label",
        default=None,
        help="Run tag for the figure title (default: trailing token of --out stem, "
        "e.g. 'run10' from host_presence_manhattan_run10.png)",
    )
    args = ap.parse_args()
    run_label = args.label or args.out.stem.split("_")[-1]

    genes = load_genes(args.genes_bed)
    csv_index = index_csvs(args.samples_csv_dir)
    panels = []
    for a in args.admix:
        meta = csv_index.get(a)
        if meta is None:
            raise SystemExit(f"ADMIX sample {a!r} not found")
        markers, pooled = build_markers(
            meta,
            a,
            args.vcf_dir,
            args.panel_suffix,
            args.min_gq,
            args.min_dp,
            args.error_rate,
            genes,
        )
        panels.append((short_label(a), markers, pooled))

    n = len(panels)
    fig, axes = plt.subplots(n, 1, figsize=(14, 2.8 * n + 1), squeeze=False, sharex=False)
    flat = [ax for row in axes for ax in row]
    for ax, (label, markers, pooled) in zip(flat, panels):
        _draw_panel(ax, label, markers, pooled)
    flat[-1].set_xlabel("genomic position (hg38)")

    handles = [
        plt.Line2D([], [], marker="o", ls="", color=ALT_COLORS[0], label="marker (blue/green just"),
        plt.Line2D(
            [], [], marker="o", ls="", color=ALT_COLORS[1], label="alternate by chromosome)"
        ),
        plt.Line2D([], [], color="black", lw=1.1, alpha=0.7, label="per-chromosome mean"),
        plt.Line2D([], [], color="#c0392b", lw=1.1, label="sample-wide pooled host fraction"),
        plt.Line2D(
            [],
            [],
            marker="o",
            ls="",
            mfc="none",
            mec="#c0392b",
            mew=1.5,
            label="upregulated (Bonferroni), gene-labelled",
        ),
        plt.Line2D(
            [],
            [],
            marker="x",
            ls="",
            color="#999999",
            mew=1.6,
            label="artifact-filtered (strand/soft-clip/read-pos bias)",
        ),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False, fontsize=8.5)
    fig.suptitle(f"Host-presence signal along the genome (CNV view) — {run_label}", fontsize=12)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    fig.savefig(args.out, dpi=150)
    print(f"Wrote {args.out} (local only; x axis shows genomic coordinates)")


if __name__ == "__main__":
    main()
