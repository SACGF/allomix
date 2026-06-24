#!/usr/bin/env python3
"""Karyogram of informative markers along the genome for host vs donor(s).

Internal diagnostic (not part of the allomix package).

Each shared marker is classified by how the host and donor genotypes relate,
which for biallelic SNPs is the identity-by-state (IBS) count:

    IBS0  opposite homozygotes (0/0 vs 1/1): the two share NO allele. This
          state is impossible wherever the pair share a haplotype, so a stretch
          of genome with no IBS0 markers is evidence of identity-by-descent.
    IBS1  homozygote vs heterozygote: exactly one allele shared (partially
          informative).
    IBS2  same genotype: both alleles compatible with sharing (non-informative).

For sibling donors the informative markers (IBS0 + IBS1) are not spread evenly:
the pair share a haplotype across long IBD blocks, so those blocks are depleted
of IBS0 markers and read as gaps in the red ticks. Plotting marker position
along the genome (rather than sorting by magnitude) makes those blocks visible,
the way a karyotype/ideogram does. With a small panel (tens of genome-wide SNPs)
the picture is coarse, but the clustering of shared vs differing markers is still
informative for a relatedness sanity check.

The x axis is genomic position, so the figure is written to a LOCAL file only
(see CLAUDE.md): open it on this machine, do not surface coordinates to stdout.

Usage:
    # Single donor from a joint-called panel VCF (host + donor as named samples):
    python scripts/plot_informative_karyogram.py \
        --vcf output/joint_called/PATIENT.panel.vcf.gz \
        --host HOST_ID --donor DONOR_ID \
        --out output/informative_karyogram.png

    # Two donors, each drawn in its own lane:
    python scripts/plot_informative_karyogram.py \
        --vcf tests/test_data/joint_multi_donor.vcf \
        --host HOST --donor DONOR1 DONOR2 \
        --out output/karyogram_multidonor.png
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

from allomix.genotype import is_sex_chrom, marker_key, parse_vcf  # noqa: E402

# hg38 chromosome lengths (bp), for a proportional genome x axis. Autosomes plus
# the sex/mito contigs (drawn only with --include-sex).
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

# IBS class colours. IBS0 is the relatedness-informative state (red, tall ticks);
# IBS1 is partial (amber); IBS2 is the shared/non-informative background (grey).
IBS_COLOR = {0: "#c0392b", 1: "#e08e0b", 2: "#bdbdbd"}
IBS_LABEL = {
    0: "IBS0 (opposite homozygotes)",
    1: "IBS1 (one shared allele)",
    2: "IBS2 (shared genotype)",
}


def genome_axis(include_sex: bool) -> tuple[list[str], dict[str, int], int]:
    """Build the concatenated-genome coordinate scaffold.

    Args:
        include_sex: Keep the sex and mitochondrial contigs in the axis.

    Returns:
        (chrom_order, offset, genome_len): the chromosomes in plotting order,
        each chromosome's start offset along the concatenated axis, and the
        total genome length.
    """
    chrom_order = [c for c in HG38 if include_sex or not is_sex_chrom(c)]
    offset: dict[str, int] = {}
    acc = 0
    for c in chrom_order:
        offset[c] = acc
        acc += HG38[c]
    return chrom_order, offset, acc


def ibs_class(host_gt: tuple[int, int], donor_gt: tuple[int, int]) -> int:
    """IBS count (0, 1, or 2) for two biallelic diploid genotypes.

    Uses allele dosage, so it does not assume which allele is REF. IBS0 is the
    only state that cannot occur where the pair share a haplotype.

    Args:
        host_gt: Host genotype as allele indices, e.g. (0, 0).
        donor_gt: Donor genotype as allele indices.

    Returns:
        2 if the genotypes are identical, 1 if they share exactly one allele,
        0 if they share no allele (opposite homozygotes).
    """
    h = host_gt[0] + host_gt[1]
    d = donor_gt[0] + donor_gt[1]
    return 2 - abs(h - d)


def load_pair(
    vcf: Path, host: str, donor: str, min_gq: int, offset: dict[str, int]
) -> list[tuple[float, int]]:
    """Classify every marker shared by host and donor into (x, ibs).

    Args:
        vcf: Joint VCF holding both samples.
        host: Host sample name (or column index as a string of digits).
        donor: Donor sample name.
        min_gq: Minimum GQ for both samples.
        offset: Chromosome start offsets from ``genome_axis``; markers on
            contigs absent from it (e.g. sex chroms when excluded) are skipped.

    Returns:
        List of (genome_x, ibs_class) for each shared, on-axis marker.
    """
    host_markers = parse_vcf(vcf, sample=host, min_gq=min_gq, gt_ad_consistency=True)
    donor_markers = parse_vcf(vcf, sample=donor, min_gq=min_gq, gt_ad_consistency=True)
    donor_idx = {marker_key(m): m for m in donor_markers}
    out: list[tuple[float, int]] = []
    for h in host_markers:
        d = donor_idx.get(marker_key(h))
        if d is None or h.chrom not in offset:
            continue
        out.append((offset[h.chrom] + h.pos, ibs_class(h.gt, d.gt)))
    return out


def _draw_lane(ax, y: float, markers: list[tuple[float, int]]) -> dict[int, int]:
    """Draw one donor's marker ticks at height ``y`` and return per-IBS counts."""
    counts = {0: 0, 1: 0, 2: 0}
    # Tick half-heights: IBS0 stands tallest so the relatedness signal pops out
    # of the shared-genotype background.
    half = {0: 0.34, 1: 0.24, 2: 0.16}
    for x, ibs in markers:
        counts[ibs] += 1
        ax.plot(
            [x, x],
            [y - half[ibs], y + half[ibs]],
            color=IBS_COLOR[ibs],
            lw=1.4 if ibs == 0 else 1.0,
            solid_capstyle="butt",
            zorder=3 if ibs == 0 else 2,
        )
    return counts


def plot(
    lanes: list[tuple[str, list[tuple[float, int]]]],
    chrom_order: list[str],
    offset: dict[str, int],
    genome_len: int,
    title: str,
    out: Path,
) -> None:
    """Draw the karyogram (one lane per donor) and write it to ``out``."""
    n = len(lanes)
    fig, ax = plt.subplots(figsize=(15, 1.1 * n + 2.2))

    # Alternating chromosome background bands so position is readable.
    for i, c in enumerate(chrom_order):
        if i % 2:
            ax.axvspan(offset[c], offset[c] + HG38[c], color="#f5f5f5", zorder=0)

    yticklabels = []
    for lane_i, (label, markers) in enumerate(lanes):
        y = n - 1 - lane_i  # first donor on top
        ax.axhline(y, color="#dddddd", lw=0.8, zorder=1)
        counts = _draw_lane(ax, y, markers)
        n_inf = counts[0] + counts[1]
        n_tot = n_inf + counts[2]
        rate = (100.0 * counts[0] / n_tot) if n_tot else 0.0
        yticklabels.append(f"{label}\nM:{n_inf}/{n_tot}  IBS0:{counts[0]} ({rate:.0f}%)")

    ax.set_ylim(-0.6, n - 0.4)
    ax.set_yticks(list(range(n)))
    ax.set_yticklabels(list(reversed(yticklabels)), fontsize=8)
    ax.set_xlim(0, genome_len)
    ticks = [offset[c] + HG38[c] / 2 for c in chrom_order]
    ax.set_xticks(ticks)
    ax.set_xticklabels([c[3:] for c in chrom_order], fontsize=7)
    ax.set_xlabel("genomic position (hg38)")
    ax.set_title(title)

    handles = [Line2D([], [], color=IBS_COLOR[k], lw=2.0, label=IBS_LABEL[k]) for k in (0, 1, 2)]
    ax.legend(
        handles=handles,
        ncol=3,
        fontsize=8,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        framealpha=0.9,
    )
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    print(f"Wrote {out}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--vcf", type=Path, required=True, help="Joint VCF with host and donor samples")
    ap.add_argument("--host", required=True, help="Host sample name (or 0-based column index)")
    ap.add_argument(
        "--donor",
        nargs="+",
        required=True,
        help="One or more donor sample names; each is drawn in its own lane",
    )
    ap.add_argument("--min-gq", type=int, default=20, help="Minimum GQ for host/donor (default 20)")
    ap.add_argument(
        "--include-sex",
        action="store_true",
        help="Keep sex/mitochondrial contigs (excluded by default, as for chimerism)",
    )
    ap.add_argument(
        "--title",
        default="Informative-marker karyogram: host vs donor(s)",
        help="Figure title",
    )
    ap.add_argument("--out", type=Path, required=True, help="Output PNG path (local file only)")
    args = ap.parse_args()

    chrom_order, offset, genome_len = genome_axis(args.include_sex)
    lanes = [
        (donor, load_pair(args.vcf, args.host, donor, args.min_gq, offset)) for donor in args.donor
    ]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    plot(lanes, chrom_order, offset, genome_len, args.title, args.out)


if __name__ == "__main__":
    main()
