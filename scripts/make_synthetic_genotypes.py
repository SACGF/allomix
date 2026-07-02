#!/usr/bin/env python3
"""Create synthetic host and donor genotype VCFs with 100 biallelic SNPs.

Generates two VCFs with realistic structure and diverse genotype combinations
to maximise the number of informative markers for chimerism testing.

Usage:
    python scripts/make_synthetic_genotypes.py --outdir tests/test_data

Output:
    tests/test_data/host.vcf
    tests/test_data/donor.vcf
"""

import argparse
import logging
import random
import sys
from pathlib import Path

import numpy as np

log = logging.getLogger(__name__)

N_MARKERS = 100

CHROMS = [f"chr{i}" for i in range(1, 23)]

TRANSITIONS = [("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")]
TRANSVERSIONS = [
    ("A", "C"),
    ("A", "T"),
    ("G", "C"),
    ("G", "T"),
    ("C", "A"),
    ("C", "G"),
    ("T", "A"),
    ("T", "G"),
]

# Fraction of SNPs that are transitions (empirically ~67% in human genome)
TRANSITION_FRACTION = 0.67

DEFAULT_DEPTH_MEAN = 2000
DEPTH_SD = 300
MIN_DEPTH = 500
HOM_SEQUENCING_ERROR_RATE = 0.002  # error rate producing noise reads in hom calls
HET_ALLELIC_IMBALANCE = 0.05  # max VAF deviation from 0.5 for realistic het calls
EXPECTED_HET_VAF = 0.5

DEFAULT_GQ = 99
DEFAULT_QUAL = 10000
CONTIG_LENGTH = 248_956_422  # hg38 chr1 length, used as placeholder for all contigs

START_POSITION = 1_000_000
POSITION_SPACING = 100_000

# (host_gt, donor_gt) pairs cycled per marker, weighted toward informative combos.
GT_COMBOS = [
    # fully informative (types 0, 1)
    ("0/0", "1/1"),
    ("1/1", "0/0"),
    ("0/0", "1/1"),
    ("1/1", "0/0"),
    # partially informative (types 10, 11, 20, 21)
    ("0/1", "0/0"),
    ("0/1", "1/1"),
    ("0/0", "0/1"),
    ("1/1", "0/1"),
    # non-informative (same genotype)
    ("0/0", "0/0"),
    ("0/1", "0/1"),
]


def _binomial(rng: random.Random, n: int, p: float) -> int:
    """Cross-version binomial draw (random.Random.binomialvariate is 3.12+ only)."""
    return int(np.random.default_rng(rng.getrandbits(32)).binomial(n, p))


def gt_to_ad(gt: str, depth: int, rng: random.Random) -> str:
    """Generate realistic AD field for a genotype at a given depth."""
    if gt == "0/0":
        alt = _binomial(rng, depth, HOM_SEQUENCING_ERROR_RATE)
        return f"{depth - alt},{alt}"
    elif gt == "1/1":
        ref = _binomial(rng, depth, HOM_SEQUENCING_ERROR_RATE)
        return f"{ref},{depth - ref}"
    else:  # 0/1
        bias = rng.uniform(-HET_ALLELIC_IMBALANCE, HET_ALLELIC_IMBALANCE)
        alt = _binomial(rng, depth, EXPECTED_HET_VAF + bias)
        return f"{depth - alt},{alt}"


def gt_to_af(_gt: str, ad: str) -> str:
    """Compute AF from AD."""
    parts = ad.split(",")
    ref_c, alt_c = int(parts[0]), int(parts[1])
    total = ref_c + alt_c
    if total == 0:
        return "0"
    return f"{alt_c / total:.4f}"


def make_vcf(
    sample_name: str,
    genotypes: list[str],
    positions: list[tuple[str, int, str, str]],
    rng: random.Random,
    depth_mean: int = DEFAULT_DEPTH_MEAN,
) -> str:
    """Build a VCF string."""
    lines = []

    lines.append("##fileformat=VCFv4.2")
    lines.append('##FILTER=<ID=PASS,Description="All filters passed">')
    lines.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    lines.append(
        "##FORMAT=<ID=AD,Number=R,Type=Integer,"
        'Description="Allelic depths for the ref and alt alleles">'
    )
    lines.append('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">')
    lines.append('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">')
    lines.append('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Variant allele frequency">')
    lines.append(
        "##FORMAT=<ID=PL,Number=G,Type=Integer,"
        'Description="Phred-scaled likelihoods for genotypes">'
    )
    lines.append('##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">')
    lines.append('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">')
    lines.append('##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles">')
    lines.append('##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth">')
    for c in CHROMS:
        lines.append(f"##contig=<ID={c},length={CONTIG_LENGTH}>")
    lines.append(f"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{sample_name}")

    for (chrom, pos, ref, alt), gt in zip(positions, genotypes):
        depth = max(MIN_DEPTH, int(rng.gauss(depth_mean, DEPTH_SD)))
        ad = gt_to_ad(gt, depth, rng)
        af = gt_to_af(gt, ad)
        dp = sum(int(x) for x in ad.split(","))
        gq = DEFAULT_GQ

        info = f"AC=1;AF=0.5;AN=2;DP={dp * 100}"
        fmt = "GT:AD:DP:GQ:AF"
        sample = f"{gt}:{ad}:{dp}:{gq}:{af}"

        lines.append(
            f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t{DEFAULT_QUAL}\tPASS\t{info}\t{fmt}\t{sample}"
        )

    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Create synthetic host + donor VCFs.")
    parser.add_argument(
        "--outdir",
        default="output/genotypes",
        help="Output directory (default: output/genotypes)",
    )
    parser.add_argument("--seed", type=int, default=12345, help="Random seed")
    parser.add_argument("--n-markers", type=int, default=N_MARKERS, help="Number of SNPs")
    args = parser.parse_args(argv)

    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    rng = random.Random(args.seed)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    n = args.n_markers

    positions = []
    for i in range(n):
        chrom = CHROMS[i % len(CHROMS)]
        pos = START_POSITION + i * POSITION_SPACING
        if rng.random() < TRANSITION_FRACTION:
            ref, alt = rng.choice(TRANSITIONS)
        else:
            ref, alt = rng.choice(TRANSVERSIONS)
        positions.append((chrom, pos, ref, alt))

    positions.sort(key=lambda x: (CHROMS.index(x[0]), x[1]))

    host_gts = []
    donor_gts = []
    for i in range(n):
        combo = GT_COMBOS[i % len(GT_COMBOS)]
        host_gts.append(combo[0])
        donor_gts.append(combo[1])

    host_vcf = make_vcf("HOST_SYNTHETIC", host_gts, positions, random.Random(rng.randint(0, 2**31)))
    donor_vcf = make_vcf(
        "DONOR_SYNTHETIC", donor_gts, positions, random.Random(rng.randint(0, 2**31))
    )

    host_path = outdir / "host.vcf"
    donor_path = outdir / "donor.vcf"
    host_path.write_text(host_vcf)
    donor_path.write_text(donor_vcf)

    n_informative = sum(1 for h, d in zip(host_gts, donor_gts) if h != d)
    log.info("Created %d markers (%d informative) in %s/", n, n_informative, outdir)
    log.info("  %s", host_path)
    log.info("  %s", donor_path)

    return 0


if __name__ == "__main__":
    sys.exit(main())
