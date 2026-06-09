#!/usr/bin/env python3
"""Recover the SRP434573 MIP capture panel as a BED from BAM coverage.

The SRP434573 thesis (Chu Xufeng, HUST 2024) describes a ~1062-autosomal-SNP
MIP capture panel but publishes no coordinates or BED, and no companion paper
is indexed (see SACGF/allomix issue #16). The panel therefore has to be rebuilt
before GATK joint calling can run. Because this is a MIP/amplicon assay, every
captured locus piles thousands of reads into a single ~95 bp footprint, so the
panel self-recovers from coverage alone: each target is a tight, high-depth
cluster, and off-target background sits far below it.

This script reproduces (and extends, using all runs rather than two pures) the
issue #16 "laptop probe" that recovered ~1053 high-depth clusters. For each
genomic position it counts how many BAMs cover it at >= --min-depth (with the
same MAPQ/BASEQ filters used downstream), keeps positions that clear that bar in
at least --min-samples BAMs, and merges adjacent kept positions (gap
<= --merge-gap) into one interval per amplicon. The result is written as a
sorted BED suitable for GATK `intervals`.

Only samtools is required (already a pipeline dependency); the interval merge is
done in Python, so bedtools is not needed.

Usage:
    python scripts/build_srp434573_panel_bed.py \
        --bam-glob 'output/bam/*.bam' \
        --out paper/public_data/SRP434573/SRP434573.bed

Re-run with a different --min-samples to trade completeness against off-target
inclusion (every recovered cluster is amplicon-shaped regardless; the threshold
mainly decides how much per-sample dropout is tolerated).
"""

import argparse
import glob
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam-glob", default="output/bam/*.bam", help="Glob for the input BAMs (default: %(default)s).")
    p.add_argument("--out", default="paper/public_data/SRP434573/SRP434573.bed", help="Output BED path (default: %(default)s).")
    p.add_argument("--samtools", default="samtools", help="samtools executable (default: %(default)s).")
    p.add_argument("--min-depth", type=int, default=100, help="Per-sample depth to count a base as covered (default: %(default)s).")
    p.add_argument("--min-mapq", type=int, default=20, help="samtools depth -Q mapping-quality filter (default: %(default)s).")
    p.add_argument("--min-baseq", type=int, default=20, help="samtools depth -q base-quality filter (default: %(default)s).")
    p.add_argument("--min-samples", type=int, default=50, help="Minimum BAMs covering a base at --min-depth to keep it (default: %(default)s).")
    p.add_argument("--merge-gap", type=int, default=10, help="Merge kept bases separated by <= this many bp into one interval (default: %(default)s).")
    return p.parse_args()


def iter_panel_positions(args: argparse.Namespace, bams: list[str]):
    """Yield (chrom, pos) 1-based positions covered at >=min_depth in >=min_samples BAMs."""
    cmd = [args.samtools, "depth", "-Q", str(args.min_mapq), "-q", str(args.min_baseq), *bams]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
    assert proc.stdout is not None
    for line in proc.stdout:
        fields = line.split("\t")
        n_covered = sum(1 for d in fields[2:] if int(d) >= args.min_depth)
        if n_covered >= args.min_samples:
            yield fields[0], int(fields[1])
    proc.stdout.close()
    if proc.wait() != 0:
        sys.exit(f"samtools depth failed (exit {proc.returncode})")


def merge_intervals(positions, merge_gap: int):
    """Merge sorted 1-based positions into 0-based half-open BED intervals."""
    cur_chrom = None
    start = end = 0
    for chrom, pos in positions:
        s, e = pos - 1, pos  # 1-based pos -> BED [s, e)
        if chrom == cur_chrom and s - end <= merge_gap:
            end = e
        else:
            if cur_chrom is not None:
                yield cur_chrom, start, end
            cur_chrom, start, end = chrom, s, e
    if cur_chrom is not None:
        yield cur_chrom, start, end


def main() -> None:
    args = parse_args()
    bams = sorted(glob.glob(args.bam_glob))
    if not bams:
        sys.exit(f"No BAMs matched {args.bam_glob!r}")

    intervals = list(merge_intervals(iter_panel_positions(args, bams), args.merge_gap))

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as fh:
        for chrom, start, end in intervals:
            fh.write(f"{chrom}\t{start}\t{end}\n")

    autosomal = sum(1 for c, _, _ in intervals if c.removeprefix("chr") not in ("X", "Y", "M", "MT"))
    print(f"BAMs used: {len(bams)}")
    print(f"Panel intervals recovered: {len(intervals)} ({autosomal} autosomal, {len(intervals) - autosomal} on X/Y/M)")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
