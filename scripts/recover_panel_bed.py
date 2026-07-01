#!/usr/bin/env python3
"""Recover an amplicon/MIP capture panel as a BED from BAM coverage.

Amplicon and MIP panels publish coordinates inconsistently: sometimes only a
kit name, sometimes a probe list in an awkward format, sometimes nothing. When
you have aligned BAMs but no usable BED, the panel can be rebuilt from coverage
alone, because these assays concentrate reads. Every captured locus piles
thousands of reads into one tight footprint (typically 80-100 bp), while
off-target background sits far below. Each target therefore shows up as a
high-depth cluster shared across samples, and the panel self-recovers.

For each genomic position this script counts how many BAMs cover it at
>= --min-depth (with the given MAPQ/BASEQ filters), keeps positions cleared in at
least --min-samples BAMs, and merges adjacent kept positions (gap <= --merge-gap)
into one interval per amplicon. The result is a sorted BED suitable for GATK
`intervals` or for `scripts/build_panel_vcf.py`.

Only samtools is required (already an allomix pipeline dependency); the interval
merge is done in Python, so bedtools is not needed.

Usage:
    python scripts/recover_panel_bed.py \
        --bam-glob 'output/bam/*.bam' \
        --out panel.bed

Set --min-samples relative to your cohort size. A base is kept only if it is
covered in at least that many BAMs, which is what separates real targets from
sporadic off-target coverage. A common choice is a clear majority of the cohort
(for 64 BAMs, 50 works well). Lower it to tolerate more per-sample dropout at the
cost of admitting more off-target positions; every recovered cluster is
amplicon-shaped regardless, so the threshold mainly trades completeness against
off-target inclusion.
"""

import argparse
import glob
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--bam-glob", required=True, help="Glob for the input BAMs.")
    p.add_argument("--out", required=True, help="Output BED path.")
    p.add_argument(
        "--samtools", default="samtools", help="samtools executable (default: %(default)s)."
    )
    p.add_argument(
        "--min-depth",
        type=int,
        default=100,
        help="Per-sample depth to count a base as covered (default: %(default)s).",
    )
    p.add_argument(
        "--min-mapq",
        type=int,
        default=20,
        help="samtools depth -Q mapping-quality filter (default: %(default)s).",
    )
    p.add_argument(
        "--min-baseq",
        type=int,
        default=20,
        help="samtools depth -q base-quality filter (default: %(default)s).",
    )
    p.add_argument(
        "--min-samples",
        type=int,
        default=50,
        help="Minimum BAMs covering a base at --min-depth to keep it; set relative "
        "to cohort size (default: %(default)s).",
    )
    p.add_argument(
        "--merge-gap",
        type=int,
        default=10,
        help="Merge kept bases separated by <= this many bp into one interval "
        "(default: %(default)s).",
    )
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

    autosomal = sum(
        1 for c, _, _ in intervals if c.removeprefix("chr") not in ("X", "Y", "M", "MT")
    )
    other = len(intervals) - autosomal
    print(f"BAMs used: {len(bams)}")
    print(f"Panel intervals recovered: {len(intervals)} ({autosomal} autosomal, {other} on X/Y/M)")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
