#!/usr/bin/env python3
"""Union two or more capture BEDs into one sorted, merged intervals BED.

The GATK two-phase pipeline (``pipeline/Snakefile``) takes a single
``intervals:`` BED. When the chimerism markers come from one kit and the
sequencing run also carries a second kit (our deployment: the IDT rhAmpSeq
Sample ID SNPs inside a haem capture run), the pipeline needs the union of
both target sets so discovery covers every captured locus. This script
builds that union. The output basename becomes the pipeline's per-patient VCF
suffix (``<patient>.<basename>.vcf.gz``), so name it deliberately.

Input BEDs may be 3+ columns; only chrom/start/end are used. ``track``,
``browser`` and ``#`` lines are skipped. Overlapping or book-ended intervals
are merged (adjacent intervals too when ``--merge-gap`` is positive). The
output is a 3-column BED, sorted by contig name then start. Pure Python, no
bedtools.

Usage:
    python scripts/build_union_bed.py \\
        /path/to/idt_rhampseq_sid_SNPsQC.bed \\
        /path/to/haem_vendor_probes.bed \\
        --out output/union_sid_haem_vendor_probes.bed
"""

import argparse
import sys
from pathlib import Path


def read_bed(path: Path) -> list[tuple[str, int, int]]:
    """Read (chrom, start, end) from a BED, skipping header/comment lines."""
    intervals: list[tuple[str, int, int]] = []
    with open(path, encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise SystemExit(f"{path}:{lineno}: expected >= 3 tab-separated columns")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if end <= start:
                raise SystemExit(f"{path}:{lineno}: end <= start ({start}, {end})")
            intervals.append((chrom, start, end))
    return intervals


def _chrom_sort_key(chrom: str) -> tuple[int, int, str]:
    """chr1..chr22, chrX, chrY, chrM, then anything else alphabetically."""
    name = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if name.isdigit():
        return (0, int(name), "")
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if name.upper() in special:
        return (0, special[name.upper()], "")
    return (1, 0, chrom)


def merge_intervals(
    intervals: list[tuple[str, int, int]], merge_gap: int = 0
) -> list[tuple[str, int, int]]:
    """Sort and merge intervals that overlap, touch, or lie within ``merge_gap``."""
    ordered = sorted(intervals, key=lambda iv: (_chrom_sort_key(iv[0]), iv[1], iv[2]))
    merged: list[tuple[str, int, int]] = []
    for chrom, start, end in ordered:
        if merged and merged[-1][0] == chrom and start <= merged[-1][2] + merge_gap:
            prev = merged[-1]
            merged[-1] = (chrom, prev[1], max(prev[2], end))
        else:
            merged.append((chrom, start, end))
    return merged


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("beds", nargs="+", type=Path, help="Two or more input BED files")
    ap.add_argument("--out", type=Path, required=True, help="Output 3-column BED")
    ap.add_argument(
        "--merge-gap",
        type=int,
        default=0,
        help="Also merge intervals separated by at most this many bases (default 0: "
        "merge only overlapping or book-ended intervals)",
    )
    args = ap.parse_args()
    if len(args.beds) < 2:
        ap.error("give at least two BED files to union")

    per_file: list[int] = []
    everything: list[tuple[str, int, int]] = []
    for bed in args.beds:
        ivs = read_bed(bed)
        per_file.append(len(ivs))
        everything.extend(ivs)

    merged = merge_intervals(everything, merge_gap=args.merge_gap)
    bases = sum(end - start for _, start, end in merged)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as out:
        for chrom, start, end in merged:
            out.write(f"{chrom}\t{start}\t{end}\n")

    for bed, n in zip(args.beds, per_file):
        print(f"{bed}: {n} intervals", file=sys.stderr)
    print(f"-> {args.out}: {len(merged)} merged intervals, {bases:,} bp", file=sys.stderr)


if __name__ == "__main__":
    main()
