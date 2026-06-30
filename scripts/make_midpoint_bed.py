#!/usr/bin/env python3
"""Derive a thin background-site BED (one position per amplicon) from a capture BED.

The host-presence detector needs a per-site sequencing-error background, and the
hardest direction to measure is ref->alt: the rate at which a true hom-ref site
emits stray ALT reads. Variant-only joint-call VCFs contain no all-hom-ref sites,
so that direction is unmeasurable from them (see doc/joint_calling.md and SACGF
issue #22).

The capture panel covers ~96 bp per amplicon but only one polymorphic SNP per
amplicon is emitted, leaving the rest of each amplicon as deep, covered,
almost-always-invariant sequence. This script picks one such position per
interval (the midpoint) so a forced raw pileup there (bcftools mpileup, NOT
GATK, which strips minority ALT reads at hom-ref blocks) measures the ref->alt
background cleanly. The midpoint is almost always hom-ref in every individual;
any that is actually polymorphic is dropped later by the estimator's VAF guard
(``allomix estimate-errors --max-vaf-homref``), so no external allele-frequency
resource is needed.

Output is a sorted, single-base BED (one ``chrom  pos0  pos0+1`` line per input
interval), suitable for ``bcftools mpileup -R``. Intervals shorter than
``--min-width`` are skipped (too small to hold a distinct background position).

Usage:
    python scripts/make_midpoint_bed.py \
        --bed paper/public_data/SRP434573/SRP434573.bed \
        --out paper/public_data/SRP434573/SRP434573.midpoints.bed
"""

import argparse
import sys
from pathlib import Path


def read_intervals(bed_path: Path) -> list[tuple[str, int, int]]:
    """Read a BED into a list of ``(chrom, start, end)`` (0-based half-open)."""
    intervals: list[tuple[str, int, int]] = []
    with open(bed_path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise SystemExit(f"Malformed BED line (need >=3 columns): {line!r}")
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if end <= start:
                raise SystemExit(f"Non-positive interval width: {line!r}")
            intervals.append((chrom, start, end))
    return intervals


def midpoints(intervals: list[tuple[str, int, int]], min_width: int) -> list[tuple[str, int]]:
    """Return one 0-based midpoint position per interval at least ``min_width`` wide."""
    out: list[tuple[str, int]] = []
    for chrom, start, end in intervals:
        if end - start < min_width:
            continue
        out.append((chrom, (start + end) // 2))
    return out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument(
        "--bed", required=True, type=Path, help="Capture-panel BED (the pipeline `intervals`)"
    )
    parser.add_argument("--out", required=True, type=Path, help="Output single-base midpoint BED")
    parser.add_argument(
        "--min-width",
        type=int,
        default=10,
        help="Skip intervals narrower than this many bp (default: 10)",
    )
    args = parser.parse_args(argv)

    intervals = read_intervals(args.bed)
    mids = midpoints(intervals, args.min_width)
    skipped = len(intervals) - len(mids)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as fh:
        for chrom, pos0 in mids:
            fh.write(f"{chrom}\t{pos0}\t{pos0 + 1}\n")

    sys.stderr.write(
        f"Wrote {len(mids)} midpoint positions from {len(intervals)} intervals "
        f"({skipped} skipped < {args.min_width} bp) -> {args.out}\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
