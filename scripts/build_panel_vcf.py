#!/usr/bin/env python3
"""Build a sites-only panel VCF from a BED of marker positions.

Reads REF/ALT alleles from a source VCF (a previous joint call covering the
target sites). The output VCF is suitable for forced genotyping:

  - GATK HaplotypeCaller ``--alleles panel.vcf.gz``
  - bcftools call ``-C alleles -T panel.targets.tsv.gz``

This means HaplotypeCaller emits a genotype at every panel site even when
neither HOST nor DONOR carries the variant in a small per-patient joint
call, recovering the marker count that a single big joint call would have
produced (because joint-calling power scales with sample count, a 2-sample
HOST+DONOR call discovers fewer variants than a 19-sample pooled call).

The BED is the authoritative list of sites. Positions in the BED but
absent from the source VCF are reported and skipped (typically sex-
chromosome or non-SNP markers that the joint call never saw an alt at).

Usage:
    python scripts/build_panel_vcf.py \\
        output/idt_rhampseq_sid_SNPsQC.bed \\
        output/joint_called/joint_called.idt_rhampseq_sid_SNPsQC.vcf.gz \\
        output/idt_rhampseq_sid_panel.vcf.gz
"""

import argparse
import subprocess
import sys
from pathlib import Path

from cyvcf2 import VCF


def _parse_bed(path: Path) -> list[tuple[str, int, str]]:
    """Return (chrom, pos_1based, rsid) for each BED row (skipping header)."""
    rows = []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom, start = parts[0], int(parts[1])
            pos_1based = start + 1
            # rsID lives in the 4th column, semicolon-delimited (rsID;...);
            # fall back to '.' if shape is unexpected.
            name = "."
            if len(parts) >= 4 and parts[3]:
                name = parts[3].split(";")[0]
            rows.append((chrom, pos_1based, name))
    return rows


def _index_source_vcf(path: Path) -> dict[tuple[str, int], tuple[str, str]]:
    """Map (chrom, pos_1based) -> (REF, ALT) for every record in path.

    Multi-allelic sites are skipped (we want a bi-allelic panel).
    """
    by_pos: dict[tuple[str, int], tuple[str, str]] = {}
    vcf = VCF(str(path))
    for v in vcf:
        if v.ALT is None or len(v.ALT) != 1:
            continue
        by_pos[(v.CHROM, v.POS)] = (v.REF, v.ALT[0])
    vcf.close()
    return by_pos


_HEADER = """\
##fileformat=VCFv4.2
##source=build_panel_vcf.py
##INFO=<ID=PANEL,Number=0,Type=Flag,Description="Site is in the chimerism panel">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    parser.add_argument("bed", type=Path, help="BED of panel marker positions")
    parser.add_argument(
        "source_vcf",
        type=Path,
        help="Existing VCF to read REF/ALT from for each panel position",
    )
    parser.add_argument(
        "out_vcf",
        type=Path,
        help="Output VCF path (.vcf or .vcf.gz; bgzipped + tabix-indexed if .vcf.gz)",
    )
    parser.add_argument(
        "--bgzip",
        default="bgzip",
        help="bgzip executable (default: bgzip on $PATH)",
    )
    parser.add_argument(
        "--tabix",
        default="tabix",
        help="tabix executable (default: tabix on $PATH)",
    )
    args = parser.parse_args()

    bed_rows = _parse_bed(args.bed)
    source = _index_source_vcf(args.source_vcf)

    matched: list[tuple[str, int, str, str, str]] = []
    missing: list[tuple[str, int, str]] = []
    for chrom, pos, name in bed_rows:
        if (chrom, pos) in source:
            ref, alt = source[(chrom, pos)]
            matched.append((chrom, pos, name, ref, alt))
        else:
            missing.append((chrom, pos, name))

    # Sort by chrom (natural-ish: numeric first, then X/Y/M) then position.
    def _chrom_key(c: str) -> tuple[int, str]:
        s = c.removeprefix("chr")
        if s.isdigit():
            return (0, f"{int(s):03d}")
        return (1, s)

    matched.sort(key=lambda r: (_chrom_key(r[0]), r[1]))

    plain_path = args.out_vcf.with_suffix("") if args.out_vcf.suffix == ".gz" else args.out_vcf
    with open(plain_path, "w") as out:
        out.write(_HEADER)
        for chrom, pos, name, ref, alt in matched:
            out.write(f"{chrom}\t{pos}\t{name}\t{ref}\t{alt}\t.\tPASS\tPANEL\n")

    if args.out_vcf.suffix == ".gz":
        subprocess.run([args.bgzip, "-f", str(plain_path)], check=True)
        subprocess.run([args.tabix, "-fp", "vcf", str(args.out_vcf)], check=True)

    print(f"Wrote {args.out_vcf}: {len(matched)} sites", file=sys.stderr)
    if missing:
        print(
            f"Skipped {len(missing)} BED sites with no REF/ALT in source VCF:",
            file=sys.stderr,
        )
        for chrom, pos, name in missing:
            print(f"  {chrom}:{pos}  {name}", file=sys.stderr)


if __name__ == "__main__":
    main()
