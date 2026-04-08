#!/usr/bin/env python3
"""Measure per-marker amplification bias from real genotyping VCFs.

For each heterozygous call at a biallelic site, computes the deviation of
observed VAF from 0.5. Reports per-marker bias statistics and an overall
summary suitable for calibrating allomix simulation parameters.

Input: a text file listing VCF paths (one per line), from e.g.:
    find /tau/data/clinical_hg38/idt_rhampseq_sid/ -name '*.vcf.gz' > vcf_list.txt

Output (to stdout): tab-separated summary statistics (no patient identifiers,
no genomic coordinates — safe to share outside /tau).

Usage:
    python scripts/measure_panel_bias.py vcf_list.txt
    python scripts/measure_panel_bias.py vcf_list.txt --min-dp 100 --min-gq 20
"""

from __future__ import annotations

import argparse
import math
import sys
from collections import defaultdict
from pathlib import Path


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Measure per-marker amplification bias from het sites in VCFs.",
    )
    p.add_argument("vcf_list", help="Text file with one VCF path per line")
    p.add_argument("--min-dp", type=int, default=100, help="Min depth (default: 100)")
    p.add_argument("--min-gq", type=int, default=20, help="Min genotype quality (default: 20)")
    p.add_argument(
        "--output", default=None,
        help="Output file (default: stdout). Two files are written: "
             "<output>_per_marker.tsv and <output>_summary.tsv",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)

    # Read VCF list
    vcf_paths = []
    with open(args.vcf_list) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                vcf_paths.append(line)

    if not vcf_paths:
        print("ERROR: no VCF paths found in input file", file=sys.stderr)
        return 1

    print(f"Processing {len(vcf_paths)} VCFs ...", file=sys.stderr)

    try:
        from cyvcf2 import VCF
    except ImportError:
        print("ERROR: cyvcf2 required. Install with: pip install cyvcf2", file=sys.stderr)
        return 1

    # Collect per-marker het VAF deviations
    # Key: marker_index (sequential, no genomic coords in output)
    # We track by (chrom, pos, ref, alt) internally but only output anonymised index
    marker_deviations: dict[tuple, list[float]] = defaultdict(list)
    marker_depths: dict[tuple, list[int]] = defaultdict(list)
    n_samples = 0
    n_het_total = 0
    n_skipped_dp = 0
    n_skipped_gq = 0
    n_skipped_multiallelic = 0

    for i, vcf_path in enumerate(vcf_paths):
        if not Path(vcf_path).exists():
            print(f"  WARNING: {vcf_path} not found, skipping", file=sys.stderr)
            continue

        try:
            vcf = VCF(vcf_path)
        except Exception as e:
            print(f"  WARNING: failed to open {vcf_path}: {e}", file=sys.stderr)
            continue

        n_samples += 1
        sample_hets = 0

        for variant in vcf:
            # Skip multiallelic
            if len(variant.ALT) != 1:
                n_skipped_multiallelic += 1
                continue

            # Get genotype
            gt = variant.genotypes[0]  # first (only) sample
            alleles = (gt[0], gt[1])

            # Skip non-het
            if alleles[0] == alleles[1]:
                continue
            if alleles[0] < 0 or alleles[1] < 0:
                continue

            # Must be 0/1 het
            if set(alleles) != {0, 1}:
                continue

            # Check depth
            dp = variant.format("DP")
            if dp is not None:
                depth = int(dp[0][0])
                if depth < args.min_dp:
                    n_skipped_dp += 1
                    continue
            else:
                n_skipped_dp += 1
                continue

            # Check GQ
            gq = variant.format("GQ")
            if gq is not None:
                gq_val = int(gq[0][0])
                if gq_val < args.min_gq:
                    n_skipped_gq += 1
                    continue

            # Get AD
            ad = variant.format("AD")
            if ad is None:
                continue
            ref_count = int(ad[0][0])
            alt_count = int(ad[0][1])
            total = ref_count + alt_count
            if total == 0:
                continue

            vaf = alt_count / total
            deviation = vaf - 0.5

            key = (variant.CHROM, variant.POS, variant.REF, variant.ALT[0])
            marker_deviations[key].append(deviation)
            marker_depths[key].append(total)
            sample_hets += 1

        n_het_total += sample_hets
        if (i + 1) % 10 == 0 or (i + 1) == len(vcf_paths):
            print(f"  Processed {i+1}/{len(vcf_paths)} VCFs, "
                  f"{n_het_total} het observations so far", file=sys.stderr)

    if not marker_deviations:
        print("ERROR: no heterozygous observations found", file=sys.stderr)
        return 1

    # Compute per-marker statistics
    marker_stats = []
    all_biases = []
    all_abs_biases = []

    for key in sorted(marker_deviations.keys()):
        devs = marker_deviations[key]
        depths = marker_depths[key]
        n = len(devs)
        if n == 0:
            continue

        median_dev = sorted(devs)[n // 2]
        mean_dev = sum(devs) / n
        abs_devs = [abs(d) for d in devs]
        mean_abs_dev = sum(abs_devs) / n
        sd_dev = math.sqrt(sum((d - mean_dev) ** 2 for d in devs) / n) if n > 1 else 0.0
        mean_depth = sum(depths) / n

        marker_stats.append({
            "n_het": n,
            "median_bias": median_dev,
            "mean_bias": mean_dev,
            "mean_abs_bias": mean_abs_dev,
            "sd_within": sd_dev,
            "mean_depth": mean_depth,
        })
        all_biases.append(median_dev)
        all_abs_biases.append(abs(median_dev))

    # Overall summary
    n_markers = len(marker_stats)
    biases_sorted = sorted(all_biases)
    abs_biases_sorted = sorted(all_abs_biases)

    overall_mean_bias = sum(all_biases) / n_markers
    overall_median_bias = biases_sorted[n_markers // 2]
    overall_sd_bias = math.sqrt(
        sum((b - overall_mean_bias) ** 2 for b in all_biases) / n_markers
    )
    overall_mean_abs_bias = sum(all_abs_biases) / n_markers
    overall_median_abs_bias = abs_biases_sorted[n_markers // 2]

    # Percentiles of absolute bias
    p25_idx = max(0, n_markers // 4)
    p75_idx = min(n_markers - 1, 3 * n_markers // 4)
    p95_idx = min(n_markers - 1, int(0.95 * n_markers))

    # Print summary to stdout (no patient identifiers, no coordinates)
    print("=" * 60)
    print("PANEL AMPLIFICATION BIAS SUMMARY")
    print("=" * 60)
    print(f"VCFs processed:              {n_samples}")
    print(f"Markers with het obs:        {n_markers}")
    print(f"Total het observations:      {n_het_total}")
    print(f"Skipped (low depth):         {n_skipped_dp}")
    print(f"Skipped (low GQ):            {n_skipped_gq}")
    print(f"Skipped (multiallelic):      {n_skipped_multiallelic}")
    print()
    print("Per-marker bias (median het VAF - 0.5):")
    print(f"  Mean bias:                 {overall_mean_bias:+.4f}")
    print(f"  Median bias:               {overall_median_bias:+.4f}")
    print(f"  SD of per-marker biases:   {overall_sd_bias:.4f}")
    print(f"  Mean |bias|:               {overall_mean_abs_bias:.4f}")
    print(f"  Median |bias|:             {overall_median_abs_bias:.4f}")
    print(f"  25th pct |bias|:           {abs_biases_sorted[p25_idx]:.4f}")
    print(f"  75th pct |bias|:           {abs_biases_sorted[p75_idx]:.4f}")
    print(f"  95th pct |bias|:           {abs_biases_sorted[p95_idx]:.4f}")
    print(f"  Max |bias|:                {max(all_abs_biases):.4f}")
    print()
    print(">>> For allomix simulation, use --bias-sd {:.3f}".format(overall_sd_bias))
    print("=" * 60)

    # Write per-marker detail (anonymised — sequential index only)
    if args.output:
        per_marker_path = f"{args.output}_per_marker.tsv"
        summary_path = f"{args.output}_summary.tsv"
    else:
        per_marker_path = None
        summary_path = None

    if per_marker_path:
        Path(per_marker_path).parent.mkdir(parents=True, exist_ok=True)
        with open(per_marker_path, "w") as f:
            f.write("marker_index\tn_het\tmedian_bias\tmean_bias\tmean_abs_bias\t"
                    "sd_within\tmean_depth\n")
            for i, s in enumerate(marker_stats):
                f.write(f"{i}\t{s['n_het']}\t{s['median_bias']:.6f}\t"
                        f"{s['mean_bias']:.6f}\t{s['mean_abs_bias']:.6f}\t"
                        f"{s['sd_within']:.6f}\t{s['mean_depth']:.0f}\n")
        print(f"\nPer-marker detail: {per_marker_path}", file=sys.stderr)

    if summary_path:
        with open(summary_path, "w") as f:
            f.write("metric\tvalue\n")
            f.write(f"n_vcfs\t{n_samples}\n")
            f.write(f"n_markers\t{n_markers}\n")
            f.write(f"n_het_total\t{n_het_total}\n")
            f.write(f"mean_bias\t{overall_mean_bias:.6f}\n")
            f.write(f"median_bias\t{overall_median_bias:.6f}\n")
            f.write(f"sd_bias\t{overall_sd_bias:.6f}\n")
            f.write(f"mean_abs_bias\t{overall_mean_abs_bias:.6f}\n")
            f.write(f"median_abs_bias\t{overall_median_abs_bias:.6f}\n")
            f.write(f"p95_abs_bias\t{abs_biases_sorted[p95_idx]:.6f}\n")
            f.write(f"max_abs_bias\t{max(all_abs_biases):.6f}\n")
        print(f"Summary: {summary_path}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
