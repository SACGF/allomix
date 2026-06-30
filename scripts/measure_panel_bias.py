#!/usr/bin/env python3
"""Measure per-marker panel characteristics from joint-called genotyping VCFs.

Collects per-marker statistics from heterozygous and all genotype calls:
  - Amplification bias (het VAF deviation from 0.5)
  - Locus dropout rate (no-call / ./.  frequency)
  - Per-marker call rate
  - Depth distribution per marker and across panel
  - Het/hom ratio vs HWE expectation (allele dropout signal)

Input: a text file listing VCF paths (one per line), from e.g.:
    find /tau/data/clinical_hg38/idt_rhampseq_sid/ -path '*/2_variants/*.gatk.hg38.vcf.gz' \
        -not -path '*/gatk_per_sample/*' > vcf_list.txt

Output (to stdout): summary statistics (no patient identifiers,
no genomic coordinates — safe to share outside /tau).

Usage:
    python scripts/measure_panel_bias.py vcf_list.txt
    python scripts/measure_panel_bias.py vcf_list.txt --min-dp 100 --min-gq 20
    python scripts/measure_panel_bias.py vcf_list.txt --output output/panel_stats
"""

import argparse
import csv
import logging
import math
import sys
from collections import defaultdict
from pathlib import Path

from cyvcf2 import VCF

log = logging.getLogger(__name__)

EXPECTED_HET_VAF = 0.5  # ideal heterozygous VAF
MIN_CALLED_FOR_HWE = 10  # min genotyped samples to compute HWE het ratio
LOW_HET_RATIO_THRESHOLD = 0.8  # het/HWE ratio below this flags allele dropout
HIGH_NOCALL_RATE = 0.05  # markers above this no-call rate are flagged
NEGLIGIBLE_ADO_THRESHOLD = 0.001  # allele dropout estimate below this is negligible


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Measure per-marker panel characteristics from genotyping VCFs.",
    )
    p.add_argument("vcf_list", help="Text file with one VCF path per line")
    p.add_argument(
        "--min-dp", type=int, default=100, help="Min depth for bias stats (default: 100)"
    )
    p.add_argument("--min-gq", type=int, default=20, help="Min GQ for bias stats (default: 20)")
    p.add_argument(
        "--output",
        default=None,
        help="Output file prefix. Writes <output>_per_marker.tsv and <output>_summary.tsv",
    )
    return p.parse_args(argv)


def _percentile(sorted_vals: list[float], pct: float) -> float:
    """Simple percentile from a pre-sorted list."""
    if not sorted_vals:
        return float("nan")
    idx = min(len(sorted_vals) - 1, int(pct / 100.0 * len(sorted_vals)))
    return sorted_vals[idx]


def _sd(vals: list[float]) -> float:
    if len(vals) < 2:
        return 0.0
    mean = sum(vals) / len(vals)
    return math.sqrt(sum((v - mean) ** 2 for v in vals) / len(vals))


def _cv(vals: list[float]) -> float:
    if len(vals) < 2:
        return 0.0
    mean = sum(vals) / len(vals)
    if mean == 0:
        return 0.0
    return _sd(vals) / mean


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    vcf_paths = []
    with open(args.vcf_list, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                vcf_paths.append(line)

    if not vcf_paths:
        log.error("no VCF paths found in input file")
        return 1

    log.info("Processing %d VCFs ...", len(vcf_paths))

    # Per-marker accumulators, keyed by (chrom, pos, ref, alt)
    marker_het_vafs: dict[tuple, list[float]] = defaultdict(list)
    marker_depths: dict[tuple, list[int]] = defaultdict(list)
    marker_gt_counts: dict[tuple, dict[str, int]] = defaultdict(
        lambda: {"hom_ref": 0, "het": 0, "hom_alt": 0, "no_call": 0, "total": 0}
    )

    sample_depth_cvs: list[float] = []  # depth CV across markers within each sample
    sample_nocall_rates: list[float] = []  # fraction of markers with no-call per sample

    n_vcfs = 0
    n_total_samples = 0
    n_het_total = 0
    n_skipped_dp = 0
    n_skipped_gq = 0
    n_skipped_multi = 0

    for vi, vcf_path in enumerate(vcf_paths):
        if not Path(vcf_path).exists():
            log.warning("%s not found, skipping", vcf_path)
            continue

        try:
            vcf = VCF(vcf_path)
        except Exception as e:  # pylint: disable=broad-exception-caught
            log.warning("failed to open %s: %s", vcf_path, e)
            continue

        n_vcfs += 1
        n_samples_in_vcf = len(vcf.samples)
        n_total_samples += n_samples_in_vcf

        sample_depths_this_vcf: list[list[int]] = [[] for _ in range(n_samples_in_vcf)]
        sample_nocalls_this_vcf: list[int] = [0] * n_samples_in_vcf
        sample_total_markers: list[int] = [0] * n_samples_in_vcf

        for variant in vcf:
            if len(variant.ALT) != 1:
                n_skipped_multi += 1
                continue

            key = (variant.CHROM, variant.POS, variant.REF, variant.ALT[0])

            for si in range(n_samples_in_vcf):
                gt = variant.genotypes[si]
                alleles = (gt[0], gt[1])
                counts = marker_gt_counts[key]
                counts["total"] += 1
                sample_total_markers[si] += 1

                if alleles[0] < 0 or alleles[1] < 0:
                    counts["no_call"] += 1
                    sample_nocalls_this_vcf[si] += 1
                    continue

                a_set = (alleles[0], alleles[1])
                is_het = a_set[0] != a_set[1]

                if is_het:
                    counts["het"] += 1
                elif alleles[0] == 0:
                    counts["hom_ref"] += 1
                else:
                    counts["hom_alt"] += 1

                # Depth for all calls, not just het
                dp = variant.format("DP")
                if dp is not None:
                    depth = int(dp[si][0])
                    marker_depths[key].append(depth)
                    sample_depths_this_vcf[si].append(depth)

                # Collect het VAF for bias measurement; filters here apply to bias only
                if is_het and set(alleles) == {0, 1}:
                    if dp is not None:
                        depth = int(dp[si][0])
                        if depth < args.min_dp:
                            n_skipped_dp += 1
                            continue
                    else:
                        n_skipped_dp += 1
                        continue

                    gq = variant.format("GQ")
                    if gq is not None:
                        gq_val = int(gq[si][0])
                        if gq_val < args.min_gq:
                            n_skipped_gq += 1
                            continue

                    ad = variant.format("AD")
                    if ad is None:
                        continue
                    ref_count = int(ad[si][0])
                    alt_count = int(ad[si][1])
                    total = ref_count + alt_count
                    if total == 0:
                        continue

                    vaf = alt_count / total
                    marker_het_vafs[key].append(vaf)
                    n_het_total += 1

        for si in range(n_samples_in_vcf):
            if sample_depths_this_vcf[si]:
                sample_depth_cvs.append(_cv(sample_depths_this_vcf[si]))
            if sample_total_markers[si] > 0:
                sample_nocall_rates.append(sample_nocalls_this_vcf[si] / sample_total_markers[si])

        if (vi + 1) % 10 == 0 or (vi + 1) == len(vcf_paths):
            log.info(
                "Processed %d/%d VCFs (%d samples, %d het obs)",
                vi + 1,
                len(vcf_paths),
                n_total_samples,
                n_het_total,
            )

    if not marker_gt_counts:
        log.error("no variant observations found")
        return 1

    # Per-marker statistics
    all_markers = sorted(marker_gt_counts.keys())
    n_markers = len(all_markers)

    marker_stats = []
    all_biases = []
    all_abs_biases = []
    all_call_rates = []
    all_nocall_rates = []
    all_mean_depths = []
    all_het_ratios = []  # observed het / expected het under HWE
    all_depth_cvs_marker = []

    for key in all_markers:
        counts = marker_gt_counts[key]
        total = counts["total"]
        n_called = counts["hom_ref"] + counts["het"] + counts["hom_alt"]
        n_nocall = counts["no_call"]

        call_rate = n_called / total if total > 0 else 0.0
        nocall_rate = n_nocall / total if total > 0 else 0.0

        depths = marker_depths.get(key, [])
        mean_depth = sum(depths) / len(depths) if depths else 0.0
        depth_cv = _cv(depths)

        het_vafs = marker_het_vafs.get(key, [])
        n_het = len(het_vafs)
        if n_het > 0:
            devs = [v - EXPECTED_HET_VAF for v in het_vafs]
            median_bias = sorted(devs)[n_het // 2]
            mean_bias = sum(devs) / n_het
            sd_within = _sd(devs)
        else:
            median_bias = mean_bias = sd_within = float("nan")

        # HWE het ratio: compare observed het rate to expected
        # Expected het = 2pq where p = (2*hom_ref + het) / (2*n_called)
        het_ratio = float("nan")
        if n_called >= MIN_CALLED_FOR_HWE:
            p = (2 * counts["hom_ref"] + counts["het"]) / (2 * n_called)
            q = 1 - p
            expected_het = 2 * p * q * n_called
            if expected_het > 0:
                het_ratio = counts["het"] / expected_het

        stat = {
            "total_obs": total,
            "n_called": n_called,
            "n_nocall": n_nocall,
            "call_rate": call_rate,
            "nocall_rate": nocall_rate,
            "n_hom_ref": counts["hom_ref"],
            "n_het": counts["het"],
            "n_hom_alt": counts["hom_alt"],
            "het_ratio_vs_hwe": het_ratio,
            "mean_depth": mean_depth,
            "depth_cv": depth_cv,
            "n_het_for_bias": n_het,
            "median_bias": median_bias,
            "mean_bias": mean_bias,
            "sd_within": sd_within,
        }
        marker_stats.append(stat)

        all_call_rates.append(call_rate)
        all_nocall_rates.append(nocall_rate)
        if depths:
            all_mean_depths.append(mean_depth)
            all_depth_cvs_marker.append(depth_cv)
        if not math.isnan(median_bias):
            all_biases.append(median_bias)
            all_abs_biases.append(abs(median_bias))
        if not math.isnan(het_ratio):
            all_het_ratios.append(het_ratio)

    # Aggregate statistics. Every panel-level number is computed once here and
    # consumed by both the printed summary and the facts CSV below, so the
    # human-readable report and the machine-readable CSV cannot drift apart.
    nocall_sorted = sorted(all_nocall_rates)
    mean_nocall = sum(all_nocall_rates) / len(all_nocall_rates) if all_nocall_rates else 0.0
    median_nocall = _percentile(nocall_sorted, 50) if all_nocall_rates else 0.0
    p95_nocall = _percentile(nocall_sorted, 95) if all_nocall_rates else 0.0
    max_nocall = max(all_nocall_rates) if all_nocall_rates else 0.0
    markers_high_nocall = sum(1 for r in all_nocall_rates if r > HIGH_NOCALL_RATE)

    n_nocall_samples = len(sample_nocall_rates)
    sample_nc_sorted = sorted(sample_nocall_rates)
    mean_sample_nocall = sum(sample_nocall_rates) / n_nocall_samples if sample_nocall_rates else 0.0
    median_sample_nocall = _percentile(sample_nc_sorted, 50) if sample_nocall_rates else 0.0
    p95_sample_nocall = _percentile(sample_nc_sorted, 95) if sample_nocall_rates else 0.0
    max_sample_nocall = max(sample_nocall_rates) if sample_nocall_rates else 0.0

    depths_sorted = sorted(all_mean_depths)
    mean_depth = sum(all_mean_depths) / len(all_mean_depths) if all_mean_depths else 0.0
    median_depth = _percentile(depths_sorted, 50) if all_mean_depths else 0.0
    p5_depth = _percentile(depths_sorted, 5) if all_mean_depths else 0.0
    min_depth = min(all_mean_depths) if all_mean_depths else 0.0
    max_depth = max(all_mean_depths) if all_mean_depths else 0.0
    if all_depth_cvs_marker:
        mean_depth_cv_marker = sum(all_depth_cvs_marker) / len(all_depth_cvs_marker)
    else:
        mean_depth_cv_marker = float("nan")

    scv_sorted = sorted(sample_depth_cvs)
    mean_sample_cv = sum(sample_depth_cvs) / len(sample_depth_cvs) if sample_depth_cvs else 0.0
    median_sample_cv = _percentile(scv_sorted, 50) if sample_depth_cvs else 0.0
    p95_sample_cv = _percentile(scv_sorted, 95) if sample_depth_cvs else 0.0

    n_bias_markers = len(all_biases)
    biases_sorted = sorted(all_biases)
    abs_biases_sorted = sorted(all_abs_biases)
    if n_bias_markers > 0:
        overall_mean_bias = sum(all_biases) / n_bias_markers
        overall_median_bias = biases_sorted[n_bias_markers // 2]
        overall_sd_bias = _sd(all_biases)
        mean_abs_bias = sum(all_abs_biases) / n_bias_markers
        median_abs_bias = _percentile(abs_biases_sorted, 50)
        p95_abs_bias = _percentile(abs_biases_sorted, 95)
        max_abs_bias = max(all_abs_biases)
    else:
        overall_mean_bias = overall_median_bias = overall_sd_bias = float("nan")
        mean_abs_bias = median_abs_bias = p95_abs_bias = max_abs_bias = float("nan")

    n_het_ratio_markers = len(all_het_ratios)
    hr_sorted = sorted(all_het_ratios)
    if all_het_ratios:
        mean_hr = sum(all_het_ratios) / n_het_ratio_markers
        median_hr = _percentile(hr_sorted, 50)
        p5_hr = _percentile(hr_sorted, 5)
        min_hr = min(all_het_ratios)
        markers_low_het = sum(1 for r in all_het_ratios if r < LOW_HET_RATIO_THRESHOLD)
        ado_estimate = max(0.0, 1.0 - mean_hr)
    else:
        mean_hr = median_hr = p5_hr = min_hr = float("nan")
        markers_low_het = 0
        ado_estimate = 0.0

    # Print summary (formatting only; all values come from the block above)
    print()
    print("=" * 65)
    print("PANEL CHARACTERISATION SUMMARY")
    print("=" * 65)

    print("\n--- DATA ---")
    print(f"VCF files processed:         {n_vcfs}")
    print(f"Total samples:               {n_total_samples}")
    print(f"Total markers (biallelic):    {n_markers}")

    print("\n--- LOCUS DROPOUT ---")
    print(f"  Mean no-call rate/marker:  {mean_nocall:.4f} ({mean_nocall * 100:.2f}%)")
    print(f"  Median no-call rate:       {median_nocall:.4f}")
    print(f"  95th pct no-call rate:     {p95_nocall:.4f}")
    print(f"  Max no-call rate:          {max_nocall:.4f}")
    print(f"  Markers with >{HIGH_NOCALL_RATE:.0%} no-call:  {markers_high_nocall}/{n_markers}")

    if sample_nocall_rates:
        print("\n  Per-sample no-call rate:")
        print(f"    Mean:                    {mean_sample_nocall:.4f}")
        print(f"    Median:                  {median_sample_nocall:.4f}")
        print(f"    95th pct:                {p95_sample_nocall:.4f}")
        print(f"    Max:                     {max_sample_nocall:.4f}")

    print("\n--- DEPTH ---")
    if all_mean_depths:
        print(f"  Mean depth/marker:         {mean_depth:.0f}x")
        print(f"  Median depth/marker:       {median_depth:.0f}x")
        print(f"  5th pct depth:             {p5_depth:.0f}x")
        print(f"  Min mean depth:            {min_depth:.0f}x")
        print(f"  Max mean depth:            {max_depth:.0f}x")
    if all_depth_cvs_marker:
        print(f"  Mean depth CV/marker:      {mean_depth_cv_marker:.3f}")
    if sample_depth_cvs:
        print("\n  Depth uniformity (per-sample CV across markers):")
        print(f"    Mean CV:                 {mean_sample_cv:.3f}")
        print(f"    Median CV:               {median_sample_cv:.3f}")
        print(f"    95th pct CV:             {p95_sample_cv:.3f}")

    print("\n--- AMPLIFICATION BIAS ---")
    print(f"  Markers with het obs:      {n_bias_markers}")
    print(f"  Total het observations:    {n_het_total}")
    print(f"  Skipped (low depth):       {n_skipped_dp}")
    print(f"  Skipped (low GQ):          {n_skipped_gq}")
    if n_bias_markers > 0:
        print(f"  Mean bias:                 {overall_mean_bias:+.4f}")
        print(f"  Median bias:               {overall_median_bias:+.4f}")
        print(f"  SD of per-marker biases:   {overall_sd_bias:.4f}")
        print(f"  Mean |bias|:               {mean_abs_bias:.4f}")
        print(f"  Median |bias|:             {median_abs_bias:.4f}")
        print(f"  95th pct |bias|:           {p95_abs_bias:.4f}")
        print(f"  Max |bias|:                {max_abs_bias:.4f}")

    print("\n--- ALLELE DROPOUT SIGNAL (het/hom ratio vs HWE) ---")
    if all_het_ratios:
        print(f"  Mean het ratio:            {mean_hr:.3f}  (1.0 = HWE, <1 = possible ADO)")
        print(f"  Median het ratio:          {median_hr:.3f}")
        print(f"  5th pct het ratio:         {p5_hr:.3f}")
        print(f"  Min het ratio:             {min_hr:.3f}")
        print(
            f"  Markers with ratio < {LOW_HET_RATIO_THRESHOLD}:  {markers_low_het}/{n_het_ratio_markers}"
        )

    print("\n--- SIMULATION PARAMETERS ---")
    if n_bias_markers > 0:
        print(f"  >>> --bias-sd {overall_sd_bias:.3f}")
    print(f"  >>> --locus-dropout-rate {mean_nocall:.4f}")
    if ado_estimate > NEGLIGIBLE_ADO_THRESHOLD:
        print(f"  >>> --allele-dropout-rate ~{ado_estimate:.3f}  (estimated from het deficit)")
    else:
        print(f"  >>> allele dropout: negligible (het ratio ~{mean_hr:.2f})")
    print("=" * 65)

    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        per_marker_path = f"{args.output}_per_marker.tsv"

        with open(per_marker_path, "w", encoding="utf-8") as f:
            header = [
                "marker_index",
                "total_obs",
                "n_called",
                "n_nocall",
                "call_rate",
                "n_hom_ref",
                "n_het",
                "n_hom_alt",
                "het_ratio_vs_hwe",
                "mean_depth",
                "depth_cv",
                "n_het_for_bias",
                "median_bias",
                "mean_bias",
                "sd_within",
            ]
            f.write("\t".join(header) + "\n")
            for i, s in enumerate(marker_stats):
                vals = [
                    str(i),
                    str(s["total_obs"]),
                    str(s["n_called"]),
                    str(s["n_nocall"]),
                    f"{s['call_rate']:.4f}",
                    str(s["n_hom_ref"]),
                    str(s["n_het"]),
                    str(s["n_hom_alt"]),
                    f"{s['het_ratio_vs_hwe']:.4f}"
                    if not math.isnan(s["het_ratio_vs_hwe"])
                    else "NA",
                    f"{s['mean_depth']:.0f}",
                    f"{s['depth_cv']:.4f}",
                    str(s["n_het_for_bias"]),
                    f"{s['median_bias']:.6f}" if not math.isnan(s["median_bias"]) else "NA",
                    f"{s['mean_bias']:.6f}" if not math.isnan(s["mean_bias"]) else "NA",
                    f"{s['sd_within']:.6f}" if not math.isnan(s["sd_within"]) else "NA",
                ]
                f.write("\t".join(vals) + "\n")
        log.info("Per-marker detail: %s", per_marker_path)

        # vibepaper facts CSV (single-row). Values come from the aggregate block
        # above; blank the fields that are undefined when a marker set is empty.
        facts_path = f"{args.output}_facts.csv"
        facts = {
            "n_vcfs": n_vcfs,
            "n_samples": n_total_samples,
            "n_markers": n_markers,
            "n_bias_markers": n_bias_markers,
            "n_het_total": n_het_total,
            "mean_nocall_rate": round(mean_nocall, 4),
            "mean_nocall_pct": round(mean_nocall * 100, 2),
            "markers_gt5pct_nocall": markers_high_nocall,
            "mean_depth": round(mean_depth),
            "mean_sample_depth_cv": round(mean_sample_cv, 3),
            "ado_estimate": round(ado_estimate, 4),
        }

        if all_mean_depths:
            facts["median_depth"] = round(median_depth)
            facts["min_depth"] = round(min_depth)
            facts["max_depth"] = round(max_depth)
        else:
            facts["median_depth"] = facts["min_depth"] = facts["max_depth"] = 0

        if n_bias_markers:
            facts["sd_bias"] = round(overall_sd_bias, 4)
            facts["mean_abs_bias"] = round(mean_abs_bias, 4)
            facts["median_abs_bias"] = round(median_abs_bias, 4)
            facts["p95_abs_bias"] = round(p95_abs_bias, 4)
            facts["max_abs_bias"] = round(max_abs_bias, 4)
        else:
            facts["sd_bias"] = facts["mean_abs_bias"] = facts["median_abs_bias"] = ""
            facts["p95_abs_bias"] = facts["max_abs_bias"] = ""

        if all_het_ratios:
            facts["mean_het_ratio"] = round(mean_hr, 3)
            facts["markers_low_het"] = markers_low_het
        else:
            facts["mean_het_ratio"] = ""
            facts["markers_low_het"] = 0
        # Explicit column order, kept stable for downstream consumers
        # regardless of the order the conditional fields are added above.
        fieldnames = [
            "n_vcfs",
            "n_samples",
            "n_markers",
            "n_bias_markers",
            "n_het_total",
            "mean_nocall_rate",
            "mean_nocall_pct",
            "markers_gt5pct_nocall",
            "mean_depth",
            "median_depth",
            "min_depth",
            "max_depth",
            "mean_sample_depth_cv",
            "sd_bias",
            "mean_abs_bias",
            "median_abs_bias",
            "p95_abs_bias",
            "max_abs_bias",
            "mean_het_ratio",
            "markers_low_het",
            "ado_estimate",
        ]
        with open(facts_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerow(facts)
        log.info("Facts CSV: %s", facts_path)
        log.info("  Copy to output/facts/panel_empirical.csv for vibepaper")

    return 0


if __name__ == "__main__":
    sys.exit(main())
