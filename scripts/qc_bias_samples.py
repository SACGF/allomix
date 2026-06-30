#!/usr/bin/env python3
"""Sample-level QC for bias-estimation training samples.

Reads a joint-called multi-sample VCF (typically restricted to the SID bed),
computes per-sample quality metrics, and outputs a list of samples that pass
QC (for use with ``allomix estimate-bias --samples ...``) plus an optional
per-sample metrics TSV.

Metrics per sample:
  - no-call rate: fraction of sites called ./. Catches panel-version
    mismatches where SID sites are absent from the capture kit.
  - het rate: fraction of called sites that are heterozygous. Unrelated
    individuals at common SNPs should land ~0.3-0.5. Far outside suggests
    a problem.
  - mean |VAF - 0.5| at het calls: allele-balance skew proxy. Samples
    with heavy CNV/LOH show systematic deviation.

Usage:
    python scripts/qc_bias_samples.py joint.vcf.gz \\
        --output-samples pass_samples.txt \\
        --output-metrics qc_metrics.tsv
"""

import argparse
import logging
import math
import sys
from collections import Counter
from pathlib import Path

from cyvcf2 import VCF

log = logging.getLogger(__name__)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Sample-level QC for bias-training joint VCFs.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("vcf", help="Joint-called multi-sample VCF")
    p.add_argument(
        "--output-samples",
        "-o",
        required=True,
        help="Output text file with pass-QC sample names, one per line",
    )
    p.add_argument(
        "--output-metrics",
        default=None,
        help="Optional TSV with per-sample metrics (all samples, pass and fail)",
    )
    p.add_argument(
        "--min-dp",
        type=int,
        default=100,
        help="Minimum depth (AD sum) for a call to count toward metrics (default: 100)",
    )
    p.add_argument(
        "--max-nocall-rate",
        type=float,
        default=0.10,
        help="Exclude samples with no-call rate above this (default: 0.10)",
    )
    p.add_argument(
        "--min-het-rate",
        type=float,
        default=0.15,
        help="Exclude samples with het rate below this (default: 0.15)",
    )
    p.add_argument(
        "--max-het-rate",
        type=float,
        default=0.60,
        help="Exclude samples with het rate above this (default: 0.60)",
    )
    p.add_argument(
        "--max-mean-vaf-dev",
        type=float,
        default=0.15,
        help="Exclude samples where mean |VAF-0.5| at het calls exceeds this (default: 0.15)",
    )
    p.add_argument(
        "--min-hets-for-vaf-dev",
        type=int,
        default=10,
        help=(
            "Minimum het calls required to evaluate the VAF-deviation check "
            "(default: 10). Samples with fewer hets skip this check."
        ),
    )
    return p.parse_args(argv)


def _metric_row(
    sample: str,
    n_sites: int,
    n_nocall: int,
    n_low_dp: int,
    n_hom_ref: int,
    n_het: int,
    n_hom_alt: int,
    sum_abs_dev: float,
    args: argparse.Namespace,
) -> dict:
    """Compute metrics and pass/fail decision for one sample."""
    called = n_hom_ref + n_het + n_hom_alt
    nocall_rate = n_nocall / n_sites if n_sites else 0.0
    het_rate = n_het / called if called else 0.0
    mean_vaf_dev = (sum_abs_dev / n_het) if n_het else float("nan")

    reasons: list[str] = []
    if nocall_rate > args.max_nocall_rate:
        reasons.append(f"nocall_rate={nocall_rate:.3f}>{args.max_nocall_rate}")
    if called > 0:
        if het_rate < args.min_het_rate:
            reasons.append(f"het_rate={het_rate:.3f}<{args.min_het_rate}")
        elif het_rate > args.max_het_rate:
            reasons.append(f"het_rate={het_rate:.3f}>{args.max_het_rate}")
    if n_het >= args.min_hets_for_vaf_dev and mean_vaf_dev > args.max_mean_vaf_dev:
        reasons.append(f"mean_vaf_dev={mean_vaf_dev:.3f}>{args.max_mean_vaf_dev}")

    return {
        "sample": sample,
        "n_sites": n_sites,
        "n_called": called,
        "n_nocall": n_nocall,
        "n_low_dp": n_low_dp,
        "n_hom_ref": n_hom_ref,
        "n_het": n_het,
        "n_hom_alt": n_hom_alt,
        "nocall_rate": nocall_rate,
        "het_rate": het_rate,
        "mean_vaf_dev": mean_vaf_dev,
        "pass": len(reasons) == 0,
        "reasons": ";".join(reasons),
    }


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    vcf = VCF(args.vcf)
    samples = list(vcf.samples)
    n_samples = len(samples)
    if n_samples == 0:
        log.error("no samples in VCF")
        return 1

    log.info("Scanning %s for %d samples", args.vcf, n_samples)

    n_sites = [0] * n_samples
    n_nocall = [0] * n_samples
    n_low_dp = [0] * n_samples
    n_hom_ref = [0] * n_samples
    n_het = [0] * n_samples
    n_hom_alt = [0] * n_samples
    sum_abs_dev = [0.0] * n_samples

    n_variants_read = 0
    for v in vcf:
        # Keep biallelic SNPs only, matching estimate-bias semantics
        if len(v.ALT) != 1:
            continue
        n_variants_read += 1
        genotypes = v.genotypes
        ad = v.format("AD")

        for si in range(n_samples):
            n_sites[si] += 1
            a0, a1 = genotypes[si][0], genotypes[si][1]
            if a0 < 0 or a1 < 0:
                n_nocall[si] += 1
                continue

            if ad is None:
                n_low_dp[si] += 1
                continue
            try:
                ref_c = int(ad[si][0])
                alt_c = int(ad[si][1])
            except (IndexError, TypeError, ValueError):
                n_low_dp[si] += 1
                continue
            dp = ref_c + alt_c
            if dp < args.min_dp:
                n_low_dp[si] += 1
                continue

            if a0 != a1:
                n_het[si] += 1
                sum_abs_dev[si] += abs(alt_c / dp - 0.5)
            elif a0 == 0:
                n_hom_ref[si] += 1
            else:
                n_hom_alt[si] += 1

    vcf.close()
    log.info("Read %d biallelic SNP sites", n_variants_read)
    if n_variants_read == 0:
        log.error("no biallelic SNP sites in VCF")
        return 1

    rows = [
        _metric_row(
            samples[si],
            n_sites[si],
            n_nocall[si],
            n_low_dp[si],
            n_hom_ref[si],
            n_het[si],
            n_hom_alt[si],
            sum_abs_dev[si],
            args,
        )
        for si in range(n_samples)
    ]
    pass_samples = [r["sample"] for r in rows if r["pass"]]

    out_samples = Path(args.output_samples)
    out_samples.parent.mkdir(parents=True, exist_ok=True)
    with open(out_samples, "w", encoding="utf-8") as f:
        for name in pass_samples:
            f.write(name + "\n")

    if args.output_metrics:
        out_metrics = Path(args.output_metrics)
        out_metrics.parent.mkdir(parents=True, exist_ok=True)
        cols = [
            "sample",
            "n_sites",
            "n_called",
            "n_nocall",
            "n_low_dp",
            "n_hom_ref",
            "n_het",
            "n_hom_alt",
            "nocall_rate",
            "het_rate",
            "mean_vaf_dev",
            "pass",
            "reasons",
        ]
        with open(out_metrics, "w", encoding="utf-8") as f:
            f.write("\t".join(cols) + "\n")
            for r in rows:
                vaf_dev = f"{r['mean_vaf_dev']:.4f}" if not math.isnan(r["mean_vaf_dev"]) else "NA"
                vals = [
                    r["sample"],
                    str(r["n_sites"]),
                    str(r["n_called"]),
                    str(r["n_nocall"]),
                    str(r["n_low_dp"]),
                    str(r["n_hom_ref"]),
                    str(r["n_het"]),
                    str(r["n_hom_alt"]),
                    f"{r['nocall_rate']:.4f}",
                    f"{r['het_rate']:.4f}",
                    vaf_dev,
                    "1" if r["pass"] else "0",
                    r["reasons"],
                ]
                f.write("\t".join(vals) + "\n")
        log.info("Metrics: %s", out_metrics)

    log.info("Pass list: %s", out_samples)
    log.info("QC complete: %d / %d samples pass", len(pass_samples), n_samples)

    fail_rows = [r for r in rows if not r["pass"]]
    if fail_rows:
        reason_counts: Counter[str] = Counter()
        for r in fail_rows:
            for reason in r["reasons"].split(";"):
                reason_counts[reason.split("=")[0]] += 1
        log.info("Failure breakdown:")
        for tag, n in reason_counts.most_common():
            log.info("  %s: %d", tag, n)

    return 0


if __name__ == "__main__":
    sys.exit(main())
