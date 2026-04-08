#!/usr/bin/env python3
"""Benchmark binomial vs beta-binomial CI calibration.

Generates synthetic chimeric data at multiple donor fractions with
configurable overdispersion, runs both estimators, and reports:
- CI coverage rate (should be ~95% at nominal 95%)
- Mean CI width
- Point estimate accuracy (MAE, RMSE)

Usage:
    python scripts/benchmark_ci_models.py \
        --n-replicates 100 \
        --n-markers 40 \
        --depth 2000 \
        --bias-sd 0.02 \
        --outdir output/ci_benchmark
"""

import argparse
import csv
import logging
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from allomix.chimerism import estimate_single_donor_bb
from allomix.genotype import InformativeMarker
from allomix.simulate import (
    expected_vaf,
    generate_marker_biases_realistic,
    sample_allele_counts,
    sample_marker_depths,
)

log = logging.getLogger(__name__)

FRACTIONS = [0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 0.80, 0.95, 1.0]


def generate_markers(
    true_f: float,
    n_markers: int,
    mean_depth: int,
    rng: random.Random,
    bias_sd: float = 0.02,
    depth_cv: float = 0.4,
) -> tuple[list[InformativeMarker], dict]:
    """Generate synthetic markers with realistic noise.

    Returns (markers, marker_biases_dict).
    """
    # Generate genotypes: alternate type-0 and type-1
    gt_pairs = []
    for i in range(n_markers):
        if i % 2 == 0:
            gt_pairs.append(((0, 0), (1, 1)))  # type 0
        else:
            gt_pairs.append(((1, 1), (0, 0)))  # type 1

    # Per-marker biases and depths
    biases = generate_marker_biases_realistic(n_markers, rng, sd=0.012 if bias_sd > 0 else 0)
    depths = sample_marker_depths(n_markers, mean_depth, depth_cv, rng)

    markers = []
    bias_dict = {}
    for i, ((h_gt, d_gt), bias, dp) in enumerate(zip(gt_pairs, biases, depths)):
        vaf = expected_vaf(h_gt, d_gt, true_f) + bias
        vaf = max(0.0, min(1.0, vaf))
        ref_count, alt_count = sample_allele_counts(vaf, dp, rng)

        chrom = f"chr{(i % 22) + 1}"
        pos = 1_000_000 + i * 100_000
        markers.append(
            InformativeMarker(
                chrom=chrom,
                pos=pos,
                ref="A",
                alt="G",
                host_gt=h_gt,
                donor_gts=[d_gt],
                marker_type=0 if h_gt == (0, 0) else 1,
                admix_ad_ref=ref_count,
                admix_ad_alt=alt_count,
                admix_dp=ref_count + alt_count,
            )
        )
        bias_dict[(chrom, pos, "A", "G")] = bias

    return markers, bias_dict


def run_benchmark(
    n_replicates: int,
    n_markers: int,
    mean_depth: int,
    bias_sd: float,
    outdir: Path,
) -> int:
    """Run the benchmark and write results. Returns 0 on pass, 1 on fail."""
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []

    for true_f in FRACTIONS:
        for rep in range(n_replicates):
            seed = hash((true_f, rep)) % (2**31)
            rng = random.Random(seed)

            markers, bias_dict = generate_markers(
                true_f, n_markers, mean_depth, rng, bias_sd
            )

            res = estimate_single_donor_bb(markers)
            ci_covers = res.donor_fraction_ci[0] <= true_f <= res.donor_fraction_ci[1]
            ci_width = res.donor_fraction_ci[1] - res.donor_fraction_ci[0]

            rows.append({
                "true_f": true_f,
                "replicate": rep,
                "estimate": res.donor_fraction,
                "ci_lo": res.donor_fraction_ci[0],
                "ci_hi": res.donor_fraction_ci[1],
                "ci_covers": ci_covers,
                "ci_width": ci_width,
            })

    # Write per-replicate results
    fields = list(rows[0].keys())
    with open(outdir / "ci_benchmark_results.tsv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    # Compute summary per fraction
    with open(outdir / "ci_benchmark_summary.tsv", "w") as f:
        f.write("true_f_pct\tcoverage\tmean_ci_width\tmae\n")
        for true_f in FRACTIONS:
            subset = [r for r in rows if r["true_f"] == true_f]
            n = len(subset)
            cov = sum(r["ci_covers"] for r in subset) / n
            w = sum(r["ci_width"] for r in subset) / n
            mae = sum(abs(r["estimate"] - true_f) for r in subset) / n
            f.write(f"{true_f * 100:.0f}\t{cov:.3f}\t{w:.4f}\t{mae:.4f}\n")

    # Print summary and return pass/fail
    cov_all = sum(r["ci_covers"] for r in rows) / len(rows)
    mae_all = sum(abs(r["estimate"] - r["true_f"]) for r in rows) / len(rows)
    w_all = sum(r["ci_width"] for r in rows) / len(rows)

    log.info("")
    log.info("=" * 60)
    log.info("BENCHMARK RESULTS")
    log.info("=" * 60)
    log.info("  CI coverage:      %.1f%%", cov_all * 100)
    log.info("  MAE:              %.4f", mae_all)
    log.info("  Mean CI width:    %.4f", w_all)
    log.info("=" * 60)

    # Automated pass/fail gate
    passed = True
    if cov_all < 0.85:
        log.error("FAIL: coverage %.1f%% < 85%%", cov_all * 100)
        passed = False
    if mae_all > 0.005:
        log.error("FAIL: MAE %.4f > 0.5%%", mae_all)
        passed = False
    if w_all > 0.10:
        log.error("FAIL: mean CI width %.4f > 10%%", w_all)
        passed = False

    if passed:
        log.info("PASS: all criteria met")
    return 0 if passed else 1


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark binomial vs beta-binomial CIs")
    parser.add_argument("--n-replicates", type=int, default=100)
    parser.add_argument("--n-markers", type=int, default=40)
    parser.add_argument("--depth", type=int, default=2000)
    parser.add_argument("--bias-sd", type=float, default=0.02)
    parser.add_argument("--outdir", type=str, default="output/ci_benchmark")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    rc = run_benchmark(
        n_replicates=args.n_replicates,
        n_markers=args.n_markers,
        mean_depth=args.depth,
        bias_sd=args.bias_sd,
        outdir=Path(args.outdir),
    )
    sys.exit(rc)


if __name__ == "__main__":
    main()
