#!/usr/bin/env python3
"""Diagnose a single admixture sample: per-marker residuals and noise model.

Internal SA Path diagnostic (not part of the allomix package). Built to chase
down two questions:

1. Why a sample fails goodness-of-fit (e.g. a CNV/LOH makes a block of markers
   disagree with a single contamination fraction). The per-marker standardized
   residuals, summarised per chromosome, localise the offending region.
2. Why a per-sample LOD is higher than a clean simulation predicts. The fitted
   beta-binomial overdispersion (rho) inflates the per-marker variance and caps
   the benefit of depth, so the script prints rho, SE(f=0), LoB and LoD and the
   implied effective-depth ceiling (rho + 1).

Data-access note (see CLAUDE.md): the full per-marker table, which contains
genomic coordinates, is written to a LOCAL file for your own CNV cross-check.
Only de-identified aggregates (fit parameters and per-chromosome counts) are
printed to stdout, so the console output is safe to share.

Usage:
    python scripts/diagnose_sample.py joint_called.vcf.gz \
        --host-sample HOST --donor-sample DONOR --sample BCOL \
        --bias-table output/haem_bias_table_donor.tsv \
        --out output/bcol_per_marker.tsv

Pass the same --bias-table / --min-dp / --min-gq / --error-rate you used for the
run you are diagnosing, or the numbers will not match that run.
"""

import argparse
import math
import statistics
from collections import defaultdict
from pathlib import Path

from allomix.bias import load_bias_table
from allomix.chimerism import (
    PanelCalibration,
    detection_limit,
    estimate_single_donor_bb,
    fraction_se,
)
from allomix.genotype import classify_markers, parse_vcf
from allomix.qc import _error_adjusted_p_alt, assess_quality


def _standardized_residual(
    residual: float, exp_vaf: float, dp: int, rho: float, error_rate: float
) -> float:
    """Residual in units of its own SD, matching the QC goodness-of-fit model.

    Uses the same beta-binomial variance the gof uses,
    ``var = ev(1-ev)(n+rho)/(n(rho+1))``, with ``ev`` the error-adjusted expected
    ALT fraction. The error adjustment is what stops homozygous markers (expected
    VAF 0 or 1) from getting a near-zero variance floor and blowing up to spurious
    outliers. ``z**2`` is then this marker's contribution to the gof chi-square.

    Args:
        residual: Observed minus expected ALT fraction (``MarkerResult.residual``).
        exp_vaf: Model-expected ALT fraction at this marker.
        dp: Marker depth.
        rho: Fitted beta-binomial concentration (inf = pure binomial).
        error_rate: Sequencing error rate used by the fit.

    Returns:
        Standardized residual (z). 0.0 if the variance is undefined.
    """
    if dp <= 0:
        return 0.0
    ev_raw = _error_adjusted_p_alt(exp_vaf, error_rate) if error_rate > 0 else exp_vaf
    ev = min(max(ev_raw, 1e-6), 1.0 - 1e-6)
    if math.isinf(rho):
        var = ev * (1.0 - ev) / dp
    else:
        var = ev * (1.0 - ev) * (dp + rho) / (dp * (rho + 1.0))
    if var <= 0.0:
        return 0.0
    return residual / math.sqrt(var)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("vcf", help="Joint-called VCF containing host, donor and sample")
    parser.add_argument("--host-sample", required=True, help="Host sample name in the VCF")
    parser.add_argument(
        "--donor-sample",
        required=True,
        action="append",
        help="Donor sample name in the VCF (repeat for multi-donor)",
    )
    parser.add_argument("--sample", required=True, help="Admixture sample to diagnose")
    parser.add_argument("--bias-table", default=None, help="Per-marker bias table TSV")
    parser.add_argument("--min-dp", type=int, default=20, help="Minimum depth (default 20)")
    parser.add_argument("--min-gq", type=int, default=20, help="Minimum GQ (default 20)")
    parser.add_argument(
        "--error-rate", type=float, default=0.01, help="Sequencing error rate (default 0.01)"
    )
    parser.add_argument(
        "--outlier-z", type=float, default=3.0, help="|z| above which a marker is an outlier"
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("output/diagnose_per_marker.tsv"),
        help="Local per-marker TSV (contains coordinates; kept on this machine)",
    )
    args = parser.parse_args()

    biases = load_bias_table(args.bias_table) if args.bias_table else None

    host = parse_vcf(args.vcf, sample=args.host_sample, min_gq=args.min_gq)
    donors = [parse_vcf(args.vcf, sample=d, min_gq=args.min_gq) for d in args.donor_sample]
    admix = parse_vcf(args.vcf, sample=args.sample, min_dp=0)
    genotypes = classify_markers(host, donors, admix, min_dp=args.min_dp, min_gq=args.min_gq)
    genotypes.sample_name = args.sample

    if len(donors) != 1:
        parser.error("This diagnostic handles single-donor samples only.")

    result = estimate_single_donor_bb(
        genotypes.informative,
        error_rate=args.error_rate,
        calibration=PanelCalibration(biases=biases or {}),
    )
    qc = assess_quality(result, genotypes)
    rho = result.rho

    # Per-marker table (local file, with coordinates).
    args.out.parent.mkdir(parents=True, exist_ok=True)
    per_chrom: dict[str, list[tuple[float, float, float]]] = defaultdict(list)
    n_outliers = 0
    with args.out.open("w", encoding="utf-8") as fh:
        fh.write("chrom\tpos\tmarker_type\tdp\texp_vaf\tobs_vaf\tresidual\tz\tincluded\toutlier\n")
        for m in result.per_marker:
            z = _standardized_residual(m.residual, m.expected_vaf, m.dp, rho, args.error_rate)
            is_outlier = abs(z) > args.outlier_z
            n_outliers += is_outlier
            per_chrom[m.chrom].append((z, m.residual, float(m.dp)))
            fh.write(
                f"{m.chrom}\t{m.pos}\t{m.marker_type}\t{m.dp}\t{m.expected_vaf:.4f}\t"
                f"{m.observed_vaf:.4f}\t{m.residual:+.4f}\t{z:+.2f}\t{m.included}\t{is_outlier}\n"
            )

    # De-identified summary to stdout.
    se0 = fraction_se(genotypes.informative, 0.0, args.error_rate, rho, biases)
    lob, lod = detection_limit(genotypes.informative, args.error_rate, rho, biases)
    eff_depth_cap = float("inf") if math.isinf(rho) else rho + 1.0
    gof = qc.goodness_of_fit_pval

    print(f"# diagnose: sample={args.sample}  status={qc.status}")
    print(f"donor_fraction   : {result.donor_fraction * 100:.3f}%")
    print(
        f"95% CI           : {result.donor_fraction_ci[0] * 100:.3f}"
        f" - {result.donor_fraction_ci[1] * 100:.3f}%"
    )
    print(f"n_informative    : {result.n_informative}  (used {result.n_markers_used})")
    print(f"gof_pval         : {gof:.4g}" if gof is not None else "gof_pval         : NA")
    print(f"rho (overdisp.)  : {rho:.1f}" if math.isfinite(rho) else "rho (overdisp.)  : inf")
    print(
        f"eff. depth cap   : {eff_depth_cap:.0f}x  (rho+1; per-marker variance stops"
        " falling with depth beyond this)"
        if math.isfinite(eff_depth_cap)
        else "eff. depth cap   : none (pure binomial)"
    )
    print(f"SE(f=0)          : {se0 * 100:.3f}%")
    print(f"LoB / LoD        : {lob * 100:.3f}% / {lod * 100:.3f}%")
    print(f"outliers (|z|>{args.outlier_z:g}): {n_outliers} of {result.n_informative}")
    print(f"\nper-marker table written to {args.out} (local only)")

    print("\n# per-chromosome residual summary (share-safe)")
    print("chrom\tn\tn_outlier\tmean|z|\tmean_signed_resid\tmean_dp")
    for chrom in sorted(per_chrom):
        zs = [z for z, _r, _d in per_chrom[chrom]]
        resids = [r for _z, r, _d in per_chrom[chrom]]
        deps = [d for _z, _r, d in per_chrom[chrom]]
        n_out = sum(1 for z in zs if abs(z) > args.outlier_z)
        print(
            f"{chrom}\t{len(zs)}\t{n_out}\t{statistics.mean(abs(z) for z in zs):.2f}\t"
            f"{statistics.mean(resids):+.4f}\t{statistics.mean(deps):.0f}"
        )


if __name__ == "__main__":
    main()
