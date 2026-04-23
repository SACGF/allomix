"""Command-line interface for allomix."""

from __future__ import annotations

import argparse
import json
import sys

from allomix import __version__
from allomix.bias import estimate_biases, load_bias_table, save_bias_table
from allomix.chimerism import estimate_multi_donor, estimate_single_donor_bb
from allomix.genotype import classify_markers, parse_vcf
from allomix.qc import assess_quality
from allomix.report import timeline_json, to_json, to_tsv


def _add_common_args(parser: argparse.ArgumentParser) -> None:
    """Add arguments shared between monitor and timeline."""
    parser.add_argument("--vcf", required=True, help="Joint-called VCF containing all samples")
    parser.add_argument("--host-sample", required=True, help="Host sample name in VCF")
    parser.add_argument(
        "--donor-sample",
        required=True,
        action="append",
        metavar="SAMPLE_NAME",
        help="Donor sample name in VCF (repeat for multi-donor)",
    )
    parser.add_argument(
        "--sample",
        required=True,
        action="append",
        metavar="SAMPLE_NAME",
        help="Admixture sample name in VCF (repeat for multiple timepoints)",
    )
    parser.add_argument("--output", "-o", default="-", help="Output file (default: stdout)")
    parser.add_argument("--min-dp", type=int, default=100, help="Minimum depth (default: 100)")
    parser.add_argument("--min-gq", type=int, default=20, help="Minimum GQ (default: 20)")
    parser.add_argument(
        "--error-rate",
        type=float,
        default=0.01,
        help="Sequencing error rate (default: 0.01)",
    )
    parser.add_argument("--verbose", action="store_true", help="Include per-marker detail")
    parser.add_argument(
        "--bias-table",
        default=None,
        help="Per-marker bias table TSV (from allomix estimate-bias or simulation)",
    )
    parser.add_argument(
        "--no-bias-correction",
        action="store_true",
        help="Disable bias correction even when a bias table is provided",
    )


def _validate_sample_names(vcf_path: str, required: list[str]) -> None:
    """Check that all required sample names exist in the VCF header.

    Raises SystemExit with a clear error listing available names if any
    are missing.
    """
    from cyvcf2 import VCF

    vcf = VCF(vcf_path)
    available = list(vcf.samples)
    vcf.close()
    missing = [s for s in required if s not in available]
    if missing:
        raise SystemExit(f"Sample(s) not found in {vcf_path}: {missing}\nAvailable: {available}")


def _run_single_sample(
    host: list,
    donors: list[list],
    vcf_path: str,
    admix_sample: str,
    min_dp: int,
    min_gq: int,
    error_rate: float,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
) -> tuple:
    """Run the chimerism pipeline for one admixture sample.

    Takes pre-parsed host and donor markers to avoid redundant VCF reads.
    Automatically uses multi-donor estimation when more than one donor
    is provided.

    Returns (ChimerismResult | MultiDonorResult, QCReport, MarkerGenotypes).
    """
    admix = parse_vcf(vcf_path, sample=admix_sample, min_dp=0)

    genotypes = classify_markers(host, donors, admix, min_dp=min_dp, min_gq=min_gq)
    genotypes.sample_name = admix_sample

    if len(donors) == 1:
        result = estimate_single_donor_bb(
            genotypes.informative,
            error_rate=error_rate,
            marker_biases=marker_biases,
        )
    else:
        result = estimate_multi_donor(
            genotypes.informative,
            n_donors=len(donors),
            error_rate=error_rate,
            marker_biases=marker_biases,
        )
    qc = assess_quality(result, genotypes)

    return result, qc, genotypes


def _open_output(path: str):
    """Open output file or return stdout."""
    if path == "-":
        return sys.stdout
    return open(path, "w", encoding="utf-8")


def _load_biases(args: argparse.Namespace) -> dict | None:
    """Load bias table if specified and not disabled."""
    if args.bias_table and not args.no_bias_correction:
        return load_bias_table(args.bias_table)
    return None


def cmd_monitor(args: argparse.Namespace) -> int:
    """Run the monitor subcommand."""
    all_names = [args.host_sample] + args.donor_sample + args.sample
    _validate_sample_names(args.vcf, all_names)

    marker_biases = _load_biases(args)

    # Parse host and donors once — they're the same for every timepoint
    host = parse_vcf(args.vcf, sample=args.host_sample, min_gq=args.min_gq)
    donors = [parse_vcf(args.vcf, sample=d, min_gq=args.min_gq) for d in args.donor_sample]

    out = _open_output(args.output)
    try:
        for sample_name in args.sample:
            result, qc, genotypes = _run_single_sample(
                host,
                donors,
                args.vcf,
                sample_name,
                args.min_dp,
                args.min_gq,
                args.error_rate,
                marker_biases=marker_biases,
            )

            if args.format == "json":
                data = to_json(result, qc, sample_name=genotypes.sample_name)
                out.write(json.dumps(data, indent=2) + "\n")
            else:
                to_tsv(result, qc, out, verbose=args.verbose, sample_name=genotypes.sample_name)
    finally:
        if out is not sys.stdout:
            out.close()

    return 0


def cmd_timeline(args: argparse.Namespace) -> int:
    """Run the timeline subcommand."""
    all_names = [args.host_sample] + args.donor_sample + args.sample
    _validate_sample_names(args.vcf, all_names)

    marker_biases = _load_biases(args)

    # Parse host and donors once
    host = parse_vcf(args.vcf, sample=args.host_sample, min_gq=args.min_gq)
    donors = [parse_vcf(args.vcf, sample=d, min_gq=args.min_gq) for d in args.donor_sample]

    results = []
    for sample_name in args.sample:
        result, qc, genotypes = _run_single_sample(
            host,
            donors,
            args.vcf,
            sample_name,
            args.min_dp,
            args.min_gq,
            args.error_rate,
            marker_biases=marker_biases,
        )
        results.append((genotypes.sample_name, result, qc))

    data = timeline_json(results)

    out = _open_output(args.output)
    try:
        out.write(json.dumps(data, indent=2) + "\n")
    finally:
        if out is not sys.stdout:
            out.close()

    return 0


def cmd_estimate_bias(args: argparse.Namespace) -> int:
    """Run the estimate-bias subcommand."""
    if args.vcfs and args.vcf:
        raise SystemExit("Use either --vcfs or --vcf/--samples, not both")
    if not args.vcfs and not args.vcf:
        raise SystemExit("One of --vcfs or --vcf is required")
    if args.vcf and not args.samples:
        raise SystemExit("--samples is required when using --vcf")

    marker_lists = []
    if args.vcfs:
        for vcf_path in args.vcfs:
            markers = parse_vcf(vcf_path, min_dp=0, min_gq=0)
            marker_lists.append(markers)
        n_source = f"{len(args.vcfs)} VCFs"
    else:
        _validate_sample_names(args.vcf, args.samples)
        for sample in args.samples:
            markers = parse_vcf(args.vcf, sample=sample, min_dp=0, min_gq=0)
            marker_lists.append(markers)
        n_source = f"{len(args.samples)} samples from {args.vcf}"

    biases = estimate_biases(marker_lists, min_het=args.min_het)
    save_bias_table(biases, args.output)
    print(
        f"Estimated bias for {len(biases)} markers from {n_source} -> {args.output}",
        file=sys.stderr,
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    """Entry point for the allomix CLI."""
    parser = argparse.ArgumentParser(
        prog="allomix",
        description="NGS-based donor chimerism monitoring for HSCT",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command")

    monitor_parser = subparsers.add_parser(
        "monitor",
        help="Calculate chimerism for one or more samples",
    )
    _add_common_args(monitor_parser)
    monitor_parser.add_argument(
        "--format",
        choices=["tsv", "json"],
        default="tsv",
        help="Output format (default: tsv)",
    )

    timeline_parser = subparsers.add_parser(
        "timeline",
        help="Generate chimerism timeline across timepoints",
    )
    _add_common_args(timeline_parser)

    bias_parser = subparsers.add_parser(
        "estimate-bias",
        help="Estimate per-marker amplification bias from VCFs",
    )
    bias_input = bias_parser.add_mutually_exclusive_group()
    bias_input.add_argument(
        "--vcfs",
        nargs="+",
        metavar="VCF",
        help="Per-sample VCFs, one per file (reads first sample from each)",
    )
    bias_input.add_argument(
        "--vcf",
        metavar="VCF",
        help="Joint-called multi-sample VCF (use with --samples)",
    )
    bias_parser.add_argument(
        "--samples",
        nargs="+",
        metavar="SAMPLE_NAME",
        help="Sample names to extract from --vcf",
    )
    bias_parser.add_argument(
        "--output",
        "-o",
        default="bias_table.tsv",
        help="Output bias table TSV (default: bias_table.tsv)",
    )
    bias_parser.add_argument(
        "--min-het",
        type=int,
        default=1,
        help="Minimum het observations per marker (default: 1)",
    )

    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 1

    if args.command == "monitor":
        return cmd_monitor(args)
    if args.command == "timeline":
        return cmd_timeline(args)
    if args.command == "estimate-bias":
        return cmd_estimate_bias(args)
    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
