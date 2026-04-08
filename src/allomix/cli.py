"""Command-line interface for allomix."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from allomix import __version__


def _add_common_args(parser: argparse.ArgumentParser) -> None:
    """Add arguments shared between monitor and timeline."""
    parser.add_argument("--host", required=True, help="Host genotype VCF")
    parser.add_argument(
        "--donor",
        required=True,
        action="append",
        help="Donor genotype VCF (repeat for multi-donor)",
    )
    parser.add_argument(
        "--sample",
        required=True,
        action="append",
        help="Post-HSCT admixture VCF (repeat for multiple timepoints)",
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
    parser.add_argument(
        "--format",
        choices=["tsv", "json"],
        default="tsv",
        help="Output format (default: tsv)",
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


def _run_single_sample(
    host_path: str,
    donor_paths: list[str],
    sample_path: str,
    min_dp: int,
    min_gq: int,
    error_rate: float,
    marker_biases: dict[tuple[str, int, str, str], float] | None = None,
) -> tuple:
    """Run the chimerism pipeline for one admixture sample.

    Automatically uses multi-donor estimation when more than one donor
    VCF is provided.

    Returns (ChimerismResult | MultiDonorResult, QCReport, MarkerGenotypes).
    """
    from allomix.chimerism import estimate_multi_donor, estimate_single_donor
    from allomix.genotype import classify_markers, parse_vcf
    from allomix.qc import assess_quality

    host = parse_vcf(host_path, min_gq=min_gq)
    donors = [parse_vcf(d, min_gq=min_gq) for d in donor_paths]
    admix = parse_vcf(sample_path, min_dp=0)  # depth filter applied in classify

    genotypes = classify_markers(host, donors, admix, min_dp=min_dp, min_gq=min_gq)
    genotypes.sample_name = Path(sample_path).stem

    if len(donor_paths) == 1:
        result = estimate_single_donor(
            genotypes.informative,
            error_rate=error_rate,
            marker_biases=marker_biases,
        )
    else:
        result = estimate_multi_donor(
            genotypes.informative,
            n_donors=len(donor_paths),
            error_rate=error_rate,
            marker_biases=marker_biases,
        )
    qc = assess_quality(result, genotypes)

    return result, qc, genotypes


def _open_output(path: str):
    """Open output file or return stdout."""
    if path == "-":
        return sys.stdout
    return open(path, "w")


def _load_biases(args: argparse.Namespace) -> dict | None:
    """Load bias table if specified and not disabled."""
    if args.bias_table and not args.no_bias_correction:
        from allomix.bias import load_bias_table

        return load_bias_table(args.bias_table)
    return None


def cmd_monitor(args: argparse.Namespace) -> int:
    """Run the monitor subcommand."""
    from allomix.report import to_json, to_tsv

    marker_biases = _load_biases(args)
    out = _open_output(args.output)
    try:
        for sample_path in args.sample:
            result, qc, genotypes = _run_single_sample(
                args.host,
                args.donor,
                sample_path,
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
    from allomix.report import timeline_json

    marker_biases = _load_biases(args)
    results = []
    for sample_path in args.sample:
        result, qc, genotypes = _run_single_sample(
            args.host,
            args.donor,
            sample_path,
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
    import sys

    from allomix.bias import estimate_biases, save_bias_table
    from allomix.genotype import parse_vcf

    marker_lists = []
    for vcf_path in args.vcfs:
        markers = parse_vcf(vcf_path, min_dp=0, min_gq=0)
        marker_lists.append(markers)

    biases = estimate_biases(marker_lists, min_het=args.min_het)

    save_bias_table(biases, args.output)
    print(
        f"Estimated bias for {len(biases)} markers from {len(args.vcfs)} VCFs -> {args.output}",
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

    timeline_parser = subparsers.add_parser(
        "timeline",
        help="Generate chimerism timeline across timepoints",
    )
    _add_common_args(timeline_parser)

    bias_parser = subparsers.add_parser(
        "estimate-bias",
        help="Estimate per-marker amplification bias from VCFs",
    )
    bias_parser.add_argument(
        "--vcfs",
        required=True,
        nargs="+",
        help="Genotyping VCFs to estimate bias from (het markers used)",
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
    elif args.command == "timeline":
        return cmd_timeline(args)
    elif args.command == "estimate-bias":
        return cmd_estimate_bias(args)
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
