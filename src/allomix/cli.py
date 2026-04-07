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
        "--donor", required=True, action="append",
        help="Donor genotype VCF (repeat for multi-donor)",
    )
    parser.add_argument(
        "--sample", required=True, action="append",
        help="Post-HSCT admixture VCF (repeat for multiple timepoints)",
    )
    parser.add_argument("--output", "-o", default="-", help="Output file (default: stdout)")
    parser.add_argument("--min-dp", type=int, default=100, help="Minimum depth (default: 100)")
    parser.add_argument("--min-gq", type=int, default=20, help="Minimum GQ (default: 20)")
    parser.add_argument(
        "--error-rate", type=float, default=0.01,
        help="Sequencing error rate (default: 0.01)",
    )
    parser.add_argument(
        "--format", choices=["tsv", "json"], default="tsv",
        help="Output format (default: tsv)",
    )
    parser.add_argument("--verbose", action="store_true", help="Include per-marker detail")


def _run_single_sample(
    host_path: str,
    donor_paths: list[str],
    sample_path: str,
    min_dp: int,
    min_gq: int,
    error_rate: float,
) -> tuple:
    """Run the chimerism pipeline for one admixture sample.

    Returns (ChimerismResult, QCReport, MarkerGenotypes).
    """
    from allomix.chimerism import estimate_single_donor
    from allomix.genotype import classify_markers, parse_vcf
    from allomix.qc import assess_quality

    host = parse_vcf(host_path, min_gq=min_gq)
    donors = [parse_vcf(d, min_gq=min_gq) for d in donor_paths]
    admix = parse_vcf(sample_path, min_dp=0)  # depth filter applied in classify

    genotypes = classify_markers(host, donors, admix, min_dp=min_dp, min_gq=min_gq)
    genotypes.sample_name = Path(sample_path).stem

    result = estimate_single_donor(genotypes.informative, error_rate=error_rate)
    qc = assess_quality(result, genotypes)

    return result, qc, genotypes


def _open_output(path: str):
    """Open output file or return stdout."""
    if path == "-":
        return sys.stdout
    return open(path, "w")


def cmd_monitor(args: argparse.Namespace) -> int:
    """Run the monitor subcommand."""
    from allomix.report import to_json, to_tsv

    out = _open_output(args.output)
    try:
        for sample_path in args.sample:
            result, qc, genotypes = _run_single_sample(
                args.host, args.donor, sample_path,
                args.min_dp, args.min_gq, args.error_rate,
            )

            if args.format == "json":
                data = to_json(result, qc, sample_name=genotypes.sample_name)
                out.write(json.dumps(data, indent=2) + "\n")
            else:
                to_tsv(result, qc, out, verbose=args.verbose)
    finally:
        if out is not sys.stdout:
            out.close()

    return 0


def cmd_timeline(args: argparse.Namespace) -> int:
    """Run the timeline subcommand."""
    from allomix.report import timeline_json

    results = []
    for sample_path in args.sample:
        result, qc, genotypes = _run_single_sample(
            args.host, args.donor, sample_path,
            args.min_dp, args.min_gq, args.error_rate,
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


def main(argv: list[str] | None = None) -> int:
    """Entry point for the allomix CLI."""
    parser = argparse.ArgumentParser(
        prog="allomix",
        description="NGS-based donor chimerism monitoring for HSCT",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command")

    monitor_parser = subparsers.add_parser(
        "monitor", help="Calculate chimerism for one or more samples",
    )
    _add_common_args(monitor_parser)

    timeline_parser = subparsers.add_parser(
        "timeline", help="Generate chimerism timeline across timepoints",
    )
    _add_common_args(timeline_parser)

    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 1

    if args.command == "monitor":
        return cmd_monitor(args)
    elif args.command == "timeline":
        return cmd_timeline(args)
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
