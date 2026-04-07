"""Command-line interface for allomix."""

import argparse
import sys

from allomix import __version__


def main(argv: list[str] | None = None) -> int:
    """Entry point for the allomix CLI."""
    parser = argparse.ArgumentParser(
        prog="allomix",
        description="NGS-based donor chimerism monitoring for HSCT",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command")

    # monitor subcommand
    monitor_parser = subparsers.add_parser(
        "monitor", help="Calculate chimerism for one or more samples"
    )
    monitor_parser.add_argument(
        "--host", required=True, help="Host genotype VCF (.vcf.gz)"
    )
    monitor_parser.add_argument(
        "--donor", required=True, action="append",
        help="Donor genotype VCF (repeat for multi-donor)",
    )
    monitor_parser.add_argument(
        "--sample", required=True, action="append",
        help="Post-HSCT admixture VCF (repeat for multiple timepoints)",
    )
    monitor_parser.add_argument(
        "--output", "-o", default="-", help="Output file (default: stdout)"
    )

    # timeline subcommand
    timeline_parser = subparsers.add_parser(
        "timeline", help="Generate chimerism timeline across timepoints"
    )
    timeline_parser.add_argument(
        "--host", required=True, help="Host genotype VCF (.vcf.gz)"
    )
    timeline_parser.add_argument(
        "--donor", required=True, action="append",
        help="Donor genotype VCF (repeat for multi-donor)",
    )
    timeline_parser.add_argument(
        "--sample", required=True, action="append",
        help="Post-HSCT admixture VCFs in chronological order",
    )
    timeline_parser.add_argument(
        "--output", "-o", default="-", help="Output file (default: stdout)"
    )

    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 1

    # Placeholder — actual implementation will come in later steps
    print(f"allomix {args.command}: not yet implemented", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main())
