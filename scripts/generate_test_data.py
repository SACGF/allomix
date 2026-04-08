#!/usr/bin/env python3
"""Generate synthetic chimeric VCFs at multiple mixture fractions.

Produces a directory of chimeric VCFs and a truth table (TSV) for validating
allomix chimerism calculations against known ground truth.

Usage:
    python scripts/generate_test_data.py \\
        --host data/host.vcf \\
        --donor data/donor.vcf \\
        --outdir test_output/simulated \\
        --depth 2000 \\
        --seed 42
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Allow running from the repo root without installing
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from script_utils import write_truth_table  # noqa: E402

from allomix.simulate import blend_vcfs, write_vcf  # noqa: E402

DEFAULT_FRACTIONS = [
    0.0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 0.80, 0.95, 0.99, 1.0,
]

SEED_HASH_MODULUS = 2**31


def _fraction_to_filename(donor_fraction: float) -> str:
    """Convert a donor fraction to a descriptive filename.

    Names show host_X_donor_Y where X and Y are integer percentages.
    E.g. donor_fraction=0.10 -> 'host_90_donor_10'
    """
    donor_pct = round(donor_fraction * 100)
    host_pct = 100 - donor_pct
    return f"host_{host_pct}_donor_{donor_pct}"


def main(argv: list[str] | None = None) -> int:
    """Entry point for the test data generator."""
    parser = argparse.ArgumentParser(
        description="Generate synthetic chimeric VCFs at multiple mixture fractions.",
    )
    parser.add_argument(
        "--host", required=True, help="Host genotype VCF (plain text .vcf)",
    )
    parser.add_argument(
        "--donor", required=True, help="Donor genotype VCF (plain text .vcf)",
    )
    parser.add_argument(
        "--outdir", default="output/simulated",
        help="Output directory for simulated VCFs (default: output/simulated)",
    )
    parser.add_argument(
        "--depth", type=int, default=None,
        help="Target depth per marker (default: use actual depths from host VCF)",
    )
    parser.add_argument(
        "--seed", type=int, default=42, help="Random seed (default: 42)",
    )
    parser.add_argument(
        "--fractions", type=float, nargs="+", default=None,
        help="Custom donor fractions (default: 0 0.001 0.005 ... 0.99 1.0)",
    )
    parser.add_argument(
        "--bias-sd", type=float, default=0.0,
        help="Per-marker capture bias SD (0=ideal, 0.02=realistic, default: 0)",
    )

    args = parser.parse_args(argv)

    fractions = args.fractions if args.fractions is not None else DEFAULT_FRACTIONS
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    truth_rows: list[dict[str, str]] = []

    for frac in fractions:
        sample_name = _fraction_to_filename(frac)
        vcf_name = f"{sample_name}.vcf"

        print(f"Generating {vcf_name} (donor fraction = {frac:.4f}) ...", file=sys.stderr)

        result = blend_vcfs(
            host_path=args.host,
            donor_path=args.donor,
            donor_fraction=frac,
            target_depth=args.depth,
            sample_name=sample_name,
            marker_bias_sd=args.bias_sd,
            seed=args.seed + hash(str(frac)) % SEED_HASH_MODULUS,
        )

        write_vcf(result, outdir / vcf_name)

        truth_rows.append({
            "sample_name": sample_name,
            "true_donor_fraction": f"{frac:.6f}",
            "num_markers": str(result.num_markers),
            "num_informative": str(result.num_informative),
        })

    # Write truth table
    truth_path = outdir / "truth_table.tsv"
    write_truth_table(
        truth_rows,
        truth_path,
        fieldnames=["sample_name", "true_donor_fraction", "num_markers", "num_informative"],
    )

    print(f"\nGenerated {len(fractions)} VCFs in {outdir}/", file=sys.stderr)
    print(f"Truth table: {truth_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
