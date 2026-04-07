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
import csv
import sys
from pathlib import Path

# Allow running from the repo root without installing
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from allomix.simulate import blend_vcfs, write_vcf  # noqa: E402

DEFAULT_FRACTIONS = [
    0.0, 0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 0.80, 0.95, 0.99, 1.0,
]


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
        "--outdir", required=True, help="Output directory for simulated VCFs",
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

    args = parser.parse_args(argv)

    fractions = args.fractions if args.fractions is not None else DEFAULT_FRACTIONS
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    truth_rows: list[dict[str, str]] = []

    for frac in fractions:
        pct_label = f"{frac * 100:.1f}".replace(".", "p")
        sample_name = f"sim_donor_{pct_label}pct"
        vcf_name = f"{sample_name}.vcf"

        print(f"Generating {vcf_name} (donor fraction = {frac:.4f}) ...", file=sys.stderr)

        result = blend_vcfs(
            host_path=args.host,
            donor_path=args.donor,
            donor_fraction=frac,
            target_depth=args.depth,
            sample_name=sample_name,
            seed=args.seed + hash(str(frac)) % (2**31),
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
    with open(truth_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["sample_name", "true_donor_fraction", "num_markers", "num_informative"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(truth_rows)

    print(f"\nGenerated {len(fractions)} VCFs in {outdir}/", file=sys.stderr)
    print(f"Truth table: {truth_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
