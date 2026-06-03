#!/usr/bin/env python3
"""Generate multi-donor test data: 3 brothers (host + 2 donors).

Creates:
    - host.vcf, donor1.vcf, donor2.vcf (sibling genotype VCFs)
    - Chimeric VCFs at a grid of (f1, f2) mixture fractions
    - truth_table.tsv with ground truth fractions

Usage:
    python scripts/generate_multidonor_test_data.py --outdir tests/test_data/multidonor
"""

import argparse
import logging
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from script_utils import write_truth_table  # noqa: E402

from allomix.simulate import (  # noqa: E402
    blend_from_genotype_dicts,
    generate_sibling_trio_genotypes,
    write_genotype_vcf,
    write_vcf,
)

log = logging.getLogger(__name__)

# Simulation parameters matching the allomix MLE error model
DEFAULT_ERROR_RATE = 0.01  # per-read sequencing error probability
DEFAULT_DEPTH_CV = 0.43  # empirical depth CV from rhAmpSeq panel characterisation
SEED_HASH_MODULUS = 2**31

# Mixture fraction grid points on the simplex f1 + f2 <= 1
MIXTURE_GRID = [
    # Pure host
    (0.00, 0.00),
    # Single-donor edges
    (0.05, 0.00),
    (0.20, 0.00),
    (0.50, 0.00),
    (0.00, 0.05),
    (0.00, 0.20),
    (0.00, 0.50),
    # Balanced two-donor
    (0.10, 0.10),
    (0.25, 0.25),
    (0.40, 0.40),
    # Asymmetric two-donor
    (0.30, 0.10),
    (0.10, 0.30),
    (0.50, 0.20),
    (0.20, 0.50),
    (0.05, 0.15),
    (0.15, 0.05),
    # High total donor
    (0.45, 0.45),
    (0.60, 0.30),
    # Low-fraction detection
    (0.02, 0.02),
    (0.01, 0.05),
    # Pure donors
    (1.00, 0.00),
    (0.00, 1.00),
]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate multi-donor test data (3 siblings).",
    )
    parser.add_argument(
        "--outdir",
        default="tests/test_data/multidonor",
        help="Output directory (default: tests/test_data/multidonor)",
    )
    parser.add_argument("--n-markers", type=int, default=100, help="Number of SNPs")
    parser.add_argument("--depth", type=int, default=1000, help="Mean sequencing depth")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    # Generate 3-sibling genotypes
    markers = generate_sibling_trio_genotypes(args.n_markers, rng)

    # Write genotype VCFs
    write_genotype_vcf(markers, outdir / "host.vcf", "HOST", key="host_gt")
    write_genotype_vcf(markers, outdir / "donor1.vcf", "DONOR1", key="donor1_gt")
    write_genotype_vcf(markers, outdir / "donor2.vcf", "DONOR2", key="donor2_gt")

    # Report informativity
    n_inf_d1 = sum(1 for m in markers if m["informative_d1"])
    n_inf_d2 = sum(1 for m in markers if m["informative_d2"])
    n_inf_any = sum(1 for m in markers if m["informative_any"])
    n_distinguishable = sum(1 for m in markers if m["donors_distinguishable"])
    log.info("Markers: %d", args.n_markers)
    log.info("  Informative for donor1: %d", n_inf_d1)
    log.info("  Informative for donor2: %d", n_inf_d2)
    log.info("  Informative for any donor: %d", n_inf_any)
    log.info("  Donors distinguishable: %d", n_distinguishable)

    # Generate chimeric VCFs at grid points
    truth_rows = []
    for f1, f2 in MIXTURE_GRID:
        name = f"host_{100 - round((f1 + f2) * 100)}_d1_{round(f1 * 100)}_d2_{round(f2 * 100)}"
        result = blend_from_genotype_dicts(
            markers,
            [f1, f2],
            target_depth=args.depth,
            seed=rng.randint(0, SEED_HASH_MODULUS),
            error_rate=DEFAULT_ERROR_RATE,
            depth_cv=DEFAULT_DEPTH_CV,
            sample_name=name,
        )
        write_vcf(result, outdir / f"{name}.vcf")
        truth_rows.append({
            "sample_name": name,
            "true_donor1_fraction": f"{f1:.4f}",
            "true_donor2_fraction": f"{f2:.4f}",
            "true_host_fraction": f"{1.0 - f1 - f2:.4f}",
            "num_markers": result.num_markers,
            "num_informative_any": result.num_informative,
        })

    # Write truth table
    write_truth_table(truth_rows, outdir / "truth_table.tsv")

    log.info("Generated %d chimeric VCFs in %s/", len(MIXTURE_GRID), outdir)

    return 0


if __name__ == "__main__":
    sys.exit(main())
