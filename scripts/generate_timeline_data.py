#!/usr/bin/env python3
"""Generate synthetic timeline data simulating a post-HSCT patient.

Creates a series of chimeric VCFs at timepoints representing a realistic
clinical scenario: initial engraftment, stable chimerism, then relapse.

Usage:
    python scripts/generate_timeline_data.py \
        --host tests/test_data/host.vcf \
        --donor tests/test_data/donor.vcf \
        --outdir tests/test_data/timeline

Output:
    tests/test_data/timeline/day030_donor_95.vcf   (early engraftment)
    tests/test_data/timeline/day060_donor_98.vcf   (stable)
    tests/test_data/timeline/day090_donor_99.vcf   (full engraftment)
    tests/test_data/timeline/day120_donor_97.vcf   (stable)
    tests/test_data/timeline/day180_donor_90.vcf   (early relapse signal)
    tests/test_data/timeline/day240_donor_70.vcf   (relapse progressing)
    tests/test_data/timeline/day300_donor_40.vcf   (significant relapse)
    tests/test_data/timeline/truth_table.tsv
"""

import argparse
import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from script_utils import write_truth_table  # noqa: E402

from allomix.simulate import blend_vcfs, write_vcf  # noqa: E402

log = logging.getLogger(__name__)

# (day, donor_fraction); simulates engraftment then relapse
TIMEPOINTS = [
    (30, 0.95),
    (60, 0.98),
    (90, 0.99),
    (120, 0.97),
    (180, 0.90),
    (240, 0.70),
    (300, 0.40),
]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate synthetic timeline data (engraftment + relapse).",
    )
    parser.add_argument("--host", required=True, help="Host genotype VCF")
    parser.add_argument("--donor", required=True, help="Donor genotype VCF")
    parser.add_argument(
        "--outdir",
        default="output/timeline",
        help="Output directory (default: output/timeline)",
    )
    parser.add_argument("--depth", type=int, default=2000, help="Target depth (default: 2000)")
    parser.add_argument("--seed", type=int, default=99, help="Random seed (default: 99)")
    parser.add_argument(
        "--bias-sd",
        type=float,
        default=0.0,
        help="Per-marker capture bias SD (0=ideal, 0.02=realistic, default: 0)",
    )
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    truth_rows: list[dict[str, str]] = []

    for day, donor_frac in TIMEPOINTS:
        donor_pct = round(donor_frac * 100)
        sample_name = f"day{day:03d}_donor_{donor_pct}"
        vcf_name = f"{sample_name}.vcf"

        log.info("Generating %s (day %d, %d%% donor) ...", vcf_name, day, donor_pct)

        result = blend_vcfs(
            host_path=args.host,
            donor_path=args.donor,
            donor_fraction=donor_frac,
            target_depth=args.depth,
            sample_name=sample_name,
            marker_bias_sd=args.bias_sd,
            seed=args.seed + day,
        )
        write_vcf(result, outdir / vcf_name)

        truth_rows.append(
            {
                "sample_name": sample_name,
                "day": str(day),
                "true_donor_fraction": f"{donor_frac:.6f}",
                "num_markers": str(result.num_markers),
                "num_informative": str(result.num_informative),
            }
        )

    truth_path = outdir / "truth_table.tsv"
    write_truth_table(
        truth_rows,
        truth_path,
        fieldnames=[
            "sample_name",
            "day",
            "true_donor_fraction",
            "num_markers",
            "num_informative",
        ],
    )

    log.info("Generated %d timepoint VCFs in %s/", len(TIMEPOINTS), outdir)
    log.info("Truth table: %s", truth_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
