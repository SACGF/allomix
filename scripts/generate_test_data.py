#!/usr/bin/env python3
"""Generate synthetic joint-called VCFs for testing allomix.

Produces joint VCFs containing host, donor, and admixture samples at multiple
mixture fractions, along with a truth table (TSV) for validation.

Usage:
    python scripts/generate_test_data.py \
        --host data/idt_rhampseq_sid_example.vcf \
        --outdir tests/test_data \
        --depth 2000 \
        --seed 42
"""

from __future__ import annotations

import argparse
import logging
import sys
from pathlib import Path

# Allow running from the repo root without installing
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from script_utils import write_truth_table  # noqa: E402

from allomix.simulate import build_joint_vcf, write_joint_vcf  # noqa: E402

log = logging.getLogger(__name__)

SEED_HASH_MODULUS = 2**31


def _create_flipped_vcf(tmp_dir: Path, source_path: Path) -> Path:
    """Create a VCF where genotypes are flipped (0/0->1/1, 1/1->0/0, 0/1 stays).

    This ensures maximum informativeness between host and donor.
    """
    lines = []
    with open(source_path) as f:
        for line in f:
            if line.startswith("#"):
                lines.append(line)
                continue
            if not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) < 10:
                lines.append(line)
                continue

            fmt_keys = fields[8].split(":")
            fmt_vals = fields[9].split(":")
            gt_idx = fmt_keys.index("GT")
            gt = fmt_vals[gt_idx]

            if gt == "0/0":
                fmt_vals[gt_idx] = "1/1"
                if "AD" in fmt_keys:
                    ad_idx = fmt_keys.index("AD")
                    ad_parts = fmt_vals[ad_idx].split(",")
                    if len(ad_parts) == 2:
                        fmt_vals[ad_idx] = f"{ad_parts[1]},{ad_parts[0]}"
                    elif len(ad_parts) == 1:
                        fmt_vals[ad_idx] = f"0,{ad_parts[0]}"
                if "AF" in fmt_keys:
                    af_idx = fmt_keys.index("AF")
                    fmt_vals[af_idx] = "1.0"
                if fields[4] == ".":
                    fields[4] = "T" if fields[3] != "T" else "A"
            elif gt == "1/1":
                fmt_vals[gt_idx] = "0/0"
                if "AD" in fmt_keys:
                    ad_idx = fmt_keys.index("AD")
                    ad_parts = fmt_vals[ad_idx].split(",")
                    if len(ad_parts) == 2:
                        fmt_vals[ad_idx] = f"{ad_parts[1]},{ad_parts[0]}"
                if "AF" in fmt_keys:
                    af_idx = fmt_keys.index("AF")
                    fmt_vals[af_idx] = "0"

            fields[9] = ":".join(fmt_vals)
            lines.append("\t".join(fields) + "\n")

    out = tmp_dir / "donor_flipped.vcf"
    with open(out, "w") as f:
        f.writelines(lines)
    return out


def main(argv: list[str] | None = None) -> int:
    """Entry point for the test data generator."""
    parser = argparse.ArgumentParser(
        description="Generate synthetic joint-called VCFs for testing allomix.",
    )
    parser.add_argument(
        "--host",
        required=True,
        help="Host genotype VCF (plain text .vcf)",
    )
    parser.add_argument(
        "--outdir",
        default="tests/test_data",
        help="Output directory (default: tests/test_data)",
    )
    parser.add_argument(
        "--depth",
        type=int,
        default=2000,
        help="Target depth per marker (default: 2000)",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed (default: 42)",
    )

    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    host_path = Path(args.host)

    # Create a flipped-genotype donor VCF for maximum informativeness
    donor_path = _create_flipped_vcf(outdir, host_path)

    # --- Single-donor joint VCF ---
    single_fracs = [0.0, 0.01, 0.10, 0.50, 1.0]
    single_names = [f"ADMIX_F{f:.2f}" for f in single_fracs]

    log.info("Building joint_single_donor.vcf with fractions %s", single_fracs)
    single_result = build_joint_vcf(
        host_path=str(host_path),
        donor_paths=[str(donor_path)],
        admix_fractions=single_fracs,
        admix_sample_names=single_names,
        host_sample_name="HOST",
        donor_sample_names=["DONOR"],
        target_depth=args.depth,
        seed=args.seed,
        error_rate=0.01,
    )
    write_joint_vcf(single_result, outdir / "joint_single_donor.vcf")

    # --- Multi-donor joint VCF (use same host, two copies of flipped donor) ---
    multi_fracs = [0.0, 0.10, 0.50]
    multi_names = [f"ADMIX_F{f:.2f}" for f in multi_fracs]

    log.info("Building joint_multi_donor.vcf with fractions %s", multi_fracs)
    multi_result = build_joint_vcf(
        host_path=str(host_path),
        donor_paths=[str(donor_path), str(donor_path)],
        admix_fractions=multi_fracs,
        admix_sample_names=multi_names,
        host_sample_name="HOST",
        donor_sample_names=["DONOR1", "DONOR2"],
        target_depth=args.depth,
        seed=args.seed + 1000,
        error_rate=0.01,
    )
    write_joint_vcf(multi_result, outdir / "joint_multi_donor.vcf")

    # --- Truth table ---
    truth_rows = []
    for frac, name in zip(single_fracs, single_names):
        truth_rows.append({
            "vcf": "joint_single_donor.vcf",
            "sample_name": name,
            "true_donor_fraction": f"{frac:.6f}",
            "n_donors": "1",
        })
    for frac, name in zip(multi_fracs, multi_names):
        truth_rows.append({
            "vcf": "joint_multi_donor.vcf",
            "sample_name": name,
            "true_donor_fraction": f"{frac:.6f}",
            "n_donors": "2",
        })

    truth_path = outdir / "truth_table.tsv"
    write_truth_table(
        truth_rows,
        truth_path,
        fieldnames=["vcf", "sample_name", "true_donor_fraction", "n_donors"],
    )

    log.info("Generated joint VCFs in %s/", outdir)
    log.info("Truth table: %s", truth_path)

    # Clean up temp donor VCF
    donor_path.unlink(missing_ok=True)

    return 0


if __name__ == "__main__":
    sys.exit(main())
