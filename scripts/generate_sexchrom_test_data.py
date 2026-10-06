#!/usr/bin/env python3
"""Generate sex-chromosome test fixtures: FF, MM and MF host/donor pairs.

Writes one joint VCF per pair (HOST, DONOR and admixture columns, the format
``allomix detect --genotype-vcf X --admix-vcf X`` consumes) plus a truth table.
Each pair carries autosomal markers, non-PAR chrX markers with sex-aware ploidy
(male calls hemizygous but GATK-style diploid-encoded, with a spurious-het
rate), two PAR chrX markers (diploid in both sexes) and one marker on an alt
contig that the analysis must exclude. The MM host is guaranteed at least one
spurious chrX het so the het-drop rule has a case to act on.

Usage:
    python scripts/generate_sexchrom_test_data.py --outdir tests/test_data/sexchrom
"""

import argparse
import logging
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from script_utils import write_truth_table  # noqa: E402

from allomix.simulate import (  # noqa: E402
    build_joint_vcf_from_genotype_dicts,
    generate_related_genotypes,
    generate_sex_chrom_genotypes,
    write_joint_vcf,
)

log = logging.getLogger(__name__)

# (label, host sex, donor sex)
PAIRS = [("FF", "F", "F"), ("MM", "M", "M"), ("MF", "M", "F")]

# Donor fractions of the admixture columns (host 5% and 20%).
ADMIX_FRACTIONS = [0.95, 0.80]

ALT_CONTIG = "chrX_KI270880v1_alt"
DEFAULT_ERROR_RATE = 0.01
SEED_HASH_MODULUS = 2**31


def _alt_contig_marker(pos: int = 100_000) -> dict:
    """One fully informative diploid marker on an alt contig.

    Fully informative so that, if it leaks into the estimate, it is at least
    blended correctly (copy number 2 everywhere) and cannot bias the fraction.
    Its job is to be counted as excluded.
    """
    return {
        "chrom": ALT_CONTIG,
        "pos": pos,
        "ref": "C",
        "alt": "T",
        "host_gt": (0, 0),
        "donor_gt": (1, 1),
        "p_alt": 0.5,
        "informative": True,
    }


def build_pair_markers(
    host_sex: str,
    donor_sex: str,
    rng: random.Random,
    n_autosomal: int,
    n_chrx: int,
    n_par: int,
    spurious_het_rate: float,
    force_host_spurious_het: bool,
) -> list[dict]:
    """Assemble one pair's marker list: autosomes, chrX (non-PAR then PAR), alt contig.

    Autosomal markers are spread over chr1-22. With ``force_host_spurious_het``
    (the MM pair) the first non-PAR chrX marker is set to a spurious ``0/1``
    host call when the random draw produced none.
    """
    autosomal = generate_related_genotypes(n_autosomal, "unrelated", rng)
    for i, m in enumerate(autosomal):
        m["chrom"] = f"chr{(i % 22) + 1}"
        m["pos"] = 1_000_000 + i * 100_000

    chrx = generate_sex_chrom_genotypes(
        n_chrx,
        host_sex,
        donor_sex,
        rng,
        n_par=n_par,
        male_chrx_spurious_het_rate=spurious_het_rate,
    )
    if force_host_spurious_het and host_sex == "M":
        nonpar = [m for m in chrx if not m["par"]]
        if nonpar and not any(m["host_gt"] != m["host_true_gt"] for m in nonpar):
            nonpar[0]["host_gt"] = (0, 1)

    return [*autosomal, *chrx, _alt_contig_marker()]


def _count_spurious(markers: list[dict], who: str) -> int:
    return sum(
        1 for m in markers if m.get(f"{who}_cn") == 1 and m[f"{who}_gt"] != m[f"{who}_true_gt"]
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate sex-chromosome test fixtures (FF, MM, MF pairs).",
    )
    parser.add_argument(
        "--outdir",
        default="tests/test_data/sexchrom",
        help="Output directory (default: tests/test_data/sexchrom)",
    )
    parser.add_argument("--n-autosomal", type=int, default=60, help="Autosomal markers")
    parser.add_argument("--n-chrx", type=int, default=20, help="Non-PAR chrX markers")
    parser.add_argument("--n-par", type=int, default=2, help="PAR chrX markers")
    parser.add_argument("--depth", type=int, default=1000, help="Depth per marker")
    parser.add_argument(
        "--spurious-het-rate",
        type=float,
        default=0.04,
        help="Spurious het rate at male non-PAR chrX (default: 0.04)",
    )
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)

    admix_names = [f"ADMIX_F{f:.2f}" for f in ADMIX_FRACTIONS]
    truth_rows = []
    for label, host_sex, donor_sex in PAIRS:
        markers = build_pair_markers(
            host_sex,
            donor_sex,
            rng,
            n_autosomal=args.n_autosomal,
            n_chrx=args.n_chrx,
            n_par=args.n_par,
            spurious_het_rate=args.spurious_het_rate,
            force_host_spurious_het=(label == "MM"),
        )
        result = build_joint_vcf_from_genotype_dicts(
            markers,
            admix_fractions=ADMIX_FRACTIONS,
            admix_sample_names=admix_names,
            target_depth=args.depth,
            seed=rng.randint(0, SEED_HASH_MODULUS),
            error_rate=DEFAULT_ERROR_RATE,
        )
        vcf_name = f"joint_{label}.vcf"
        write_joint_vcf(result, outdir / vcf_name)

        n_autosomal_inf = sum(
            1
            for m in markers
            if m["chrom"].startswith("chr")
            and m["chrom"] != "chrX"
            and m["chrom"] != ALT_CONTIG
            and m["informative"]
        )
        n_chrx_inf = sum(
            1 for m in markers if m["chrom"] == "chrX" and not m.get("par") and m["informative"]
        )
        host_spurious = _count_spurious(markers, "host")
        donor_spurious = _count_spurious(markers, "donor")
        log.info(
            "%s: %d markers, %d informative (autosomal %d, chrX non-PAR %d); "
            "spurious male chrX hets host=%d donor=%d",
            vcf_name,
            result.num_markers,
            result.num_informative,
            n_autosomal_inf,
            n_chrx_inf,
            host_spurious,
            donor_spurious,
        )

        for frac, name in zip(ADMIX_FRACTIONS, admix_names):
            truth_rows.append(
                {
                    "vcf": vcf_name,
                    "sample_name": name,
                    "true_donor_fraction": f"{frac:.6f}",
                    "host_sex": host_sex,
                    "donor_sex": donor_sex,
                    "n_autosomal": args.n_autosomal,
                    "n_autosomal_informative": n_autosomal_inf,
                    "n_chrx_nonpar": args.n_chrx,
                    "n_chrx_nonpar_informative": n_chrx_inf,
                    "n_chrx_par": args.n_par,
                    "n_other_contig": 1,
                    "n_chrx_host_spurious_het": host_spurious,
                    "n_chrx_donor_spurious_het": donor_spurious,
                }
            )

    truth_path = outdir / "truth_table.tsv"
    write_truth_table(truth_rows, truth_path)
    log.info("Generated %d joint VCFs in %s/", len(PAIRS), outdir)
    log.info("Truth table: %s", truth_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
