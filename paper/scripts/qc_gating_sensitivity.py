"""Quantify the QC effect-size gating on the SRP434573 real mixtures for the supplement.

Reruns allomix on the committed SRP434573 two-person mixtures and compares the
clinical-gated QC (the default: promote a reliability flag to REVIEW only when the
misfit is large) against the legacy p-value-only rule. Writes
``output/facts/qc_gating.csv`` (template variables for the supplementary QC section).

The point the numbers make: at panel depth (>1000x) a chi-squared goodness-of-fit
test and a discordant-site count are significant for a clinically trivial misfit, so
a bare p-value flags most of the series even though the beta-binomial model fits
(post-trim reduced chi-squared ~1). Gating on effect size confines REVIEW to the
samples whose interpretation actually turns on it. No sequencing or joint calling;
reads the committed genotype VCFs only.
"""

import csv
import statistics
import sys
from pathlib import Path

from srp434573_common import resolve_srp434573_genotypes_dir

from allomix.analysis import analyse_sample
from allomix.calibration.error_rates import load_error_table
from allomix.estimate.likelihood import PanelCalibration
from allomix.genotype import parse_vcf
from allomix.qc.qc import (
    GOF_PRETRIM_LOD_MULTIPLE,
    GOF_REVIEW_REDUCED_CHISQ,
    SWAP_REVIEW_FRACTION,
    _compute_gof,
    assess_quality,
)
from allomix.qc.relatedness import Relatedness

OUT = Path("output")
FACTS_DIR = OUT / "facts"
TWO_TSV = OUT / "srp434573_two_person.tsv"
GEN = resolve_srp434573_genotypes_dir()

# Match the run_srp434573_allomix.py / CLI defaults so the facts line up with the
# main SRP434573 results.
ERROR_RATE = 0.01
MIN_DP = 100


def _dilution_rows() -> dict[str, list[dict]]:
    """Group the two-person dilution timepoints (known fraction present) by mixture."""
    if not TWO_TSV.exists():
        sys.exit(f"Missing {TWO_TSV}; run run_srp434573_allomix.py first.")
    by_mix: dict[str, list[dict]] = {}
    with open(TWO_TSV, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if not r["known_pct"].strip():
                continue  # skip the pure host/donor endpoints
            by_mix.setdefault(r["mixture"], []).append(r)
    return by_mix


def main() -> None:
    by_mix = _dilution_rows()

    n_review_clinical = 0
    n_review_legacy = 0
    n_pass_clinical = 0
    reduced_chisq: list[float] = []
    n_gof_p_significant = 0  # post-trim GoF p < 0.01 (would-be legacy GoF flag)

    for mix, rows in sorted(by_mix.items()):
        host, donor = rows[0]["host"], rows[0]["donor"]
        gvcf = GEN / f"{mix}.SRP434573.vcf.gz"
        avcf = GEN / f"{mix}.admix.vcf.gz"
        if not gvcf.exists():
            continue
        h = parse_vcf(gvcf, sample=host, gt_ad_consistency=True)
        d = parse_vcf(gvcf, sample=donor, gt_ad_consistency=True)
        etab = GEN / f"{mix}.error_table.tsv"
        cal = PanelCalibration(errors=load_error_table(etab) if etab.exists() else {})

        for row in rows:
            a = parse_vcf(avcf, sample=row["sample"])
            if not a:
                continue
            an = analyse_sample(
                h,
                [d],
                a,
                min_dp=MIN_DP,
                min_gq=0,
                error_rate=ERROR_RATE,
                calibration=cal,
                robust="auto",
                marker_type_overdispersion=True,
                expected_relatedness=[Relatedness.UNRELATED],
                sample_name=row["sample"],
            )
            # Legacy status: re-assess the same fitted result under the p-only rules.
            legacy = assess_quality(
                an.result,
                an.genotypes,
                expected_relatedness=[Relatedness.UNRELATED],
                clinical_gating=False,
            )
            if an.qc.status == "REVIEW":
                n_review_clinical += 1
            if an.qc.status == "PASS":
                n_pass_clinical += 1
            if legacy.status == "REVIEW":
                n_review_legacy += 1

            rho = getattr(an.result, "rho", float("inf"))
            gof = _compute_gof(
                an.result.per_marker,
                rho=rho,
                n_fitted_params=2,
                error_rate=an.result.error_rate,
            )
            if gof is not None:
                reduced_chisq.append(gof.reduced_chisq)
                if gof.pval < 0.01:
                    n_gof_p_significant += 1

    n_timepoints = sum(len(v) for v in by_mix.values())
    facts = {
        "n_timepoints": str(n_timepoints),
        "n_pass_clinical": str(n_pass_clinical),
        "n_review_clinical": str(n_review_clinical),
        "n_review_legacy": str(n_review_legacy),
        "gof_p_significant": str(n_gof_p_significant),
        "gof_reduced_median": f"{statistics.median(reduced_chisq):.2f}",
        "gof_reduced_max": f"{max(reduced_chisq):.2f}",
        "gof_reduced_threshold": f"{GOF_REVIEW_REDUCED_CHISQ:g}",
        "pretrim_lod_multiple": f"{GOF_PRETRIM_LOD_MULTIPLE:g}",
        "swap_review_fraction_pct": f"{SWAP_REVIEW_FRACTION * 100:g}",
    }

    FACTS_DIR.mkdir(parents=True, exist_ok=True)
    out = FACTS_DIR / "qc_gating.csv"
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(list(facts))
        w.writerow([facts[k] for k in facts])
    sys.stderr.write(
        f"Wrote {out}: clinical REVIEW {n_review_clinical}/{n_timepoints}, "
        f"legacy {n_review_legacy}/{n_timepoints}, "
        f"reduced chi-sq median {facts['gof_reduced_median']} max {facts['gof_reduced_max']}, "
        f"{n_gof_p_significant} significant by p-value.\n"
    )


if __name__ == "__main__":
    main()
