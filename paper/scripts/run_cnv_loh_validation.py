#!/usr/bin/env python3
"""Measure how host-genome CN-LoH degrades chimerism estimates (issue #13).

The HSCT recipient is usually a haematological malignancy patient, so the
residual or relapsing host clone routinely carries somatic copy-number
changes. Copy-neutral loss of heterozygosity (CN-LoH, acquired uniparental
disomy) is common in AML/MDS: the clone retains two copies of one germline
homolog, turning a host heterozygous marker into an effective homozygote
without any copy-number change. The host genotype VCF is taken from a clean
germline reference, so classification is unaffected, but the affected markers
in the admixture sample carry a biased ALT VAF, which can pull the donor
fraction estimate.

This sweep blends synthetic chimeric samples with a controllable burden of
host CN-LoH and reports the effect on point-estimate accuracy (MAE, signed
bias), CI coverage, and how many rogue markers the estimator's 3-SD outlier
flag catches. The fraction_affected=0 cells are the no-aberration baseline.

Design follows run_lod_validation.py: genotypes and per-marker capture biases
are fixed per (relatedness, replicate) "pair" and reused across every
aberration/fraction cell, so the only thing changing within a pair is the
CN-LoH burden and the sequencing draw.

Outputs:
  output/facts/cnv_loh_raw.csv       # one row per replicate per cell
  output/facts/cnv_loh_summary.csv   # one row per (rel, clonal, burden, true_frac)
  output/facts/cnv_loh_headline.csv  # single-row headline snapshot

Usage:
    python paper/scripts/run_cnv_loh_validation.py
    python paper/scripts/run_cnv_loh_validation.py --n-reps 10 --n-workers 8
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import random
import statistics
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from allomix.chimerism import estimate_single_donor_bb  # noqa: E402
from allomix.genotype import classify_markers, parse_vcf  # noqa: E402
from allomix.simulate import (  # noqa: E402
    assign_cnv_aberrations,
    blend_vcfs,
    generate_marker_biases_realistic,
    generate_related_genotypes,
    write_genotype_vcf,
    write_vcf,
)

# --- Sweep grid --------------------------------------------------------------

RELATEDNESS_LEVELS = ["unrelated", "sibling"]
# Copy-number aberration kinds applied to the host clone. CN-LoH is copy-neutral
# (allele-balance effect at het markers); deletion (CN1) and gain (CN3) also
# change the locus DNA mass, so they shift the local mixing fraction even at
# homozygous markers.
KINDS = ["cnloh", "deletion", "gain"]
# Fraction of eligible host markers carrying the aberration. For cnloh, only het
# markers are eligible; for deletion/gain every marker is. 0.0 is the baseline,
# run once per pair and shared across kinds.
BURDEN_LEVELS = [0.0, 0.1, 0.25, 0.5]
# Fraction of host cells that are the aberrant clone (1.0 = pure clone, e.g.
# diagnosis/relapse; 0.5 = clone is half the residual host).
CLONAL_FRACTIONS = [0.5, 1.0]
# True donor fractions spanning host-dominant to donor-dominant chimerism.
TRUE_FRACTIONS = [0.2, 0.5, 0.8, 0.9, 0.95, 0.99]

N_MARKERS = 100
DEPTH = 1000
MAF_RANGE = (0.2, 0.5)
ERROR_RATE = 0.01
LOCUS_DROPOUT_RATE = 0.016
DEPTH_CV = 0.43
ESTIMATOR_GRID_STEPS = 201

DEFAULT_N_REPS = 20

FACTS_DIR = Path("output/facts")
WORK_DIR = Path("output/cnv_loh_validation")


def derive_seed(*parts: object) -> int:
    """Deterministic seed from arbitrary parts (stable across processes)."""
    digest = hashlib.sha256(repr(parts).encode("utf-8")).digest()[:4]
    return int.from_bytes(digest, "big")


def _eval_cell(
    host_vcf: Path,
    donor_vcf: Path,
    host_md: list,
    donor_md: list,
    biases: list[float],
    aberrations: list | None,
    true_frac: float,
    blend_seed: int,
    admix_path: Path,
) -> dict:
    """Blend one admix sample and estimate chimerism; return result fields."""
    blend = blend_vcfs(
        host_path=str(host_vcf),
        donor_path=str(donor_vcf),
        donor_fraction=true_frac,
        target_depth=DEPTH,
        sample_name="admix",
        seed=blend_seed,
        fixed_biases=biases,
        error_rate=ERROR_RATE,
        locus_dropout_rate=LOCUS_DROPOUT_RATE,
        depth_cv=DEPTH_CV,
        host_aberrations=aberrations,
    )
    bias_dict = (
        {(c, p, r, a): b for c, p, r, a, b in blend.marker_biases}
        if blend.marker_biases is not None
        else None
    )
    write_vcf(blend, admix_path)
    admix_md = parse_vcf(str(admix_path), min_dp=0, min_gq=0)

    genos = classify_markers(host_md, [donor_md], admix_md, min_dp=0, min_gq=0, pass_only=False)
    if len(genos.informative) < 1:
        return dict(
            est_frac=float("nan"), ci_lo=float("nan"), ci_hi=float("nan"),
            n_informative=0, n_flagged=0,
            est_frac_robust=float("nan"), ci_lo_robust=float("nan"),
            ci_hi_robust=float("nan"), n_robust_excluded=0,
        )

    result = estimate_single_donor_bb(
        genos.informative, error_rate=ERROR_RATE,
        grid_steps=ESTIMATOR_GRID_STEPS, marker_biases=bias_dict,
    )
    # Robust refit (default policy) to show the mitigation effect in the paper.
    robust = estimate_single_donor_bb(
        genos.informative, error_rate=ERROR_RATE,
        grid_steps=ESTIMATOR_GRID_STEPS, marker_biases=bias_dict, robust="auto",
    )
    n_flagged = sum(1 for mr in result.per_marker if not mr.included)
    return dict(
        est_frac=result.donor_fraction,
        ci_lo=result.donor_fraction_ci[0],
        ci_hi=result.donor_fraction_ci[1],
        n_informative=result.n_informative,
        n_flagged=n_flagged,
        est_frac_robust=robust.donor_fraction,
        ci_lo_robust=robust.donor_fraction_ci[0],
        ci_hi_robust=robust.donor_fraction_ci[1],
        n_robust_excluded=robust.n_robust_excluded,
    )


def run_pair(relatedness: str, rep: int, base_seed: int) -> list[dict]:
    """Evaluate every aberration cell for one fixed genotype pair.

    Genotypes and capture biases are fixed for this (relatedness, rep). Within
    the pair we vary the aberration kind, burden, clonal fraction, and donor
    fraction. The no-aberration baseline (burden 0) is run once and shared
    across kinds.
    """
    pair_dir = WORK_DIR / f"{relatedness}_rep{rep}"
    pair_dir.mkdir(parents=True, exist_ok=True)

    gt_rng = random.Random(derive_seed("gt", relatedness, rep, base_seed))
    markers = generate_related_genotypes(N_MARKERS, relatedness, gt_rng, maf_range=MAF_RANGE)

    host_vcf = pair_dir / "host.vcf"
    donor_vcf = pair_dir / "donor.vcf"
    write_genotype_vcf(markers, host_vcf, "host", key="host_gt")
    write_genotype_vcf(markers, donor_vcf, "donor", key="donor_gt")
    host_md = parse_vcf(str(host_vcf), min_dp=0, min_gq=0)
    donor_md = parse_vcf(str(donor_vcf), min_dp=0, min_gq=0)

    bias_rng = random.Random(derive_seed("bias", relatedness, rep, base_seed))
    biases = generate_marker_biases_realistic(N_MARKERS, bias_rng)

    admix_path = pair_dir / "admix.vcf"
    rows: list[dict] = []

    # Baseline (no aberration), shared by every kind.
    for true_frac in TRUE_FRACTIONS:
        blend_seed = derive_seed("blend", relatedness, rep, "baseline", true_frac)
        row = {
            "relatedness": relatedness, "rep": rep, "kind": "baseline",
            "burden": 0.0, "clonal_fraction": 0.0, "true_frac": true_frac,
            "n_affected": 0, "seed": blend_seed,
        }
        row.update(_eval_cell(host_vcf, donor_vcf, host_md, donor_md, biases,
                              None, true_frac, blend_seed, admix_path))
        rows.append(row)

    for kind in KINDS:
        for burden in BURDEN_LEVELS:
            if burden == 0.0:
                continue  # baseline already covered above
            for clonal in CLONAL_FRACTIONS:
                aberr_rng = random.Random(
                    derive_seed("aberr", relatedness, rep, kind, burden, clonal)
                )
                aberrations = assign_cnv_aberrations(
                    markers, burden, clonal, aberr_rng, kind=kind
                )
                n_affected = sum(1 for a in aberrations if a is not None)

                for true_frac in TRUE_FRACTIONS:
                    blend_seed = derive_seed(
                        "blend", relatedness, rep, kind, burden, clonal, true_frac
                    )
                    row = {
                        "relatedness": relatedness, "rep": rep, "kind": kind,
                        "burden": burden, "clonal_fraction": clonal, "true_frac": true_frac,
                        "n_affected": n_affected, "seed": blend_seed,
                    }
                    row.update(_eval_cell(host_vcf, donor_vcf, host_md, donor_md, biases,
                                          aberrations, true_frac, blend_seed, admix_path))
                    rows.append(row)

    return rows


def summarise(raw: list[dict]) -> list[dict]:
    """Aggregate raw rows into per-cell summary statistics."""
    cells: dict[tuple, list[dict]] = {}
    for r in raw:
        key = (r["relatedness"], r["kind"], r["clonal_fraction"], r["burden"], r["true_frac"])
        cells.setdefault(key, []).append(r)

    out: list[dict] = []
    for (relatedness, kind, clonal, burden, true_frac), rs in sorted(cells.items()):
        valid = [r for r in rs if r["est_frac"] == r["est_frac"]]  # drop NaN
        if not valid:
            continue
        errs = [r["est_frac"] - true_frac for r in valid]
        covered = [1 if r["ci_lo"] <= true_frac <= r["ci_hi"] else 0 for r in valid]
        rob_valid = [r for r in valid if r["est_frac_robust"] == r["est_frac_robust"]]
        rob_errs = [r["est_frac_robust"] - true_frac for r in rob_valid]
        rob_covered = [
            1 if r["ci_lo_robust"] <= true_frac <= r["ci_hi_robust"] else 0 for r in rob_valid
        ]
        out.append(
            {
                "relatedness": relatedness,
                "kind": kind,
                "clonal_fraction": clonal,
                "burden": burden,
                "true_frac": true_frac,
                "n_reps": len(valid),
                "mean_est": statistics.fmean(r["est_frac"] for r in valid),
                "bias": statistics.fmean(errs),
                "mae": statistics.fmean(abs(e) for e in errs),
                "rmse": (statistics.fmean(e * e for e in errs)) ** 0.5,
                "ci_coverage": statistics.fmean(covered),
                "mae_robust": statistics.fmean(abs(e) for e in rob_errs) if rob_errs else float("nan"),
                "bias_robust": statistics.fmean(rob_errs) if rob_errs else float("nan"),
                "ci_coverage_robust": statistics.fmean(rob_covered) if rob_covered else float("nan"),
                "mean_n_affected": statistics.fmean(r["n_affected"] for r in valid),
                "mean_n_flagged": statistics.fmean(r["n_flagged"] for r in valid),
                "mean_n_robust_excluded": statistics.fmean(r["n_robust_excluded"] for r in valid),
                "mean_n_informative": statistics.fmean(r["n_informative"] for r in valid),
            }
        )
    return out


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in rows:
            writer.writerow(r)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-reps", type=int, default=DEFAULT_N_REPS)
    parser.add_argument("--n-workers", type=int, default=4)
    parser.add_argument("--base-seed", type=int, default=2026)
    args = parser.parse_args()

    WORK_DIR.mkdir(parents=True, exist_ok=True)

    jobs = [(rel, rep) for rel in RELATEDNESS_LEVELS for rep in range(args.n_reps)]
    raw: list[dict] = []

    if args.n_workers <= 1:
        for rel, rep in jobs:
            raw.extend(run_pair(rel, rep, args.base_seed))
    else:
        with ProcessPoolExecutor(max_workers=args.n_workers) as ex:
            futures = {ex.submit(run_pair, rel, rep, args.base_seed): (rel, rep) for rel, rep in jobs}
            for fut in as_completed(futures):
                raw.extend(fut.result())

    raw_fields = [
        "relatedness", "rep", "kind", "burden", "clonal_fraction", "true_frac", "n_affected",
        "seed", "est_frac", "ci_lo", "ci_hi", "n_informative", "n_flagged",
        "est_frac_robust", "ci_lo_robust", "ci_hi_robust", "n_robust_excluded",
    ]
    write_csv(FACTS_DIR / "cnv_loh_raw.csv", raw, raw_fields)

    summary = summarise(raw)
    summary_fields = [
        "relatedness", "kind", "clonal_fraction", "burden", "true_frac", "n_reps",
        "mean_est", "bias", "mae", "rmse", "ci_coverage",
        "mae_robust", "bias_robust", "ci_coverage_robust",
        "mean_n_affected", "mean_n_flagged", "mean_n_robust_excluded", "mean_n_informative",
    ]
    write_csv(FACTS_DIR / "cnv_loh_summary.csv", summary, summary_fields)

    headline = build_headline(summary)
    write_csv(FACTS_DIR / "cnv_loh_headline.csv", headline, ["metric", "value"])

    print(f"Wrote {len(raw)} raw rows, {len(summary)} summary cells to {FACTS_DIR}")
    for h in headline:
        print(f"  {h['metric']}: {h['value']}")


def build_headline(summary: list[dict]) -> list[dict]:
    """Pull interpretable headline numbers per aberration kind."""
    max_burden = max(BURDEN_LEVELS)

    low_burden = min(b for b in BURDEN_LEVELS if b > 0)

    def mae_mean(rows, field="mae"):
        vals = [s[field] for s in rows if s[field] == s[field]]
        return statistics.fmean(vals) if vals else float("nan")

    headline: list[dict] = []
    for rel in RELATEDNESS_LEVELS:
        base = mae_mean([s for s in summary if s["relatedness"] == rel and s["kind"] == "baseline"])
        headline.append({"metric": f"mae_baseline_{rel}", "value": round(base, 5)})
        for kind in KINDS:
            worst_rows = [
                s for s in summary
                if s["relatedness"] == rel and s["kind"] == kind
                and abs(s["burden"] - max_burden) < 1e-9
                and abs(s["clonal_fraction"] - 1.0) < 1e-9
            ]
            worst = mae_mean(worst_rows)
            headline.append({"metric": f"mae_{kind}_b{max_burden}_pureclone_{rel}", "value": round(worst, 5)})
            headline.append(
                {"metric": f"mae_inflation_x_{kind}_{rel}",
                 "value": round(worst / base, 2) if base > 0 else float("nan")}
            )
            # Robust mitigation at a realistic low burden, pure clone.
            low_rows = [
                s for s in summary
                if s["relatedness"] == rel and s["kind"] == kind
                and abs(s["burden"] - low_burden) < 1e-9
                and abs(s["clonal_fraction"] - 1.0) < 1e-9
            ]
            headline.append({"metric": f"mae_{kind}_b{low_burden}_std_{rel}",
                             "value": round(mae_mean(low_rows), 5)})
            headline.append({"metric": f"mae_{kind}_b{low_burden}_robust_{rel}",
                             "value": round(mae_mean(low_rows, "mae_robust"), 5)})
    # Outlier flagging at the highest burden, pure clone (across kinds/relatedness).
    flagged_rows = [s for s in summary if s["kind"] != "baseline"
                    and abs(s["burden"] - max_burden) < 1e-9
                    and abs(s["clonal_fraction"] - 1.0) < 1e-9]
    if flagged_rows:
        headline.append({"metric": "mean_markers_affected_highburden",
                         "value": round(statistics.fmean(s["mean_n_affected"] for s in flagged_rows), 2)})
        headline.append({"metric": "mean_markers_flagged_highburden",
                         "value": round(statistics.fmean(s["mean_n_flagged"] for s in flagged_rows), 2)})
    return headline


if __name__ == "__main__":
    main()
