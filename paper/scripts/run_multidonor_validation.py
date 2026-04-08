#!/usr/bin/env python3
"""Run allomix multi-donor validation on the 3-brothers test data.

Loads sibling donor test data, runs estimate_multi_donor() on each chimeric VCF,
compares against truth, and outputs facts for the paper.

Usage:
    python paper/scripts/run_multidonor_validation.py
"""

from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from allomix.chimerism import estimate_multi_donor
from allomix.genotype import classify_markers, parse_vcf

DATA_DIR = Path("tests/test_data/multidonor")
FACTS_DIR = Path("output/facts")


def load_truth_table(path: Path) -> list[dict]:
    """Load truth_table.tsv into list of dicts."""
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(
                {
                    "sample_name": row["sample_name"],
                    "true_f1": float(row["true_donor1_fraction"]),
                    "true_f2": float(row["true_donor2_fraction"]),
                    "true_host": float(row["true_host_fraction"]),
                    "num_markers": int(row["num_markers"]),
                    "num_informative_any": int(row["num_informative_any"]),
                }
            )
    return rows


def run_validation() -> list[dict]:
    """Run multi-donor estimation on all chimeric VCFs and return results."""
    host_vcf = str(DATA_DIR / "host.vcf")
    donor1_vcf = str(DATA_DIR / "donor1.vcf")
    donor2_vcf = str(DATA_DIR / "donor2.vcf")

    truth = load_truth_table(DATA_DIR / "truth_table.tsv")

    host = parse_vcf(host_vcf, min_dp=0, min_gq=0)
    donor1 = parse_vcf(donor1_vcf, min_dp=0, min_gq=0)
    donor2 = parse_vcf(donor2_vcf, min_dp=0, min_gq=0)

    results = []
    for t in truth:
        sample_path = str(DATA_DIR / f"{t['sample_name']}.vcf")
        admix = parse_vcf(sample_path, min_dp=0, min_gq=0)

        genotypes = classify_markers(
            host, [donor1, donor2], admix, min_dp=0, min_gq=0, pass_only=False
        )
        result = estimate_multi_donor(genotypes.informative, n_donors=2, error_rate=0.01)

        est_f1, est_f2 = result.donor_fractions
        err_f1 = est_f1 - t["true_f1"]
        err_f2 = est_f2 - t["true_f2"]
        err_total = (est_f1 + est_f2) - (t["true_f1"] + t["true_f2"])

        ci1_lo, ci1_hi = result.donor_fraction_cis[0]
        ci2_lo, ci2_hi = result.donor_fraction_cis[1]
        ci1_covers = ci1_lo <= t["true_f1"] <= ci1_hi
        ci2_covers = ci2_lo <= t["true_f2"] <= ci2_hi

        # Check asymmetric ranking: if truth says d1 > d2, does estimate agree?
        rank_correct = True
        if t["true_f1"] > t["true_f2"] + 0.01:
            rank_correct = est_f1 > est_f2
        elif t["true_f2"] > t["true_f1"] + 0.01:
            rank_correct = est_f2 > est_f1

        per_donor_inf = result.per_donor_n_informative or [0, 0]

        results.append(
            {
                "sample_name": t["sample_name"],
                "true_f1": t["true_f1"],
                "true_f2": t["true_f2"],
                "true_host": t["true_host"],
                "est_f1": est_f1,
                "est_f2": est_f2,
                "est_host": result.host_fraction,
                "err_f1": err_f1,
                "err_f2": err_f2,
                "err_total": err_total,
                "ci1_lo": ci1_lo,
                "ci1_hi": ci1_hi,
                "ci2_lo": ci2_lo,
                "ci2_hi": ci2_hi,
                "ci1_covers": ci1_covers,
                "ci2_covers": ci2_covers,
                "rank_correct": rank_correct,
                "n_informative": result.n_informative,
                "n_inf_d1": per_donor_inf[0],
                "n_inf_d2": per_donor_inf[1],
                "log_likelihood": result.log_likelihood,
            }
        )

    return results


def compute_metrics(rows: list[dict], donor_key: str = "f1") -> dict:
    """Compute aggregate metrics for one donor, excluding boundary fractions."""
    true_key = f"true_{donor_key}"
    err_key = f"err_{donor_key}"
    ci_lo_key = f"ci{donor_key[-1]}_lo"
    ci_hi_key = f"ci{donor_key[-1]}_hi"
    ci_covers_key = f"ci{donor_key[-1]}_covers"

    # Interior: exclude rows where this donor's true fraction is exactly 0 or 1
    interior = [r for r in rows if 0.0 < r[true_key] < 1.0]
    n = len(interior)
    if n == 0:
        return {}

    errors = [r[err_key] for r in interior]
    abs_errors = [abs(e) for e in errors]
    sq_errors = [e * e for e in errors]

    ci_covers = sum(1 for r in rows if r[ci_covers_key])
    ci_widths = [r[ci_hi_key] - r[ci_lo_key] for r in rows]

    return {
        "n_interior": n,
        "n_all": len(rows),
        "mae": sum(abs_errors) / n,
        "rmse": math.sqrt(sum(sq_errors) / n),
        "max_error": max(abs_errors),
        "ci_coverage": ci_covers / len(rows),
        "mean_ci_width": sum(ci_widths) / len(rows),
    }


def compute_total_metrics(rows: list[dict]) -> dict:
    """Compute metrics on total donor fraction (f1 + f2)."""
    # Interior: exclude pure host (f1+f2=0) and pure donor (f1+f2=1)
    interior = [r for r in rows if 0.0 < (r["true_f1"] + r["true_f2"]) < 1.0]
    n = len(interior)
    if n == 0:
        return {}

    errors = [r["err_total"] for r in interior]
    abs_errors = [abs(e) for e in errors]
    sq_errors = [e * e for e in errors]

    return {
        "n_interior": n,
        "mae": sum(abs_errors) / n,
        "rmse": math.sqrt(sum(sq_errors) / n),
        "max_error": max(abs_errors),
    }


def write_fact(name: str, data: dict) -> None:
    """Write a single-row facts CSV."""
    path = FACTS_DIR / f"{name}.csv"
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(data.keys()))
        writer.writeheader()
        writer.writerow(data)
    print(f"  Wrote {path}", file=sys.stderr)


def main() -> int:
    FACTS_DIR.mkdir(parents=True, exist_ok=True)

    print("Running multi-donor validation on 3-brothers test data...", file=sys.stderr)
    results = run_validation()

    # Print summary table to stderr
    print(
        f"\n{'Sample':<30} {'True d1':>8} {'True d2':>8} {'Est d1':>8} {'Est d2':>8} "
        f"{'Err d1':>8} {'Err d2':>8} {'CI1':>5} {'CI2':>5} {'Rank':>5}",
        file=sys.stderr,
    )
    print("-" * 110, file=sys.stderr)
    for r in results:
        print(
            f"{r['sample_name']:<30} {r['true_f1'] * 100:>7.1f}% {r['true_f2'] * 100:>7.1f}% "
            f"{r['est_f1'] * 100:>7.2f}% {r['est_f2'] * 100:>7.2f}% "
            f"{r['err_f1'] * 100:>+7.2f}% {r['err_f2'] * 100:>+7.2f}% "
            f"{'Y' if r['ci1_covers'] else 'N':>5} {'Y' if r['ci2_covers'] else 'N':>5} "
            f"{'Y' if r['rank_correct'] else 'N':>5}",
            file=sys.stderr,
        )

    # Compute per-donor metrics
    m_d1 = compute_metrics(results, "f1")
    m_d2 = compute_metrics(results, "f2")
    m_total = compute_total_metrics(results)

    print(
        f"\nDonor 1 — MAE: {m_d1['mae'] * 100:.2f}%  RMSE: {m_d1['rmse'] * 100:.2f}%  "
        f"Max: {m_d1['max_error'] * 100:.2f}%  CI cov: {m_d1['ci_coverage'] * 100:.1f}%",
        file=sys.stderr,
    )
    print(
        f"Donor 2 — MAE: {m_d2['mae'] * 100:.2f}%  RMSE: {m_d2['rmse'] * 100:.2f}%  "
        f"Max: {m_d2['max_error'] * 100:.2f}%  CI cov: {m_d2['ci_coverage'] * 100:.1f}%",
        file=sys.stderr,
    )
    print(
        f"Total   — MAE: {m_total['mae'] * 100:.2f}%  RMSE: {m_total['rmse'] * 100:.2f}%  "
        f"Max: {m_total['max_error'] * 100:.2f}%",
        file=sys.stderr,
    )

    # Ranking accuracy (only count asymmetric mixes)
    asymmetric = [r for r in results if abs(r["true_f1"] - r["true_f2"]) > 0.01]
    n_asymmetric = len(asymmetric)
    n_rank_correct = sum(1 for r in asymmetric if r["rank_correct"])
    print(
        f"\nRanking: {n_rank_correct}/{n_asymmetric} asymmetric mixes correctly ranked",
        file=sys.stderr,
    )

    # Informativity stats (same for all samples since same genotypes)
    n_inf_any = results[0]["n_informative"]
    n_inf_d1 = results[0]["n_inf_d1"]
    n_inf_d2 = results[0]["n_inf_d2"]

    # Write facts CSV for vibepaper
    write_fact(
        "multidonor",
        {
            "n_markers": 100,
            "n_informative_any": n_inf_any,
            "n_informative_d1": n_inf_d1,
            "n_informative_d2": n_inf_d2,
            "depth": 1000,
            "n_samples": len(results),
            "mae_d1_pct": round(m_d1["mae"] * 100, 2),
            "rmse_d1_pct": round(m_d1["rmse"] * 100, 2),
            "max_error_d1_pct": round(m_d1["max_error"] * 100, 2),
            "ci_coverage_d1_pct": round(m_d1["ci_coverage"] * 100, 1),
            "mae_d2_pct": round(m_d2["mae"] * 100, 2),
            "rmse_d2_pct": round(m_d2["rmse"] * 100, 2),
            "max_error_d2_pct": round(m_d2["max_error"] * 100, 2),
            "ci_coverage_d2_pct": round(m_d2["ci_coverage"] * 100, 1),
            "mae_total_pct": round(m_total["mae"] * 100, 2),
            "rmse_total_pct": round(m_total["rmse"] * 100, 2),
            "max_error_total_pct": round(m_total["max_error"] * 100, 2),
            "n_asymmetric": n_asymmetric,
            "n_rank_correct": n_rank_correct,
        },
    )

    # Write per-sample results CSV
    results_path = FACTS_DIR / "multidonor_results.csv"
    fields = [
        "sample_name",
        "true_f1_pct",
        "true_f2_pct",
        "est_f1_pct",
        "est_f2_pct",
        "err_f1_pct",
        "err_f2_pct",
        "ci1_lo_pct",
        "ci1_hi_pct",
        "ci2_lo_pct",
        "ci2_hi_pct",
        "ci1_covers",
        "ci2_covers",
        "rank_correct",
    ]
    with open(results_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in results:
            w.writerow(
                {
                    "sample_name": r["sample_name"],
                    "true_f1_pct": f"{r['true_f1'] * 100:.1f}",
                    "true_f2_pct": f"{r['true_f2'] * 100:.1f}",
                    "est_f1_pct": f"{r['est_f1'] * 100:.2f}",
                    "est_f2_pct": f"{r['est_f2'] * 100:.2f}",
                    "err_f1_pct": f"{r['err_f1'] * 100:.4f}",
                    "err_f2_pct": f"{r['err_f2'] * 100:.4f}",
                    "ci1_lo_pct": f"{r['ci1_lo'] * 100:.2f}",
                    "ci1_hi_pct": f"{r['ci1_hi'] * 100:.2f}",
                    "ci2_lo_pct": f"{r['ci2_lo'] * 100:.2f}",
                    "ci2_hi_pct": f"{r['ci2_hi'] * 100:.2f}",
                    "ci1_covers": r["ci1_covers"],
                    "ci2_covers": r["ci2_covers"],
                    "rank_correct": r["rank_correct"],
                }
            )
    print(f"  Wrote {results_path}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
