#!/usr/bin/env python3
"""Generate simulated engraftment timeline figure for the allomix paper.

Simulates a 6-timepoint post-HSCT engraftment trajectory (day +14 to +365)
with a clinically interesting dip at day +180. Runs N=5 independent replicates
with different random seeds.

Usage:
    python paper/scripts/generate_timeline_figure.py
"""

import csv
import math
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import paper_quick  # noqa: E402, F401  -- quick-build watermark (import for side effect)

from allomix.chimerism import estimate_single_donor_bb  # noqa: E402
from allomix.genotype import classify_markers, parse_vcf  # noqa: E402
from allomix.simulate import (  # noqa: E402
    blend_vcfs,
    generate_marker_biases_realistic,
    write_vcf,
)
from allomix.simulate import (  # noqa: E402
    parse_text_vcf as sim_parse_vcf,
)

FACTS_DIR = Path("output/facts")

TIMEPOINTS = [
    {"day": 14, "donor_frac": 0.15, "label": "Day +14"},
    {"day": 28, "donor_frac": 0.55, "label": "Day +28"},
    {"day": 60, "donor_frac": 0.85, "label": "Day +60"},
    {"day": 100, "donor_frac": 0.95, "label": "Day +100"},
    {"day": 180, "donor_frac": 0.92, "label": "Day +180"},
    {"day": 365, "donor_frac": 0.97, "label": "Day +365"},
]

# Simulation parameters (match main validation)
DEPTH = 500
N_REPLICATES = 5
DEPTH_CV = 0.43
LOCUS_DROPOUT_RATE = 0.016
ERROR_RATE = 0.01
BASE_SEED = 42


def run_timeline(
    host_vcf: str,
    donor_vcf: str,
    outdir: Path,
) -> list[list[dict]]:
    """Run timeline simulation with N replicates.

    Returns a list of replicates, each a list of per-timepoint result dicts.
    """
    _, host_records = sim_parse_vcf(host_vcf)
    _, donor_records = sim_parse_vcf(donor_vcf)
    donor_loci = {r.locus for r in donor_records}
    n_shared = sum(1 for r in host_records if r.locus in donor_loci)

    all_replicates = []

    for rep in range(N_REPLICATES):
        rep_seed = BASE_SEED + rep * 1000
        vcf_dir = outdir / f"seed_{rep_seed}"
        vcf_dir.mkdir(parents=True, exist_ok=True)

        bias_rng = random.Random(rep_seed)
        fixed_biases = generate_marker_biases_realistic(n_shared, bias_rng)

        rep_results = []
        for tp in TIMEPOINTS:
            sample_name = f"day_{tp['day']}"
            sample_seed = rep_seed + tp["day"]

            result = blend_vcfs(
                host_path=host_vcf,
                donor_path=donor_vcf,
                donor_fraction=tp["donor_frac"],
                target_depth=DEPTH,
                sample_name=sample_name,
                seed=sample_seed,
                fixed_biases=fixed_biases,
                locus_dropout_rate=LOCUS_DROPOUT_RATE,
                depth_cv=DEPTH_CV,
            )
            write_vcf(result, vcf_dir / f"{sample_name}.vcf")

            host = parse_vcf(host_vcf, min_dp=0, min_gq=0)
            donor = parse_vcf(donor_vcf, min_dp=0, min_gq=0)
            admix = parse_vcf(str(vcf_dir / f"{sample_name}.vcf"), min_dp=0, min_gq=0)
            genotypes = classify_markers(
                host, [donor], admix, min_dp=0, min_gq=0, pass_only=False
            )
            est = estimate_single_donor_bb(genotypes.informative, error_rate=ERROR_RATE)

            error = est.donor_fraction - tp["donor_frac"]
            ci_covers = est.donor_fraction_ci[0] <= tp["donor_frac"] <= est.donor_fraction_ci[1]

            rep_results.append({
                "day": tp["day"],
                "label": tp["label"],
                "true_frac": tp["donor_frac"],
                "est_frac": est.donor_fraction,
                "error": error,
                "ci_lo": est.donor_fraction_ci[0],
                "ci_hi": est.donor_fraction_ci[1],
                "ci_covers": ci_covers,
                "n_informative": est.n_informative,
            })

        all_replicates.append(rep_results)

        max_err = max(abs(r["error"]) for r in rep_results)
        mae = sum(abs(r["error"]) for r in rep_results) / len(rep_results)
        print(
            f"  Rep {rep}: MAE={mae*100:.3f}%  Max={max_err*100:.3f}%",
            file=sys.stderr,
        )

    return all_replicates


def write_facts(all_replicates: list[list[dict]]) -> None:
    """Write timeline facts CSV for vibepaper."""
    FACTS_DIR.mkdir(parents=True, exist_ok=True)

    n_reps = len(all_replicates)
    n_timepoints = len(TIMEPOINTS)

    rep_maes = []
    rep_max_errs = []
    for rep in all_replicates:
        abs_errors = [abs(r["error"]) for r in rep]
        rep_maes.append(sum(abs_errors) / len(abs_errors))
        rep_max_errs.append(max(abs_errors))

    mae_mean = sum(rep_maes) / n_reps
    mae_sd = math.sqrt(sum((v - mae_mean) ** 2 for v in rep_maes) / (n_reps - 1))
    max_err_mean = sum(rep_max_errs) / n_reps
    max_err_sd = math.sqrt(
        sum((v - max_err_mean) ** 2 for v in rep_max_errs) / (n_reps - 1)
    )

    dip_errors = [rep[4]["error"] for rep in all_replicates]  # index 4 = day +180 dip
    dip_abs_mean = sum(abs(e) for e in dip_errors) / n_reps

    total_ci = sum(1 for rep in all_replicates for r in rep if r["ci_covers"])
    ci_coverage = total_ci / (n_reps * n_timepoints)

    facts = {
        "n_replicates": n_reps,
        "n_timepoints": n_timepoints,
        "depth": DEPTH,
        "mae_pct": round(mae_mean * 100, 2),
        "mae_sd_pct": round(mae_sd * 100, 2),
        "max_error_pct": round(max_err_mean * 100, 2),
        "max_error_sd_pct": round(max_err_sd * 100, 2),
        "dip_abs_error_pct": round(dip_abs_mean * 100, 2),
        "ci_coverage_pct": round(ci_coverage * 100, 1),
    }

    path = FACTS_DIR / "timeline.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(facts.keys()))
        writer.writeheader()
        writer.writerow(facts)
    print(f"  Wrote {path}", file=sys.stderr)


def plot_figure(all_replicates: list[list[dict]]) -> None:
    """Generate the timeline figure."""
    fig, ax = plt.subplots(figsize=(7, 4.5))

    days = [tp["day"] for tp in TIMEPOINTS]
    true_fracs = [tp["donor_frac"] * 100 for tp in TIMEPOINTS]

    for rep in all_replicates:
        est_fracs = [r["est_frac"] * 100 for r in rep]
        ax.plot(days, est_fracs, "o-", color="steelblue", alpha=0.2, linewidth=1, markersize=4)

    n_tp = len(TIMEPOINTS)
    mean_est = []
    sd_est = []
    mean_ci_lo = []
    mean_ci_hi = []
    for i in range(n_tp):
        ests = [rep[i]["est_frac"] * 100 for rep in all_replicates]
        m = sum(ests) / len(ests)
        s = math.sqrt(sum((v - m) ** 2 for v in ests) / (len(ests) - 1))
        mean_est.append(m)
        sd_est.append(s)
        mean_ci_lo.append(sum(rep[i]["ci_lo"] * 100 for rep in all_replicates) / len(all_replicates))
        mean_ci_hi.append(sum(rep[i]["ci_hi"] * 100 for rep in all_replicates) / len(all_replicates))

    ax.plot(days, true_fracs, "s--", color="0.4", linewidth=1.5, markersize=5, label="True", zorder=4)

    ax.fill_between(days, mean_ci_lo, mean_ci_hi, color="steelblue", alpha=0.15, label="Mean 95% CI")
    ax.plot(
        days, mean_est, "o-",
        color="steelblue", linewidth=2.5, markersize=7,
        markeredgecolor="white", markeredgewidth=0.8,
        label="Mean estimate", zorder=5,
    )

    phase_labels = [
        (14, "Early\nengraftment"),
        (100, "Stable"),
        (180, "Dip"),
        (365, "Recovery"),
    ]
    for day, phase in phase_labels:
        idx = next(i for i, tp in enumerate(TIMEPOINTS) if tp["day"] == day)
        y_offset = -8 if day != 180 else -8
        ax.annotate(
            phase,
            xy=(day, true_fracs[idx]),
            xytext=(0, y_offset),
            textcoords="offset points",
            fontsize=8,
            color="0.4",
            ha="center",
            va="top",
        )

    ax.set_xlabel("Day post-HSCT", fontsize=11)
    ax.set_ylabel("Donor chimerism (%)", fontsize=11)
    ax.set_title(
        "Simulated post-HSCT engraftment monitoring",
        fontsize=12, fontweight="bold", loc="left",
    )
    ax.set_xlim(-10, 390)
    ax.set_ylim(0, 105)
    ax.set_xticks(days)
    ax.set_xticklabels([f"+{d}" for d in days])
    ax.grid(True, alpha=0.2)
    ax.legend(fontsize=9, loc="lower right")

    fig.tight_layout()

    FACTS_DIR.mkdir(parents=True, exist_ok=True)
    fig_path = FACTS_DIR / "fig_timeline.png"
    fig.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved {fig_path}", file=sys.stderr)


def main() -> int:
    print("Generating timeline figure...", file=sys.stderr)

    host_vcf = "tests/test_data/host.vcf"
    donor_vcf = "tests/test_data/donor.vcf"
    outdir = Path("output/timeline")
    outdir.mkdir(parents=True, exist_ok=True)

    all_replicates = run_timeline(host_vcf, donor_vcf, outdir)
    write_facts(all_replicates)
    plot_figure(all_replicates)

    return 0


if __name__ == "__main__":
    sys.exit(main())
