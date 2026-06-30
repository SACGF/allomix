#!/usr/bin/env python3
"""Validate allomix across donor-host relatedness levels.

Generates synthetic host-donor pairs at different relatedness levels
(unrelated, cousin, half-sibling, sibling), blends chimeric samples,
and shows how accuracy and informative marker count vary.

Usage:
    python scripts/run_relatedness_validation.py
"""

import csv
import math
import random
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from paper_quick import qval  # noqa: E402  (also patches savefig for the watermark)

from allomix.estimate.chimerism import estimate_single_donor_bb  # noqa: E402
from allomix.genotype import classify_markers, parse_vcf  # noqa: E402
from allomix.simulate import (  # noqa: E402
    blend_vcfs,
    generate_paired_related_genotypes,
    write_genotype_vcf,
    write_vcf,
)

FRACTIONS = [0.0, 0.01, 0.02, 0.05, 0.10, 0.20, 0.50, 0.80, 0.95, 1.0]
RELATEDNESS_LEVELS = ["unrelated", "cousin", "half-sibling", "sibling"]
N_MARKERS = 100
DEPTH = 500
N_REPLICATES = qval(20, 4)  # paired donor-host panels (quick-build mode uses few)
FACTS_DIR = Path("output/facts")


def fraction_to_name(f: float) -> str:
    d = round(f * 100)
    h = 100 - d
    return f"host_{h}_donor_{d}"


def run_one_replicate(
    relatedness: str,
    rep: int,
    markers: list[dict],
    outdir: Path,
    seed: int,
) -> dict:
    """Run validation across fractions for one (pre-generated) donor-host pair.

    ``markers`` is the panel for this (relatedness, rep) from the paired design,
    so the host genotypes are shared across relatedness levels and only the donor
    differs by IBD sharing.
    """
    rep_dir = outdir / relatedness / f"rep_{rep}"
    rep_dir.mkdir(parents=True, exist_ok=True)

    n_informative = sum(1 for m in markers if m["informative"])

    host_vcf = rep_dir / "host.vcf"
    donor_vcf = rep_dir / "donor.vcf"
    write_genotype_vcf(markers, host_vcf, "host", key="host_gt")
    write_genotype_vcf(markers, donor_vcf, "donor", key="donor_gt")

    frac_results = []
    for frac in FRACTIONS:
        name = fraction_to_name(frac)
        blend_result = blend_vcfs(
            host_path=str(host_vcf),
            donor_path=str(donor_vcf),
            donor_fraction=frac,
            target_depth=DEPTH,
            sample_name=name,
            seed=seed + hash(str(frac)) % (2**31),
            realistic_biases=True,
            error_rate=0.01,
            locus_dropout_rate=0.016,
            depth_cv=0.43,
        )
        sample_vcf = rep_dir / f"{name}.vcf"
        write_vcf(blend_result, sample_vcf)

        host = parse_vcf(str(host_vcf), min_dp=0, min_gq=0)
        donor = parse_vcf(str(donor_vcf), min_dp=0, min_gq=0)
        admix = parse_vcf(str(sample_vcf), min_dp=0, min_gq=0)

        genotypes = classify_markers(host, [donor], admix, min_dp=0, min_gq=0, pass_only=False)

        if len(genotypes.informative) < 1:
            frac_results.append({
                "true_frac": frac,
                "est_frac": float("nan"),
                "error": float("nan"),
                "ci_covers": False,
                "n_informative": 0,
            })
            continue

        result = estimate_single_donor_bb(genotypes.informative, error_rate=0.01)
        error = result.donor_fraction - frac
        ci_covers = result.donor_fraction_ci[0] <= frac <= result.donor_fraction_ci[1]

        frac_results.append({
            "true_frac": frac,
            "est_frac": result.donor_fraction,
            "error": error,
            "ci_covers": ci_covers,
            "n_informative": result.n_informative,
        })

    # Aggregate metrics, interior fractions only.
    interior = [r for r in frac_results if 0.0 < r["true_frac"] < 1.0 and not math.isnan(r["error"])]
    if interior:
        errors = [r["error"] for r in interior]
        mae = sum(abs(e) for e in errors) / len(errors)
        rmse = math.sqrt(sum(e * e for e in errors) / len(errors))
    else:
        mae = rmse = float("nan")

    return {
        "relatedness": relatedness,
        "rep": rep,
        "n_markers": N_MARKERS,
        "n_informative": n_informative,
        "n_informative_allomix": frac_results[0]["n_informative"] if frac_results else 0,
        "mae": mae,
        "rmse": rmse,
        "frac_results": frac_results,
    }


def plot_results(all_results: dict[str, list[dict]], outdir: Path) -> None:
    """Generate relatedness comparison figures."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    colors = {
        "unrelated": "#2196F3",
        "cousin": "#4CAF50",
        "half-sibling": "#FF9800",
        "sibling": "#F44336",
    }
    labels = {
        "unrelated": "Unrelated",
        "cousin": "1st cousin",
        "half-sibling": "Half-sibling",
        "sibling": "Full sibling",
    }

    # Panel 1: Informative markers by relatedness.
    ax = axes[0]
    for rel in RELATEDNESS_LEVELS:
        reps = all_results[rel]
        n_inf = [r["n_informative"] for r in reps]
        x_jitter = [RELATEDNESS_LEVELS.index(rel) + random.gauss(0, 0.05) for _ in reps]
        ax.scatter(x_jitter, n_inf, c=colors[rel], s=40, alpha=0.6, zorder=3)
        mean_inf = sum(n_inf) / len(n_inf)
        ax.plot([RELATEDNESS_LEVELS.index(rel) - 0.2, RELATEDNESS_LEVELS.index(rel) + 0.2],
                [mean_inf, mean_inf], c=colors[rel], linewidth=3, zorder=4)

    ax.set_xticks(range(len(RELATEDNESS_LEVELS)))
    ax.set_xticklabels([labels[r] for r in RELATEDNESS_LEVELS], fontsize=10)
    ax.set_ylabel("Informative markers", fontsize=11)
    ax.set_title("Informative markers vs relatedness", fontsize=12)
    ax.grid(True, alpha=0.2, axis="y")

    # Panel 2: MAE by relatedness.
    ax = axes[1]
    for rel in RELATEDNESS_LEVELS:
        reps = all_results[rel]
        maes = [r["mae"] * 100 for r in reps if not math.isnan(r["mae"])]
        x_jitter = [RELATEDNESS_LEVELS.index(rel) + random.gauss(0, 0.05) for _ in maes]
        ax.scatter(x_jitter, maes, c=colors[rel], s=40, alpha=0.6, zorder=3)
        if maes:
            mean_mae = sum(maes) / len(maes)
            ax.plot([RELATEDNESS_LEVELS.index(rel) - 0.2, RELATEDNESS_LEVELS.index(rel) + 0.2],
                    [mean_mae, mean_mae], c=colors[rel], linewidth=3, zorder=4)

    ax.set_xticks(range(len(RELATEDNESS_LEVELS)))
    ax.set_xticklabels([labels[r] for r in RELATEDNESS_LEVELS], fontsize=10)
    ax.set_ylabel("Mean absolute error (%)", fontsize=11)
    ax.set_title("Accuracy vs relatedness", fontsize=12)
    ax.grid(True, alpha=0.2, axis="y")

    # Panel 3: Truth vs estimated across relatedness.
    ax = axes[2]
    for rel in RELATEDNESS_LEVELS:
        reps = all_results[rel]
        for rep in reps:
            truths = [r["true_frac"] * 100 for r in rep["frac_results"]
                      if not math.isnan(r.get("est_frac", float("nan")))]
            ests = [r["est_frac"] * 100 for r in rep["frac_results"]
                    if not math.isnan(r.get("est_frac", float("nan")))]
            ax.scatter(truths, ests, c=colors[rel], s=15, alpha=0.3, zorder=3)

    ax.plot([0, 100], [0, 100], "k--", alpha=0.4, linewidth=1)
    legend_elements = [Line2D([0], [0], marker="o", color="w",
                              markerfacecolor=colors[r], markersize=8, label=labels[r])
                       for r in RELATEDNESS_LEVELS]
    ax.legend(handles=legend_elements, fontsize=9, loc="upper left")
    ax.set_xlabel("True donor %", fontsize=11)
    ax.set_ylabel("Estimated donor %", fontsize=11)
    ax.set_title("Estimation across relatedness", fontsize=12)
    ax.set_aspect("equal")
    ax.set_xlim(-2, 102)
    ax.set_ylim(-2, 102)
    ax.grid(True, alpha=0.2)

    fig.suptitle(
        f"allomix performance by donor-host relatedness "
        f"({N_MARKERS} markers, {DEPTH}x depth, n={N_REPLICATES} replicates)",
        fontsize=13, y=1.02,
    )
    fig.tight_layout()
    fig.savefig(outdir / "fig4_relatedness.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved: {outdir / 'fig4_relatedness.png'}", file=sys.stderr)


def main() -> int:
    outdir = Path("output/relatedness_validation")
    outdir.mkdir(parents=True, exist_ok=True)
    FACTS_DIR.mkdir(parents=True, exist_ok=True)

    all_results: dict[str, list[dict]] = {}

    # Paired design: one shared host panel per replicate, each relatedness level's donor
    # derived from that same host through its IBD probabilities. Only the IBD sharing
    # varies across levels, so the informative-marker count is monotone non-increasing
    # with relatedness instead of being swamped by independent host draws (which
    # previously let 1st cousin show more informative markers than unrelated by chance).
    # Same approach as the nested LoD sweep.
    paired_panels = [
        generate_paired_related_genotypes(
            N_MARKERS, RELATEDNESS_LEVELS, random.Random(20240601 + rep)
        )
        for rep in range(N_REPLICATES)
    ]

    for rel in RELATEDNESS_LEVELS:
        print(f"\n{'='*50}", file=sys.stderr)
        print(f"Relatedness: {rel}", file=sys.stderr)
        print(f"{'='*50}", file=sys.stderr)

        reps = []
        for rep in range(N_REPLICATES):
            seed = 42 + rep * 1000 + RELATEDNESS_LEVELS.index(rel) * 7919
            markers = paired_panels[rep][rel]
            result = run_one_replicate(rel, rep, markers, outdir, seed)
            reps.append(result)
            print(
                f"  Rep {rep}: {result['n_informative']} informative, "
                f"MAE={result['mae']*100:.2f}%",
                file=sys.stderr,
            )

        all_results[rel] = reps

        n_infs = [r["n_informative"] for r in reps]
        maes = [r["mae"] * 100 for r in reps if not math.isnan(r["mae"])]
        rmses = [r["rmse"] * 100 for r in reps if not math.isnan(r["rmse"])]

        mean_inf = sum(n_infs) / len(n_infs)
        mean_mae = sum(maes) / len(maes) if maes else float("nan")
        mean_rmse = sum(rmses) / len(rmses) if rmses else float("nan")

        print(f"  Mean informative: {mean_inf:.1f}", file=sys.stderr)
        print(f"  Mean MAE: {mean_mae:.2f}%", file=sys.stderr)
        print(f"  Mean RMSE: {mean_rmse:.2f}%", file=sys.stderr)

        csv_path = FACTS_DIR / f"rel_{rel.replace('-', '_')}.csv"
        with open(csv_path, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=[
                "relatedness", "n_replicates", "n_markers",
                "mean_informative", "min_informative", "max_informative",
                "mean_mae_pct", "mean_rmse_pct",
            ])
            w.writeheader()
            w.writerow({
                "relatedness": rel,
                "n_replicates": N_REPLICATES,
                "n_markers": N_MARKERS,
                "mean_informative": round(mean_inf, 1),
                "min_informative": min(n_infs),
                "max_informative": max(n_infs),
                "mean_mae_pct": round(mean_mae, 2),
                "mean_rmse_pct": round(mean_rmse, 2),
            })

    summary_path = outdir / "relatedness_summary.tsv"
    with open(summary_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["relatedness", "mean_informative", "min", "max", "mean_mae_pct", "mean_rmse_pct"])
        for rel in RELATEDNESS_LEVELS:
            reps = all_results[rel]
            n_infs = [r["n_informative"] for r in reps]
            maes = [r["mae"] * 100 for r in reps if not math.isnan(r["mae"])]
            rmses = [r["rmse"] * 100 for r in reps if not math.isnan(r["rmse"])]
            w.writerow([
                rel,
                f"{sum(n_infs)/len(n_infs):.1f}",
                min(n_infs),
                max(n_infs),
                f"{sum(maes)/len(maes):.2f}" if maes else "NA",
                f"{sum(rmses)/len(rmses):.2f}" if rmses else "NA",
            ])
    print(f"\nSummary: {summary_path}", file=sys.stderr)

    plot_results(all_results, outdir)

    src = outdir / "fig4_relatedness.png"
    if src.exists():
        shutil.copy2(src, FACTS_DIR / "fig4_relatedness.png")
        print(f"Copied to {FACTS_DIR / 'fig4_relatedness.png'}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
