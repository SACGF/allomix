#!/usr/bin/env python3
"""Compare the symmetric vs one-sided robust trim at low host fraction.

This is the in-silico check for the direction-aware robust trim (Obs 1 in
``claude/further_improvements.md``). The motivating real-data observation on
SRP434573 was that ``allomix``'s robust refit drops a runaway fraction of markers
as the host fraction falls and the MLE host estimate collapses toward zero (in
the worst case reading exactly 0% when the truth is ~1%), because the symmetric
median/MAD cut trims the host-carrying markers, which at low host fraction sit
off the donor-dominated fit and read as outliers.

The fix makes the trim one-sided: a marker whose residual deviates toward host
presence is never trimmed (``chimerism.ROBUST_ONE_SIDED``). This script blends
synthetic declining-chimerism mixtures (host is the minor component swept from
0.5% to 10%) over many seeds and, per blended sample, runs the robust estimator
under both trim modes plus the host-presence detector. It reports, per known
host fraction, the mean host estimate, the robust drop fraction, and the rate at
which the MLE collapses to ~0, for each mode.

Ground truth is known here, which is the point: the SRP434573 numbers motivated
the change, but the threshold-free claim ("one-sided recovers low-fraction host
signal the symmetric cut destroys") has to hold on independent simulated data
before we trust it.

Usage:
    python scripts/validate_onesided_trim.py
    python scripts/validate_onesided_trim.py --n-reps 30 --n-markers 200
"""

import argparse
import random
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from allomix import chimerism  # noqa: E402
from allomix.chimerism import estimate_single_donor_bb  # noqa: E402
from allomix.genotype import classify_markers, parse_vcf  # noqa: E402
from allomix.host_presence import host_presence_test  # noqa: E402
from allomix.likelihood import PanelCalibration  # noqa: E402
from allomix.simulate import (  # noqa: E402
    blend_vcfs,
    generate_marker_biases_realistic,
    generate_related_genotypes,
    write_genotype_vcf,
    write_vcf,
)

# Declining-chimerism / relapse series: host is the minor (titrated) component,
# matching the SRP434573 mapping. Includes 0 for the blank.
HOST_FRACTIONS = [0.0, 0.005, 0.01, 0.0125, 0.025, 0.05, 0.10]
MAF_RANGE = (0.2, 0.5)
ERROR_RATE = 0.01
LOCUS_DROPOUT_RATE = 0.016
DEPTH_CV = 0.43
GRID_STEPS = 201
# An MLE host estimate at or below this is treated as a "collapse" (the clinical
# failure mode: reads as no host when there is real residual host).
COLLAPSE_HOST_FRAC = 0.001


def blend_and_estimate(
    host_vcf: Path,
    donor_vcf: Path,
    host_md: list,
    donor_md: list,
    biases: list[float],
    host_frac: float,
    blend_seed: int,
    admix_path: Path,
    n_markers: int,
    depth: int,
    rho: float,
) -> dict:
    """Blend one sample and estimate host fraction under both trim modes."""
    blend = blend_vcfs(
        host_path=str(host_vcf),
        donor_path=str(donor_vcf),
        donor_fraction=1.0 - host_frac,
        target_depth=depth,
        sample_name="admix",
        seed=blend_seed,
        fixed_biases=biases,
        error_rate=ERROR_RATE,
        locus_dropout_rate=LOCUS_DROPOUT_RATE,
        depth_cv=DEPTH_CV,
        rho=rho,
    )
    bias_dict = (
        {(c, p, r, a): b for c, p, r, a, b in blend.marker_biases}
        if blend.marker_biases is not None
        else None
    )
    write_vcf(blend, admix_path)
    admix_md = parse_vcf(str(admix_path), min_dp=0, min_gq=0)

    genos = classify_markers(
        host_md[:n_markers],
        [donor_md[:n_markers]],
        admix_md[:n_markers],
        min_dp=0,
        min_gq=0,
        pass_only=False,
    )
    cal = PanelCalibration(biases=bias_dict)

    def host_est(one_sided: bool) -> tuple[float, float]:
        chimerism.ROBUST_ONE_SIDED = one_sided
        res = estimate_single_donor_bb(
            genos.informative,
            error_rate=ERROR_RATE,
            grid_steps=GRID_STEPS,
            calibration=cal,
            robust="auto",
        )
        return 1.0 - res.donor_fraction, res.robust_drop_fraction

    host_sym, drop_sym = host_est(False)
    host_one, drop_one = host_est(True)

    hp = host_presence_test(genos.informative, error_rate=ERROR_RATE, artifact_filter=False)

    return {
        "host_frac": host_frac,
        "host_sym": host_sym,
        "host_one": host_one,
        "drop_sym": drop_sym,
        "drop_one": drop_one,
        "pres_f": hp.f_host_mle,
        "pres_p": hp.lrt_pval,
    }


def run_rep(
    rep: int, base_seed: int, n_markers: int, depth: int, rho: float, work_dir: Path
) -> list[dict]:
    """One genotype pair (unrelated), all host fractions."""
    pair_dir = work_dir / f"rep{rep}"
    pair_dir.mkdir(parents=True, exist_ok=True)

    gt_rng = random.Random(base_seed * 1000 + rep)
    markers = generate_related_genotypes(n_markers, "unrelated", gt_rng, maf_range=MAF_RANGE)
    host_vcf = pair_dir / "host.vcf"
    donor_vcf = pair_dir / "donor.vcf"
    write_genotype_vcf(markers, host_vcf, "host", key="host_gt")
    write_genotype_vcf(markers, donor_vcf, "donor", key="donor_gt")
    host_md = parse_vcf(str(host_vcf), min_dp=0, min_gq=0)
    donor_md = parse_vcf(str(donor_vcf), min_dp=0, min_gq=0)

    bias_rng = random.Random(base_seed * 7919 + rep)
    biases = generate_marker_biases_realistic(n_markers, bias_rng)

    admix_path = pair_dir / "admix.vcf"
    rows = []
    for host_frac in HOST_FRACTIONS:
        blend_seed = base_seed * 31 + rep * 101 + int(host_frac * 1e5)
        rows.append(
            blend_and_estimate(
                host_vcf,
                donor_vcf,
                host_md,
                donor_md,
                biases,
                host_frac,
                blend_seed,
                admix_path,
                n_markers,
                depth,
                rho,
            )
        )
    return rows


def summarise(rows: list[dict]) -> None:
    """Print per-host-fraction means for both trim modes."""
    by_frac: dict[float, list[dict]] = {}
    for r in rows:
        by_frac.setdefault(r["host_frac"], []).append(r)

    def mean(rs, key):
        return statistics.mean(r[key] for r in rs)

    def collapse_rate(rs, key):
        return sum(1 for r in rs if r[key] <= COLLAPSE_HOST_FRAC) / len(rs)

    print(
        f"\n{'known':>7} | {'sym MLE':>8} {'1-sided':>8} {'presF':>7} | "
        f"{'drop sym':>8} {'drop 1s':>8} | {'collapse sym':>12} {'collapse 1s':>11}"
    )
    print("-" * 92)
    for frac in sorted(by_frac):
        rs = by_frac[frac]
        print(
            f"{frac * 100:>6.2f}% | "
            f"{mean(rs, 'host_sym') * 100:>7.3f}% {mean(rs, 'host_one') * 100:>7.3f}% "
            f"{mean(rs, 'pres_f') * 100:>6.3f}% | "
            f"{mean(rs, 'drop_sym') * 100:>7.1f}% {mean(rs, 'drop_one') * 100:>7.1f}% | "
            f"{collapse_rate(rs, 'host_sym') * 100:>11.0f}% "
            f"{collapse_rate(rs, 'host_one') * 100:>10.0f}%"
        )

    nonblank = [r for r in rows if r["host_frac"] > 0]
    mae_sym = statistics.mean(abs(r["host_sym"] - r["host_frac"]) for r in nonblank)
    mae_one = statistics.mean(abs(r["host_one"] - r["host_frac"]) for r in nonblank)
    print(
        f"\nMAE over non-blank fractions: symmetric {mae_sym * 100:.3f}%, "
        f"one-sided {mae_one * 100:.3f}%"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-reps", type=int, default=20)
    parser.add_argument("--n-markers", type=int, default=200)
    parser.add_argument("--base-seed", type=int, default=2026)
    parser.add_argument(
        "--depth", type=int, default=2000, help="Mean depth (lower = harder; v1-library regime)"
    )
    parser.add_argument(
        "--rho",
        type=float,
        default=float("inf"),
        help="Beta-binomial overdispersion (finite = noisier library)",
    )
    args = parser.parse_args()

    work_dir = Path("output/validate_onesided_trim")
    work_dir.mkdir(parents=True, exist_ok=True)

    print(
        f"depth={args.depth} rho={args.rho} n_markers={args.n_markers} n_reps={args.n_reps}",
        file=sys.stderr,
    )
    rows: list[dict] = []
    for rep in range(args.n_reps):
        rows.extend(run_rep(rep, args.base_seed, args.n_markers, args.depth, args.rho, work_dir))
        print(f"rep {rep + 1}/{args.n_reps} done", file=sys.stderr)

    summarise(rows)


if __name__ == "__main__":
    main()
