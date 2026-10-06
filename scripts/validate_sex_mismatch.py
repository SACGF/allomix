#!/usr/bin/env python3
"""In-silico check of the sex-mismatch cross-check (#48).

Blends synthetic host/donor/admixture trios with the simulator's sex-aware
chrX model (hemizygous male calls written GATK-style as diploid homozygotes
with a spurious-het rate; copy-number-weighted blending) for a sex-mismatched
pair in both directions (host M into donor F, host F into donor M) over N
independent seeds, runs ``analyse_sample`` with the sexes declared, and
compares the chrX copy-number-weighted estimate of the donor fraction
(``result.sex_mismatch``) with the truth and with the autosomal MLE at 1%, 5%
and 10% host.

Per arm and host fraction the table reports the mean usable chrX markers and
male het drops, the bias and SD of the chrX estimate against truth (host
percentage points), its maximum absolute error, its mean CI width, the same
for the autosomal MLE, the mean chrX-minus-MLE gap, and the fraction of seeds
flagged concordant. The sexes are declared so that a male with two spurious
hets in a short chrX marker list does not come out ambiguous and disable the
check (the inference still runs and would still fail QC on a confident
conflict). The chrY depth readout has no simulator support and is not
exercised here.

Usage:
    python scripts/validate_sex_mismatch.py
    python scripts/validate_sex_mismatch.py --n-seeds 10 --n-chrx 30 --out output/sex_mismatch.csv
"""

import argparse
import csv
import random
import statistics
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from allomix.analysis import analyse_sample  # noqa: E402
from allomix.constants import DEFAULT_ERROR_RATE, DEFAULT_MIN_DP, DEFAULT_MIN_GQ  # noqa: E402
from allomix.genotype import parse_vcf  # noqa: E402
from allomix.sex_types import Sex  # noqa: E402
from allomix.simulate import (  # noqa: E402
    build_joint_vcf_from_genotype_dicts,
    generate_related_genotypes,
    generate_sex_chrom_genotypes,
    write_joint_vcf,
)

# Empirical panel noise (see CLAUDE.md / paper/empirical_results).
DEPTH_CV = 0.43
MARKER_BIAS_SD = 0.018
SEED_MODULUS = 2**31

# (label, host sex, donor sex)
ARMS = [("MF", "M", "F"), ("FM", "F", "M")]
SEX_OF = {"F": Sex.FEMALE, "M": Sex.MALE}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--n-seeds", type=int, default=5, help="Independent seeds (default: 5)")
    parser.add_argument("--n-autosomal", type=int, default=60, help="Autosomal markers (60)")
    parser.add_argument("--n-chrx", type=int, default=20, help="Non-PAR chrX markers (20)")
    parser.add_argument("--depth", type=int, default=1000, help="Mean depth per marker (1000)")
    parser.add_argument(
        "--host-fractions",
        type=float,
        nargs="+",
        default=[0.01, 0.05, 0.10],
        help="Host fractions of the admixture (default: 0.01 0.05 0.10)",
    )
    parser.add_argument(
        "--spurious-het-rate",
        type=float,
        default=0.04,
        help="Spurious het rate at male non-PAR chrX sites (default: 0.04)",
    )
    parser.add_argument("--seed", type=int, default=11, help="Base seed (default: 11)")
    parser.add_argument("--out", type=Path, default=None, help="Optional CSV of the table")
    return parser.parse_args(argv)


def build_markers(
    host_sex: str, donor_sex: str, rng: random.Random, n_autosomal: int, n_chrx: int, het: float
) -> list[dict]:
    autosomal = generate_related_genotypes(n_autosomal, "unrelated", rng)
    for i, m in enumerate(autosomal):
        m["chrom"] = f"chr{(i % 22) + 1}"
        m["pos"] = 1_000_000 + i * 100_000
    chrx = generate_sex_chrom_genotypes(
        n_chrx, host_sex, donor_sex, rng, n_par=0, male_chrx_spurious_het_rate=het
    )
    return [*autosomal, *chrx]


def run_one(
    markers: list[dict],
    donor_fractions: list[float],
    depth: int,
    seed: int,
    host_sex: str,
    donor_sex: str,
    tmpdir: Path,
) -> list[dict]:
    """Blend one trio, write it, and analyse every admixture column."""
    names = [f"ADMIX_F{f:.2f}" for f in donor_fractions]
    joint = build_joint_vcf_from_genotype_dicts(
        markers,
        admix_fractions=donor_fractions,
        admix_sample_names=names,
        target_depth=depth,
        seed=seed,
        error_rate=DEFAULT_ERROR_RATE,
        depth_cv=DEPTH_CV,
        marker_bias_sd=MARKER_BIAS_SD,
    )
    path = tmpdir / f"trio_{host_sex}{donor_sex}_{seed}.vcf"
    write_joint_vcf(joint, path)
    host = parse_vcf(path, sample="HOST", min_gq=DEFAULT_MIN_GQ, gt_ad_consistency=True)
    donor = parse_vcf(path, sample="DONOR", min_gq=DEFAULT_MIN_GQ, gt_ad_consistency=True)
    rows = []
    for f, name in zip(donor_fractions, names):
        admix = parse_vcf(path, sample=name, min_dp=0)
        a = analyse_sample(
            host,
            [donor],
            admix,
            min_dp=DEFAULT_MIN_DP,
            min_gq=DEFAULT_MIN_GQ,
            error_rate=DEFAULT_ERROR_RATE,
            declared_host_sex=SEX_OF[host_sex],
            declared_donor_sexes=[SEX_OF[donor_sex]],
            robust="auto",
        )
        r = a.result
        sm = r.sex_mismatch
        if sm is None or sm.basis is None:
            sys.stderr.write(
                f"seed {seed} {name}: cross-check did not run "
                f"(pair {r.sex.pair.value}, chrx_n={getattr(sm, 'chrx_n', 'NA')})\n"
            )
            continue
        lo, hi = r.donor_fraction_ci
        rows.append(
            {
                "host_frac": round(1.0 - f, 4),
                "chrx_n": sm.chrx_n,
                "chrx_het_dropped": sm.chrx_n_male_het_dropped,
                "sex_err_pp": 100 * ((1.0 - sm.frac_donor) - (1.0 - f)),
                "sex_ci_width_pp": 100 * (sm.ci_high - sm.ci_low),
                "mle_err_pp": 100 * ((1.0 - r.donor_fraction) - (1.0 - f)),
                "mle_ci_width_pp": 100 * (hi - lo),
                "gap_pp": 100 * (r.donor_fraction - sm.frac_donor),
                "concordant": bool(sm.concordant),
                "qc": a.qc.status,
            }
        )
    return rows


def summarise(label: str, rows: list[dict]) -> list[dict]:
    out = []
    for hf in sorted({r["host_frac"] for r in rows}):
        sub = [r for r in rows if r["host_frac"] == hf]
        se = [r["sex_err_pp"] for r in sub]
        me = [r["mle_err_pp"] for r in sub]
        out.append(
            {
                "arm": label,
                "host_pct": 100 * hf,
                "n_seeds": len(sub),
                "mean_chrx_n": statistics.mean(r["chrx_n"] for r in sub),
                "mean_het_dropped": statistics.mean(r["chrx_het_dropped"] for r in sub),
                "sex_bias_pp": statistics.mean(se),
                "sex_sd_pp": statistics.stdev(se) if len(se) > 1 else 0.0,
                "sex_max_abs_err_pp": max(abs(e) for e in se),
                "sex_ci_width_pp": statistics.mean(r["sex_ci_width_pp"] for r in sub),
                "mle_bias_pp": statistics.mean(me),
                "mle_sd_pp": statistics.stdev(me) if len(me) > 1 else 0.0,
                "mle_ci_width_pp": statistics.mean(r["mle_ci_width_pp"] for r in sub),
                "mean_gap_pp": statistics.mean(r["gap_pp"] for r in sub),
                "concordant_frac": statistics.mean(1.0 if r["concordant"] else 0.0 for r in sub),
                "n_review": sum(1 for r in sub if r["qc"] == "REVIEW"),
            }
        )
    return out


def _fmt(v) -> str:
    if isinstance(v, float):
        return f"{v:.3f}"
    return str(v)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    base = random.Random(args.seed)
    donor_fractions = [1.0 - h for h in args.host_fractions]
    table: list[dict] = []
    with tempfile.TemporaryDirectory() as td:
        tmpdir = Path(td)
        for label, host_sex, donor_sex in ARMS:
            rows: list[dict] = []
            for _ in range(args.n_seeds):
                seed = base.randint(0, SEED_MODULUS)
                markers = build_markers(
                    host_sex,
                    donor_sex,
                    random.Random(seed),
                    args.n_autosomal,
                    args.n_chrx,
                    args.spurious_het_rate,
                )
                rows.extend(
                    run_one(markers, donor_fractions, args.depth, seed, host_sex, donor_sex, tmpdir)
                )
            table.extend(summarise(label, rows))

    cols = list(table[0].keys())
    widths = {c: max(len(c), *(len(_fmt(r[c])) for r in table)) for c in cols}
    print("  ".join(c.ljust(widths[c]) for c in cols))
    for r in table:
        print("  ".join(_fmt(r[c]).ljust(widths[c]) for c in cols))
    print(
        f"\n{args.n_seeds} seeds, {args.n_autosomal} autosomal + {args.n_chrx} non-PAR chrX "
        f"markers, depth {args.depth} (CV {DEPTH_CV}), marker bias SD {MARKER_BIAS_SD}, "
        f"male spurious het rate {args.spurious_het_rate}. Errors in host percentage points "
        "(estimate minus truth); gap is autosomal MLE minus chrX estimate (donor pp)."
    )
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with open(args.out, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(table)
        print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
