#!/usr/bin/env python3
"""In-silico check of the sex-aware chrX routing (#46).

Blends synthetic host/donor/admixture trios with the simulator's sex-aware chrX
model (hemizygous male calls written GATK-style as diploid homozygotes, with a
spurious-het rate; copy-number-weighted blending) and runs ``analyse_sample``
under the contig policies, over N independent seeds. Three arms, each at a
matched autosomal marker count:

1. **FF**: a female/female pair with chrX routed in (``SEX_AWARE``) against
   autosomes only (``AUTOSOMES_ONLY``). Expectation: more markers used, no bias
   shift, a tighter interval and a lower LoD in proportion to the added
   informative markers.
2. **MF**: a male/female pair with chrX forced in (``ALL_PRIMARY``) against
   the default gate (``SEX_AWARE``, which excludes chrX for a mismatched pair).
   This documents the bias the gate prevents: the blend follows the true copy
   numbers (male 1, female 2) while the estimator assumes two copies in each.
3. **MM**: a male/male pair with spurious chrX hets, with the het drop
   (``SEX_AWARE``: hom/hom sites only) against every chrX site kept
   (``ALL_PRIMARY``). A het call on a hemizygous chromosome puts a marker with
   the wrong expected VAF into the fit.

Host and donor sex are declared to ``analyse_sample`` in every arm, as a
clinical run would, so the routing demonstration is not confounded by an
ambiguous inference on a small chrX marker count (the inference still runs
and still fails QC on a confident conflict). The per-arm table reports, per
donor fraction and policy, the mean markers used, mean chrX markers used and
dropped, the mean estimate error and its SD across seeds (host percentage
points), the mean CI width, the mean LoD, and the mean robust-refit exclusions.

Usage:
    python scripts/validate_sexchrom_routing.py
    python scripts/validate_sexchrom_routing.py --n-seeds 10 --n-chrx 30 --out output/sexchrom.csv
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
from allomix.genotype import ContigPolicy, parse_vcf  # noqa: E402
from allomix.sex_types import Sex  # noqa: E402
from allomix.simulate import (  # noqa: E402
    build_joint_vcf_from_genotype_dicts,
    generate_related_genotypes,
    generate_sex_chrom_genotypes,
    write_joint_vcf,
)

# Empirical panel noise (see CLAUDE.md / paper/empirical_results): per-marker
# depth CV and bias SD measured on the deployment panel.
DEPTH_CV = 0.43
MARKER_BIAS_SD = 0.018
SEED_MODULUS = 2**31

ARMS = [
    # (label, host sex, donor sex, policy under test, reference policy)
    ("FF", "F", "F", ContigPolicy.SEX_AWARE, ContigPolicy.AUTOSOMES_ONLY),
    ("MF", "M", "F", ContigPolicy.ALL_PRIMARY, ContigPolicy.SEX_AWARE),
    ("MM", "M", "M", ContigPolicy.SEX_AWARE, ContigPolicy.ALL_PRIMARY),
]
SEX_OF = {"F": Sex.FEMALE, "M": Sex.MALE}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--n-seeds", type=int, default=5, help="Independent seeds (default: 5)")
    parser.add_argument("--n-autosomal", type=int, default=60, help="Autosomal markers (60)")
    parser.add_argument("--n-chrx", type=int, default=20, help="Non-PAR chrX markers (20)")
    parser.add_argument("--depth", type=int, default=1000, help="Mean depth per marker (1000)")
    parser.add_argument(
        "--fractions",
        type=float,
        nargs="+",
        default=[0.90, 0.99],
        help="Donor fractions of the admixture (default: 0.90 0.99)",
    )
    parser.add_argument(
        "--spurious-het-rate",
        type=float,
        default=0.04,
        help="Spurious het rate at male non-PAR chrX sites (default: 0.04)",
    )
    parser.add_argument("--seed", type=int, default=7, help="Base seed (default: 7)")
    parser.add_argument("--out", type=Path, default=None, help="Optional CSV of the table")
    return parser.parse_args(argv)


def build_markers(
    host_sex: str,
    donor_sex: str,
    rng: random.Random,
    n_autosomal: int,
    n_chrx: int,
    het_rate: float,
) -> list[dict]:
    autosomal = generate_related_genotypes(n_autosomal, "unrelated", rng)
    for i, m in enumerate(autosomal):
        m["chrom"] = f"chr{(i % 22) + 1}"
        m["pos"] = 1_000_000 + i * 100_000
    chrx = generate_sex_chrom_genotypes(
        n_chrx, host_sex, donor_sex, rng, n_par=0, male_chrx_spurious_het_rate=het_rate
    )
    return [*autosomal, *chrx]


def run_one(
    markers: list[dict],
    fractions: list[float],
    depth: int,
    seed: int,
    policy: ContigPolicy,
    host_sex: str,
    donor_sex: str,
    tmpdir: Path,
) -> list[dict]:
    """Blend one trio, write it, and analyse every admixture under ``policy``."""
    names = [f"ADMIX_F{f:.2f}" for f in fractions]
    joint = build_joint_vcf_from_genotype_dicts(
        markers,
        admix_fractions=fractions,
        admix_sample_names=names,
        target_depth=depth,
        seed=seed,
        error_rate=DEFAULT_ERROR_RATE,
        depth_cv=DEPTH_CV,
        marker_bias_sd=MARKER_BIAS_SD,
    )
    path = tmpdir / f"trio_{seed}_{policy.value}.vcf"
    write_joint_vcf(joint, path)
    host = parse_vcf(path, sample="HOST", min_gq=DEFAULT_MIN_GQ, gt_ad_consistency=True)
    donor = parse_vcf(path, sample="DONOR", min_gq=DEFAULT_MIN_GQ, gt_ad_consistency=True)
    rows = []
    for f, name in zip(fractions, names):
        admix = parse_vcf(path, sample=name, min_dp=0)
        a = analyse_sample(
            host,
            [donor],
            admix,
            min_dp=DEFAULT_MIN_DP,
            min_gq=DEFAULT_MIN_GQ,
            error_rate=DEFAULT_ERROR_RATE,
            contig_policy=policy,
            declared_host_sex=SEX_OF[host_sex],
            declared_donor_sexes=[SEX_OF[donor_sex]],
            robust="auto",
        )
        r = a.result
        lo, hi = r.donor_fraction_ci
        rows.append(
            {
                "fraction": f,
                "est_err_pp": 100 * (r.donor_fraction - f),
                "ci_width_pp": 100 * (hi - lo),
                "lod_pp": 100 * r.lod_fraction,
                "n_used": r.n_markers_used,
                "n_chrx_used": a.genotypes.n_chrx_used,
                "n_het_dropped": a.genotypes.n_chrx_male_het_dropped,
                "n_robust_excluded": r.n_robust_excluded,
                "pair": a.result.sex.pair.value,
            }
        )
    return rows


def summarise(label: str, policy: ContigPolicy, rows: list[dict]) -> list[dict]:
    out = []
    for f in sorted({r["fraction"] for r in rows}):
        sub = [r for r in rows if r["fraction"] == f]
        errs = [r["est_err_pp"] for r in sub]
        out.append(
            {
                "arm": label,
                "policy": policy.value,
                "donor_frac": f,
                "n_seeds": len(sub),
                "mean_n_used": statistics.mean(r["n_used"] for r in sub),
                "mean_n_chrx_used": statistics.mean(r["n_chrx_used"] for r in sub),
                "mean_n_het_dropped": statistics.mean(r["n_het_dropped"] for r in sub),
                "bias_pp": statistics.mean(errs),
                "sd_pp": statistics.stdev(errs) if len(errs) > 1 else 0.0,
                "max_abs_err_pp": max(abs(e) for e in errs),
                "mean_ci_width_pp": statistics.mean(r["ci_width_pp"] for r in sub),
                "mean_lod_pp": statistics.mean(r["lod_pp"] for r in sub),
                "mean_robust_excluded": statistics.mean(r["n_robust_excluded"] for r in sub),
            }
        )
    return out


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    base = random.Random(args.seed)
    table: list[dict] = []
    with tempfile.TemporaryDirectory() as td:
        tmpdir = Path(td)
        for label, host_sex, donor_sex, test_policy, ref_policy in ARMS:
            per_policy: dict[ContigPolicy, list[dict]] = {test_policy: [], ref_policy: []}
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
                # The same blend (same seed) is analysed under both policies, so
                # the two rows differ only in the routing.
                for policy in (test_policy, ref_policy):
                    per_policy[policy].extend(
                        run_one(
                            markers,
                            args.fractions,
                            args.depth,
                            seed,
                            policy,
                            host_sex,
                            donor_sex,
                            tmpdir,
                        )
                    )
            for policy in (test_policy, ref_policy):
                table.extend(summarise(label, policy, per_policy[policy]))

    cols = list(table[0].keys())
    widths = {c: max(len(c), *(len(_fmt(r[c])) for r in table)) for c in cols}
    print("  ".join(c.ljust(widths[c]) for c in cols))
    for r in table:
        print("  ".join(_fmt(r[c]).ljust(widths[c]) for c in cols))
    print(
        f"\n{args.n_seeds} seeds, {args.n_autosomal} autosomal + {args.n_chrx} non-PAR chrX "
        f"markers, depth {args.depth} (CV {DEPTH_CV}), marker bias SD {MARKER_BIAS_SD}, "
        f"male spurious het rate {args.spurious_het_rate}. Errors in host percentage points "
        "(estimate minus truth)."
    )
    if args.out is not None:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        with open(args.out, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(table)
        print(f"Wrote {args.out}")
    return 0


def _fmt(v) -> str:
    if isinstance(v, float):
        return f"{v:.3f}"
    return str(v)


if __name__ == "__main__":
    sys.exit(main())
