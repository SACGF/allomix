#!/usr/bin/env python3
"""Simulated presence-test LoD sweep: companion to the MLE Figure 1 (issue #25).

``run_lod_validation.py`` sweeps the magnitude (MLE) estimator's limit of
detection across panel size, depth, and relatedness and produces Figure 1
(``fig5_lod_curves.png``). That figure is MLE-only. This script produces the
parallel curve for the host-presence detector (``allomix.host_presence.host_presence_test``)
so the simulated side is symmetric with the real-data subsample figures (issue
#24), which already treat the MLE magnitude estimate and the presence test as
separate, complementary tests.

It reuses the same two-level pair design and helpers as ``run_lod_validation.py``:

  - K donor/host PAIRS per relatedness (genotypes fixed per pair, markers nested
    so the n_markers=50 panel is a strict prefix of n_markers=400).
  - M SEQUENCING replicates per pair (only the blend seed varies).
  - Per-pair LoD from a logistic fit of detection rate vs log10(fraction); the
    curve is the MEDIAN LoD across pairs and the band is the 10th-90th percentile.

Two deliberate differences from the MLE sweep:

  1. Detection rule. The presence test controls its own type-I error against the
     sequencing-error background, so EP17's "limit of blank" role is filled by the
     test's null: a cell is detected when ``lrt_pval < 0.05``. No frac=0 blanks or
     LoB are needed (this is the blank-free per-sample rule issue #24 also uses),
     so the fraction grid here is positive-only.
  2. Per-marker bias is OFF (``realistic_biases=False, marker_bias_sd=0.0``). The
     simulator adds bias to the expected VAF before drawing reads, so at a
     donor-homozygous marker (where vaf is 0 or 1) bias injects a real per-marker
     ALT/REF excess on top of the e/3 error background and miscalibrates the
     presence null (see ``claude/20_host_presence_detection_plan.md``). The MLE
     figure keeps realistic bias; the two tests therefore need different blends
     and cannot share one pass.
     depth_cv and locus_dropout are kept (they do not add background).

Overdispersion: like Figure 1 this runs binomial (``rho=inf``) by default. A
finite ``--rho`` must be paired with ``--rho-marker-type het_only`` so the
donor-absent background stays binomial; a uniform ``rho`` over all markers would
overdisperse the clean background and miscalibrate the null (same objection as the
bias one above).

Outputs (under output/facts/):
  presence_lod_grid_raw.csv         # one row per (pair, seq_rep) per cell
  presence_lod_curve_summary.csv    # per (relatedness, depth, n_markers): median + band
  presence_lod_curve_headline.csv   # named facts for later prose

Usage:
    python paper/scripts/run_presence_lod_validation.py
    python paper/scripts/run_presence_lod_validation.py --n-pairs 5 --n-seq-reps 10 -j 8
"""

import argparse
import csv
import math
import random
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from allomix.genotype import classify_markers, parse_vcf  # noqa: E402
from allomix.host_presence import host_presence_test  # noqa: E402
from allomix.simulate import (  # noqa: E402
    blend_vcfs,
    generate_related_genotypes,
    write_genotype_vcf,
    write_vcf,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paper_quick import qval  # noqa: E402

# Reuse the MLE sweep's grids, two-level defaults, and fit/seed helpers so the
# presence curve is drawn on the same axes as Figure 1 (do not re-derive).
from run_lod_validation import (  # noqa: E402
    BAND_HI_Q,
    BAND_LO_Q,
    DEFAULT_N_PAIRS,
    DEFAULT_N_SEQ_REPS,
    DEPTH_CV,
    DEPTHS,
    LOCUS_DROPOUT_RATE,
    MAF_RANGE,
    N_MARKERS_GRID,
    RELATEDNESS_LEVELS,
    derive_seed,
    fit_lod,
)

# Positive-only fraction grid (no 0.0 blank: the presence null is the test's own
# error background). Matches the positive fractions of the MLE sweep.
POSITIVE_FRACTIONS = qval(
    [0.001, 0.002, 0.005, 0.01, 0.02, 0.05],
    [0.005, 0.01, 0.02, 0.05],
)
ERROR_RATE = 0.01  # symmetric sequencing error; per-direction floor e/3 at donor-hom markers
ALPHA = 0.05  # presence-test LRT significance threshold (detection rule)

# Sentinels for per-pair LoD edge cases (mirror run_lod_validation.py).
LOD_BELOW_RANGE = -1.0
LOD_ABOVE_RANGE = float("inf")
MIN_POS_FRACTION = min(POSITIVE_FRACTIONS)
MAX_FRACTION = max(POSITIVE_FRACTIONS)

FACTS_DIR = Path("output/facts")
WORK_DIR = Path("output/presence_lod_validation")


def run_pair(
    relatedness: str,
    pair_idx: int,
    n_seq_reps: int,
    base_seed: int,
    work_root: Path,
    rho: float,
    rho_marker_type: str,
) -> list[dict]:
    """One fixed host/donor pair, then ``n_seq_reps`` sequencing replicates across
    every (n_markers, depth, fraction) cell, run through the presence test.

    Genotypes are drawn once per pair (seeded by ``pair_idx``). The admix VCF for
    a given (seq_rep, depth, frac) is generated once at max_markers and sliced per
    panel size so the n_markers grid is strictly nested (a smaller panel is a
    prefix of a larger one). Returns one row per (n_markers, depth, frac, seq_rep).
    """
    pair_dir = work_root / relatedness / f"pair_{pair_idx}"
    pair_dir.mkdir(parents=True, exist_ok=True)

    max_markers = max(N_MARKERS_GRID)
    gt_seed = derive_seed("pgt", relatedness, pair_idx, base_seed)
    rng = random.Random(gt_seed)
    all_markers = generate_related_genotypes(max_markers, relatedness, rng, maf_range=MAF_RANGE)

    host_vcf = pair_dir / "host.vcf"
    donor_vcf = pair_dir / "donor.vcf"
    write_genotype_vcf(all_markers, host_vcf, "host", key="host_gt")
    write_genotype_vcf(all_markers, donor_vcf, "donor", key="donor_gt")
    host_md_full = parse_vcf(str(host_vcf), min_dp=0, min_gq=0)
    donor_md_full = parse_vcf(str(donor_vcf), min_dp=0, min_gq=0)

    rows: list[dict] = []
    admix_path = pair_dir / "admix.vcf"

    for seq_rep in range(n_seq_reps):
        for depth in DEPTHS:
            for frac in POSITIVE_FRACTIONS:
                # Host is the minor (residual) contributor: donor_fraction = 1 - f.
                blend_seed = derive_seed(
                    "pblend", relatedness, pair_idx, seq_rep, depth, frac, base_seed
                )
                blend = blend_vcfs(
                    host_path=str(host_vcf),
                    donor_path=str(donor_vcf),
                    donor_fraction=1.0 - frac,
                    target_depth=depth,
                    sample_name="admix",
                    seed=blend_seed,
                    error_rate=ERROR_RATE,
                    locus_dropout_rate=LOCUS_DROPOUT_RATE,
                    depth_cv=DEPTH_CV,
                    realistic_biases=False,  # bias OFF: keep the presence null calibrated
                    marker_bias_sd=0.0,
                    rho=rho,
                    rho_marker_type=rho_marker_type,
                )
                write_vcf(blend, admix_path)
                admix_md_full = parse_vcf(str(admix_path), min_dp=0, min_gq=0)

                for n_markers in N_MARKERS_GRID:
                    host_md = host_md_full[:n_markers]
                    donor_md = donor_md_full[:n_markers]
                    admix_md = admix_md_full[:n_markers]

                    genos = classify_markers(
                        host_md,
                        [donor_md],
                        admix_md,
                        min_dp=0,
                        min_gq=0,
                        pass_only=False,
                    )
                    pres = host_presence_test(genos.informative, error_rate=ERROR_RATE)
                    rows.append(
                        {
                            "relatedness": relatedness,
                            "depth": depth,
                            "n_markers": n_markers,
                            "frac": frac,
                            "pair": pair_idx,
                            "seq_rep": seq_rep,
                            "seed": blend_seed,
                            "lrt_pval": pres.lrt_pval,
                            "f_host_mle": pres.f_host_mle,
                            "n_presence_markers": pres.n_markers,
                            "detected": int(pres.lrt_pval < ALPHA),
                        }
                    )

    return rows


def compute_pair_lod(cell_rows: list[dict]) -> dict:
    """LoD for ONE pair at one (relatedness, depth, n_markers) cell.

    Detection rate per fraction over the pair's sequencing replicates, fed to the
    shared logistic fit. Returns ``lod`` (fraction or sentinel), a ``note``, and
    the mean number of donor-homozygous markers the test ran on.
    """
    by_frac: dict[float, list[int]] = defaultdict(list)
    n_pres: list[int] = []
    for r in cell_rows:
        by_frac[r["frac"]].append(r["detected"])
        n_pres.append(r["n_presence_markers"])

    fractions = sorted(by_frac)
    rates = [float(np.mean(by_frac[f])) for f in fractions]
    weights = [len(by_frac[f]) for f in fractions]

    fit = fit_lod(fractions, rates, weights)
    note = ""
    if fit is None:
        if all(r >= 0.95 for r in rates):
            lod, note = LOD_BELOW_RANGE, "all_detected"
        elif all(r < 0.05 for r in rates):
            lod, note = LOD_ABOVE_RANGE, "none_detected"
        else:
            lod, note = float("nan"), "fit_failed"
    else:
        lod = fit[0]

    return {
        "lod": lod,
        "note": note,
        "mean_n_presence_markers": float(np.mean(n_pres)) if n_pres else 0.0,
    }


def _lod_for_aggregation(lod: float) -> float | None:
    """Map a per-pair LoD (possibly a sentinel) to a finite fraction, or None."""
    if lod == LOD_BELOW_RANGE:
        return MIN_POS_FRACTION
    if lod == LOD_ABOVE_RANGE:
        return MAX_FRACTION
    if lod is None or math.isnan(lod):
        return None
    return lod


def _to_pct(x: float) -> float:
    """Fraction -> percent, preserving sentinels and NaN."""
    if x == LOD_BELOW_RANGE:
        return -1.0
    if x == LOD_ABOVE_RANGE:
        return float("inf")
    if x is None or math.isnan(x):
        return float("nan")
    return round(x * 100, 4)


def summarise_cell(
    relatedness: str,
    depth: int,
    n_markers: int,
    pair_summaries: list[dict],
    n_seq_reps: int,
) -> dict:
    """Aggregate per-pair LoDs into a median curve point and a 10-90% band.

    Columns mirror ``lod_summary.csv`` (relatedness, depth, n_markers, lod_pct,
    band) so the figure reads on the same axes as Figure 1.
    """
    lods = [v for ps in pair_summaries if (v := _lod_for_aggregation(ps["lod"])) is not None]
    n_pres = [ps["mean_n_presence_markers"] for ps in pair_summaries]
    n_dropped = len(pair_summaries) - len(lods)

    if lods:
        arr = np.asarray(lods, dtype=float)
        lod_med = float(np.median(arr))
        lod_lo = float(np.quantile(arr, BAND_LO_Q))
        lod_hi = float(np.quantile(arr, BAND_HI_Q))
    else:
        lod_med = lod_lo = lod_hi = float("nan")

    note = ""
    if not lods:
        note = "no_pairs_fit"
    elif n_dropped:
        note = f"{n_dropped}_pairs_dropped"

    return {
        "relatedness": relatedness,
        "depth": depth,
        "n_markers": n_markers,
        "lod_pct": _to_pct(lod_med),
        "lod_pct_ci_lo": _to_pct(lod_lo),
        "lod_pct_ci_hi": _to_pct(lod_hi),
        "mean_n_presence_markers": round(float(np.mean(n_pres)), 1) if n_pres else 0.0,
        "median_n_presence_markers": round(float(np.median(n_pres)), 1) if n_pres else 0.0,
        "n_pairs": len(pair_summaries),
        "n_pairs_used": len(lods),
        "n_pairs_dropped": n_dropped,
        "n_seq_reps": n_seq_reps,
        "note": note,
    }


def write_grid_raw(rows: list[dict], path: Path) -> None:
    fields = [
        "relatedness",
        "depth",
        "n_markers",
        "frac",
        "pair",
        "seq_rep",
        "seed",
        "lrt_pval",
        "f_host_mle",
        "n_presence_markers",
        "detected",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})


def write_summary(summaries: list[dict], path: Path) -> None:
    fields = [
        "relatedness",
        "depth",
        "n_markers",
        "lod_pct",
        "lod_pct_ci_lo",
        "lod_pct_ci_hi",
        "mean_n_presence_markers",
        "median_n_presence_markers",
        "n_pairs",
        "n_pairs_used",
        "n_pairs_dropped",
        "n_seq_reps",
        "note",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for s in summaries:
            w.writerow({k: s.get(k, "") for k in fields})


def write_headline(
    summaries: list[dict], path: Path, n_pairs: dict[str, int], n_seq_reps: int
) -> None:
    """One-row CSV of cells later prose may quote (parallels lod_headline.csv)."""

    def lookup(rel: str, d: int, nm: int) -> float:
        for s in summaries:
            if s["relatedness"] == rel and s["depth"] == d and s["n_markers"] == nm:
                return s["lod_pct"]
        return float("nan")

    headline = {
        "n_pairs_unrelated": n_pairs.get("unrelated", 0),
        "n_pairs_sibling": n_pairs.get("sibling", 0),
        "n_seq_reps": n_seq_reps,
        "presence_unrelated_lod_1000x_100markers_pct": lookup("unrelated", 1000, 100),
        "presence_sibling_lod_1000x_100markers_pct": lookup("sibling", 1000, 100),
        "presence_unrelated_lod_2000x_400markers_pct": lookup("unrelated", 2000, 400),
        "presence_sibling_lod_2000x_400markers_pct": lookup("sibling", 2000, 400),
    }
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(headline.keys()))
        w.writeheader()
        w.writerow(
            {
                k: (round(v, 4) if isinstance(v, float) and math.isfinite(v) else v)
                for k, v in headline.items()
            }
        )


def _worker(arg: tuple) -> list[dict]:
    rel, pair_idx, n_seq_reps, base_seed, work_root, rho, rho_marker_type = arg
    return run_pair(rel, pair_idx, n_seq_reps, base_seed, work_root, rho, rho_marker_type)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--n-pairs",
        type=int,
        default=None,
        help=f"Donor/host pairs per relatedness. Overrides the defaults {DEFAULT_N_PAIRS}.",
    )
    parser.add_argument(
        "--n-seq-reps", type=int, default=DEFAULT_N_SEQ_REPS, help="Sequencing replicates per pair."
    )
    parser.add_argument("--n-workers", "-j", type=int, default=1)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--relatedness",
        nargs="+",
        default=RELATEDNESS_LEVELS,
        help=f"Subset of {RELATEDNESS_LEVELS}.",
    )
    parser.add_argument("--workdir", default=str(WORK_DIR))
    parser.add_argument(
        "--rho",
        type=float,
        default=float("inf"),
        help="Beta-binomial overdispersion (inf = binomial, default, matches Figure 1).",
    )
    parser.add_argument(
        "--rho-marker-type",
        choices=["all", "het_only"],
        default="all",
        help="Where to apply rho. With a finite --rho use 'het_only' to keep the "
        "donor-absent background binomial; 'all' would miscalibrate the null.",
    )
    args = parser.parse_args(argv)

    if math.isfinite(args.rho) and args.rho_marker_type == "all":
        print(
            "WARNING: finite --rho with --rho-marker-type all overdisperses the "
            "clean donor-absent background and miscalibrates the presence null. "
            "Use --rho-marker-type het_only.",
            file=sys.stderr,
        )

    FACTS_DIR.mkdir(parents=True, exist_ok=True)
    work_root = Path(args.workdir)
    work_root.mkdir(parents=True, exist_ok=True)

    n_pairs = {
        rel: (args.n_pairs if args.n_pairs is not None else DEFAULT_N_PAIRS.get(rel, 10))
        for rel in args.relatedness
    }

    tasks = [
        (rel, pair_idx, args.n_seq_reps, args.seed, work_root, args.rho, args.rho_marker_type)
        for rel in args.relatedness
        for pair_idx in range(n_pairs[rel])
    ]
    print(
        f"Presence sweep: relatedness={args.relatedness}, depths={DEPTHS}, "
        f"n_markers={N_MARKERS_GRID}, fractions={POSITIVE_FRACTIONS}, "
        f"pairs={n_pairs}, seq_reps={args.n_seq_reps}, rho={args.rho}",
        file=sys.stderr,
    )

    all_rows: list[dict] = []
    if args.n_workers <= 1:
        for i, t in enumerate(tasks, 1):
            all_rows.extend(_worker(t))
            if i % 5 == 0 or i == len(tasks):
                print(f"  [{i}/{len(tasks)}] tasks done", file=sys.stderr)
    else:
        with ProcessPoolExecutor(max_workers=args.n_workers) as pool:
            futures = [pool.submit(_worker, t) for t in tasks]
            done = 0
            for fut in as_completed(futures):
                all_rows.extend(fut.result())
                done += 1
                if done % 5 == 0 or done == len(futures):
                    print(f"  [{done}/{len(futures)}] tasks done", file=sys.stderr)

    all_rows.sort(
        key=lambda r: (
            r["relatedness"],
            r["n_markers"],
            r["depth"],
            r["pair"],
            r["seq_rep"],
            r["frac"],
        )
    )

    write_grid_raw(all_rows, FACTS_DIR / "presence_lod_grid_raw.csv")
    print(
        f"Wrote {FACTS_DIR / 'presence_lod_grid_raw.csv'} ({len(all_rows)} rows)", file=sys.stderr
    )

    by_pair_cell: dict[tuple, list[dict]] = defaultdict(list)
    for r in all_rows:
        by_pair_cell[(r["relatedness"], r["depth"], r["n_markers"], r["pair"])].append(r)

    pairs_by_cell: dict[tuple, list[dict]] = defaultdict(list)
    for (rel, depth, nm, _pair), rows in by_pair_cell.items():
        pairs_by_cell[(rel, depth, nm)].append(compute_pair_lod(rows))

    summaries: list[dict] = []
    for rel in args.relatedness:
        for depth in DEPTHS:
            for nm in N_MARKERS_GRID:
                ps_list = pairs_by_cell.get((rel, depth, nm), [])
                if not ps_list:
                    continue
                summaries.append(summarise_cell(rel, depth, nm, ps_list, args.n_seq_reps))

    write_summary(summaries, FACTS_DIR / "presence_lod_curve_summary.csv")
    print(
        f"Wrote {FACTS_DIR / 'presence_lod_curve_summary.csv'} ({len(summaries)} rows)",
        file=sys.stderr,
    )

    write_headline(
        summaries, FACTS_DIR / "presence_lod_curve_headline.csv", n_pairs, args.n_seq_reps
    )
    print(f"Wrote {FACTS_DIR / 'presence_lod_curve_headline.csv'}", file=sys.stderr)

    flagged = [s for s in summaries if s["note"]]
    if flagged:
        print("\nCells flagged for review (LoD edge cases):", file=sys.stderr)
        for s in flagged:
            print(
                f"  {s['relatedness']:9s} depth={s['depth']:5d} nm={s['n_markers']:4d} "
                f"-> {s['note']} (lod={s['lod_pct']})",
                file=sys.stderr,
            )

    return 0


if __name__ == "__main__":
    sys.exit(main())
