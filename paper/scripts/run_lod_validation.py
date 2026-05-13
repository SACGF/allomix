#!/usr/bin/env python3
"""Sweep limit-of-detection (LoD) across panel size, depth, and relatedness.

Follows CLSI EP17-A2:
  - LoB = 95th percentile of estimated donor fraction at true fraction = 0
  - LoD = lowest true fraction at which >=95% of replicates have est_frac > LoB,
          read from a 2-parameter logistic fit P(detected | f) = sigmoid(a + b*log10(f)).

For each (relatedness, n_markers, rep) the host/donor genotypes are generated
once and reused across every (depth, true_frac) cell. The estimator is then run
on a freshly blended admixture VCF per cell.

Outputs:
  output/facts/lod_grid_raw.csv     # one row per replicate per cell
  output/facts/lod_summary.csv      # one row per (relatedness, depth, n_markers)
  output/facts/lod_headline.csv     # single-row snapshot of headline numbers

Usage:
    python paper/scripts/run_lod_validation.py
    python paper/scripts/run_lod_validation.py --n-replicates 10 --n-workers 8
"""

from __future__ import annotations

import argparse
import csv
import math
import random
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from allomix.chimerism import estimate_single_donor_bb  # noqa: E402
from allomix.genotype import classify_markers, parse_vcf  # noqa: E402
from allomix.simulate import (  # noqa: E402
    blend_vcfs,
    generate_marker_biases_realistic,
    generate_related_genotypes,
    write_genotype_vcf,
    write_vcf,
)

# --- Sweep grid (plan section "Sweep design") -------------------------------

RELATEDNESS_LEVELS = ["unrelated", "sibling"]
DEPTHS = [100, 250, 500, 1000, 2000]
N_MARKERS_GRID = [25, 50, 75, 100, 200, 400]
TRUE_FRACTIONS = [0.0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05]
DEFAULT_N_REPLICATES = 60
DEFAULT_N_BOOTSTRAP = 200

MAF_RANGE = (0.2, 0.5)
ERROR_RATE = 0.01
LOCUS_DROPOUT_RATE = 0.016
DEPTH_CV = 0.43
# The estimator uses grid_steps=1001 by default. For LoD characterisation we
# only need detection-rule precision (est_frac > LoB), not MLE precision, so a
# coarser grid is fine and gives a ~5x speedup. Nelder-Mead still refines from
# the grid maximum, so f-estimates remain accurate to <1e-3.
ESTIMATOR_GRID_STEPS = 201

FACTS_DIR = Path("output/facts")
WORK_DIR = Path("output/lod_validation")

# Numerical constants for logistic / probit work.
LOGIT_95 = math.log(0.95 / 0.05)  # 2.944...
# Sentinels for the two edge cells the plan asks us to stop on.
LOD_BELOW_RANGE = -1.0
LOD_ABOVE_RANGE = float("inf")


def derive_seed(*parts: object) -> int:
    """Deterministic seed from a tuple of arbitrary parts.

    Hashing through repr to avoid Python's randomised string hashing.
    """
    h = hash(repr(parts)) & 0xFFFFFFFF
    return h


def detection_rate(est_fracs: list[float], lob: float) -> float:
    """Fraction of replicates whose estimate exceeds LoB."""
    if not est_fracs:
        return 0.0
    return sum(1 for e in est_fracs if e > lob) / len(est_fracs)


def compute_lob(est_fracs_at_zero: list[float]) -> float:
    """LoB = 95th percentile of est_frac across blank replicates."""
    if not est_fracs_at_zero:
        return float("nan")
    return float(np.quantile(np.asarray(est_fracs_at_zero), 0.95))


def _logistic(log10_f: np.ndarray, a: float, b: float) -> np.ndarray:
    """P(detected | f) = 1 / (1 + exp(-(a + b * log10(f))))."""
    return 1.0 / (1.0 + np.exp(-(a + b * log10_f)))


def _interp_lod(
    fractions: list[float], rates: list[float], target: float = 0.95,
) -> float | None:
    """Linear interpolation in log10(f) between the two fractions that bracket
    ``target`` detection rate. Returns None if the rates never cross ``target``.

    This is the fallback when the logistic fit fails (typically when detection
    jumps too sharply from 0 to 1 to identify the logistic's slope).
    """
    pairs = sorted([(f, r) for f, r in zip(fractions, rates) if f > 0])
    for (f_lo, r_lo), (f_hi, r_hi) in zip(pairs, pairs[1:]):
        if r_lo <= target <= r_hi:
            if r_hi == r_lo:
                return f_lo
            log_lo, log_hi = math.log10(f_lo), math.log10(f_hi)
            frac = (target - r_lo) / (r_hi - r_lo)
            return 10.0 ** (log_lo + frac * (log_hi - log_lo))
    return None


def fit_lod(
    fractions: list[float],
    detection_rates: list[float],
    weights: list[int] | None = None,
) -> tuple[float, float, float] | None:
    """Fit logistic in log10(f), solve for f at P = 0.95.

    Args:
        fractions: Non-zero true fractions (parallel with detection_rates).
        detection_rates: Empirical P(detected) at each fraction.
        weights: Optional per-point weights (number of replicates), passed as
            inverse-variance-like sigma to curve_fit.

    Returns:
        (f95, a, b) on success, or None if the fit failed.
    """
    pos = [(f, r) for f, r in zip(fractions, detection_rates) if f > 0]
    if len(pos) < 2:
        return None
    log10_f = np.array([math.log10(f) for f, _ in pos], dtype=float)
    rates = np.array([r for _, r in pos], dtype=float)

    if weights is not None and len(weights) == len(pos):
        # Use binomial SE as sigma: sqrt(p(1-p)/N), floored.
        n = np.array(weights, dtype=float)
        var = np.clip(rates * (1.0 - rates), 1e-3, None) / np.clip(n, 1.0, None)
        sigma = np.sqrt(var)
    else:
        sigma = None

    try:
        popt, _ = curve_fit(
            _logistic,
            log10_f,
            rates,
            p0=[2.0, 2.0],
            sigma=sigma,
            absolute_sigma=False,
            maxfev=10000,
        )
        a, b = float(popt[0]), float(popt[1])
    except (RuntimeError, ValueError):
        a = b = float("nan")

    f95: float | None = None
    if math.isfinite(a) and math.isfinite(b) and abs(b) > 1e-9:
        log10_f95 = (LOGIT_95 - a) / b
        try:
            cand = 10.0**log10_f95
        except OverflowError:
            cand = float("inf")
        if math.isfinite(cand) and cand > 0:
            f95 = cand

    if f95 is None:
        # Fallback: linear interpolation in log10(f) between bracketing points.
        # The logistic cannot identify its slope when detection jumps too
        # sharply from 0 to 1; the interpolated crossing is still meaningful.
        interp = _interp_lod([f for f, _ in pos], [r for _, r in pos])
        if interp is None:
            return None
        f95 = interp
        if not math.isfinite(a) or not math.isfinite(b):
            b = 5.0
            a = LOGIT_95 - b * math.log10(f95)

    return f95, a, b


def bootstrap_lod_ci(
    per_frac_booleans: dict[float, list[bool]],
    n_bootstrap: int,
    rng: random.Random,
) -> tuple[float, float]:
    """Bootstrap CI on LoD by resampling replicate-level detection booleans.

    Resamples the per-replicate booleans at each fraction (with replacement)
    and refits the logistic each time.
    """
    fractions = sorted(f for f in per_frac_booleans if f > 0)
    if len(fractions) < 2:
        return float("nan"), float("nan")

    boots: list[float] = []
    for _ in range(n_bootstrap):
        rates = []
        weights = []
        for f in fractions:
            booleans = per_frac_booleans[f]
            if not booleans:
                rates.append(0.0)
                weights.append(0)
                continue
            sample = [booleans[rng.randrange(len(booleans))] for _ in booleans]
            rates.append(sum(sample) / len(sample))
            weights.append(len(sample))
        fit = fit_lod(fractions, rates, weights)
        if fit is None:
            continue
        boots.append(fit[0])

    if len(boots) < 10:
        return float("nan"), float("nan")
    arr = np.asarray(boots)
    return float(np.quantile(arr, 0.025)), float(np.quantile(arr, 0.975))


# --- Per-replicate work ------------------------------------------------------


def run_replicate(
    relatedness: str,
    n_markers: int,
    rep: int,
    base_seed: int,
    work_root: Path,
) -> list[dict]:
    """Generate one host/donor pair, then run estimator across depth x fraction.

    Returns one row per (depth, true_frac).
    """
    rep_dir = work_root / relatedness / f"nm_{n_markers}" / f"rep_{rep}"
    rep_dir.mkdir(parents=True, exist_ok=True)

    # Genotype RNG is deterministic in (relatedness, n_markers, rep) but not in
    # depth/fraction, so the same pair is reused across all inner cells.
    gt_seed = derive_seed("gt", relatedness, n_markers, rep, base_seed)
    rng = random.Random(gt_seed)
    markers = generate_related_genotypes(n_markers, relatedness, rng, maf_range=MAF_RANGE)
    n_informative_truth = sum(1 for m in markers if m["informative"])

    host_vcf = rep_dir / "host.vcf"
    donor_vcf = rep_dir / "donor.vcf"
    write_genotype_vcf(markers, host_vcf, "host", key="host_gt")
    write_genotype_vcf(markers, donor_vcf, "donor", key="donor_gt")

    host_md = parse_vcf(str(host_vcf), min_dp=0, min_gq=0)
    donor_md = parse_vcf(str(donor_vcf), min_dp=0, min_gq=0)

    # Pre-generate per-marker amplification biases for this replicate. Reusing
    # the same biases across all (depth, fraction) cells in this rep means we
    # can hand them to the estimator as known bias correction values, isolating
    # LoD from bias-calibration uncertainty (the "panel calibrated" ceiling).
    bias_rng = random.Random(derive_seed("bias", relatedness, n_markers, rep, base_seed))
    fixed_biases = generate_marker_biases_realistic(n_markers, bias_rng)

    rows: list[dict] = []
    admix_path = rep_dir / "admix.vcf"  # reused across cells, saves disk
    bias_dict: dict[tuple[str, int, str, str], float] | None = None

    for depth in DEPTHS:
        for frac in TRUE_FRACTIONS:
            blend_seed = derive_seed("blend", relatedness, n_markers, rep, depth, frac, base_seed)
            blend = blend_vcfs(
                host_path=str(host_vcf),
                donor_path=str(donor_vcf),
                donor_fraction=frac,
                target_depth=depth,
                sample_name="admix",
                seed=blend_seed,
                fixed_biases=fixed_biases,
                error_rate=ERROR_RATE,
                locus_dropout_rate=LOCUS_DROPOUT_RATE,
                depth_cv=DEPTH_CV,
            )
            if bias_dict is None and blend.marker_biases is not None:
                bias_dict = {(c, p, r, a): b for c, p, r, a, b in blend.marker_biases}
            write_vcf(blend, admix_path)
            admix_md = parse_vcf(str(admix_path), min_dp=0, min_gq=0)

            genos = classify_markers(host_md, [donor_md], admix_md, min_dp=0, min_gq=0,
                                     pass_only=False)
            if len(genos.informative) < 1:
                rows.append({
                    "relatedness": relatedness,
                    "depth": depth,
                    "n_markers": n_markers,
                    "true_frac": frac,
                    "rep": rep,
                    "seed": blend_seed,
                    "est_frac": float("nan"),
                    "ci_lo": float("nan"),
                    "ci_hi": float("nan"),
                    "n_informative": 0,
                    "n_informative_truth": n_informative_truth,
                })
                continue

            result = estimate_single_donor_bb(
                genos.informative, error_rate=ERROR_RATE,
                grid_steps=ESTIMATOR_GRID_STEPS,
                marker_biases=bias_dict,
            )
            rows.append({
                "relatedness": relatedness,
                "depth": depth,
                "n_markers": n_markers,
                "true_frac": frac,
                "rep": rep,
                "seed": blend_seed,
                "est_frac": result.donor_fraction,
                "ci_lo": result.donor_fraction_ci[0],
                "ci_hi": result.donor_fraction_ci[1],
                "n_informative": result.n_informative,
                "n_informative_truth": n_informative_truth,
            })

    return rows


# --- Cell-level LoB / LoD summarisation -------------------------------------


def summarise_cell(
    relatedness: str,
    depth: int,
    n_markers: int,
    cell_rows: list[dict],
    n_bootstrap: int,
    boot_rng: random.Random,
) -> dict:
    """Compute LoB, LoD, bootstrap CI for one (relatedness, depth, n_markers) cell."""
    by_frac: dict[float, list[float]] = {f: [] for f in TRUE_FRACTIONS}
    n_inf_by_frac: dict[float, list[int]] = {f: [] for f in TRUE_FRACTIONS}
    for r in cell_rows:
        if math.isnan(r["est_frac"]):
            continue
        by_frac[r["true_frac"]].append(r["est_frac"])
        n_inf_by_frac[r["true_frac"]].append(r["n_informative"])

    blanks = by_frac.get(0.0, [])
    lob = compute_lob(blanks)

    fractions = [f for f in TRUE_FRACTIONS if f > 0]
    rates = [detection_rate(by_frac[f], lob) for f in fractions]
    weights = [len(by_frac[f]) for f in fractions]
    booleans_by_frac = {f: [e > lob for e in by_frac[f]] for f in fractions}

    fit = fit_lod(fractions, rates, weights)
    note = ""
    if fit is None:
        # Two pathological cases the plan asks us to stop on.
        if all(r >= 0.95 for r in rates):
            lod = LOD_BELOW_RANGE
            note = "all_detected"
        elif all(r < 0.05 for r in rates):
            lod = LOD_ABOVE_RANGE
            note = "none_detected"
        else:
            lod = float("nan")
            note = "fit_failed"
        lod_ci_lo = lod_ci_hi = float("nan")
        a_fit = b_fit = float("nan")
    else:
        lod, a_fit, b_fit = fit
        lod_ci_lo, lod_ci_hi = bootstrap_lod_ci(booleans_by_frac, n_bootstrap, boot_rng)

    all_n_inf = [n for vals in n_inf_by_frac.values() for n in vals]
    mean_inf = float(np.mean(all_n_inf)) if all_n_inf else 0.0
    median_inf = float(np.median(all_n_inf)) if all_n_inf else 0.0

    return {
        "relatedness": relatedness,
        "depth": depth,
        "n_markers": n_markers,
        "lob_pct": _to_pct(lob),
        "lod_pct": _to_pct(lod),
        "lod_pct_ci_lo": _to_pct(lod_ci_lo),
        "lod_pct_ci_hi": _to_pct(lod_ci_hi),
        "logistic_a": a_fit,
        "logistic_b": b_fit,
        "mean_n_informative": round(mean_inf, 1),
        "median_n_informative": round(median_inf, 1),
        "n_replicates": max(len(v) for v in by_frac.values()) if by_frac else 0,
        "note": note,
        # Per-fraction detection rates as columns for traceability.
        **{f"det_rate_f{int(f * 10000):05d}": round(r, 3) for f, r in zip(fractions, rates)},
    }


def _to_pct(x: float) -> float:
    """Convert fraction to percent; preserve sentinels and NaN."""
    if x == LOD_BELOW_RANGE:
        return -1.0
    if x == LOD_ABOVE_RANGE:
        return float("inf")
    if math.isnan(x):
        return float("nan")
    return round(x * 100, 4)


# --- Output writers ----------------------------------------------------------


def write_grid_raw(rows: list[dict], path: Path) -> None:
    fields = [
        "relatedness", "depth", "n_markers", "true_frac", "rep", "seed",
        "est_frac", "ci_lo", "ci_hi", "n_informative", "n_informative_truth",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})


def write_summary(summaries: list[dict], path: Path) -> None:
    # Compose fieldnames from the first row (deterministic order, det_rate_* at end).
    base = [
        "relatedness", "depth", "n_markers",
        "lob_pct", "lod_pct", "lod_pct_ci_lo", "lod_pct_ci_hi",
        "logistic_a", "logistic_b",
        "mean_n_informative", "median_n_informative",
        "n_replicates", "note",
    ]
    extra = sorted({k for s in summaries for k in s if k.startswith("det_rate_")})
    fields = base + extra
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for s in summaries:
            w.writerow({k: s.get(k, "") for k in fields})


def write_headline(summaries: list[dict], path: Path, n_replicates: int) -> None:
    """One-row CSV of cells the Results prose quotes."""
    def lookup(rel: str, d: int, nm: int) -> float:
        for s in summaries:
            if s["relatedness"] == rel and s["depth"] == d and s["n_markers"] == nm:
                return s["lod_pct"]
        return float("nan")

    headline = {
        "n_replicates": n_replicates,
        "unrelated_lod_1000x_75markers_pct": lookup("unrelated", 1000, 75),
        "unrelated_lod_1000x_100markers_pct": lookup("unrelated", 1000, 100),
        "sibling_lod_1000x_75markers_pct": lookup("sibling", 1000, 75),
        "sibling_lod_1000x_100markers_pct": lookup("sibling", 1000, 100),
        "unrelated_lod_500x_200markers_pct": lookup("unrelated", 500, 200),
        "unrelated_lod_2000x_400markers_pct": lookup("unrelated", 2000, 400),
        "sibling_lod_2000x_400markers_pct": lookup("sibling", 2000, 400),
    }
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(headline.keys()))
        w.writeheader()
        w.writerow({k: (round(v, 4) if isinstance(v, float) and math.isfinite(v) else v)
                    for k, v in headline.items()})


# --- Driver ------------------------------------------------------------------


def _worker(args: tuple) -> list[dict]:
    relatedness, n_markers, rep, base_seed, work_root = args
    return run_replicate(relatedness, n_markers, rep, base_seed, work_root)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--n-replicates", "-n", type=int, default=DEFAULT_N_REPLICATES)
    parser.add_argument("--n-bootstrap", type=int, default=DEFAULT_N_BOOTSTRAP)
    parser.add_argument("--n-workers", "-j", type=int, default=1,
                        help="Process pool size for replicate-level parallelism.")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--relatedness", nargs="+", default=RELATEDNESS_LEVELS,
        help=f"Subset of {RELATEDNESS_LEVELS}.",
    )
    parser.add_argument(
        "--n-markers", type=int, nargs="+", default=N_MARKERS_GRID,
        help=f"Subset of {N_MARKERS_GRID}.",
    )
    parser.add_argument(
        "--depths", type=int, nargs="+", default=None,
        help="Override the depth grid (default uses the full sweep).",
    )
    parser.add_argument("--workdir", default=str(WORK_DIR))
    args = parser.parse_args(argv)

    global DEPTHS
    if args.depths is not None:
        DEPTHS = sorted(args.depths)

    FACTS_DIR.mkdir(parents=True, exist_ok=True)
    work_root = Path(args.workdir)
    work_root.mkdir(parents=True, exist_ok=True)

    tasks = [
        (rel, nm, rep, args.seed, work_root)
        for rel in args.relatedness
        for nm in args.n_markers
        for rep in range(args.n_replicates)
    ]
    print(
        f"Sweep: relatedness={args.relatedness}, depths={DEPTHS}, "
        f"n_markers={args.n_markers}, fractions={TRUE_FRACTIONS}, "
        f"reps={args.n_replicates}",
        file=sys.stderr,
    )
    print(
        f"Total replicate-runs: {len(tasks)} "
        f"(each produces {len(DEPTHS) * len(TRUE_FRACTIONS)} estimator calls)",
        file=sys.stderr,
    )

    all_rows: list[dict] = []
    if args.n_workers <= 1:
        for i, t in enumerate(tasks, 1):
            all_rows.extend(_worker(t))
            if i % 10 == 0 or i == len(tasks):
                print(f"  [{i}/{len(tasks)}] reps done", file=sys.stderr)
    else:
        with ProcessPoolExecutor(max_workers=args.n_workers) as pool:
            futures = [pool.submit(_worker, t) for t in tasks]
            done = 0
            for fut in as_completed(futures):
                all_rows.extend(fut.result())
                done += 1
                if done % 10 == 0 or done == len(futures):
                    print(f"  [{done}/{len(futures)}] reps done", file=sys.stderr)

    # Sort rows for stable output regardless of worker order.
    all_rows.sort(key=lambda r: (
        r["relatedness"], r["n_markers"], r["depth"], r["rep"], r["true_frac"],
    ))

    write_grid_raw(all_rows, FACTS_DIR / "lod_grid_raw.csv")
    print(f"Wrote {FACTS_DIR / 'lod_grid_raw.csv'} ({len(all_rows)} rows)", file=sys.stderr)

    # Cell summaries.
    summaries: list[dict] = []
    boot_rng = random.Random(args.seed + 7919)
    for rel in args.relatedness:
        for depth in DEPTHS:
            for nm in args.n_markers:
                cell_rows = [r for r in all_rows
                             if r["relatedness"] == rel
                             and r["depth"] == depth
                             and r["n_markers"] == nm]
                summaries.append(summarise_cell(rel, depth, nm, cell_rows,
                                                args.n_bootstrap, boot_rng))

    write_summary(summaries, FACTS_DIR / "lod_summary.csv")
    print(f"Wrote {FACTS_DIR / 'lod_summary.csv'} ({len(summaries)} rows)", file=sys.stderr)

    write_headline(summaries, FACTS_DIR / "lod_headline.csv", args.n_replicates)
    print(f"Wrote {FACTS_DIR / 'lod_headline.csv'}", file=sys.stderr)

    # Flag pathological cells for the user to look at.
    flagged = [s for s in summaries if s["note"]]
    if flagged:
        print("\nCells flagged for review (LoD edge cases):", file=sys.stderr)
        for s in flagged:
            print(
                f"  {s['relatedness']:9s} depth={s['depth']:5d} nm={s['n_markers']:4d} "
                f"-> {s['note']} (lob={s['lob_pct']}, lod={s['lod_pct']})",
                file=sys.stderr,
            )

    return 0


if __name__ == "__main__":
    sys.exit(main())
