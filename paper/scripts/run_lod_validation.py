#!/usr/bin/env python3
"""Sweep limit-of-detection (LoD) across panel size, depth, and relatedness.

Follows CLSI EP17-A2:
  - LoB = mean + 1.645 * SD of estimated donor fraction at true fraction = 0
          (parametric form for approximately-normal blanks; lower SE than the
          empirical 95th percentile at fixed n).
  - LoD = lowest true fraction at which >=95% of replicates have est_frac > LoB,
          read from a 2-parameter logistic fit P(detected | f) = sigmoid(a + b*log10(f)).

Two-level replicate design. Pooling fresh-pair + fresh-draw replicates per cell
let the sibling pair-to-pair variation in IBD sharing (how many markers are
informative) leak into the pooled LoB and logistic fit, making the
LoD-vs-panel-size curve non-monotone. We separate the two variation sources:

  - K donor/host PAIRS per relatedness (genotypes + per-marker capture biases
    fixed per pair, reused across every depth/panel/fraction cell). Markers are
    nested, so n_markers=50 is a strict prefix of n_markers=400: adding markers
    to a pair can only add informative ones.
  - M SEQUENCING replicates per pair (only the blend seed varies), so within a
    pair the sole variation is sequencing/sampling noise.

LoB + LoD are computed per pair (over its M reps); the curve is the MEDIAN LoD
across the K pairs and the band is the 10th-90th percentile across pairs. Each
pair's LoD curve is monotone in panel size, so the median is too; the band is
the IBD-driven spread (narrowing as markers are added). Unrelated pairs vary
little pair-to-pair, so they run with fewer pairs.

Outputs:
  output/facts/lod_grid_raw.csv     # one row per (pair, seq_rep) per cell
  output/facts/lod_per_pair.csv     # one row per (relatedness, depth, n_markers, pair)
  output/facts/lod_summary.csv      # one row per (relatedness, depth, n_markers)
  output/facts/lod_headline.csv     # single-row snapshot of headline numbers

Usage:
    python paper/scripts/run_lod_validation.py
    python paper/scripts/run_lod_validation.py --n-pairs 5 --n-seq-reps 10 --n-workers 8

    # Opt in to the fast vectorized grid estimator (~6.5x faster, max lod_pct
    # deviation 0.0011 pp vs the exact default; see estimate_single_donor_bb_grid).
    # Use for fast iteration; OMIT for the final publication run.
    python paper/scripts/run_lod_validation.py --fast-grid
    ALLOMIX_FAST_GRID=1 python paper/scripts/run_lod_validation.py

    # Composes with quick mode (ALLOMIX_PAPER_QUICK / --config quick=1): quick
    # shrinks the grid, fast-grid swaps the estimator.
    ALLOMIX_PAPER_QUICK=1 ALLOMIX_FAST_GRID=1 python paper/scripts/run_lod_validation.py
"""

import argparse
import csv
import hashlib
import math
import os
import random
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from allomix.estimate.chimerism import estimate_single_donor_bb  # noqa: E402
from allomix.genotype import classify_markers, parse_vcf  # noqa: E402
from allomix.estimate.likelihood import PanelCalibration  # noqa: E402
from allomix.simulate import (  # noqa: E402
    blend_vcfs,
    generate_marker_biases_realistic,
    generate_related_genotypes,
    write_genotype_vcf,
    write_vcf,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from fast_grid_estimator import estimate_single_donor_bb_grid  # noqa: E402
from paper_quick import qval  # noqa: E402

# Sweep grid (plan section "Sweep design").
# Quick-build mode (ALLOMIX_PAPER_QUICK=1) shrinks every loop, turning the
# slowest rule from ~tens of minutes into ~a minute at the cost of precision
# (the resulting figures are watermarked, not for publication).

RELATEDNESS_LEVELS = ["unrelated", "sibling"]
DEPTHS = qval([100, 250, 500, 1000, 2000], [100, 1000])
N_MARKERS_GRID = qval([25, 50, 75, 100, 200, 400], [50, 100, 200])
TRUE_FRACTIONS = qval(
    [0.0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05],
    [0.0, 0.005, 0.01, 0.02, 0.05],
)
# Pairs per relatedness. Siblings need many to characterise the IBD-driven
# spread; unrelated pairs barely vary, so the budget goes to siblings.
DEFAULT_N_PAIRS = qval({"unrelated": 10, "sibling": 40}, {"unrelated": 2, "sibling": 3})
# Sequencing replicates per pair: enough blanks for a stable parametric LoB and
# resolution on the 95% detection point (Wald SE ~3% at this count).
DEFAULT_N_SEQ_REPS = qval(30, 5)
# Across-pair band reported around the median LoD curve.
BAND_LO_Q = 0.10
BAND_HI_Q = 0.90

MAF_RANGE = (0.2, 0.5)
ERROR_RATE = 0.01
LOCUS_DROPOUT_RATE = 0.016
DEPTH_CV = 0.43
# Estimator default is grid_steps=1001. LoD only needs detection-rule precision
# (est_frac > LoB), not MLE precision, so a coarser grid is fine and ~5x faster;
# Nelder-Mead still refines from the grid maximum (f accurate to <1e-3).
ESTIMATOR_GRID_STEPS = 201

# Fast grid estimator params. Routed in via ALLOMIX_FAST_GRID=1 / --fast-grid;
# matches the exact estimator to well under 0.01 pp on donor_fraction (all the
# LoD rule consumes). Equivalence pinned in tests/test_fast_grid_estimator.py.
FAST_GRID_N_F = 201
FAST_GRID_N_RHO = 32


def _fast_grid_enabled() -> bool:
    """True if the opt-in fast grid estimator is requested via env var."""
    return os.environ.get("ALLOMIX_FAST_GRID", "") not in ("", "0", "false", "False")


def estimate_fraction(genos_informative, bias_dict):
    """Run the configured single-donor estimator.

    Routes through the fast vectorized grid estimator when ALLOMIX_FAST_GRID is
    set, otherwise the exact (default) estimator.
    """
    cal = PanelCalibration(biases=bias_dict)
    if _fast_grid_enabled():
        return estimate_single_donor_bb_grid(
            genos_informative,
            error_rate=ERROR_RATE,
            calibration=cal,
            n_f=FAST_GRID_N_F,
            n_rho=FAST_GRID_N_RHO,
        )
    return estimate_single_donor_bb(
        genos_informative,
        error_rate=ERROR_RATE,
        grid_steps=ESTIMATOR_GRID_STEPS,
        calibration=cal,
    )

FACTS_DIR = Path("output/facts")
WORK_DIR = Path("output/lod_validation")

LOGIT_95 = math.log(0.95 / 0.05)  # 2.944...
# Sentinels for the two edge cells: LoD below or above the probed range.
LOD_BELOW_RANGE = -1.0
LOD_ABOVE_RANGE = float("inf")
# Probed-fraction bounds clamping per-pair sentinel LoDs before the across-pair
# median/percentiles (detects-everything -> smallest, detects-nothing -> largest).
MIN_POS_FRACTION = min(f for f in TRUE_FRACTIONS if f > 0)
MAX_FRACTION = max(TRUE_FRACTIONS)


def derive_seed(*parts: object) -> int:
    """Deterministic seed from a tuple of arbitrary parts.

    Uses SHA-256 of the repr so seeds are stable across processes. Python's
    built-in hash() is per-process randomised for str/bytes (PEP 456), which
    previously made each run produce different "deterministic" output.
    """
    digest = hashlib.sha256(repr(parts).encode("utf-8")).digest()[:4]
    return int.from_bytes(digest, "big")


def detection_rate(est_fracs: list[float], lob: float) -> float:
    """Fraction of replicates whose estimate exceeds LoB."""
    if not est_fracs:
        return 0.0
    return sum(1 for e in est_fracs if e > lob) / len(est_fracs)


def compute_lob(est_fracs_at_zero: list[float]) -> float:
    """LoB = mean + 1.645 * SD across blank replicates (parametric, CLSI EP17-A2).

    EP17-A2 allows the parametric form (more efficient than the empirical 95th
    percentile, ~1/sqrt(n) lower SE) when blanks pass a normality test. Our
    blank est_fracs are sums of many independent marker contributions (CLT) and
    pass Anderson-Darling in practice. 1.645 is the standard-normal 95th pctile.
    """
    if not est_fracs_at_zero:
        return float("nan")
    arr = np.asarray(est_fracs_at_zero, dtype=float)
    return float(arr.mean() + 1.645 * arr.std(ddof=1))


def _logistic(log10_f: np.ndarray, a: float, b: float) -> np.ndarray:
    """P(detected | f) = 1 / (1 + exp(-(a + b * log10(f))))."""
    return 1.0 / (1.0 + np.exp(-(a + b * log10_f)))


def _interp_lod(
    fractions: list[float], rates: list[float], target: float = 0.95,
) -> float | None:
    """Linear interpolation in log10(f) between the fractions bracketing ``target``
    detection rate; None if the rates never cross it.

    Fallback when the logistic fit fails (detection jumping too sharply from 0
    to 1 makes the logistic slope unidentifiable).
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

    ``weights`` are per-point replicate counts, turned into binomial-SE sigma
    for curve_fit. Returns (f95, a, b), or None if the fit failed.
    """
    pos = [(f, r) for f, r in zip(fractions, detection_rates) if f > 0]
    if len(pos) < 2:
        return None
    log10_f = np.array([math.log10(f) for f, _ in pos], dtype=float)
    rates = np.array([r for _, r in pos], dtype=float)

    if weights is not None and len(weights) == len(pos):
        # binomial SE as sigma: sqrt(p(1-p)/N), floored
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
    # Require positive slope (detection rises with fraction). A negative b is a
    # degenerate curve_fit solution (the "ultra-easy" corner where detection is
    # ~1.0 everywhere and the slope is unidentifiable); use the interp fallback.
    if math.isfinite(a) and math.isfinite(b) and b > 1e-9:
        log10_f95 = (LOGIT_95 - a) / b
        try:
            cand = 10.0**log10_f95
        except OverflowError:
            cand = float("inf")
        if math.isfinite(cand) and cand > 0:
            f95 = cand

    if f95 is None:
        # Fallback: interpolate the crossing when the logistic slope is
        # unidentifiable (sharp 0->1 jump). The interpolated crossing is valid.
        interp = _interp_lod([f for f, _ in pos], [r for _, r in pos])
        if interp is None:
            return None
        f95 = interp
        if not math.isfinite(a) or not math.isfinite(b):
            b = 5.0
            a = LOGIT_95 - b * math.log10(f95)

    return f95, a, b


def run_pair(
    relatedness: str,
    pair_idx: int,
    n_seq_reps: int,
    base_seed: int,
    work_root: Path,
) -> list[dict]:
    """Generate one fixed host/donor pair, then run ``n_seq_reps`` sequencing
    replicates across every (n_markers, depth, true_frac) cell.

    Genotypes and per-marker biases are drawn once per pair (seeded by
    ``pair_idx``, not the seq replicate) and reused across replicates and panel
    sizes; only the blend seed varies per replicate, so within a pair the sole
    variation is sequencing noise. A caller can then estimate an LoD per pair
    and characterise the across-pair spread separately.

    Strict nesting invariant: the n_markers=50 result must be a bit-identical
    prefix of n_markers=400. The admix VCF for a given (seq_rep, depth, frac) is
    generated ONCE at max_markers and sliced per panel size, sharing the rng
    consumption order (depth, locus-dropout, allele-count draws) so the read
    counts at marker i never depend on the downstream panel size.

    Returns one row per (n_markers, depth, true_frac, seq_rep).
    """
    pair_dir = work_root / relatedness / f"pair_{pair_idx}"
    pair_dir.mkdir(parents=True, exist_ok=True)

    max_markers = max(N_MARKERS_GRID)
    gt_seed = derive_seed("gt", relatedness, pair_idx, base_seed)
    rng = random.Random(gt_seed)
    all_markers = generate_related_genotypes(
        max_markers, relatedness, rng, maf_range=MAF_RANGE,
    )

    host_vcf = pair_dir / "host.vcf"
    donor_vcf = pair_dir / "donor.vcf"
    write_genotype_vcf(all_markers, host_vcf, "host", key="host_gt")
    write_genotype_vcf(all_markers, donor_vcf, "donor", key="donor_gt")
    host_md_full = parse_vcf(str(host_vcf), min_dp=0, min_gq=0)
    donor_md_full = parse_vcf(str(donor_vcf), min_dp=0, min_gq=0)

    bias_rng = random.Random(derive_seed("bias", relatedness, pair_idx, base_seed))
    all_biases = generate_marker_biases_realistic(max_markers, bias_rng)

    rows: list[dict] = []
    admix_path = pair_dir / "admix.vcf"

    for seq_rep in range(n_seq_reps):
        for depth in DEPTHS:
            for frac in TRUE_FRACTIONS:
                blend_seed = derive_seed(
                    "blend", relatedness, pair_idx, seq_rep, depth, frac, base_seed,
                )
                blend = blend_vcfs(
                    host_path=str(host_vcf),
                    donor_path=str(donor_vcf),
                    donor_fraction=frac,
                    target_depth=depth,
                    sample_name="admix",
                    seed=blend_seed,
                    fixed_biases=all_biases,
                    error_rate=ERROR_RATE,
                    locus_dropout_rate=LOCUS_DROPOUT_RATE,
                    depth_cv=DEPTH_CV,
                )
                bias_dict = (
                    {(c, p, r, a): b for c, p, r, a, b in blend.marker_biases}
                    if blend.marker_biases is not None
                    else None
                )
                # Do not "optimise" this write_vcf + parse_vcf round-trip away:
                # it is ~0.2% of wall time. ~99% is in estimate_single_donor_bb
                # (rho profiling + Nelder-Mead). To speed the sweep up, optimise
                # the estimator hot path (chimerism._ll_from_p_alt / grid+NM),
                # not the simulation or IO.
                write_vcf(blend, admix_path)
                admix_md_full = parse_vcf(str(admix_path), min_dp=0, min_gq=0)

                for n_markers in N_MARKERS_GRID:
                    host_md = host_md_full[:n_markers]
                    donor_md = donor_md_full[:n_markers]
                    admix_md = admix_md_full[:n_markers]
                    n_informative_truth = sum(
                        1 for m in all_markers[:n_markers] if m["informative"]
                    )

                    genos = classify_markers(host_md, [donor_md], admix_md, min_dp=0,
                                             min_gq=0, pass_only=False)
                    if len(genos.informative) < 1:
                        rows.append({
                            "relatedness": relatedness,
                            "depth": depth,
                            "n_markers": n_markers,
                            "true_frac": frac,
                            "pair": pair_idx,
                            "seq_rep": seq_rep,
                            "seed": blend_seed,
                            "est_frac": float("nan"),
                            "ci_lo": float("nan"),
                            "ci_hi": float("nan"),
                            "n_informative": 0,
                            "n_informative_truth": n_informative_truth,
                        })
                        continue

                    result = estimate_fraction(genos.informative, bias_dict)
                    rows.append({
                        "relatedness": relatedness,
                        "depth": depth,
                        "n_markers": n_markers,
                        "true_frac": frac,
                        "pair": pair_idx,
                        "seq_rep": seq_rep,
                        "seed": blend_seed,
                        "est_frac": result.donor_fraction,
                        "ci_lo": result.donor_fraction_ci[0],
                        "ci_hi": result.donor_fraction_ci[1],
                        "n_informative": result.n_informative,
                        "n_informative_truth": n_informative_truth,
                    })

    return rows


def compute_pair_lod(cell_rows: list[dict]) -> dict:
    """LoB + LoD for ONE pair at one (relatedness, depth, n_markers) cell.

    LoB comes from the blank (frac=0) replicates; LoD from a logistic fit of
    detection rate vs log10(fraction). ``lod`` is in fraction units, or a
    sentinel (LOD_BELOW_RANGE / LOD_ABOVE_RANGE) when detection is saturated, or
    NaN when the fit fails for a non-saturated cell.
    """
    by_frac: dict[float, list[float]] = {f: [] for f in TRUE_FRACTIONS}
    n_inf: list[int] = []
    for r in cell_rows:
        if math.isnan(r["est_frac"]):
            continue
        by_frac[r["true_frac"]].append(r["est_frac"])
        n_inf.append(r["n_informative"])

    lob = compute_lob(by_frac.get(0.0, []))

    fractions = [f for f in TRUE_FRACTIONS if f > 0]
    rates = [detection_rate(by_frac[f], lob) for f in fractions]
    weights = [len(by_frac[f]) for f in fractions]

    fit = fit_lod(fractions, rates, weights)
    note = ""
    if fit is None:
        # Two saturated cases: every cell detected, or none detected.
        if all(r >= 0.95 for r in rates):
            lod = LOD_BELOW_RANGE
            note = "all_detected"
        elif all(r < 0.05 for r in rates):
            lod = LOD_ABOVE_RANGE
            note = "none_detected"
        else:
            lod = float("nan")
            note = "fit_failed"
    else:
        lod = fit[0]

    return {
        "lod": lod,
        "lob": lob,
        "note": note,
        "mean_n_informative": float(np.mean(n_inf)) if n_inf else 0.0,
    }


def _lod_for_aggregation(lod: float) -> float | None:
    """Map a per-pair LoD (possibly a sentinel) to a finite fraction for the
    across-pair median / percentile aggregation, or None to drop it.

    - all_detected (LoD below the probed range): clamp to the smallest probed
      fraction; we cannot resolve finer, so this is a conservative floor.
    - none_detected (LoD above the probed range): clamp to the largest probed
      fraction; a conservative ceiling.
    - fit_failed / NaN: drop (the pair contributes nothing to the percentiles).
    """
    if lod == LOD_BELOW_RANGE:
        return MIN_POS_FRACTION
    if lod == LOD_ABOVE_RANGE:
        return MAX_FRACTION
    if math.isnan(lod):
        return None
    return lod


def summarise_cell(
    relatedness: str,
    depth: int,
    n_markers: int,
    pair_summaries: list[dict],
    n_seq_reps: int,
) -> dict:
    """Aggregate per-pair LoB/LoD into a median curve point and a 10-90% band.

    The LoD reported per (relatedness, depth, n_markers) is the median across
    pairs; the band is the BAND_LO_Q..BAND_HI_Q percentile across pairs (the
    IBD-driven spread, replacing the old bootstrap CI). LoB is the across-pair
    median. Pairs whose per-pair fit failed are excluded from the percentiles
    and counted in ``n_pairs_dropped``.
    """
    lods = [v for ps in pair_summaries
            if (v := _lod_for_aggregation(ps["lod"])) is not None]
    lobs = [ps["lob"] for ps in pair_summaries if not math.isnan(ps["lob"])]
    n_inf = [ps["mean_n_informative"] for ps in pair_summaries]
    n_dropped = len(pair_summaries) - len(lods)

    if lods:
        arr = np.asarray(lods, dtype=float)
        lod_med = float(np.median(arr))
        lod_lo = float(np.quantile(arr, BAND_LO_Q))
        lod_hi = float(np.quantile(arr, BAND_HI_Q))
    else:
        lod_med = lod_lo = lod_hi = float("nan")

    lob_med = float(np.median(lobs)) if lobs else float("nan")

    note = ""
    if not lods:
        note = "no_pairs_fit"
    elif n_dropped:
        note = f"{n_dropped}_pairs_dropped"

    return {
        "relatedness": relatedness,
        "depth": depth,
        "n_markers": n_markers,
        "lob_pct": _to_pct(lob_med),
        "lod_pct": _to_pct(lod_med),
        "lod_pct_ci_lo": _to_pct(lod_lo),
        "lod_pct_ci_hi": _to_pct(lod_hi),
        "mean_n_informative": round(float(np.mean(n_inf)), 1) if n_inf else 0.0,
        "median_n_informative": round(float(np.median(n_inf)), 1) if n_inf else 0.0,
        "n_pairs": len(pair_summaries),
        "n_pairs_used": len(lods),
        "n_pairs_dropped": n_dropped,
        "n_seq_reps": n_seq_reps,
        "note": note,
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


def write_grid_raw(rows: list[dict], path: Path) -> None:
    fields = [
        "relatedness", "depth", "n_markers", "true_frac", "pair", "seq_rep", "seed",
        "est_frac", "ci_lo", "ci_hi", "n_informative", "n_informative_truth",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})


def write_per_pair(pair_rows: list[dict], path: Path) -> None:
    """Per-pair LoD/LoB, one row per (relatedness, depth, n_markers, pair).

    Input to the across-pair median/band, kept for tracing which sibling pairs
    drive the band width.
    """
    fields = [
        "relatedness", "depth", "n_markers", "pair",
        "lod_pct", "lob_pct", "note", "mean_n_informative",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in pair_rows:
            w.writerow({k: r.get(k, "") for k in fields})


def write_summary(summaries: list[dict], path: Path) -> None:
    fields = [
        "relatedness", "depth", "n_markers",
        "lob_pct", "lod_pct", "lod_pct_ci_lo", "lod_pct_ci_hi",
        "mean_n_informative", "median_n_informative",
        "n_pairs", "n_pairs_used", "n_pairs_dropped", "n_seq_reps", "note",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for s in summaries:
            w.writerow({k: s.get(k, "") for k in fields})


def write_headline(
    summaries: list[dict], path: Path, n_pairs: dict[str, int], n_seq_reps: int,
) -> None:
    """One-row CSV of cells the Results prose quotes."""
    def lookup(rel: str, d: int, nm: int) -> float:
        for s in summaries:
            if s["relatedness"] == rel and s["depth"] == d and s["n_markers"] == nm:
                return s["lod_pct"]
        return float("nan")

    headline = {
        "n_pairs_unrelated": n_pairs.get("unrelated", 0),
        "n_pairs_sibling": n_pairs.get("sibling", 0),
        "n_seq_reps": n_seq_reps,
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


def _worker(args: tuple) -> list[dict]:
    relatedness, pair_idx, n_seq_reps, base_seed, work_root, n_markers_subset = args
    rows = run_pair(relatedness, pair_idx, n_seq_reps, base_seed, work_root)
    if n_markers_subset is not None:
        keep = set(n_markers_subset)
        rows = [r for r in rows if r["n_markers"] in keep]
    return rows


def main(argv: list[str] | None = None) -> int:
    global DEPTHS, ERROR_RATE
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--n-pairs", type=int, default=None,
        help="Donor/host pairs per relatedness. Overrides the per-relatedness "
             f"defaults {DEFAULT_N_PAIRS} with a single value for all levels.",
    )
    parser.add_argument("--n-seq-reps", type=int, default=DEFAULT_N_SEQ_REPS,
                        help="Sequencing replicates per pair.")
    parser.add_argument("--n-workers", "-j", type=int, default=1,
                        help="Process pool size for pair-level parallelism.")
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
    parser.add_argument(
        "--error-rate", type=float, default=ERROR_RATE,
        help=f"Symmetric sequencing error rate (default {ERROR_RATE}). For a fair "
             "overlay against the presence LoD (plot_lod_grid --presence-summary), "
             "run BOTH sweeps at the same error rate: the presence test is far more "
             "error-sensitive, so comparing them at mismatched rates is misleading.",
    )
    parser.add_argument(
        "--fast-grid", action="store_true",
        help="Opt in to the fast vectorized grid estimator (sets ALLOMIX_FAST_GRID=1 "
             "for the worker processes). The exact estimator is the default; use the "
             "fast path for quick sweep iteration, not the published full-accuracy pass.",
    )
    args = parser.parse_args(argv)

    if args.fast_grid:
        os.environ["ALLOMIX_FAST_GRID"] = "1"

    if args.depths is not None:
        DEPTHS = sorted(args.depths)
    ERROR_RATE = args.error_rate

    FACTS_DIR.mkdir(parents=True, exist_ok=True)
    work_root = Path(args.workdir)
    work_root.mkdir(parents=True, exist_ok=True)

    n_pairs = {
        rel: (args.n_pairs if args.n_pairs is not None
              else DEFAULT_N_PAIRS.get(rel, 10))
        for rel in args.relatedness
    }

    # One task per (relatedness, pair); the worker iterates internally over seq
    # replicates and n_markers, sharing the max-size admix VCF across panel sizes
    # so the cells are strictly nested.
    n_markers_subset = (
        None if args.n_markers == N_MARKERS_GRID else list(args.n_markers)
    )
    tasks = [
        (rel, pair_idx, args.n_seq_reps, args.seed, work_root, n_markers_subset)
        for rel in args.relatedness
        for pair_idx in range(n_pairs[rel])
    ]
    print(
        f"Sweep: relatedness={args.relatedness}, depths={DEPTHS}, "
        f"n_markers={args.n_markers}, fractions={TRUE_FRACTIONS}, "
        f"pairs={n_pairs}, seq_reps={args.n_seq_reps}",
        file=sys.stderr,
    )
    print(
        f"Total (relatedness, pair) tasks: {len(tasks)} "
        f"(each produces "
        f"{args.n_seq_reps * len(args.n_markers) * len(DEPTHS) * len(TRUE_FRACTIONS)} "
        f"estimator calls)",
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

    # Sort for output stable across worker completion order.
    all_rows.sort(key=lambda r: (
        r["relatedness"], r["n_markers"], r["depth"], r["pair"], r["seq_rep"],
        r["true_frac"],
    ))

    write_grid_raw(all_rows, FACTS_DIR / "lod_grid_raw.csv")
    print(f"Wrote {FACTS_DIR / 'lod_grid_raw.csv'} ({len(all_rows)} rows)", file=sys.stderr)

    # Per-pair LoD: one fit per (rel, depth, n_markers, pair) over its seq reps.
    by_pair_cell: dict[tuple, list[dict]] = defaultdict(list)
    for r in all_rows:
        by_pair_cell[(r["relatedness"], r["depth"], r["n_markers"], r["pair"])].append(r)

    pair_rows: list[dict] = []
    pairs_by_cell: dict[tuple, list[dict]] = defaultdict(list)
    for (rel, depth, nm, pair), rows in by_pair_cell.items():
        ps = compute_pair_lod(rows)
        pairs_by_cell[(rel, depth, nm)].append(ps)
        pair_rows.append({
            "relatedness": rel,
            "depth": depth,
            "n_markers": nm,
            "pair": pair,
            "lod_pct": _to_pct(ps["lod"]),
            "lob_pct": _to_pct(ps["lob"]),
            "note": ps["note"],
            "mean_n_informative": round(ps["mean_n_informative"], 1),
        })

    pair_rows.sort(key=lambda r: (
        r["relatedness"], r["n_markers"], r["depth"], r["pair"],
    ))
    write_per_pair(pair_rows, FACTS_DIR / "lod_per_pair.csv")
    print(f"Wrote {FACTS_DIR / 'lod_per_pair.csv'} ({len(pair_rows)} rows)", file=sys.stderr)

    # Across-pair summaries: median curve + 10-90% band.
    summaries: list[dict] = []
    for rel in args.relatedness:
        for depth in DEPTHS:
            for nm in args.n_markers:
                ps_list = pairs_by_cell.get((rel, depth, nm), [])
                summaries.append(
                    summarise_cell(rel, depth, nm, ps_list, args.n_seq_reps)
                )

    write_summary(summaries, FACTS_DIR / "lod_summary.csv")
    print(f"Wrote {FACTS_DIR / 'lod_summary.csv'} ({len(summaries)} rows)", file=sys.stderr)

    write_headline(summaries, FACTS_DIR / "lod_headline.csv", n_pairs, args.n_seq_reps)
    print(f"Wrote {FACTS_DIR / 'lod_headline.csv'}", file=sys.stderr)

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
