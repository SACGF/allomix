#!/usr/bin/env python3
"""Measure how host-genome CNV / CN-LoH degrades the donor limit of detection (issue #13).

The HSCT recipient is usually a haematological malignancy patient, so the
residual host clone routinely carries somatic copy-number changes: copy-neutral
LoH (CN-LoH, acquired uniparental disomy; het -> hom in the clone), deletions
(CN1) and gains (CN3). The diploid VAF model does not account for these, so the
affected markers carry a biased ALT VAF.

The host (major, background component) carries the aberration and the donor is
the minor component swept for the **limit of detection** (LoD), the metric the
rest of the paper reports (run_lod_validation.py, CLSI EP17-A2). This is the
mixed-chimerism / substantial-recipient regime, where the recipient's CN-LoH/CNV
background degrades the donor LoD. The donor fraction is swept at low log-spaced
values including 0:

  - LoB = mean + 1.645 * SD of the estimated donor fraction over blank
          replicates (true donor = 0, i.e. a pure host carrying the aberration).
  - LoD = lowest true donor fraction at which >=95% of replicates have
          est_frac > LoB, from a logistic fit of P(detected) vs log10(f).

The LoD is plotted on a log donor-% axis (0.3, 0.5, 1, 2, 5, 10, 20 %), matching
the depth x markers LoD curves (plot_lod_curves.py). Both standard and robust
(``--robust auto``) are computed per cell, so the figure shows the LoD inflation
from host CNV/LoH and the robust recovery. Genotypes and biases are fixed per
(relatedness, replicate).

Outputs:
  output/facts/cnv_loh_raw.csv       # one row per replicate per cell
  output/facts/cnv_loh_summary.csv   # donor LoB/LoD per (rel, kind, burden, clonal)
  output/facts/cnv_loh_headline.csv  # headline donor-LoD snapshot

Usage:
    python paper/scripts/run_cnv_loh_validation.py
    python paper/scripts/run_cnv_loh_validation.py --n-reps 40 --n-workers 8
"""

import argparse
import csv
import hashlib
import math
import random
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

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
# Detection direction (which low-fraction component the LoD is for; the recipient
# always carries the aberration). "donor": detect donor against a CN-LoH host
# background (mixed-chimerism). "host": detect the recipient relapse clone (it
# carries the aberration) against a clean donor background (early warning).
MODES = ["donor", "host"]
# Copy-number aberration kinds applied to the host clone. CN-LoH is copy-neutral
# (allele-balance effect at het markers); deletion (CN1) and gain (CN3) also
# change the locus DNA mass, so they shift the local mixing fraction even at
# homozygous markers.
KINDS = ["cnloh", "deletion", "gain"]
# Fraction of eligible host markers carrying the aberration. For cnloh, only het
# markers are eligible; for deletion/gain every marker is. 0.0 is the baseline,
# run once per replicate and shared across kinds.
BURDEN_LEVELS = [0.0, 0.1, 0.25, 0.5]
# Pure clone (whole host is the aberrant clone), the worst case for the host
# component carrying the aberration.
CLONAL_FRACTIONS = [1.0]
# Low, log-spaced true donor fractions (incl. 0 for LoB). Extends above the usual
# 5% ceiling because host CNV/LoH inflates the donor LoD well past it; a LoD
# beyond MAX_PROBED is reported as "above range" (donor undetectable here).
TRUE_FRACTIONS = [0.0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]
MAX_PROBED = 0.2

N_MARKERS = 100
DEPTH = 1000
MAF_RANGE = (0.2, 0.5)
ERROR_RATE = 0.01
LOCUS_DROPOUT_RATE = 0.016
DEPTH_CV = 0.43
ESTIMATOR_GRID_STEPS = 201

DEFAULT_N_REPS = 40

# CLSI EP17-A2 LoD helpers (mirrors run_lod_validation.py).
LOGIT_95 = math.log(0.95 / 0.05)

FACTS_DIR = Path("output/facts")
WORK_DIR = Path("output/cnv_loh_validation")


def derive_seed(*parts: object) -> int:
    """Deterministic seed from arbitrary parts (stable across processes)."""
    digest = hashlib.sha256(repr(parts).encode("utf-8")).digest()[:4]
    return int.from_bytes(digest, "big")


# --- CLSI EP17-A2 LoB / LoD (copied from run_lod_validation.py for independence) ---


def compute_lob(est_fracs_at_zero: list[float]) -> float:
    """LoB = mean + 1.645 * SD across blank (true=0) replicates."""
    if len(est_fracs_at_zero) < 2:
        return float("nan")
    arr = np.asarray(est_fracs_at_zero, dtype=float)
    return float(arr.mean() + 1.645 * arr.std(ddof=1))


def detection_rate(est_fracs: list[float], lob: float) -> float:
    """Fraction of replicates whose estimate exceeds LoB."""
    if not est_fracs:
        return 0.0
    return sum(1 for e in est_fracs if e > lob) / len(est_fracs)


def _logistic(log10_f: np.ndarray, a: float, b: float) -> np.ndarray:
    z = np.clip(a + b * log10_f, -500.0, 500.0)
    return 1.0 / (1.0 + np.exp(-z))


def _interp_lod(fractions: list[float], rates: list[float], target: float = 0.95) -> float | None:
    pairs = sorted([(f, r) for f, r in zip(fractions, rates) if f > 0])
    for (f_lo, r_lo), (f_hi, r_hi) in zip(pairs, pairs[1:]):
        if r_lo <= target <= r_hi:
            if r_hi == r_lo:
                return f_lo
            log_lo, log_hi = math.log10(f_lo), math.log10(f_hi)
            frac = (target - r_lo) / (r_hi - r_lo)
            return 10.0 ** (log_lo + frac * (log_hi - log_lo))
    return None


def fit_lod(fractions: list[float], detection_rates: list[float]) -> float | None:
    """Logistic fit in log10(f); solve for f at P=0.95, with interp fallback."""
    pos = [(f, r) for f, r in zip(fractions, detection_rates) if f > 0]
    if len(pos) < 2:
        return None
    log10_f = np.array([math.log10(f) for f, _ in pos], dtype=float)
    rates = np.array([r for _, r in pos], dtype=float)
    a = b = float("nan")
    try:
        popt, _ = curve_fit(_logistic, log10_f, rates, p0=[2.0, 2.0], maxfev=10000)
        a, b = float(popt[0]), float(popt[1])
    except (RuntimeError, ValueError):
        pass

    f95: float | None = None
    if math.isfinite(a) and math.isfinite(b) and b > 1e-9:
        try:
            cand = 10.0 ** ((LOGIT_95 - a) / b)
        except OverflowError:
            cand = float("inf")
        if math.isfinite(cand) and cand > 0:
            f95 = cand
    if f95 is None:
        f95 = _interp_lod([f for f, _ in pos], [r for _, r in pos])
    return f95


def _eval_cell(
    mode: str,
    host_vcf: Path,
    donor_vcf: Path,
    host_md: list,
    donor_md: list,
    biases: list[float],
    aberrations: list | None,
    minor_frac: float,
    blend_seed: int,
    admix_path: Path,
) -> dict:
    """Blend one sample and estimate the minor component for the LoD.

    ``minor_frac`` is the true fraction of the minor (detected) component. The
    recipient/host always carries the aberration; ``mode`` sets which role it
    plays:

      - ``"donor"``: detect the donor (minor); host is the major background
        carrying the aberration. Mixed-chimerism / substantial-recipient regime.
      - ``"host"``: detect the recipient relapse (minor) carrying the aberration;
        donor is the clean major background. Early-warning relapse regime.

    The returned ``est`` is the estimated minor-component fraction (donor f, or
    1 - f for host mode), the detection statistic for the LoD.
    """
    donor_frac = minor_frac if mode == "donor" else 1.0 - minor_frac
    blend = blend_vcfs(
        host_path=str(host_vcf),
        donor_path=str(donor_vcf),
        donor_fraction=donor_frac,
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
        return dict(est=float("nan"), est_robust=float("nan"),
                    n_informative=0, n_robust_excluded=0)

    result = estimate_single_donor_bb(
        genos.informative, error_rate=ERROR_RATE,
        grid_steps=ESTIMATOR_GRID_STEPS, marker_biases=bias_dict,
    )
    # Robust refit (default policy) to show the mitigation effect in the paper.
    robust = estimate_single_donor_bb(
        genos.informative, error_rate=ERROR_RATE,
        grid_steps=ESTIMATOR_GRID_STEPS, marker_biases=bias_dict, robust="auto",
    )

    def minor(f: float) -> float:
        return f if mode == "donor" else 1.0 - f

    return dict(
        est=minor(result.donor_fraction),
        est_robust=minor(robust.donor_fraction),
        n_informative=result.n_informative,
        n_robust_excluded=robust.n_robust_excluded,
    )


def run_pair(mode: str, relatedness: str, rep: int, base_seed: int) -> list[dict]:
    """Evaluate every aberration cell for one fixed genotype pair, one mode.

    Genotypes and capture biases are fixed for this (relatedness, rep). Within
    the pair we vary the aberration kind, burden, and minor-component fraction.
    The no-aberration baseline (burden 0) is run once and shared across kinds.
    ``mode`` is "donor" or "host" (see ``_eval_cell``).
    """
    pair_dir = WORK_DIR / f"{mode}_{relatedness}_rep{rep}"
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
        blend_seed = derive_seed("blend", mode, relatedness, rep, "baseline", true_frac)
        row = {
            "mode": mode, "relatedness": relatedness, "rep": rep, "kind": "baseline",
            "burden": 0.0, "clonal_fraction": 0.0, "true_frac": true_frac,
            "n_affected": 0, "seed": blend_seed,
        }
        row.update(_eval_cell(mode, host_vcf, donor_vcf, host_md, donor_md, biases,
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
                        "blend", mode, relatedness, rep, kind, burden, clonal, true_frac
                    )
                    row = {
                        "mode": mode, "relatedness": relatedness, "rep": rep, "kind": kind,
                        "burden": burden, "clonal_fraction": clonal, "true_frac": true_frac,
                        "n_affected": n_affected, "seed": blend_seed,
                    }
                    row.update(_eval_cell(mode, host_vcf, donor_vcf, host_md, donor_md, biases,
                                          aberrations, true_frac, blend_seed, admix_path))
                    rows.append(row)

    return rows


def _cell_lod(rows: list[dict], est_key: str) -> tuple[float, float]:
    """Compute (LoB, LoD) for one (rel, kind, burden, clonal) cell.

    Groups the cell's replicate donor-fraction estimates by true donor fraction,
    forms the LoB from the blanks (true donor = 0, i.e. pure host carrying the
    aberration), then fits the >=95% detection point. LoD is a donor fraction.
    """
    by_frac: dict[float, list[float]] = defaultdict(list)
    for r in rows:
        v = r[est_key]
        if v == v:  # not NaN
            by_frac[r["true_frac"]].append(v)
    blanks = by_frac.get(0.0, [])
    lob = compute_lob(blanks)
    if not math.isfinite(lob):
        return float("nan"), float("nan")
    fracs = sorted(f for f in by_frac if f > 0)
    rates = [detection_rate(by_frac[f], lob) for f in fracs]
    lod = fit_lod(fracs, rates)
    # A LoD past the probed ceiling means the donor is not detectable here;
    # report it as above-range (inf) rather than an extrapolated number.
    if lod is None or not math.isfinite(lod) or lod > MAX_PROBED:
        lod = float("inf")
    return lob, lod


def summarise(raw: list[dict]) -> list[dict]:
    """Aggregate raw rows into per-cell minor-component LoB/LoD (std and robust)."""
    cells: dict[tuple, list[dict]] = defaultdict(list)
    for r in raw:
        cells[(r["mode"], r["relatedness"], r["kind"], r["clonal_fraction"], r["burden"])].append(r)

    out: list[dict] = []
    for (mode, relatedness, kind, clonal, burden), rs in sorted(cells.items()):
        lob_std, lod_std = _cell_lod(rs, "est")
        lob_rob, lod_rob = _cell_lod(rs, "est_robust")
        n_reps = max((sum(1 for r in rs if r["true_frac"] == f) for f in TRUE_FRACTIONS), default=0)
        out.append(
            {
                "mode": mode,
                "relatedness": relatedness,
                "kind": kind,
                "clonal_fraction": clonal,
                "burden": burden,
                "n_reps_per_frac": n_reps,
                "lob_std": lob_std,
                "lod_std": lod_std,
                "lob_robust": lob_rob,
                "lod_robust": lod_rob,
                "mean_n_affected": sum(r["n_affected"] for r in rs) / len(rs),
                "mean_n_robust_excluded": sum(r["n_robust_excluded"] for r in rs) / len(rs),
                "mean_n_informative": sum(r["n_informative"] for r in rs) / len(rs),
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
    parser.add_argument("--modes", nargs="+", default=list(MODES), choices=MODES)
    args = parser.parse_args()

    WORK_DIR.mkdir(parents=True, exist_ok=True)

    jobs = [
        (mode, rel, rep)
        for mode in args.modes
        for rel in RELATEDNESS_LEVELS
        for rep in range(args.n_reps)
    ]
    raw: list[dict] = []
    if args.n_workers <= 1:
        for mode, rel, rep in jobs:
            raw.extend(run_pair(mode, rel, rep, args.base_seed))
    else:
        with ProcessPoolExecutor(max_workers=args.n_workers) as ex:
            futures = [ex.submit(run_pair, mode, rel, rep, args.base_seed) for mode, rel, rep in jobs]
            for fut in as_completed(futures):
                raw.extend(fut.result())

    raw_fields = [
        "mode", "relatedness", "rep", "kind", "burden", "clonal_fraction", "true_frac",
        "n_affected", "seed", "est", "est_robust", "n_informative", "n_robust_excluded",
    ]
    write_csv(FACTS_DIR / "cnv_loh_raw.csv", raw, raw_fields)

    summary = summarise(raw)
    summary_fields = [
        "mode", "relatedness", "kind", "clonal_fraction", "burden", "n_reps_per_frac",
        "lob_std", "lod_std", "lob_robust", "lod_robust",
        "mean_n_affected", "mean_n_robust_excluded", "mean_n_informative",
    ]
    write_csv(FACTS_DIR / "cnv_loh_summary.csv", summary, summary_fields)

    headline = build_headline(summary)
    write_csv(FACTS_DIR / "cnv_loh_headline.csv", headline, ["metric", "value"])

    print(f"Wrote {len(raw)} raw rows, {len(summary)} summary cells to {FACTS_DIR}")
    for h in headline:
        print(f"  {h['metric']}: {h['value']}")


def build_headline(summary: list[dict]) -> list[dict]:
    """Pull interpretable headline LoD numbers per mode / aberration kind."""
    max_burden = max(BURDEN_LEVELS)

    def cell(mode, rel, kind, burden):
        for s in summary:
            if (s["mode"] == mode and s["relatedness"] == rel and s["kind"] == kind
                    and abs(s["burden"] - burden) < 1e-9):
                return s
        return None

    def lod_pct(lod):
        """Minor-component LoD as a %; '>ceiling' when above the probed range."""
        if lod != lod:  # NaN
            return float("nan")
        if not math.isfinite(lod):  # above probed ceiling = undetectable
            return f">{MAX_PROBED * 100:g}"
        return round(lod * 100, 4)

    headline: list[dict] = []
    for mode in MODES:
        for rel in RELATEDNESS_LEVELS:
            base = cell(mode, rel, "baseline", 0.0)
            if base:
                headline.append(
                    {"metric": f"{mode}_lod_pct_baseline_{rel}", "value": lod_pct(base["lod_std"])}
                )
            for kind in KINDS:
                worst = cell(mode, rel, kind, max_burden)
                if not worst:
                    continue
                headline.append(
                    {"metric": f"{mode}_lod_pct_{kind}_b{max_burden}_std_{rel}",
                     "value": lod_pct(worst["lod_std"])}
                )
                headline.append(
                    {"metric": f"{mode}_lod_pct_{kind}_b{max_burden}_robust_{rel}",
                     "value": lod_pct(worst["lod_robust"])}
                )
                if base and math.isfinite(base["lod_std"]) and base["lod_std"] > 0:
                    w = worst["lod_std"]
                    infl = round(w / base["lod_std"], 2) if math.isfinite(w) else "above_range"
                    headline.append(
                        {"metric": f"{mode}_lod_inflation_x_{kind}_b{max_burden}_{rel}", "value": infl}
                    )
    return headline


if __name__ == "__main__":
    main()
