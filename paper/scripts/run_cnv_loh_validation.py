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

from allomix.chimerism import PanelCalibration, estimate_single_donor_bb  # noqa: E402
from allomix.genotype import classify_markers, parse_vcf  # noqa: E402
from allomix.simulate import (  # noqa: E402
    assign_cnv_aberrations,
    blend_vcfs,
    generate_marker_biases_realistic,
    generate_related_genotypes,
    write_genotype_vcf,
    write_vcf,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paper_quick import qval  # noqa: E402

# --- Sweep grid --------------------------------------------------------------

# Quick-build mode (ALLOMIX_PAPER_QUICK=1) only cuts the replicate count (see
# DEFAULT_N_REPS below), not the sweep grid: the headline snapshot needs every
# relatedness/kind/burden cell and the full true-fraction grid for the LoD
# interpolation, so trimming the grid leaves nan/"above range" headline values
# that break the paper render.
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

# Depth and panel-size grids, swept like run_lod_validation.py / plot_lod_curves.py.
# Defaults are a single operating point (100 markers, 1000x) so the Snakefile / CI
# build stays cheap and reproducible. To produce the full depth x markers curve
# (the heavy run; do it on a dedicated machine), pass e.g.
#   --depths 100 250 500 1000 2000 --n-markers 25 50 75 100 200 400
# Markers are nested: each cell blends the largest panel once per depth and
# evaluates every smaller panel as a strict prefix (as in run_lod_validation).
DEPTHS = [1000]
N_MARKERS_GRID = [100]
FULL_DEPTHS = [100, 250, 500, 1000, 2000]
FULL_N_MARKERS_GRID = [25, 50, 75, 100, 200, 400]

MAF_RANGE = (0.2, 0.5)
ERROR_RATE = 0.01
LOCUS_DROPOUT_RATE = 0.016
DEPTH_CV = 0.43
ESTIMATOR_GRID_STEPS = 201

# Quick-build mode (ALLOMIX_PAPER_QUICK=1) cuts the replicate count; the LoD
# estimates get noisier but the rule finishes in a fraction of the time. The
# resulting figure is watermarked, not for publication.
DEFAULT_N_REPS = qval(40, 8)

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


def _blend_and_estimate(
    mode: str,
    host_vcf: Path,
    donor_vcf: Path,
    host_md_full: list,
    donor_md_full: list,
    biases: list[float],
    aberrations: list | None,
    minor_frac: float,
    depth: int,
    blend_seed: int,
    admix_path: Path,
    n_markers_grid: list[int],
) -> dict[int, dict]:
    """Blend one sample at ``depth`` and estimate the minor component per panel.

    ``minor_frac`` is the true fraction of the minor (detected) component. The
    recipient/host always carries the aberration; ``mode`` sets which role it
    plays:

      - ``"donor"``: detect the donor (minor); host is the major background
        carrying the aberration. Mixed-chimerism / substantial-recipient regime.
      - ``"host"``: detect the recipient relapse (minor) carrying the aberration;
        donor is the clean major background. Early-warning relapse regime.

    The full panel is blended once; each panel size in ``n_markers_grid`` is then
    evaluated as a strict prefix (nested, as in run_lod_validation). Returns a
    dict keyed by panel size; ``est`` is the estimated minor-component fraction.
    """
    donor_frac = minor_frac if mode == "donor" else 1.0 - minor_frac
    blend = blend_vcfs(
        host_path=str(host_vcf),
        donor_path=str(donor_vcf),
        donor_fraction=donor_frac,
        target_depth=depth,
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
    admix_md_full = parse_vcf(str(admix_path), min_dp=0, min_gq=0)

    def minor(f: float) -> float:
        return f if mode == "donor" else 1.0 - f

    out: dict[int, dict] = {}
    for n_markers in n_markers_grid:
        genos = classify_markers(
            host_md_full[:n_markers], [donor_md_full[:n_markers]], admix_md_full[:n_markers],
            min_dp=0, min_gq=0, pass_only=False,
        )
        if len(genos.informative) < 1:
            out[n_markers] = dict(est=float("nan"), est_robust=float("nan"),
                                  n_informative=0, n_robust_excluded=0)
            continue
        calibration = PanelCalibration(biases=bias_dict)
        result = estimate_single_donor_bb(
            genos.informative, error_rate=ERROR_RATE,
            grid_steps=ESTIMATOR_GRID_STEPS, calibration=calibration,
        )
        robust = estimate_single_donor_bb(
            genos.informative, error_rate=ERROR_RATE,
            grid_steps=ESTIMATOR_GRID_STEPS, calibration=calibration, robust="auto",
        )
        out[n_markers] = dict(
            est=minor(result.donor_fraction),
            est_robust=minor(robust.donor_fraction),
            n_informative=result.n_informative,
            n_robust_excluded=robust.n_robust_excluded,
        )
    return out


def run_pair(
    mode: str,
    relatedness: str,
    rep: int,
    base_seed: int,
    depths: list[int],
    n_markers_grid: list[int],
) -> list[dict]:
    """Evaluate every aberration / depth / panel cell for one genotype pair, one mode.

    Genotypes and capture biases are fixed for this (relatedness, rep). Within
    the pair we vary the aberration kind, burden, minor-component fraction,
    depth, and panel size (markers nested as prefixes). The no-aberration
    baseline (burden 0) is run once and shared across kinds.
    """
    max_markers = max(n_markers_grid)
    pair_dir = WORK_DIR / f"{mode}_{relatedness}_rep{rep}"
    pair_dir.mkdir(parents=True, exist_ok=True)

    gt_rng = random.Random(derive_seed("gt", relatedness, rep, base_seed))
    markers = generate_related_genotypes(max_markers, relatedness, gt_rng, maf_range=MAF_RANGE)

    host_vcf = pair_dir / "host.vcf"
    donor_vcf = pair_dir / "donor.vcf"
    write_genotype_vcf(markers, host_vcf, "host", key="host_gt")
    write_genotype_vcf(markers, donor_vcf, "donor", key="donor_gt")
    host_md = parse_vcf(str(host_vcf), min_dp=0, min_gq=0)
    donor_md = parse_vcf(str(donor_vcf), min_dp=0, min_gq=0)

    bias_rng = random.Random(derive_seed("bias", relatedness, rep, base_seed))
    biases = generate_marker_biases_realistic(max_markers, bias_rng)

    admix_path = pair_dir / "admix.vcf"
    rows: list[dict] = []

    def emit(kind, burden, clonal, aberrations, true_frac, depth):
        blend_seed = derive_seed("blend", mode, relatedness, rep, kind, burden, clonal,
                                 true_frac, depth)
        per_panel = _blend_and_estimate(
            mode, host_vcf, donor_vcf, host_md, donor_md, biases, aberrations,
            true_frac, depth, blend_seed, admix_path, n_markers_grid,
        )
        for n_markers, res in per_panel.items():
            n_aff = 0 if aberrations is None else sum(1 for a in aberrations[:n_markers] if a)
            rows.append({
                "mode": mode, "relatedness": relatedness, "rep": rep, "kind": kind,
                "burden": burden, "clonal_fraction": clonal, "true_frac": true_frac,
                "depth": depth, "n_markers": n_markers, "n_affected": n_aff, "seed": blend_seed,
                **res,
            })

    # Baseline (no aberration), shared by every kind.
    for true_frac in TRUE_FRACTIONS:
        for depth in depths:
            emit("baseline", 0.0, 0.0, None, true_frac, depth)

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
                for true_frac in TRUE_FRACTIONS:
                    for depth in depths:
                        emit(kind, burden, clonal, aberrations, true_frac, depth)

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
        cells[(r["mode"], r["relatedness"], r["kind"], r["clonal_fraction"],
               r["burden"], r["depth"], r["n_markers"])].append(r)

    out: list[dict] = []
    for (mode, relatedness, kind, clonal, burden, depth, n_markers), rs in sorted(cells.items()):
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
                "depth": depth,
                "n_markers": n_markers,
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
    parser.add_argument("--depths", nargs="+", type=int, default=DEPTHS,
                        help=f"Mean depths to sweep (default {DEPTHS}; full curve: {FULL_DEPTHS})")
    parser.add_argument("--n-markers", nargs="+", type=int, default=N_MARKERS_GRID,
                        dest="n_markers_grid",
                        help=f"Panel sizes to sweep, nested (default {N_MARKERS_GRID}; "
                             f"full curve: {FULL_N_MARKERS_GRID})")
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
            raw.extend(run_pair(mode, rel, rep, args.base_seed, args.depths, args.n_markers_grid))
    else:
        with ProcessPoolExecutor(max_workers=args.n_workers) as ex:
            futures = [
                ex.submit(run_pair, mode, rel, rep, args.base_seed, args.depths, args.n_markers_grid)
                for mode, rel, rep in jobs
            ]
            for fut in as_completed(futures):
                raw.extend(fut.result())

    raw_fields = [
        "mode", "relatedness", "rep", "kind", "burden", "clonal_fraction", "true_frac",
        "depth", "n_markers", "n_affected", "seed", "est", "est_robust",
        "n_informative", "n_robust_excluded",
    ]
    write_csv(FACTS_DIR / "cnv_loh_raw.csv", raw, raw_fields)

    summary = summarise(raw)
    summary_fields = [
        "mode", "relatedness", "kind", "clonal_fraction", "burden", "depth", "n_markers",
        "n_reps_per_frac", "lob_std", "lod_std", "lob_robust", "lod_robust",
        "mean_n_affected", "mean_n_robust_excluded", "mean_n_informative",
    ]
    write_csv(FACTS_DIR / "cnv_loh_summary.csv", summary, summary_fields)

    headline = build_headline(summary)
    # Single wide row, like lod_headline.csv, so the paper can template
    # `{{ cnv_loh_headline.<field> }}`.
    write_csv(FACTS_DIR / "cnv_loh_headline.csv", [headline], list(headline.keys()))

    print(f"Wrote {len(raw)} raw rows, {len(summary)} summary cells to {FACTS_DIR}")
    for k, v in headline.items():
        print(f"  {k}: {v}")


def build_headline(summary: list[dict]) -> list[dict]:
    """Pull interpretable headline LoD numbers per mode / aberration kind.

    Headline numbers are reported at a reference operating point (the largest
    swept depth and panel size present), so they are well-defined whether the
    sweep ran a single cell or the full depth x markers grid.
    """
    max_burden = max(BURDEN_LEVELS)
    ref_depth = max(s["depth"] for s in summary) if summary else None
    ref_markers = max(s["n_markers"] for s in summary) if summary else None

    def cell(mode, rel, kind, burden):
        for s in summary:
            if (s["mode"] == mode and s["relatedness"] == rel and s["kind"] == kind
                    and abs(s["burden"] - burden) < 1e-9
                    and s["depth"] == ref_depth and s["n_markers"] == ref_markers):
                return s
        return None

    low_burden = min(b for b in BURDEN_LEVELS if b > 0)

    def lod_pct(lod):
        """Minor-component LoD as a %; '>ceiling' string when above probed range."""
        if lod != lod:
            return float("nan")
        if not math.isfinite(lod):
            return f">{MAX_PROBED * 100:g}"
        return round(lod * 100, 4)

    def lp(mode, rel, kind, burden, field):
        c = cell(mode, rel, kind, burden)
        return lod_pct(c[field]) if c else float("nan")

    # Single wide row of named facts (paper templating consumes
    # `{{ cnv_loh_headline.<field> }}`); reference operating point = largest
    # swept depth and panel size.
    h: dict = {
        "n_reps_per_frac": summary[0]["n_reps_per_frac"] if summary else 0,
        "ref_depth": ref_depth,
        "ref_markers": ref_markers,
        # Relapse (host) detection: baseline and the worst aberration cell.
        "relapse_lod_baseline_unrel_pct": lp("host", "unrelated", "baseline", 0.0, "lod_std"),
        "relapse_lod_baseline_sib_pct": lp("host", "sibling", "baseline", 0.0, "lod_std"),
        # Donor detection (mixed chimerism): baselines and key cells.
        "donor_lod_baseline_unrel_pct": lp("donor", "unrelated", "baseline", 0.0, "lod_std"),
        "donor_lod_baseline_sib_pct": lp("donor", "sibling", "baseline", 0.0, "lod_std"),
        "donor_lod_cnloh_low_unrel_std_pct": lp("donor", "unrelated", "cnloh", low_burden, "lod_std"),
        "donor_lod_deletion_low_unrel_std_pct": lp("donor", "unrelated", "deletion", low_burden, "lod_std"),
        "donor_lod_deletion_low_unrel_robust_pct": lp("donor", "unrelated", "deletion", low_burden, "lod_robust"),
        "donor_lod_gain_high_unrel_std_pct": lp("donor", "unrelated", "gain", max_burden, "lod_std"),
        "donor_lod_gain_high_unrel_robust_pct": lp("donor", "unrelated", "gain", max_burden, "lod_robust"),
    }
    # Worst (largest) relapse LoD across all aberration cells at the reference
    # point, to bound the "relapse detection unaffected" claim.
    relapse_cells = [
        s for s in summary
        if s["mode"] == "host" and s["kind"] != "baseline"
        and s["depth"] == ref_depth and s["n_markers"] == ref_markers
        and math.isfinite(s["lod_std"])
    ]
    h["relapse_lod_max_pct"] = round(max(s["lod_std"] for s in relapse_cells) * 100, 4) if relapse_cells else float("nan")
    return h


if __name__ == "__main__":
    main()
