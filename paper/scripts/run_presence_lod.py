#!/usr/bin/env python3
"""Calibration + LoD harness for the host-presence detection test.

Production code now lives in ``src/allomix/detect.py``
(``host_presence_test``, ``HostPresenceResult``). This script remains the
calibration / gate-evidence harness — it sweeps the synthetic grid used to
clear the three acceptance gates in the plan, so the LRT and pooled-Poisson
helpers below are kept as their inline-tested form. New downstream code
should call ``allomix.detect.host_presence_test`` instead.

See ``claude/20_host_presence_detection_plan.md``, section "Prototype spec".

We simulate post-HSCT admixture samples at very low host fractions, restrict to
donor-homozygous markers where the host carries the donor-absent allele (Vynck
types 0, 1, 10, 11), and run two presence statistics on the donor-absent allele
read counts:

  - Pooled Poisson:  Y = sum y_i, Lam = sum n_i e_i, p_pois = P(Poisson(Lam) >= Y).
  - LRT:             q_i(f_h) = e_i + (h_i/2) f_h; loglik over a binomial per marker,
                     bounded at f_h >= 0; chi-bar-square one-sided p-value.

The simulator's symmetric global error rate ``e`` puts an effective per-direction
floor of ``e/3`` at a donor-homozygous marker, so we give the detector ``e_i = e/3``
for every marker and sweep ``e``. The presence test is then calibrated by
construction; sweeping ``e`` substitutes for the (not-yet-built) per-site error
table from Step 14.

We also call ``chimerism.estimate_single_donor_bb`` on the same admixture sample
so we can compare presence-LoD to the MLE-LoB LoD at matched cells, and check
that ``f_h_hat`` from the LRT tracks ``1 - donor_fraction`` from the MLE.

Outputs:
  output/facts/presence_lod_raw.csv      # one row per replicate
  output/facts/presence_lod_summary.csv  # per-cell summary (FP rate, LoDs)

Usage:
    # Pilot (fast): 10 reps per cell, reduced grid
    python paper/scripts/run_presence_lod.py --pilot

    # Full run as specified in the plan
    python paper/scripts/run_presence_lod.py --n-blanks 200 --n-positives 60 \
        --n-workers 8

Note on rho: the simulator gained a ``rho`` argument (Step 21). By default this
script runs the binomial case (``rho=inf, rho_marker_type='all'``) for backward
compatibility with the original prototype run. Passing ``--rho 100
--rho-marker-type het_only`` runs the marker-type-aware regime: overdispersion
applies to het/intermediate markers only, leaving the donor-absent allele
background binomial. That is the configuration that honestly exposes acceptance
gate #3: the MLE pays the overdispersion tax while the presence detector,
restricted to donor-homozygous markers at the binomial error floor, does not.
A uniform ``--rho 100 --rho-marker-type all`` would miscalibrate the presence
null (overdispersing the clean background); do not use it for the gate check.
"""

import argparse
import csv
import hashlib
import math
import random
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import chi2, poisson

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from allomix.chimerism import PanelCalibration, estimate_single_donor_bb  # noqa: E402
from allomix.genotype import classify_markers, parse_vcf  # noqa: E402
from allomix.simulate import (  # noqa: E402
    blend_vcfs,
    generate_related_genotypes,
    write_genotype_vcf,
    write_vcf,
)

# Import LoD-fit helpers from the issue #8 sweep (reuse, do not re-derive).
sys.path.insert(0, str(Path(__file__).resolve().parent))
from run_lod_validation import compute_lob, fit_lod  # noqa: E402

# --- Sweep grid (plan section "Prototype spec") -----------------------------

RELATEDNESS_LEVELS = ["unrelated", "sibling"]
N_MARKERS_GRID = [76]
DEPTHS = [1000, 2000, 5000]
ERROR_RATES = [1e-2, 3e-3, 1e-3, 3e-4]
# f_h = 0 is the negative control; the rest are positives.
HOST_FRACTIONS = [0.0, 1e-5, 2e-5, 5e-5, 1e-4, 2e-4, 5e-4, 1e-3]
POSITIVE_FRACTIONS = [f for f in HOST_FRACTIONS if f > 0]

MAF_RANGE = (0.2, 0.5)
LOCUS_DROPOUT_RATE = 0.016
DEPTH_CV = 0.43
# Cheaper grid for the MLE — we only use it for the LoB/LoD comparison, not
# fine f-estimates.
MLE_GRID_STEPS = 201

ALPHA = 0.05

FACTS_DIR = Path("output/facts")


# --- Helpers ----------------------------------------------------------------


def derive_seed(*parts: object) -> int:
    """Deterministic, process-stable seed derived from a tuple of arbitrary parts."""
    digest = hashlib.sha256(repr(parts).encode("utf-8")).digest()[:4]
    return int.from_bytes(digest, "big")


def select_donor_hom_markers(informative_markers: list) -> list[tuple[int, int, int]]:
    """Pull (y_i, n_i, h_i) for the markers the presence test uses.

    Restricts to Vynck types 0, 1, 10, 11 (donor homozygous, host carries the
    donor-absent allele). y_i is the donor-absent-allele read count, n_i is the
    admixture depth, h_i is the host dose of the donor-absent allele (1 or 2).
    """
    out: list[tuple[int, int, int]] = []
    for m in informative_markers:
        mt = m.marker_type
        if mt == 0:
            # host 0/0, donor 1/1 -> donor-absent = REF
            out.append((m.admix_ad_ref, m.admix_dp, 2))
        elif mt == 1:
            # host 1/1, donor 0/0 -> donor-absent = ALT
            out.append((m.admix_ad_alt, m.admix_dp, 2))
        elif mt == 10:
            # host 0/1, donor 0/0 -> donor-absent = ALT, host dose 1
            out.append((m.admix_ad_alt, m.admix_dp, 1))
        elif mt == 11:
            # host 0/1, donor 1/1 -> donor-absent = REF, host dose 1
            out.append((m.admix_ad_ref, m.admix_dp, 1))
        # types 20, 21 (donor het) excluded: no donor-absent allele.
    return out


def presence_pooled_poisson(
    rows: list[tuple[int, int, int]], e_per_dir: float,
) -> tuple[int, float, float]:
    """Pooled count test under H0: y_i ~ Binomial(n_i, e_per_dir).

    Returns (Y, Lam, p_pois) where p_pois = P(Poisson(Lam) >= Y). Poisson is the
    natural approximation since e_per_dir is tiny and n_i large; with Lam small
    the survival function at Y is the right tail probability we want.
    """
    Y = sum(y for y, _, _ in rows)
    Lam = sum(n * e_per_dir for _, n, _ in rows)
    if Lam <= 0:
        # Degenerate (no usable markers); treat as non-significant.
        return Y, Lam, 1.0
    p = float(poisson.sf(Y - 1, Lam)) if Y > 0 else 1.0
    return Y, Lam, p


def presence_lrt(
    rows: list[tuple[int, int, int]], e_per_dir: float,
) -> tuple[float, float, float]:
    """Bounded-MLE LRT for f_h >= 0.

    Returns (f_h_hat, D, p_lrt) with p_lrt the chi-bar-square one-sided p-value
    (0.5 * P(chi2_1 >= D) for D > 0, else 1).
    """
    if not rows:
        return 0.0, 0.0, 1.0

    ys = np.asarray([y for y, _, _ in rows], dtype=float)
    ns = np.asarray([n for _, n, _ in rows], dtype=float)
    hs = np.asarray([h for _, _, h in rows], dtype=float)
    # Marker contribution to q_i is (h_i/2). Multiplying by f_h and adding the
    # constant background e_per_dir gives q_i(f_h).
    coef = hs / 2.0

    def loglik(f_h: float) -> float:
        q = e_per_dir + coef * f_h
        # Clip to a safe range so log() is finite even when q is at the boundary.
        q = np.clip(q, 1e-15, 1.0 - 1e-12)
        return float(np.sum(ys * np.log(q) + (ns - ys) * np.log1p(-q)))

    ll0 = loglik(0.0)

    # Bracket: f_h in [0, 1]. The LRT is concave in f_h here; use a bounded
    # minimiser of -loglik.
    res = minimize_scalar(
        lambda fh: -loglik(fh),
        bounds=(0.0, 1.0),
        method="bounded",
        options={"xatol": 1e-8},
    )
    f_hat = float(max(0.0, min(1.0, res.x)))
    ll_hat = -float(res.fun)

    # Bound at f_h = 0: if the unconstrained maximiser is essentially 0, the LRT
    # collapses to 0.
    if ll_hat <= ll0 + 1e-9:
        return 0.0, 0.0, 1.0

    D = 2.0 * (ll_hat - ll0)
    if D <= 0:
        return f_hat, 0.0, 1.0
    p_lrt = 0.5 * float(chi2.sf(D, 1))
    return f_hat, D, p_lrt


# --- Per-replicate ----------------------------------------------------------


def _run_replicate(
    relatedness: str,
    n_markers: int,
    depth: int,
    error_rate: float,
    f_h: float,
    rep_idx: int,
    seed_base: int,
    workdir: Path,
    rho: float,
    rho_marker_type: str,
) -> dict:
    """One synthetic replicate.

    Builds fresh host + donor genotypes, blends to a donor_fraction of (1 - f_h)
    so the host is the minor contributor (the relapse scenario), runs the
    presence test on donor-homozygous markers, and also runs the full MLE on the
    same admixture sample for the MLE-LoD comparison.
    """
    cell_seed = derive_seed(
        "presence", relatedness, n_markers, depth, error_rate, f_h, rep_idx, seed_base,
    )

    rng = random.Random(derive_seed("gt", cell_seed))
    markers = generate_related_genotypes(
        n_markers, relatedness, rng, maf_range=MAF_RANGE,
    )

    rep_dir = workdir / f"rep_{cell_seed}"
    rep_dir.mkdir(parents=True, exist_ok=True)
    host_vcf = rep_dir / "host.vcf"
    donor_vcf = rep_dir / "donor.vcf"
    admix_vcf = rep_dir / "admix.vcf"

    write_genotype_vcf(markers, host_vcf, "host", key="host_gt")
    write_genotype_vcf(markers, donor_vcf, "donor", key="donor_gt")

    # NB: per-marker bias is intentionally OFF here. The simulator adds bias to
    # the *expected VAF* before drawing reads (see simulate.blend_vcfs ~line
    # 845): at a donor-hom-ref marker, vaf=0 + bias clamps to a real per-marker
    # ALT contribution on top of the e/3 error background, breaking the
    # "e_i = e/3" self-consistency the prototype relies on. This is the same
    # class of objection the plan raises for a uniform rho (see the
    # "Update (2026-05-28)" note in claude/20_host_presence_detection_plan.md):
    # a knob calibrated for het / intermediate markers contaminates the clean
    # near-zero background. depth_cv and locus_dropout do not add background,
    # so we keep them at the issue #8 sweep values.
    blend = blend_vcfs(
        host_path=str(host_vcf),
        donor_path=str(donor_vcf),
        donor_fraction=1.0 - f_h,
        target_depth=depth,
        sample_name="admix",
        seed=derive_seed("blend", cell_seed),
        error_rate=error_rate,
        locus_dropout_rate=LOCUS_DROPOUT_RATE,
        depth_cv=DEPTH_CV,
        realistic_biases=False,
        marker_bias_sd=0.0,
        rho=rho,
        rho_marker_type=rho_marker_type,
    )
    write_vcf(blend, admix_vcf)

    host_md = parse_vcf(str(host_vcf), min_dp=0, min_gq=0)
    donor_md = parse_vcf(str(donor_vcf), min_dp=0, min_gq=0)
    admix_md = parse_vcf(str(admix_vcf), min_dp=0, min_gq=0)

    genos = classify_markers(
        host_md, [donor_md], admix_md, min_dp=0, min_gq=0, pass_only=False,
    )

    # Presence statistic on donor-homozygous markers.
    rows_presence = select_donor_hom_markers(genos.informative)
    e_per_dir = error_rate / 3.0
    Y, Lam, p_pois = presence_pooled_poisson(rows_presence, e_per_dir)
    f_hat, D, p_lrt = presence_lrt(rows_presence, e_per_dir)

    # MLE on all informative markers (the existing fraction estimator).
    bias_dict = (
        {(c, p, r, a): b for c, p, r, a, b in blend.marker_biases}
        if blend.marker_biases is not None
        else None
    )
    if len(genos.informative) >= 1:
        mle = estimate_single_donor_bb(
            genos.informative,
            error_rate=error_rate,
            grid_steps=MLE_GRID_STEPS,
            calibration=PanelCalibration(biases=bias_dict or {}),
        )
        donor_fraction = mle.donor_fraction
        mle_host_est = 1.0 - mle.donor_fraction
    else:
        donor_fraction = float("nan")
        mle_host_est = float("nan")

    # Cleanup the per-replicate VCFs (each replicate is independent; we never
    # come back to these files).
    for p in (host_vcf, donor_vcf, admix_vcf):
        try:
            p.unlink()
        except OSError:
            pass
    try:
        rep_dir.rmdir()
    except OSError:
        pass

    return {
        "relatedness": relatedness,
        "n_markers": n_markers,
        "depth": depth,
        "error_rate": error_rate,
        "f_h": f_h,
        "rep": rep_idx,
        "seed": cell_seed,
        "n_presence_markers": len(rows_presence),
        "Y": Y,
        "Lam": Lam,
        "p_pois": p_pois,
        "f_h_hat": f_hat,
        "D": D,
        "p_lrt": p_lrt,
        "donor_fraction_mle": donor_fraction,
        "mle_host_est": mle_host_est,
    }


# --- Worker entry (pickled for ProcessPoolExecutor) -------------------------


def _worker(args: tuple) -> dict:
    return _run_replicate(*args)


# --- Cell aggregation -------------------------------------------------------


def summarise_cell(rows: list[dict]) -> dict:
    """One row per (relatedness, n_markers, depth, error_rate)."""
    rel = rows[0]["relatedness"]
    nm = rows[0]["n_markers"]
    depth = rows[0]["depth"]
    e = rows[0]["error_rate"]

    blanks = [r for r in rows if r["f_h"] == 0.0]

    # False-positive rate at the LRT's alpha=0.05.
    if blanks:
        fp_rate_lrt = sum(1 for r in blanks if r["p_lrt"] < ALPHA) / len(blanks)
        fp_rate_pois = sum(1 for r in blanks if r["p_pois"] < ALPHA) / len(blanks)
        # 95th-percentile Y/Lam (the empirical LoB on the count axis).
        y_lam_ratios = [
            (r["Y"] / r["Lam"]) for r in blanks if r["Lam"] > 0
        ]
        lob_y_lam = float(np.quantile(y_lam_ratios, 0.95)) if y_lam_ratios else float("nan")
    else:
        fp_rate_lrt = fp_rate_pois = float("nan")
        lob_y_lam = float("nan")

    # Detection rate per positive fraction (presence test, alpha=0.05 LRT).
    by_pos = defaultdict(list)
    for r in rows:
        if r["f_h"] > 0:
            by_pos[r["f_h"]].append(r)
    pos_fracs = sorted(by_pos.keys())
    det_rates = [
        (sum(1 for r in by_pos[f] if r["p_lrt"] < ALPHA) / len(by_pos[f]))
        for f in pos_fracs
    ]
    n_per_pos = [len(by_pos[f]) for f in pos_fracs]
    fit = fit_lod(pos_fracs, det_rates, weights=n_per_pos)
    presence_lod = fit[0] if fit is not None else float("nan")

    # MLE-LoB / LoD comparison on the SAME replicates.
    # LoB on the MLE's donor_fraction estimate at f_h = 0: at f_h = 0 the host
    # estimate is 1 - donor_fraction ~ 0, so we use mle_host_est at f_h = 0 as
    # blanks (estimate of host present). The LoD is the host fraction at which
    # the MLE host estimate exceeds the LoB in >= 95% of reps.
    blank_mle = [r["mle_host_est"] for r in blanks if math.isfinite(r["mle_host_est"])]
    mle_lob = compute_lob(blank_mle) if len(blank_mle) >= 2 else float("nan")
    if math.isfinite(mle_lob):
        mle_det_rates = [
            (
                sum(
                    1
                    for r in by_pos[f]
                    if math.isfinite(r["mle_host_est"]) and r["mle_host_est"] > mle_lob
                )
                / len(by_pos[f])
            )
            for f in pos_fracs
        ]
        mle_fit = fit_lod(pos_fracs, mle_det_rates, weights=n_per_pos)
        mle_lod = mle_fit[0] if mle_fit is not None else float("nan")
    else:
        mle_det_rates = [float("nan")] * len(pos_fracs)
        mle_lod = float("nan")

    return {
        "relatedness": rel,
        "n_markers": nm,
        "depth": depth,
        "error_rate": e,
        "n_blanks": len(blanks),
        "fp_rate_lrt": fp_rate_lrt,
        "fp_rate_pois": fp_rate_pois,
        "lob_y_over_lam_q95": lob_y_lam,
        "presence_lod": presence_lod,
        "mle_lob_host": mle_lob,
        "mle_lod": mle_lod,
        "pos_fractions": ";".join(f"{f:.0e}" for f in pos_fracs),
        "presence_det_rates": ";".join(f"{r:.3f}" for r in det_rates),
        "mle_det_rates": ";".join(
            f"{r:.3f}" if math.isfinite(r) else "nan" for r in mle_det_rates
        ),
    }


# --- Output writers ---------------------------------------------------------


def write_raw(rows: list[dict], path: Path) -> None:
    fields = [
        "relatedness", "n_markers", "depth", "error_rate", "f_h", "rep", "seed",
        "n_presence_markers", "Y", "Lam", "p_pois",
        "f_h_hat", "D", "p_lrt",
        "donor_fraction_mle", "mle_host_est",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})


def write_summary(summaries: list[dict], path: Path) -> None:
    fields = [
        "relatedness", "n_markers", "depth", "error_rate",
        "n_blanks", "fp_rate_lrt", "fp_rate_pois", "lob_y_over_lam_q95",
        "presence_lod", "mle_lob_host", "mle_lod",
        "pos_fractions", "presence_det_rates", "mle_det_rates",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for s in summaries:
            w.writerow({k: s.get(k, "") for k in fields})


# --- Console summary --------------------------------------------------------


def print_console_summary(summaries: list[dict]) -> None:
    print()
    print("=" * 96)
    print("Presence-test calibration + LoD summary")
    print("=" * 96)
    print(
        f"{'rel':<10} {'nm':>4} {'depth':>6} {'e':>8} "
        f"{'fp_lrt':>7} {'fp_pois':>8} "
        f"{'pLoD':>10} {'mleLoD':>10}"
    )
    for s in summaries:
        plod = s["presence_lod"]
        mlod = s["mle_lod"]
        print(
            f"{s['relatedness']:<10} {s['n_markers']:>4} {s['depth']:>6} "
            f"{s['error_rate']:>8.1e} "
            f"{s['fp_rate_lrt']:>7.3f} {s['fp_rate_pois']:>8.3f} "
            f"{plod:>10.2e} {mlod:>10.2e}"
        )
    print("=" * 96)


# --- Driver -----------------------------------------------------------------


def build_tasks(
    relatedness_levels: list[str],
    n_markers_grid: list[int],
    depths: list[int],
    error_rates: list[float],
    fractions: list[float],
    n_blanks: int,
    n_positives: int,
    seed: int,
    workdir: Path,
    rho: float,
    rho_marker_type: str,
) -> list[tuple]:
    """Build a flat task list: one tuple per replicate."""
    tasks: list[tuple] = []
    for rel in relatedness_levels:
        for nm in n_markers_grid:
            for depth in depths:
                for e in error_rates:
                    for f_h in fractions:
                        n_reps = n_blanks if f_h == 0.0 else n_positives
                        for rep_idx in range(n_reps):
                            tasks.append(
                                (rel, nm, depth, e, f_h, rep_idx, seed, workdir,
                                 rho, rho_marker_type)
                            )
    return tasks


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--pilot", action="store_true",
        help="Pilot: 10 reps everywhere, reduced grid (depth=2000, e in {1e-2, 3e-4}).",
    )
    parser.add_argument("--n-blanks", type=int, default=200)
    parser.add_argument("--n-positives", type=int, default=60)
    parser.add_argument("--n-workers", "-j", type=int, default=1)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--relatedness", nargs="+", default=RELATEDNESS_LEVELS,
        help=f"Subset of {RELATEDNESS_LEVELS}.",
    )
    parser.add_argument(
        "--n-markers", type=int, nargs="+", default=N_MARKERS_GRID,
    )
    parser.add_argument("--depths", type=int, nargs="+", default=None)
    parser.add_argument(
        "--error-rates", type=float, nargs="+", default=None,
    )
    parser.add_argument(
        "--host-fractions", type=float, nargs="+", default=None,
        help="Override the host-fraction grid (must include 0.0 as the blank). "
             "Default is the very-low grid used for the gate evidence; pass a "
             "wider grid (e.g. up to 0.05) to resolve the LoD across the "
             "panel-size / low-depth cells that match fig5_lod_curves.",
    )
    parser.add_argument(
        "--out-raw", default=str(FACTS_DIR / "presence_lod_raw.csv"),
    )
    parser.add_argument(
        "--out-summary", default=str(FACTS_DIR / "presence_lod_summary.csv"),
    )
    parser.add_argument("--workdir", default=None,
                        help="Per-replicate temp dir (default: tempfile.mkdtemp).")
    parser.add_argument(
        "--rho", type=float, default=float("inf"),
        help="Beta-binomial overdispersion (inf = binomial, default). Use 100 "
             "for the realistic value identified in claude/21_*.",
    )
    parser.add_argument(
        "--rho-marker-type", choices=["all", "het_only"], default="all",
        help="Where to apply rho. 'all' (default) overdisperses every marker; "
             "'het_only' applies rho only at intermediate VAF, keeping the "
             "donor-absent allele background binomial. Use 'het_only' with "
             "--rho 100 for the honest gate-3 check.",
    )
    args = parser.parse_args(argv)

    if args.pilot:
        depths = [2000]
        error_rates = [1e-2, 3e-4]
        n_blanks = 10
        n_positives = 10
    else:
        depths = args.depths if args.depths is not None else DEPTHS
        error_rates = (
            args.error_rates if args.error_rates is not None else ERROR_RATES
        )
        n_blanks = args.n_blanks
        n_positives = args.n_positives

    host_fractions = (
        args.host_fractions if args.host_fractions is not None else HOST_FRACTIONS
    )
    if 0.0 not in host_fractions:
        parser.error("--host-fractions must include 0.0 (the negative-control blank)")
    positive_fractions = [f for f in host_fractions if f > 0]

    FACTS_DIR.mkdir(parents=True, exist_ok=True)
    if args.workdir is not None:
        workdir = Path(args.workdir)
        workdir.mkdir(parents=True, exist_ok=True)
    else:
        workdir = Path(tempfile.mkdtemp(prefix="presence_lod_"))

    tasks = build_tasks(
        args.relatedness, args.n_markers, depths, error_rates,
        host_fractions, n_blanks, n_positives, args.seed, workdir,
        args.rho, args.rho_marker_type,
    )

    print(
        f"Cells: rel={args.relatedness} nm={args.n_markers} depths={depths} "
        f"e={error_rates}",
        file=sys.stderr,
    )
    print(
        f"Overdispersion: rho={args.rho} rho_marker_type={args.rho_marker_type}",
        file=sys.stderr,
    )
    print(
        f"Replicates per cell: blanks={n_blanks}, positives={n_positives} "
        f"(at each of {len(positive_fractions)} positive fractions)",
        file=sys.stderr,
    )
    print(f"Total replicate tasks: {len(tasks)}", file=sys.stderr)
    print(f"Per-replicate workdir: {workdir}", file=sys.stderr)

    rows: list[dict] = []
    if args.n_workers <= 1:
        for i, t in enumerate(tasks, 1):
            rows.append(_worker(t))
            if i % 200 == 0 or i == len(tasks):
                print(f"  [{i}/{len(tasks)}] done", file=sys.stderr)
    else:
        with ProcessPoolExecutor(max_workers=args.n_workers) as pool:
            futures = [pool.submit(_worker, t) for t in tasks]
            done = 0
            for fut in as_completed(futures):
                rows.append(fut.result())
                done += 1
                if done % 200 == 0 or done == len(futures):
                    print(f"  [{done}/{len(futures)}] done", file=sys.stderr)

    rows.sort(key=lambda r: (
        r["relatedness"], r["n_markers"], r["depth"], r["error_rate"],
        r["f_h"], r["rep"],
    ))
    write_raw(rows, Path(args.out_raw))
    print(f"Wrote {args.out_raw} ({len(rows)} rows)", file=sys.stderr)

    # Aggregate per cell.
    by_cell: dict[tuple, list[dict]] = defaultdict(list)
    for r in rows:
        by_cell[
            (r["relatedness"], r["n_markers"], r["depth"], r["error_rate"])
        ].append(r)
    summaries = [summarise_cell(cell_rows) for cell_rows in by_cell.values()]
    summaries.sort(key=lambda s: (
        s["relatedness"], s["n_markers"], s["depth"], s["error_rate"],
    ))
    write_summary(summaries, Path(args.out_summary))
    print(f"Wrote {args.out_summary} ({len(summaries)} rows)", file=sys.stderr)

    print_console_summary(summaries)

    return 0


if __name__ == "__main__":
    sys.exit(main())
