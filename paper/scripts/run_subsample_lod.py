#!/usr/bin/env python3
"""Real-data LoD curves by sub-sampling the SRP434573 public mixtures (issue #24).

The simulated LoD figure (``run_lod_validation.py`` -> Figure 1/5) is a best-case
analytical number on near-binomial simulated data. This script reruns the same
panel-size / depth LoD sweep on real reads: it throws away markers and reads from
the high-depth SRP434573 titrated mixtures until the LoD rises into the
measurable 0.5-10% window, keeping what a simulator cannot fully reproduce (true
per-marker capture bias, real between-marker overdispersion, and this dataset's
known co-pooled contamination floor).

Two complementary readouts run in the same pass over every cell; neither is
primary:

  - Magnitude (MLE) LoD: the beta-binomial donor/host-fraction estimate
    (``estimate_single_donor_bb``). A cell is detected when the 95% CI for the
    host (minor) fraction excludes 0 (``donor_fraction_ci[1] < 1``).
  - Presence-test LoD: the donor-homozygous residual-host detector
    (``host_presence_test``). A cell is detected when ``lrt_pval < 0.05``.

Both use blank-free per-sample detection rules, so neither needs the EP17/blanks
construction (real data has no clean zero-analyte sample). For each test the LoD
is the lowest titration fraction with >=95% detection, read from the same
2-parameter logistic fit ``run_lod_validation.py`` uses; the median across
mixtures is the curve and the 10th-90th percentile across mixtures is the band.

A bonus information-theoretic curve (``test=mle_analytical``) aggregates the
per-sample analytical ``lod_fraction`` (Fisher information on the real markers'
real bias and overdispersion) that ``ChimerismResult`` already exposes; it is
floor-independent and is overlaid faintly on the MLE figure as a consistency
check.

Sub-sampling mechanics (plan section 4):
  - Depth: one global binomial keep-rate per sample,
    ``rate = min(1, target_mean_depth / observed_mean_depth)``, applied uniformly
    to every marker (``allomix.simulate.thin_informative_markers``). This is the
    exact statistical analog of FASTQ/BAM read subsampling for this dataset and
    preserves the real locus-to-locus depth CV. ``min_dp`` is re-applied AFTER
    thinning (low-depth locus dropout).
  - Panel size: nested random subsets (prefixes of one permutation per
    (mixture, seed)) of the INFORMATIVE markers, so each curve is monotone.

Genotype VCFs come from the committed snapshot
(``paper/public_data/SRP434573/genotypes``) unless a freshly joint-called
``output/genotypes/SRP434573`` is present. Writes nothing to /tau.

Outputs (under output/facts/), a ``test`` column keeps the analyses separable:
  output/facts/subsample_lod_raw.csv          # one row per (mixture, sample,
                                              #   fraction, depth, n_markers, seed)
  output/facts/subsample_lod_per_mixture.csv  # per (test, mixture, depth, n_markers) LoD
  output/facts/subsample_lod_summary.csv      # per (test, depth, n_markers): median + 10/90 band
  output/facts/subsample_lod_headline.csv     # named facts for later prose

Usage:
    python paper/scripts/run_subsample_lod.py
    python paper/scripts/run_subsample_lod.py --n-seeds 10 --n-workers 8
"""

import argparse
import csv
import math
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from cyvcf2 import VCF

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from allomix.estimate.chimerism import estimate_single_donor_bb  # noqa: E402
from allomix.genotype import DEFAULT_MIN_DP, classify_markers, parse_vcf  # noqa: E402
from allomix.qc.host_presence import host_presence_test  # noqa: E402
from allomix.simulate import thin_informative_markers  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from paper_quick import qval  # noqa: E402
from run_lod_validation import derive_seed, fit_lod  # noqa: E402
from run_srp434573_allomix import MIXES, known_host_pct  # noqa: E402
from srp434573_common import resolve_srp434573_genotypes_dir  # noqa: E402

# Sweep grid (plan section 6b).
# Quick-build mode (ALLOMIX_PAPER_QUICK=1 / --config quick=1) shrinks the depth,
# panel-size, and seed grids (tens of minutes -> ~a minute; watermarked figures,
# not for publication). Grids are sized like run_lod_validation.py's (the sweep
# runs both estimators in one pass, so the quick budget is not doubled).
DEPTHS = qval([100, 250, 500, 1000, 2000], [250, 1000])
N_MARKERS_GRID = qval([25, 50, 75, 100, 200, 400], [50, 200])
# Titration fractions are fixed by the SRP434573 dilution series (minor = host),
# expressed in percent. RATIO_PCT in run_srp434573_allomix maps admix aliases to
# these same values.
FRACTIONS_PCT = [10.0, 5.0, 2.5, 1.25, 1.0, 0.5]
# Resampling replicates per (mixture, fraction-sample): each seed redraws the
# panel permutation and the binomial thinning, giving pseudo-replicates from one
# real read draw (plan section 4c).
DEFAULT_N_SEEDS = qval(20, 3)

# Across-mixture band reported around the median LoD curve.
BAND_LO_Q = 0.10
BAND_HI_Q = 0.90

FACTS_DIR = Path("output/facts")

# Two-person mixtures only: the LoD sweep is the real analog of Figure 1's
# unrelated left panel. The three-person mix (two donors) is excluded.
TWO_PERSON_MIXES = {
    name: (host, donors[0]) for name, (host, donors) in MIXES.items() if len(donors) == 1
}

# Sentinels for per-mixture LoD edge cases (mirrors run_lod_validation.py).
LOD_BELOW_RANGE = -1.0
LOD_ABOVE_RANGE = float("inf")
MIN_POS_FRACTION = min(FRACTIONS_PCT) / 100.0
MAX_FRACTION = max(FRACTIONS_PCT) / 100.0

# The three tests written to the facts files.
TESTS = ("mle", "presence", "mle_analytical")


def rate_for_mean_depth(markers: list, target_mean_depth: float) -> float:
    """Global binomial keep-rate to hit ``target_mean_depth`` for one sample.

    Computed from the FULL informative set's observed mean admix depth (a stable
    per-sample quantity), then applied to whichever panel subset is drawn.
    Returns ``min(1.0, target / observed)``; 1.0 (no thinning) when the target
    is at or above the observed mean.
    """
    d_obs = float(np.mean([m.admix_dp for m in markers]))
    if d_obs <= 0:
        return 1.0
    return min(1.0, target_mean_depth / d_obs)


def load_mixture(name: str, host: str, donor: str) -> dict[float, list[tuple[str, list]]]:
    """Parse one two-person mixture into informative markers per titration sample.

    ``host`` is the minor (titrated) contributor, ``donor`` the major
    (background) contributor. Returns ``fraction_pct -> [(admix_sample_name,
    informative_markers), ...]``; each admix sample (including v1/v2 repeats)
    contributes one entry under its known minor (= host) percent.
    """
    gen = resolve_srp434573_genotypes_dir()
    panel = gen / f"{name}.SRP434573.vcf.gz"
    admix = gen / f"{name}.admix.vcf.gz"
    host_md = parse_vcf(panel, sample=host, min_gq=0, gt_ad_consistency=True)
    donor_md = parse_vcf(panel, sample=donor, min_gq=0, gt_ad_consistency=True)

    by_fraction: dict[float, list[tuple[str, list]]] = defaultdict(list)
    for sample in VCF(str(admix)).samples:
        pct = known_host_pct(sample)
        if pct is None:
            continue
        admix_md = parse_vcf(admix, sample=sample, min_dp=0)
        # min_dp=0: the depth filter is applied after thinning, so the rate is
        # computed from the complete informative set.
        genos = classify_markers(host_md, [donor_md], admix_md, min_dp=0)
        if genos.informative:
            by_fraction[pct].append((sample, genos.informative))
    return dict(by_fraction)


def cell_results(markers_thinned_full: list, n_markers: int) -> dict:
    """Run both estimators on one (already-thinned) panel prefix.

    ``n_markers`` is the nominal panel size (before low-depth dropout); a
    length-``n_markers`` prefix of the globally-thinned permutation is taken.
    Returns a dict with both readouts and ``n_used`` (markers surviving min_dp).
    """
    subset = markers_thinned_full[:n_markers]
    thinned = [m for m in subset if m.admix_dp >= DEFAULT_MIN_DP]
    if not thinned:
        return {
            "n_used": 0,
            "mle_host_frac": float("nan"),
            "mle_detected": False,
            "mle_analytical_lod": float("nan"),
            "presence_f": float("nan"),
            "presence_detected": False,
        }
    mle = estimate_single_donor_bb(thinned)
    host_ci_lo = 1.0 - mle.donor_fraction_ci[1]
    pres = host_presence_test(thinned)
    return {
        "n_used": len(thinned),
        "mle_host_frac": mle.host_fraction,
        "mle_detected": host_ci_lo > 0.0,
        "mle_analytical_lod": mle.lod_fraction,
        "presence_f": pres.f_host_mle,
        "presence_detected": pres.lrt_pval < 0.05,
    }


def run_mixture(name: str, host: str, donor: str, n_seeds: int) -> list[dict]:
    """Sweep one mixture over (sample, fraction, depth, n_markers, seed).

    One panel permutation per (sample, seed), shared across depths so panel
    composition is stable; per depth the ordered list is thinned once at that
    depth's global rate, so the n_markers grid is strict prefixes (monotone).

    Returns one raw row per (mixture, sample, fraction, depth, n_markers, seed).
    """
    by_fraction = load_mixture(name, host, donor)
    rows: list[dict] = []
    for pct, samples in by_fraction.items():
        for sample, markers_full in samples:
            n_inf = len(markers_full)
            panel_sizes = [n for n in N_MARKERS_GRID if n <= n_inf] or [n_inf]
            for seed in range(n_seeds):
                perm_rng = np.random.default_rng(derive_seed("perm", name, sample, seed))
                order = perm_rng.permutation(n_inf)
                ordered = [markers_full[i] for i in order]
                for depth in DEPTHS:
                    rate = rate_for_mean_depth(markers_full, depth)
                    thin_rng = np.random.default_rng(derive_seed("thin", name, sample, depth, seed))
                    thinned_full = thin_informative_markers(ordered, rate, thin_rng)
                    for n_markers in panel_sizes:
                        res = cell_results(thinned_full, n_markers)
                        rows.append(
                            {
                                "mixture": name,
                                "sample": sample,
                                "fraction_pct": pct,
                                "depth": depth,
                                "n_markers": n_markers,
                                "seed": seed,
                                "n_used": res["n_used"],
                                "mle_host_frac": res["mle_host_frac"],
                                "mle_detected": int(res["mle_detected"]),
                                "presence_f": res["presence_f"],
                                "presence_detected": int(res["presence_detected"]),
                                "mle_analytical_lod": res["mle_analytical_lod"],
                            }
                        )
    return rows


def _empirical_pair_lod(rows: list[dict], detected_key: str) -> dict:
    """LoD for ONE mixture at one (depth, n_markers) cell, for one empirical test.

    Detection rate per titration fraction over the cell's (sample x seed)
    replicates, fed to the shared logistic fit. Returns ``lod`` in fraction units
    (or a sentinel), a ``note``, and the mean ``n_used``.
    """
    by_frac: dict[float, list[int]] = defaultdict(list)
    n_used: list[int] = []
    for r in rows:
        by_frac[r["fraction_pct"]].append(int(r[detected_key]))
        n_used.append(r["n_used"])

    fractions = sorted(by_frac)  # ascending percent
    fr_units = [f / 100.0 for f in fractions]
    rates = [float(np.mean(by_frac[f])) for f in fractions]
    weights = [len(by_frac[f]) for f in fractions]

    fit = fit_lod(fr_units, rates, weights)
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
        if lod > MAX_FRACTION:
            # 95% crossing extrapolates beyond the largest probed fraction
            # (detection never reached 95% in range). Treat as above-range
            # (conservative ceiling), not a far-extrapolated value.
            lod, note = LOD_ABOVE_RANGE, "above_range_extrap"

    return {
        "lod": lod,
        "note": note,
        "mean_n_used": float(np.mean(n_used)) if n_used else 0.0,
        # This mixture's own lowest probed fraction. Used to clamp an
        # all-detected sentinel to a bound we actually tested, rather than the
        # global grid minimum (which a mixture missing the lowest titration
        # point never reached).
        "min_fraction": min(fr_units) if fr_units else MIN_POS_FRACTION,
    }


def _analytical_pair_lod(rows: list[dict]) -> dict:
    """Per-mixture analytical LoD at one (depth, n_markers) cell.

    The analytical ``lod_fraction`` is a property of the sample's markers (depth,
    panel, overdispersion), not of the titration fraction, so it is pooled across
    all fractions/samples/seeds in the cell and the median is taken.
    """
    vals = [
        r["mle_analytical_lod"]
        for r in rows
        if r["mle_analytical_lod"] is not None and math.isfinite(r["mle_analytical_lod"])
    ]
    n_used = [r["n_used"] for r in rows]
    lod = float(np.median(vals)) if vals else float("nan")
    return {
        "lod": lod,
        "note": "" if vals else "no_finite_analytical",
        "mean_n_used": float(np.mean(n_used)) if n_used else 0.0,
        # The analytical LoD is a real number from Fisher information, never an
        # all-detected sentinel, so it carries no own-floor clamp.
        "min_fraction": None,
    }


def _lod_for_aggregation(ps: dict) -> tuple[float, bool] | None:
    """Map a per-mixture LoD (possibly a sentinel) to a (fraction, censored) pair.

    ``all_detected`` is an upper bound: every probed fraction was detected, so
    the true LoD is at or below this mixture's OWN lowest probed fraction (not
    the global grid minimum, which a mixture missing the lowest titration point
    never reached). It is returned as that fraction, flagged censored.
    ``none_detected`` clamps to the largest probed fraction (a conservative
    ceiling). NaN is dropped.
    """
    lod = ps["lod"]
    if lod == LOD_BELOW_RANGE:
        return (ps.get("min_fraction") or MIN_POS_FRACTION, True)
    if lod == LOD_ABOVE_RANGE:
        return (MAX_FRACTION, False)
    if lod is None or math.isnan(lod):
        return None
    return (lod, False)


def _to_pct(x: float) -> float:
    """Fraction -> percent, preserving sentinels and NaN."""
    if x == LOD_BELOW_RANGE:
        return -1.0
    if x == LOD_ABOVE_RANGE:
        return float("inf")
    if x is None or math.isnan(x):
        return float("nan")
    return round(x * 100, 4)


def summarise_cell(test: str, depth: int, n_markers: int, pair_summaries: list[dict]) -> dict:
    """Aggregate per-mixture LoDs into a median curve point and a 10-90% band.

    The cell is flagged ``censored`` when its median lands on a mixture that
    only yielded an upper bound (all_detected). The real LoD is then below the
    dilution grid at that point, so the plotted value must be read as
    "<= lod_pct", not a resolved estimate.
    """
    pairs = [p for ps in pair_summaries if (p := _lod_for_aggregation(ps)) is not None]
    n_used = [ps["mean_n_used"] for ps in pair_summaries]
    n_dropped = len(pair_summaries) - len(pairs)

    if pairs:
        vals = np.asarray([v for v, _ in pairs], dtype=float)
        cens = np.asarray([c for _, c in pairs], dtype=bool)
        lod_med = float(np.median(vals))
        lod_lo = float(np.quantile(vals, BAND_LO_Q))
        lod_hi = float(np.quantile(vals, BAND_HI_Q))
        # The median is itself an upper bound when the mixture(s) at the median
        # position only gave "<= own floor" sentinels.
        cens_sorted = cens[np.argsort(vals)]
        n = len(vals)
        cell_censored = (
            bool(cens_sorted[n // 2])
            if n % 2
            else bool(cens_sorted[n // 2 - 1] or cens_sorted[n // 2])
        )
    else:
        lod_med = lod_lo = lod_hi = float("nan")
        cell_censored = False

    note = ""
    if not pairs:
        note = "no_mixtures_fit"
    elif n_dropped:
        note = f"{n_dropped}_mixtures_dropped"

    return {
        "test": test,
        "depth": depth,
        "n_markers": n_markers,
        "lod_pct": _to_pct(lod_med),
        "lod_pct_ci_lo": _to_pct(lod_lo),
        "lod_pct_ci_hi": _to_pct(lod_hi),
        "censored": cell_censored,
        "median_n_used": round(float(np.median(n_used)), 1) if n_used else 0.0,
        "n_mixtures": len(pair_summaries),
        "n_mixtures_used": len(pairs),
        "n_mixtures_dropped": n_dropped,
        "note": note,
    }


def _monotonize_over_panels(by_nm: dict[int, dict]) -> None:
    """Enforce non-increasing per-mixture LoD across ascending panel size, in place.

    Panels are nested (larger is a superset, one shared permutation), so true
    LoD cannot rise as markers are added; a fitted LoD above the running minimum
    is estimation noise and is clamped to it. Sentinels order below-range (best,
    dilution-grid floor) < finite LoD < above-range (worst); NaN cells (fit
    failed) are left untouched and do not update the running minimum.
    """
    best = math.inf
    for nm in sorted(by_nm):
        ps = by_nm[nm]
        lod = ps["lod"]
        if lod == LOD_BELOW_RANGE:
            v = -math.inf
        elif lod == LOD_ABOVE_RANGE:
            v = math.inf
        elif lod is None or (isinstance(lod, float) and math.isnan(lod)):
            continue
        else:
            v = lod
        best = min(best, v)
        if best == -math.inf:
            ps["lod"] = LOD_BELOW_RANGE
        elif best == math.inf:
            ps["lod"] = LOD_ABOVE_RANGE
        else:
            ps["lod"] = best


def write_raw(rows: list[dict], path: Path) -> None:
    fields = [
        "mixture",
        "sample",
        "fraction_pct",
        "depth",
        "n_markers",
        "seed",
        "n_used",
        "mle_host_frac",
        "mle_detected",
        "presence_f",
        "presence_detected",
        "mle_analytical_lod",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})


def write_per_mixture(rows: list[dict], path: Path) -> None:
    fields = ["test", "mixture", "depth", "n_markers", "lod_pct", "note", "mean_n_used"]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in fields})


def write_summary(rows: list[dict], path: Path) -> None:
    fields = [
        "test",
        "mixture_set",
        "depth",
        "n_markers",
        "lod_pct",
        "lod_pct_ci_lo",
        "lod_pct_ci_hi",
        "censored",
        "median_n_used",
        "n_mixtures",
        "n_mixtures_used",
        "n_mixtures_dropped",
        "note",
    ]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for s in rows:
            w.writerow({k: s.get(k, "") for k in fields})


def _contamination_floor_pct(raw_rows: list[dict]) -> float:
    """Empirical level the MLE host-fraction estimate bottoms out at.

    Median MLE host-fraction estimate at the lowest titration fraction in the
    richest cell (deepest depth, largest panel) across mixtures. This is the
    concrete data version of the co-pooled contamination floor that makes the
    MLE "CI excludes 0" rule read low at the bottom of the titration (plan
    section 5). Reported descriptively; not a claim that this is allomix's
    intrinsic limit.
    """
    min_frac = min(FRACTIONS_PCT)
    best_depth = max(DEPTHS)
    best_nm = max(N_MARKERS_GRID)
    vals = [
        r["mle_host_frac"]
        for r in raw_rows
        if r["fraction_pct"] == min_frac
        and r["depth"] == best_depth
        and r["n_markers"] == best_nm
        and r["mle_host_frac"] is not None
        and math.isfinite(r["mle_host_frac"])
    ]
    if not vals:
        return float("nan")
    return round(float(np.median(vals)) * 100.0, 4)


def write_headline(
    summaries: list[dict], raw_rows: list[dict], path: Path, n_mixtures: int, n_seeds: int
) -> None:
    """One-row CSV of cells/values later prose may quote."""

    def lookup(test: str, d: int, nm: int) -> float:
        for s in summaries:
            if (
                s["test"] == test
                and s["depth"] == d
                and s["n_markers"] == nm
                and s.get("mixture_set", "all") == "all"
            ):
                return s["lod_pct"]
        return float("nan")

    headline = {
        "n_mixtures": n_mixtures,
        "n_seeds": n_seeds,
        "n_fractions": len(FRACTIONS_PCT),
        "min_titration_pct": min(FRACTIONS_PCT),
        "mle_lod_1000x_100markers_pct": lookup("mle", 1000, 100),
        "presence_lod_1000x_100markers_pct": lookup("presence", 1000, 100),
        "mle_analytical_lod_1000x_100markers_pct": lookup("mle_analytical", 1000, 100),
        "mle_lod_2000x_400markers_pct": lookup("mle", 2000, 400),
        "presence_lod_2000x_400markers_pct": lookup("presence", 2000, 400),
        "mle_lod_250x_50markers_pct": lookup("mle", 250, 50),
        "presence_lod_250x_50markers_pct": lookup("presence", 250, 50),
        "contamination_floor_pct": _contamination_floor_pct(raw_rows),
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


def _read_raw_csv(path: Path) -> list[dict]:
    """Load a previously written subsample_lod_raw.csv back into typed rows.

    Lets aggregation/summary/plotting be recomputed from an existing sweep
    without repeating the expensive read sub-sampling. Column types mirror
    ``run_mixture``'s output.
    """

    def _f(x: str) -> float:
        return float(x) if x not in ("", "nan", "None") else float("nan")

    rows: list[dict] = []
    with open(path, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rows.append(
                {
                    "mixture": r["mixture"],
                    "sample": r["sample"],
                    "fraction_pct": float(r["fraction_pct"]),
                    "depth": int(r["depth"]),
                    "n_markers": int(r["n_markers"]),
                    "seed": int(r["seed"]),
                    "n_used": int(r["n_used"]),
                    "mle_host_frac": _f(r["mle_host_frac"]),
                    "mle_detected": int(r["mle_detected"]),
                    "presence_f": _f(r["presence_f"]),
                    "presence_detected": int(r["presence_detected"]),
                    "mle_analytical_lod": _f(r["mle_analytical_lod"]),
                }
            )
    return rows


def _worker(arg: tuple) -> list[dict]:
    name, host, donor, n_seeds = arg
    return run_mixture(name, host, donor, n_seeds)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--n-seeds",
        type=int,
        default=DEFAULT_N_SEEDS,
        help="Resampling replicates per (mixture, fraction-sample).",
    )
    parser.add_argument(
        "--n-workers", "-j", type=int, default=1, help="Process pool size (one task per mixture)."
    )
    parser.add_argument(
        "--from-raw",
        type=str,
        default=None,
        help="Recompute per-mixture, summary, and headline CSVs from an existing "
        "subsample_lod_raw.csv, skipping the (expensive) read sub-sampling sweep.",
    )
    args = parser.parse_args(argv)

    FACTS_DIR.mkdir(parents=True, exist_ok=True)

    if args.from_raw:
        all_rows = _read_raw_csv(Path(args.from_raw))
        n_seeds_eff = len({r["seed"] for r in all_rows})
        print(
            f"Loaded {len(all_rows)} raw rows from {args.from_raw}; skipping sweep "
            f"(n_seeds={n_seeds_eff})",
            file=sys.stderr,
        )
    else:
        n_seeds_eff = args.n_seeds
        tasks = [
            (name, host, donor, args.n_seeds) for name, (host, donor) in TWO_PERSON_MIXES.items()
        ]
        print(
            f"Sweep: mixtures={len(tasks)}, depths={DEPTHS}, n_markers={N_MARKERS_GRID}, "
            f"fractions={FRACTIONS_PCT}, n_seeds={args.n_seeds}",
            file=sys.stderr,
        )

        all_rows = []
        if args.n_workers <= 1:
            for i, t in enumerate(tasks, 1):
                all_rows.extend(_worker(t))
                print(f"  [{i}/{len(tasks)}] {t[0]} done", file=sys.stderr)
        else:
            with ProcessPoolExecutor(max_workers=args.n_workers) as pool:
                futures = {pool.submit(_worker, t): t[0] for t in tasks}
                done = 0
                for fut in as_completed(futures):
                    all_rows.extend(fut.result())
                    done += 1
                    print(f"  [{done}/{len(futures)}] {futures[fut]} done", file=sys.stderr)

        all_rows.sort(
            key=lambda r: (
                r["mixture"],
                r["sample"],
                r["n_markers"],
                r["depth"],
                r["fraction_pct"],
                r["seed"],
            )
        )

        write_raw(all_rows, FACTS_DIR / "subsample_lod_raw.csv")
        print(
            f"Wrote {FACTS_DIR / 'subsample_lod_raw.csv'} ({len(all_rows)} rows)", file=sys.stderr
        )

    # Per-mixture LoD: group rows by (mixture, depth, n_markers), fit one LoD per
    # mixture per empirical test; the analytical test pools the cell directly.
    by_cell: dict[tuple, list[dict]] = defaultdict(list)
    for r in all_rows:
        by_cell[(r["mixture"], r["depth"], r["n_markers"])].append(r)

    # One LoD per (test, mixture, depth, n_markers) cell, grouped so the panel
    # axis can be monotonized per mixture before aggregating across mixtures.
    ps_by_group: dict[tuple, dict[int, dict]] = defaultdict(dict)
    mixtures_seen = set()
    for (mixture, depth, nm), rows in by_cell.items():
        mixtures_seen.add(mixture)
        for test in TESTS:
            if test == "mle":
                ps = _empirical_pair_lod(rows, "mle_detected")
            elif test == "presence":
                ps = _empirical_pair_lod(rows, "presence_detected")
            else:
                ps = _analytical_pair_lod(rows)
            ps["mixture"] = mixture
            ps_by_group[(test, mixture, depth)][nm] = ps

    # Monotonize each per-mixture curve (nested panels: LoD cannot rise with
    # more markers). Median/quantiles of monotone curves are themselves monotone,
    # so the summary curve and band stay smooth without aggregate enforcement.
    for by_nm in ps_by_group.values():
        _monotonize_over_panels(by_nm)

    per_mixture: list[dict] = []
    pairs_by_summary: dict[tuple, list[dict]] = defaultdict(list)
    for (test, mixture, depth), by_nm in ps_by_group.items():
        for nm, ps in by_nm.items():
            pairs_by_summary[(test, depth, nm)].append(ps)
            per_mixture.append(
                {
                    "test": test,
                    "mixture": mixture,
                    "depth": depth,
                    "n_markers": nm,
                    "lod_pct": _to_pct(ps["lod"]),
                    "note": ps["note"],
                    "mean_n_used": round(ps["mean_n_used"], 1),
                }
            )

    per_mixture.sort(key=lambda r: (r["test"], r["n_markers"], r["depth"], r["mixture"]))
    write_per_mixture(per_mixture, FACTS_DIR / "subsample_lod_per_mixture.csv")
    print(
        f"Wrote {FACTS_DIR / 'subsample_lod_per_mixture.csv'} ({len(per_mixture)} rows)",
        file=sys.stderr,
    )

    # Mixtures whose titration series reaches the lowest probed fraction. The
    # all-mixture median is pinned above that fraction by the mixtures lacking
    # it, so a second curve over just these shows the LoD the data can resolve at
    # the bottom of the grid (overlaid in the figure, not a replacement).
    min_frac_pct = min(FRACTIONS_PCT)
    mix_fracs: dict[str, set] = defaultdict(set)
    for r in all_rows:
        mix_fracs[r["mixture"]].add(r["fraction_pct"])
    lowdf_mixtures = {m for m, fs in mix_fracs.items() if min_frac_pct in fs}
    lowdf_set_name = f"to_{min_frac_pct:g}pct"
    # Disjoint top-row set: the mixtures NOT titrated to the lowest fraction, named
    # by their own lowest reached fraction (e.g. 1%). The figure draws these two
    # subsets as its two rows so they are disjoint (top = stops at 1%, bottom =
    # reaches 0.5%), instead of top=all overlapping bottom=subset.
    highdf_mixtures = set(mix_fracs) - lowdf_mixtures
    highdf_min_pct = min((min(mix_fracs[m]) for m in highdf_mixtures), default=min_frac_pct)
    highdf_set_name = f"to_{highdf_min_pct:g}pct"
    print(
        f"Mixture subset reaching {min_frac_pct:g}%: {len(lowdf_mixtures)} of {len(mix_fracs)} "
        f"({', '.join(sorted(lowdf_mixtures))}); "
        f"top subset stopping at {highdf_min_pct:g}%: {len(highdf_mixtures)} "
        f"({', '.join(sorted(highdf_mixtures))})",
        file=sys.stderr,
    )

    summaries: list[dict] = []
    for test in TESTS:
        for depth in DEPTHS:
            for nm in N_MARKERS_GRID:
                ps_list = pairs_by_summary.get((test, depth, nm), [])
                if not ps_list:
                    continue
                # "all" retained for the headline facts (read at mixture_set=="all"),
                # not drawn; the figure uses the two disjoint subsets below.
                summaries.append({**summarise_cell(test, depth, nm, ps_list), "mixture_set": "all"})
                highset = [ps for ps in ps_list if ps.get("mixture") in highdf_mixtures]
                if highset:
                    summaries.append(
                        {**summarise_cell(test, depth, nm, highset), "mixture_set": highdf_set_name}
                    )
                subset = [ps for ps in ps_list if ps.get("mixture") in lowdf_mixtures]
                if subset:
                    summaries.append(
                        {**summarise_cell(test, depth, nm, subset), "mixture_set": lowdf_set_name}
                    )

    write_summary(summaries, FACTS_DIR / "subsample_lod_summary.csv")
    print(
        f"Wrote {FACTS_DIR / 'subsample_lod_summary.csv'} ({len(summaries)} rows)", file=sys.stderr
    )

    write_headline(
        summaries,
        all_rows,
        FACTS_DIR / "subsample_lod_headline.csv",
        n_mixtures=len(mixtures_seen),
        n_seeds=n_seeds_eff,
    )
    print(f"Wrote {FACTS_DIR / 'subsample_lod_headline.csv'}", file=sys.stderr)

    flagged = [s for s in summaries if s["note"]]
    if flagged:
        print("\nCells flagged for review (LoD edge cases):", file=sys.stderr)
        for s in flagged:
            print(
                f"  {s['test']:14s} depth={s['depth']:5d} nm={s['n_markers']:4d} "
                f"-> {s['note']} (lod={s['lod_pct']})",
                file=sys.stderr,
            )

    return 0


if __name__ == "__main__":
    sys.exit(main())
