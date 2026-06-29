#!/usr/bin/env python3
"""Validate the fast grid single-donor estimator against the exact one.

Generates a battery of synthetic informative-marker sets spanning the LoD
parameter space (panel size, true fraction, depth, with/without per-marker bias,
with/without an asymmetric error table), runs BOTH the exact estimator
(``estimate_single_donor_bb``, the default) and the fast grid estimator
(``estimate_single_donor_bb_grid``) on each, and reports the max/median absolute
donor-fraction deviation in percentage points plus the per-call wall-clock
speedup.

The marker sets are built the same way ``run_lod_validation.run_pair`` builds
them: ``generate_related_genotypes`` -> ``blend_vcfs`` -> ``classify_markers``,
sliced to nested panel prefixes. This exercises the real noise model (per-marker
capture bias, depth CV, locus dropout) rather than a hand-rolled binomial.

Usage:
    python paper/scripts/validate_grid_estimator.py
    python paper/scripts/validate_grid_estimator.py --quick
"""

import argparse
import hashlib
import random
import statistics
import sys
import tempfile
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from fast_grid_estimator import estimate_single_donor_bb_grid  # noqa: E402

from allomix.chimerism import estimate_single_donor_bb  # noqa: E402
from allomix.error_rates import MarkerErrorRates  # noqa: E402
from allomix.genotype import classify_markers, parse_vcf  # noqa: E402
from allomix.likelihood import PanelCalibration  # noqa: E402
from allomix.simulate import (  # noqa: E402
    blend_vcfs,
    generate_marker_biases_realistic,
    generate_related_genotypes,
    write_genotype_vcf,
    write_vcf,
)

# Battery grid (spans the LoD parameter space).
PANEL_SIZES = [25, 50, 100, 200, 400]
TRUE_FRACTIONS = [0.0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.5]
DEPTHS = [100, 250, 1000, 2000]
RELATEDNESS = ["unrelated", "sibling"]
ERROR_RATE = 0.01
LOCUS_DROPOUT_RATE = 0.016
DEPTH_CV = 0.43
MAF_RANGE = (0.2, 0.5)
ESTIMATOR_GRID_STEPS = 201  # matches run_lod_validation


def _seed(*parts: object) -> int:
    """Deterministic 32-bit seed from arbitrary parts (process-stable).

    Python's built-in ``hash()`` is salted per process for str/bytes (PEP 456),
    so it cannot be used for reproducible seeds across runs. SHA-256 is.
    """
    return int.from_bytes(hashlib.sha256(repr(parts).encode()).digest()[:4], "big")


def _build_error_table(markers, rng: random.Random):
    """Build an asymmetric per-marker error table for the calibration test.

    Draws small per-direction substitution rates around the global rate so the
    asymmetric likelihood path is exercised in both estimators.
    """
    errors = {}
    for m in markers:
        e_ra = max(1e-4, rng.gauss(ERROR_RATE, ERROR_RATE / 3))
        e_ar = max(1e-4, rng.gauss(ERROR_RATE, ERROR_RATE / 3))
        errors[(m.chrom, m.pos, m.ref, m.alt)] = MarkerErrorRates(e_refalt=e_ra, e_altref=e_ar)
    return errors


def run_battery(quick: bool, workdir: Path, marker_type_overdispersion: bool = False) -> None:
    depths = DEPTHS[:2] if quick else DEPTHS
    fractions = TRUE_FRACTIONS if not quick else [0.0, 0.005, 0.01, 0.05, 0.1]
    panels = PANEL_SIZES if not quick else [50, 100, 200]
    relatedness_levels = RELATEDNESS if not quick else ["sibling"]
    max_markers = max(panels)

    devs_pp: list[float] = []
    devs_pp_sub05: list[float] = []  # sub-0.5% regime, where the grid is coarsest
    exact_time = 0.0
    grid_time = 0.0
    n_calls = 0
    n_two_rho_engaged = 0
    worst = (0.0, None)

    for rel in relatedness_levels:
        for bias_on in (False, True):
            for error_on in (False, True):
                pair_seed = _seed("pair", rel, bias_on, error_on)
                gt_rng = random.Random(pair_seed)
                all_markers = generate_related_genotypes(
                    max_markers, rel, gt_rng, maf_range=MAF_RANGE
                )
                host_vcf = workdir / "host.vcf"
                donor_vcf = workdir / "donor.vcf"
                write_genotype_vcf(all_markers, host_vcf, "host", key="host_gt")
                write_genotype_vcf(all_markers, donor_vcf, "donor", key="donor_gt")
                host_md_full = parse_vcf(str(host_vcf), min_dp=0, min_gq=0)
                donor_md_full = parse_vcf(str(donor_vcf), min_dp=0, min_gq=0)

                bias_rng = random.Random(pair_seed + 1)
                all_biases = (
                    generate_marker_biases_realistic(max_markers, bias_rng)
                    if bias_on
                    else None
                )

                for depth in depths:
                    for frac in fractions:
                        blend_seed = _seed("blend", pair_seed, depth, frac)
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
                        admix_path = workdir / "admix.vcf"
                        write_vcf(blend, admix_path)
                        admix_md_full = parse_vcf(str(admix_path), min_dp=0, min_gq=0)

                        for n_markers in panels:
                            host_md = host_md_full[:n_markers]
                            donor_md = donor_md_full[:n_markers]
                            admix_md = admix_md_full[:n_markers]
                            genos = classify_markers(
                                host_md, [donor_md], admix_md,
                                min_dp=0, min_gq=0, pass_only=False,
                            )
                            if len(genos.informative) < 1:
                                continue

                            err_table = None
                            if error_on:
                                err_table = _build_error_table(
                                    genos.informative, random.Random(_seed("err", blend_seed))
                                )
                            cal = PanelCalibration(biases=bias_dict, errors=err_table)

                            t0 = time.perf_counter()
                            re = estimate_single_donor_bb(
                                genos.informative, error_rate=ERROR_RATE,
                                grid_steps=ESTIMATOR_GRID_STEPS, calibration=cal,
                                marker_type_overdispersion=marker_type_overdispersion,
                            )
                            t1 = time.perf_counter()
                            rg = estimate_single_donor_bb_grid(
                                genos.informative, error_rate=ERROR_RATE,
                                calibration=cal,
                                marker_type_overdispersion=marker_type_overdispersion,
                            )
                            t2 = time.perf_counter()

                            exact_time += t1 - t0
                            grid_time += t2 - t1
                            n_calls += 1
                            if getattr(re, "rho_hom", None) is not None:
                                n_two_rho_engaged += 1

                            dev = abs(re.donor_fraction - rg.donor_fraction) * 100.0
                            devs_pp.append(dev)
                            if frac < 0.005:
                                devs_pp_sub05.append(dev)
                            if dev > worst[0]:
                                worst = (dev, {
                                    "rel": rel, "bias": bias_on, "error": error_on,
                                    "depth": depth, "frac": frac, "n_markers": n_markers,
                                    "n_inf": len(genos.informative),
                                    "exact": re.donor_fraction, "grid": rg.donor_fraction,
                                })

    mode = "two-rho (per-marker-type)" if marker_type_overdispersion else "single-rho"
    print(f"\nBattery: {n_calls} estimator calls ({mode})")
    if marker_type_overdispersion:
        print(f"  two-rho engaged in {n_two_rho_engaged}/{n_calls} cells "
              f"(rest fell back to shared rho, both estimators)")
    print(f"  max  |f_grid - f_exact|  = {max(devs_pp):.6f} pp")
    print(f"  median |f_grid - f_exact| = {statistics.median(devs_pp):.6f} pp")
    if devs_pp_sub05:
        print(f"  max  |f_grid - f_exact|  (sub-0.5%) = {max(devs_pp_sub05):.6f} pp "
              f"({len(devs_pp_sub05)} cells)")
    print(f"  exact total time = {exact_time:.2f}s  ({exact_time / n_calls * 1000:.2f} ms/call)")
    print(f"  grid  total time = {grid_time:.2f}s  ({grid_time / n_calls * 1000:.2f} ms/call)")
    print(f"  per-call speedup = {exact_time / grid_time:.1f}x")
    # Single-rho holds the tighter 0.01 pp contract it has always met. The two-rho
    # path's target is the plan's fraction < 1e-3 (0.1 pp) over the LoD/ladder
    # regime (issue #33, step 21 plan), since the het-class rho sits near the
    # grid's coarsest bound; the sub-0.5% max above is the cell that actually
    # matters and it stays far under that.
    tol = 0.1 if marker_type_overdispersion else 0.01
    print(f"  tolerance {tol} pp: {'PASS' if max(devs_pp) < tol else 'FAIL'}")
    if worst[1] is not None:
        print(f"  worst case: {worst[0]:.6f} pp at {worst[1]}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--quick", action="store_true", help="Smaller battery.")
    parser.add_argument(
        "--marker-type-overdispersion",
        action="store_true",
        help="Validate the two-rho path (issue #33) instead of single-rho: run "
             "both estimators with per-marker-type overdispersion on.",
    )
    args = parser.parse_args(argv)

    with tempfile.TemporaryDirectory(prefix="grid_validate_") as d:
        run_battery(args.quick, Path(d), marker_type_overdispersion=args.marker_type_overdispersion)
    return 0


if __name__ == "__main__":
    sys.exit(main())
