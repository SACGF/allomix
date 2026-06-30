#!/usr/bin/env python3
"""In-silico validation of the contamination floor feeding the LoD and presence test.

This is the check for Observation 2 in ``claude/further_improvements.md``: the
in-data contamination floor allomix estimates per sample should (a) floor the
reported limit of detection, and (b) raise the per-marker background the
host-presence test compares donor-absent reads against.

``simulate.blend_vcfs`` has no contamination knob, so the two parts use different
ground truth:

  Part A (LoD floor): validated against a *synthetic contamination scalar*. We
    blend a real 0.5% host sample, take its analytical LoB/LoD from the fitted
    noise model, then apply the production flooring rule
    (``analysis._floor_detection_limits``) with a 0.2% scalar and confirm the
    reported LoD lands at or above the floor. The floor logic is what is under
    test, not the in-data estimator (that is covered in ``test_sample_contamination``).

  Part B (raised presence background, the Observation 6 gate): run on clean,
    uncontaminated high-depth blends through the real pipeline. The in-data
    contamination floor on clean data should be ~0, so the raised background must
    not suppress genuine low-fraction host signal at deployment depth (>1000x).
    We measure the in-data floor distribution, the presence-test false-positive
    rate on true blanks, and the detection rate at a true 0.5% host with the
    floor fed in vs held at zero, and confirm they agree.

Multiple independent replicates (N>=5, different seeds), per the project rule.

Usage:
    python scripts/validate_contamination_lod_floor.py
    python scripts/validate_contamination_lod_floor.py --n-reps 30 --depth 2000
"""

import argparse
import random
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from allomix.analysis import _floor_detection_limits  # noqa: E402
from allomix.chimerism import estimate_single_donor_bb  # noqa: E402
from allomix.genotype import classify_markers, parse_vcf  # noqa: E402
from allomix.host_presence import host_presence_test  # noqa: E402
from allomix.likelihood import PanelCalibration  # noqa: E402
from allomix.sample_contamination import estimate_contamination  # noqa: E402
from allomix.simulate import (  # noqa: E402
    blend_vcfs,
    generate_marker_biases_realistic,
    generate_related_genotypes,
    write_genotype_vcf,
    write_vcf,
)

MAF_RANGE = (0.2, 0.5)
ERROR_RATE = 0.01
LOCUS_DROPOUT_RATE = 0.016
DEPTH_CV = 0.43
GRID_STEPS = 201

# Synthetic contamination floor for Part A (the ~0.2% co-pooled floor measured on
# SRP434573). Part A validates the flooring rule against this scalar; the in-data
# estimator that produces it in deployment is validated separately.
SYNTHETIC_FLOOR = 0.002
# Host fraction used for the "0.5% call against a 0.2% floor" case.
LOW_HOST_FRAC = 0.005
ALPHA = 0.05


def _build_pair(rep: int, base_seed: int, n_markers: int, work_dir: Path):
    """Generate one unrelated host/donor genotype pair and parse the markers."""
    pair_dir = work_dir / f"rep{rep}"
    pair_dir.mkdir(parents=True, exist_ok=True)
    gt_rng = random.Random(base_seed * 1000 + rep)
    markers = generate_related_genotypes(n_markers, "unrelated", gt_rng, maf_range=MAF_RANGE)
    host_vcf = pair_dir / "host.vcf"
    donor_vcf = pair_dir / "donor.vcf"
    write_genotype_vcf(markers, host_vcf, "host", key="host_gt")
    write_genotype_vcf(markers, donor_vcf, "donor", key="donor_gt")
    host_md = parse_vcf(str(host_vcf), min_dp=0, min_gq=0)
    donor_md = parse_vcf(str(donor_vcf), min_dp=0, min_gq=0)
    bias_rng = random.Random(base_seed * 7919 + rep)
    biases = generate_marker_biases_realistic(n_markers, bias_rng)
    return pair_dir, host_vcf, donor_vcf, host_md, donor_md, biases


def _blend(host_vcf, donor_vcf, biases, host_frac, depth, seed, admix_path):
    """Blend one admixture sample and return parsed markers."""
    blend = blend_vcfs(
        host_path=str(host_vcf),
        donor_path=str(donor_vcf),
        donor_fraction=1.0 - host_frac,
        target_depth=depth,
        sample_name="admix",
        seed=seed,
        fixed_biases=biases,
        error_rate=ERROR_RATE,
        locus_dropout_rate=LOCUS_DROPOUT_RATE,
        depth_cv=DEPTH_CV,
    )
    write_vcf(blend, admix_path)
    return parse_vcf(str(admix_path), min_dp=0, min_gq=0)


def run_rep(rep: int, base_seed: int, n_markers: int, depth: int, work_dir: Path) -> dict:
    """One genotype pair: Part A (LoD floor) and Part B (presence background)."""
    pair_dir, host_vcf, donor_vcf, host_md, donor_md, biases = _build_pair(
        rep, base_seed, n_markers, work_dir
    )
    cal = PanelCalibration()
    admix_path = pair_dir / "admix.vcf"

    def classify(admix_md):
        return classify_markers(
            host_md[:n_markers],
            [donor_md[:n_markers]],
            admix_md[:n_markers],
            min_dp=0,
            min_gq=0,
            pass_only=False,
        )

    # Part A: a 0.5% host call, analytical LoD floored at a synthetic 0.2%.
    admix_md = _blend(
        host_vcf, donor_vcf, biases, LOW_HOST_FRAC, depth, base_seed * 31 + rep, admix_path
    )
    genos = classify(admix_md)
    res = estimate_single_donor_bb(
        genos.informative, error_rate=ERROR_RATE, grid_steps=GRID_STEPS, calibration=cal
    )
    analytical_lod = res.lod_fraction
    _floor_detection_limits(res, SYNTHETIC_FLOOR)
    floored_lod = res.lod_fraction
    floored_lob = res.lob_fraction

    # Part B: clean blank and clean 0.5% host through the real pipeline.
    # Blank (no host): in-data floor and presence FP rate.
    blank_md = _blend(host_vcf, donor_vcf, biases, 0.0, depth, base_seed * 53 + rep, admix_path)
    blank_genos = classify(blank_md)
    blank_contam = estimate_contamination(
        host_md[:n_markers],
        [donor_md[:n_markers]],
        blank_md[:n_markers],
        error_rate=ERROR_RATE,
        min_dp=0,
    )
    blank_floor = blank_contam.contamination_fraction
    blank_pres = host_presence_test(
        blank_genos.informative,
        error_rate=ERROR_RATE,
        contamination_floor=blank_floor,
        artifact_filter=False,
    )

    # True 0.5% host: power with the in-data floor fed in vs held at zero.
    host_contam = estimate_contamination(
        host_md[:n_markers],
        [donor_md[:n_markers]],
        admix_md[:n_markers],
        error_rate=ERROR_RATE,
        min_dp=0,
    )
    host_floor = host_contam.contamination_fraction
    pres_floor = host_presence_test(
        genos.informative,
        error_rate=ERROR_RATE,
        contamination_floor=host_floor,
        artifact_filter=False,
    )
    pres_nofloor = host_presence_test(
        genos.informative,
        error_rate=ERROR_RATE,
        contamination_floor=0.0,
        artifact_filter=False,
    )

    return {
        "analytical_lod": analytical_lod,
        "floored_lod": floored_lod,
        "floored_lob": floored_lob,
        "blank_floor": blank_floor,
        "blank_pres_p": blank_pres.lrt_pval,
        "host_floor": host_floor,
        "pres_floor_p": pres_floor.lrt_pval,
        "pres_nofloor_p": pres_nofloor.lrt_pval,
        "pres_floor_f": pres_floor.f_host_mle,
        "pres_nofloor_f": pres_nofloor.f_host_mle,
    }


def summarise(rows: list[dict], depth: int) -> int:
    """Print the validation tables and return a process exit code (0 = pass)."""
    n = len(rows)

    def mean(key):
        return statistics.mean(r[key] for r in rows)

    def rate(pred):
        return sum(1 for r in rows if pred(r)) / n

    print(
        f"\n=== Part A: LoD floored at synthetic {SYNTHETIC_FLOOR:.3%} "
        f"({LOW_HOST_FRAC:.1%} host, depth {depth}x) ==="
    )
    print(f"analytical LoD (mean):  {mean('analytical_lod'):.4%}")
    print(f"floored LoD   (mean):   {mean('floored_lod'):.4%}")
    print(f"floored LoB   (mean):   {mean('floored_lob'):.4%}")
    floor_holds = all(r["floored_lod"] >= SYNTHETIC_FLOOR - 1e-12 for r in rows)
    floor_binds = rate(lambda r: r["analytical_lod"] < SYNTHETIC_FLOOR)
    print(f"all reported LoD >= floor:         {floor_holds}")
    print(f"floor binds (analytical < floor):  {floor_binds:.0%} of reps")

    print("\n=== Part B: raised presence background on clean data (Obs 6 gate) ===")
    print(f"in-data floor on blanks (mean):    {mean('blank_floor'):.4%}")
    print(f"in-data floor at 0.5% host (mean): {mean('host_floor'):.4%}")
    blank_fp = rate(lambda r: r["blank_pres_p"] < ALPHA)
    print(f"presence FP rate on blanks (a={ALPHA}): {blank_fp:.2%} ({n} reps)")
    det_floor = rate(lambda r: r["pres_floor_p"] < ALPHA)
    det_nofloor = rate(lambda r: r["pres_nofloor_p"] < ALPHA)
    print(f"0.5% host detection, floor fed:    {det_floor:.0%}")
    print(f"0.5% host detection, floor=0:      {det_nofloor:.0%}")
    print(
        f"mean f_host (floor / no-floor):    "
        f"{mean('pres_floor_f'):.4%} / {mean('pres_nofloor_f'):.4%}"
    )

    ok = True
    if not floor_holds:
        print("FAIL: a reported LoD fell below the contamination floor.")
        ok = False
    # On clean high-depth data the in-data floor is ~0, so it must not destroy
    # the genuine 0.5% signal: detection with the floor should track the no-floor
    # detection closely.
    if det_nofloor - det_floor > 0.10:
        print(
            "FAIL: the in-data floor suppressed genuine 0.5% host signal "
            f"(detection {det_nofloor:.0%} -> {det_floor:.0%})."
        )
        ok = False
    print(f"\n{'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-reps", type=int, default=20)
    parser.add_argument("--n-markers", type=int, default=200)
    parser.add_argument("--base-seed", type=int, default=2026)
    parser.add_argument(
        "--depth", type=int, default=2000, help="Mean depth (deployment regime is >1000x)"
    )
    args = parser.parse_args()

    work_dir = Path("output/validate_contamination_lod_floor")
    work_dir.mkdir(parents=True, exist_ok=True)
    print(f"depth={args.depth} n_markers={args.n_markers} n_reps={args.n_reps}", file=sys.stderr)

    rows = []
    for rep in range(args.n_reps):
        rows.append(run_rep(rep, args.base_seed, args.n_markers, args.depth, work_dir))
        print(f"rep {rep + 1}/{args.n_reps} done", file=sys.stderr)

    sys.exit(summarise(rows, args.depth))


if __name__ == "__main__":
    main()
