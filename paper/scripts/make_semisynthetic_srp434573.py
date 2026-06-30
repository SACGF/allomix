#!/usr/bin/env python3
"""Generate semi-synthetic sub-0.5% mixtures for the SRP434573 paper section.

The public SRP434573 titration series bottoms out at a 0.5% minor (host)
fraction. This driver builds *lower* points (host fractions 0.1-0.5% by default)
by blending pure reference BAMs with ``samtools view --subsample`` (wrapped by
``scripts/mix_bams.sh``, which depth-normalizes on on-target reads). The resulting
reads are real (real panel noise, real GATK/bcftools path) but the mixing ratio
is artificial, so the points are *semi-synthetic* and must be labelled as such
downstream.

Two kinds of mixture are produced:

  - Two-person dilution series (one donor): host titrated along the low ladder in
    a single donor background, the sub-0.5% counterpart of the real titration.
  - Three-person host + 2 donor mixtures (the real F2/M1/M2 trio): host titrated
    along the same ladder while the two donors split the remaining background by
    each ``--donor-splits`` ratio. This exercises allomix's 2-donor capability on
    real-noise reads (issue #5, requested for double-graft monitoring).

This step runs ONCE, TAU-side, because the source BAMs live on ``/tau``. It does
two things:

  1. For every mixture, subsample+merge the pure reference BAMs at each
     (fraction, [split,] seed) into a mixed BAM under ``--bam-dir``. Each
     component is passed to ``mix_bams.sh`` with its explicit target fraction, so
     there is no minor/major argument-order trap to get wrong.
  2. Write one synthetic per-mixture CSV (``<name>.synthetic.csv``) whose
     HOST/DONOR rows reuse the pure reference BAMs (so phase-1 genotypes match the
     real run) and whose ADMIX rows point at the mixed BAMs.

It then prints the follow-on commands to run by hand (it does NOT launch GATK
itself): the existing ``pipeline/Snakefile`` over the synthetic CSVs, then a copy
of the genotype + admix VCFs into the committed snapshot at
``paper/public_data/SRP434573/genotypes_synthetic``. The paper build consumes
that snapshot (see ``paper/scripts/run_srp434573_allomix.py``).

Usage (TAU-side, where the BAMs are):

    python paper/scripts/make_semisynthetic_srp434573.py \
        --bam-dir /tau/data/chimerism/SRP434573/synthetic_bam

Inspect without touching any BAM:

    python paper/scripts/make_semisynthetic_srp434573.py --dry-run
"""

import argparse
import csv
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent
DEFAULT_CSV_DIR = REPO / "paper/public_data/SRP434573/sample_csvs"
DEFAULT_OUT_CSV_DIR = REPO / "output/semisynthetic_csv"
DEFAULT_BAM_DIR = REPO / "output/semisynthetic_bam"
DEFAULT_MIX_SCRIPT = REPO / "scripts/mix_bams.sh"
DEFAULT_PANEL_BED = REPO / "paper/public_data/SRP434573/SRP434573.bed"
SNAPSHOT_DIR = REPO / "paper/public_data/SRP434573/genotypes_synthetic"
PIPELINE_OUTPUT_DIR = "output/genotypes/SRP434573_synthetic"

# Host (minor) fractions to synthesise, as percentages. Brackets the real 0.5%
# point (which the anchored pairs also have for a synthetic-vs-real cross-check).
DEFAULT_FRACTIONS_PCT = [0.1, 0.2, 0.3, 0.4, 0.5]
DEFAULT_REPS = 5

# Three-person (host + 2 donor) mixtures reuse the real F2/M1/M2 trio. The host is
# titrated at the same low ladder; the two donors split the remaining background by
# these ratios (donor1 : donor2). "eq" is a realistic balanced double graft; the
# unequal split also stresses allomix resolving two donors at different levels. The
# realised donor percents are written into each sample name (the ground truth the
# paper runner decodes), so these ratios live in one place.
DEFAULT_DONOR_SPLITS = {
    "eq": (1.0, 1.0),
    "2to1": (2.0, 1.0),
}


def read_mix_csv(
    path: Path,
) -> tuple[tuple[str, str], list[tuple[str, str]]] | None:
    """Return ``((host_id, host_bam), [(donor_id, donor_bam), ...])`` for a CSV.

    The host is the minor (titrated) contributor, the donor(s) the majority
    background. Handles both the two-person pairs (one donor) and the three-person
    host + 2 donor mixture. Returns ``None`` for CSVs with no host or with neither
    one nor two donors.
    """
    host: tuple[str, str] | None = None
    donors: list[tuple[str, str]] = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            stype = row["sample_type"].strip().upper()
            entry = (row["sample_id"], row["bam_filename"])
            if stype == "HOST":
                host = entry
            elif stype == "DONOR":
                donors.append(entry)
    if host is None or len(donors) not in (1, 2):
        return None
    return host, donors


def sample_name(minor: str, major: str, pct: float, rep: int) -> str:
    """Two-person synthetic admix sample id, e.g. ``syn_F2-M1_f0.1_rep3``.

    The known host fraction and replicate read straight off the name, mirroring
    how the real ``1_N_X-Y`` aliases encode their fraction.
    """
    return f"syn_{minor}-{major}_f{pct:g}_rep{rep}"


def sample_name3(
    host: str, d1: str, d2: str, hpct: float, d1pct: float, d2pct: float, rep: int
) -> str:
    """Three-person synthetic admix sample id.

    e.g. ``syn3_F2-M1-M2_h0.5_d66.47-33.23_rep1``. The known host percent and the
    two donor percents (the ground truth ``run_srp434573_allomix.py`` decodes)
    read straight off the name, so no separate split table has to be kept in sync.
    """
    return f"syn3_{host}-{d1}-{d2}_h{hpct:g}_d{d1pct:g}-{d2pct:g}_rep{rep}"


def mixed_bam_path(bam_dir: Path, minor: str, major: str, pct: float, rep: int) -> Path:
    return bam_dir / f"{minor}-{major}_f{pct:g}_rep{rep}.bam"


def mix_cmd(
    mix_script: Path,
    out_bam: Path,
    panel_bed: Path,
    sid: str,
    seed: int,
    components: list[tuple[str, float]],
) -> list[str]:
    """Build a mix_bams.sh command: fixed args then (BAM, fraction) pairs."""
    cmd = [str(mix_script), str(out_bam), str(panel_bed), sid, str(seed)]
    for bam, frac in components:
        cmd += [bam, f"{frac:g}"]
    return cmd


def _emit(cmd: list[str], dry_run: bool) -> None:
    if dry_run:
        print("  mix: " + " ".join(cmd))
    else:
        subprocess.run(cmd, check=True)


def _write_csv(
    csv_path: Path,
    rows: list[tuple[str, str, str]],
    n_ref: int,
    dry_run: bool,
) -> None:
    if dry_run:
        print(f"  would write {csv_path} ({len(rows)} rows: {n_ref} reference, "
              f"{len(rows) - n_ref} ADMIX)")
        return
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample_id", "bam_filename", "sample_type"])
        w.writerows(rows)
    print(f"Wrote {csv_path} ({len(rows) - n_ref} synthetic ADMIX rows)")


def build_mix(
    patient: str,
    host: tuple[str, str],
    donors: list[tuple[str, str]],
    fractions_pct: list[float],
    reps: int,
    donor_splits: dict[str, tuple[float, float]],
    bam_dir: Path,
    out_csv_dir: Path,
    mix_script: Path,
    panel_bed: Path,
    dry_run: bool,
) -> int:
    """Mix one host/donor(s) group across the grid and write its synthetic CSV.

    One donor -> two-person dilution series (host minor fraction ladder). Two
    donors -> three-person host + 2 donor mixtures: the host is titrated at the
    same ladder while the two donors split the remaining background by each ratio
    in ``donor_splits``. Returns the number of ADMIX rows produced.
    """
    host_id, host_bam = host
    rows: list[tuple[str, str, str]] = [(host_id, host_bam, "HOST")]
    rows += [(d_id, d_bam, "DONOR") for d_id, d_bam in donors]
    n_ref = len(rows)

    for pct in fractions_pct:
        t = pct / 100.0
        for rep in range(1, reps + 1):
            # Spaced so a component's seed (seed + i) never collides across reps.
            seed = 100 * rep
            if len(donors) == 1:
                d_id, d_bam = donors[0]
                sid = sample_name(host_id, d_id, pct, rep)
                out_bam = mixed_bam_path(bam_dir, host_id, d_id, pct, rep)
                # Donor (background) first so it keeps the seed the old 2-person
                # runs used (seed); host (minor) gets seed+1 as before.
                components = [(d_bam, 1.0 - t), (host_bam, t)]
                _emit(mix_cmd(mix_script, out_bam, panel_bed, sid, seed, components),
                      dry_run)
                rows.append((sid, str(out_bam), "ADMIX"))
            else:
                (d1_id, d1_bam), (d2_id, d2_bam) = donors
                for split_name, (w1, w2) in donor_splits.items():
                    bg = 1.0 - t
                    d1f = bg * w1 / (w1 + w2)
                    d2f = bg * w2 / (w1 + w2)
                    sid = sample_name3(
                        host_id, d1_id, d2_id, pct, d1f * 100, d2f * 100, rep
                    )
                    out_bam = bam_dir / (
                        f"{host_id}-{d1_id}-{d2_id}_h{pct:g}_{split_name}_rep{rep}.bam"
                    )
                    components = [(d1_bam, d1f), (d2_bam, d2f), (host_bam, t)]
                    _emit(
                        mix_cmd(mix_script, out_bam, panel_bed, sid, seed, components),
                        dry_run,
                    )
                    rows.append((sid, str(out_bam), "ADMIX"))

    _write_csv(out_csv_dir / f"{patient}.synthetic.csv", rows, n_ref, dry_run)
    return len(rows) - n_ref
    return len(rows) - 2


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv-dir", type=Path, default=DEFAULT_CSV_DIR,
                    help="Directory of real per-pair sample CSVs to source HOST/DONOR BAMs from.")
    ap.add_argument("--out-csv-dir", type=Path, default=DEFAULT_OUT_CSV_DIR,
                    help="Where to write the synthetic per-pair CSVs.")
    ap.add_argument("--bam-dir", type=Path, default=DEFAULT_BAM_DIR,
                    help="Where to write the mixed BAMs (TAU-side, large).")
    ap.add_argument("--mix-script", type=Path, default=DEFAULT_MIX_SCRIPT,
                    help="Path to scripts/mix_bams.sh.")
    ap.add_argument("--panel-bed", type=Path, default=DEFAULT_PANEL_BED,
                    help="Capture-panel BED used to depth-normalize the subsampling "
                         "(on-target reads). Must match the pipeline `intervals` BED.")
    ap.add_argument("--fractions", type=float, nargs="+", default=DEFAULT_FRACTIONS_PCT,
                    metavar="PCT", help="Host (minor) fractions to synthesise, as percentages.")
    ap.add_argument("--reps", type=int, default=DEFAULT_REPS,
                    help="Independent subsample seeds per (mixture, fraction).")
    ap.add_argument("--donor-splits", nargs="+", default=list(DEFAULT_DONOR_SPLITS),
                    metavar="NAME", choices=list(DEFAULT_DONOR_SPLITS),
                    help="Donor1:donor2 background splits for the three-person (host + "
                         f"2 donor) mixtures. Choices: {', '.join(DEFAULT_DONOR_SPLITS)}.")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print the mix_bams commands and follow-on steps; touch nothing.")
    args = ap.parse_args()

    donor_splits = {k: DEFAULT_DONOR_SPLITS[k] for k in args.donor_splits}

    csv_paths = [p for p in sorted(args.csv_dir.glob("*.csv"))
                 if not p.name.endswith(".synthetic.csv")]
    if not csv_paths:
        sys.exit(f"No sample CSVs found in {args.csv_dir}")

    n_two = 0
    n_three = 0
    n_skipped = 0
    total_rows = 0
    for path in csv_paths:
        patient = path.stem
        mix = read_mix_csv(path)
        if mix is None:
            print(f"Skipping {patient} (no host, or not 1-2 donors)")
            n_skipped += 1
            continue
        host, donors = mix
        if len(donors) == 1:
            print(f"Pair {patient}: host(minor)={host[0]} donor(major)={donors[0][0]}")
            n_two += 1
        else:
            print(f"Trio {patient}: host(minor)={host[0]} "
                  f"donors={donors[0][0]},{donors[1][0]} splits={list(donor_splits)}")
            n_three += 1
        total_rows += build_mix(
            patient, host, donors, args.fractions, args.reps, donor_splits,
            args.bam_dir, args.out_csv_dir, args.mix_script, args.panel_bed,
            args.dry_run,
        )

    print(
        f"\n{n_two} pairs + {n_three} trios processed, {n_skipped} skipped, "
        f"{total_rows} synthetic mixtures."
    )

    csv_glob = f"{args.out_csv_dir}/*.synthetic.csv"
    print("\nFollow-on steps (run TAU-side, where the BAMs and GATK are):")
    print("\n  # 1. Joint-call + pileup the synthetic mixtures with the normal pipeline")
    print("  snakemake -s pipeline/Snakefile \\")
    print("      --configfile paper/public_data/SRP434573/config.yaml \\")
    print(f"      --config samples_csv_dir={args.out_csv_dir} "
          f"output_dir={PIPELINE_OUTPUT_DIR} \\")
    print("      --cores 16")
    print("\n  # 2. Copy the per-pair genotype + admix VCFs into the committed snapshot")
    print(f"  mkdir -p {SNAPSHOT_DIR}")
    print(f"  cp {PIPELINE_OUTPUT_DIR}/*.synthetic.*.vcf.gz* {SNAPSHOT_DIR}/")
    print(f"\n  (then rebuild the paper; synthetic CSV glob: {csv_glob})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
