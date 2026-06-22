#!/usr/bin/env python3
"""Generate semi-synthetic sub-0.5% mixtures for the SRP434573 paper section.

The public SRP434573 titration series bottoms out at a 0.5% minor (host)
fraction. This driver builds *lower* points (host fractions 0.1-0.5% by default)
by blending each two-person mixture's two pure reference BAMs with
``samtools view --subsample`` (wrapped by ``scripts/mix_bams.sh``). The resulting
reads are real (real panel noise, real GATK/bcftools path) but the mixing ratio
is artificial, so the points are *semi-synthetic* and must be labelled as such
downstream.

This step runs ONCE, TAU-side, because the source BAMs live on ``/tau``. It does
two things:

  1. For every two-person pair, subsample+merge the pure HOST and DONOR BAMs at
     each (fraction, seed) into a mixed BAM under ``--bam-dir``.
  2. Write one synthetic per-pair CSV (``<pair>.synthetic.csv``) whose HOST/DONOR
     rows reuse the pure reference BAMs (so phase-1 genotypes match the real run)
     and whose ADMIX rows point at the mixed BAMs.

It then prints the two follow-on commands to run by hand (it does NOT launch GATK
itself): the existing ``pipeline/Snakefile`` over the synthetic CSVs, then a copy
of the per-pair genotype + admix VCFs into the committed snapshot at
``paper/public_data/SRP434573/genotypes_synthetic``. The paper build consumes
that snapshot (see ``paper/scripts/run_srp434573_allomix.py``).

Convention trap: ``mix_bams.sh HOST_BAM DONOR_BAM DONOR_FRACTION`` treats its
DONOR_BAM as the *minor* (titrated) contributor. In SRP434573/allomix the *host*
is the minor monitored fraction. So the allomix-HOST individual is passed as
mix_bams' DONOR_BAM, and the allomix-DONOR (background) as mix_bams' HOST_BAM.

Usage (TAU-side, where the BAMs are):

    python scripts/make_semisynthetic_srp434573.py \
        --bam-dir /tau/data/chimerism/SRP434573/synthetic_bam

Inspect without touching any BAM:

    python scripts/make_semisynthetic_srp434573.py --dry-run
"""

import argparse
import csv
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
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


def read_pair_csv(path: Path) -> tuple[tuple[str, str], tuple[str, str]] | None:
    """Return ((host_id, host_bam), (donor_id, donor_bam)) for a two-person CSV.

    The host is the minor (titrated) contributor, the donor the majority
    background. Returns ``None`` for CSVs that are not a clean one-host/one-donor
    pair (e.g. the three-person mixture), which are skipped.
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
    if host is None or len(donors) != 1:
        return None
    return host, donors[0]


def sample_name(minor: str, major: str, pct: float, rep: int) -> str:
    """Synthetic admix sample id, e.g. ``syn_F2-M1_f0.1_rep3``.

    The known host fraction and replicate read straight off the name, mirroring
    how the real ``1_N_X-Y`` aliases encode their fraction.
    """
    return f"syn_{minor}-{major}_f{pct:g}_rep{rep}"


def mixed_bam_path(bam_dir: Path, minor: str, major: str, pct: float, rep: int) -> Path:
    return bam_dir / f"{minor}-{major}_f{pct:g}_rep{rep}.bam"


def build_pair(
    patient: str,
    host: tuple[str, str],
    donor: tuple[str, str],
    fractions_pct: list[float],
    reps: int,
    bam_dir: Path,
    out_csv_dir: Path,
    mix_script: Path,
    panel_bed: Path,
    dry_run: bool,
) -> int:
    """Mix one pair across the grid and write its synthetic CSV.

    Returns the number of ADMIX (synthetic mixture) rows produced.
    """
    host_id, host_bam = host          # allomix host = minor (titrated)
    donor_id, donor_bam = donor       # allomix donor = majority background

    rows: list[tuple[str, str, str]] = [
        (host_id, host_bam, "HOST"),
        (donor_id, donor_bam, "DONOR"),
    ]
    for pct in fractions_pct:
        frac = pct / 100.0
        for rep in range(1, reps + 1):
            sid = sample_name(host_id, donor_id, pct, rep)
            out_bam = mixed_bam_path(bam_dir, host_id, donor_id, pct, rep)
            # Spaced so host (seed) and donor (seed+1) never collide across reps.
            seed = 100 * rep
            cmd = [
                str(mix_script),
                donor_bam,   # mix_bams HOST_BAM = majority background
                host_bam,    # mix_bams DONOR_BAM = minor (titrated) = allomix host
                f"{frac:g}",
                str(out_bam),
                str(panel_bed),  # on-target depth normalization (see mix_bams.sh)
                sid,
                str(seed),
            ]
            if dry_run:
                print("  mix: " + " ".join(cmd))
            else:
                subprocess.run(cmd, check=True)
            rows.append((sid, str(out_bam), "ADMIX"))

    csv_path = out_csv_dir / f"{patient}.synthetic.csv"
    if dry_run:
        print(f"  would write {csv_path} ({len(rows)} rows: 1 HOST, 1 DONOR, "
              f"{len(rows) - 2} ADMIX)")
    else:
        out_csv_dir.mkdir(parents=True, exist_ok=True)
        with open(csv_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["sample_id", "bam_filename", "sample_type"])
            w.writerows(rows)
        print(f"Wrote {csv_path} ({len(rows) - 2} synthetic ADMIX rows)")
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
                    help="Independent subsample seeds per (pair, fraction).")
    ap.add_argument("--dry-run", action="store_true",
                    help="Print the mix_bams commands and follow-on steps; touch nothing.")
    args = ap.parse_args()

    csv_paths = [p for p in sorted(args.csv_dir.glob("*.csv"))
                 if not p.name.endswith(".synthetic.csv")]
    if not csv_paths:
        sys.exit(f"No sample CSVs found in {args.csv_dir}")

    n_pairs = 0
    n_skipped = 0
    total_rows = 0
    for path in csv_paths:
        patient = path.stem
        pair = read_pair_csv(path)
        if pair is None:
            print(f"Skipping {patient} (not a one-host/one-donor pair)")
            n_skipped += 1
            continue
        host, donor = pair
        print(f"Pair {patient}: host(minor)={host[0]} donor(major)={donor[0]}")
        total_rows += build_pair(
            patient, host, donor, args.fractions, args.reps,
            args.bam_dir, args.out_csv_dir, args.mix_script, args.panel_bed,
            args.dry_run,
        )
        n_pairs += 1

    print(
        f"\n{n_pairs} pairs processed, {n_skipped} skipped, "
        f"{total_rows} synthetic mixtures "
        f"({len(args.fractions)} fractions x {args.reps} reps x {n_pairs} pairs)."
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
