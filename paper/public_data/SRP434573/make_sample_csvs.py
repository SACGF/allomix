#!/usr/bin/env python3
"""Generate allomix joint-calling CSVs + a ground-truth manifest for SRP434573.

SRP434573 / PRJNA960854 is the public deposit behind Chu Xufeng's HUST PhD
thesis (2024, DOI 10.27157/d.cnki.ghzku.2024.000769): a 1062-autosomal-SNP
MIP capture panel run on two-person (and one three-person) DNA mixtures of
seven unrelated individuals (F1-F3, M1-M4) at known ratios. See README.md and
SACGF/allomix issue #16 for the full provenance.

This script reads the ENA run->alias table (``ena_runs.tsv``) and emits one
per-patient CSV per contributor mixture, in the format the joint-calling
pipeline expects (``sample_id,bam_filename,sample_type``).

Mixture aliases are ``1_<N>_<X>-<Y>`` = a major:minor = 1:N mixture, so the
minor contributor's fraction is ``1/(1+N)``. The thesis confirms the 1:N
convention; which of ``X-Y`` is the minor is NOT stated in the thesis text, so
we take the FIRST-listed individual as the minor. The data structure supports
this (each first-position individual is the one titrated across several ratios
and backgrounds), but it is an inference, not author-confirmed.

Role mapping: the minor (titrated) contributor is labelled HOST and the major
(background) contributor DONOR. This mirrors our common clinical case, where the
residual / recurring patient (host) is the small fraction detected against a
donor-dominated graft, so the titration series exercises allomix exactly as a
relapse / declining-chimerism series would. The labels do not affect the
joint-calling genotypes (HOST and DONOR are genotyped identically); they only
set which contributor allomix treats as host vs donor downstream. Flip
``MINOR_IS_FIRST`` if the authors confirm the minor is the second-listed name.

The ENA files are pear-merged single-end FASTQs. The pipeline consumes aligned
BAMs, so these CSVs reference ``<bam_dir>/<run>.<bam_suffix>`` which an upstream
alignment step must produce (see README.md).
"""

from __future__ import annotations

import argparse
import csv
import os
import re
from collections import defaultdict

MINOR_IS_FIRST = True  # alias "1_N_X-Y": X is the minor (titrated) contributor

PURE_RE = re.compile(r"^([FM][1-9])$")
SINGLE_RE = re.compile(r"^single\d+$")
TWO_PERSON_RE = re.compile(r"^1_(\d+)_([FM][1-9])-([FM][1-9])(?:_v(\d))?$")
DEGRADED_RE = re.compile(r"^1_(\d+)_([FM][1-9])-([FM][1-9])-degraded$")
THREE_PERSON_RE = re.compile(r"^1_(\d+)_(\d+)_([FM][1-9])-([FM][1-9])-([FM][1-9])$")


def minor_fraction(n: int) -> float:
    """Minor-contributor fraction for a major:minor = 1:N mixture."""
    return 1.0 / (1.0 + n)


def parse_runs(tsv_path: str) -> list[dict[str, str]]:
    with open(tsv_path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def bam_path(bam_dir: str, run: str, suffix: str) -> str:
    return os.path.join(bam_dir, f"{run}.{suffix}")


def main() -> None:
    here = os.path.dirname(os.path.abspath(__file__))
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--ena-tsv", default=os.path.join(here, "ena_runs.tsv"))
    ap.add_argument("--out-dir", default=os.path.join(here, "sample_csvs"))
    ap.add_argument("--manifest", default=os.path.join(here, "manifest.tsv"))
    ap.add_argument(
        "--bam-dir",
        default="/tau/data/chimerism/SRP434573/bam",
        help="Directory holding the aligned BAMs referenced by the CSVs.",
    )
    ap.add_argument(
        "--bam-suffix",
        default="hg38.bam",
        help="Suffix appended to each run accession to form the BAM filename.",
    )
    args = ap.parse_args()

    rows = parse_runs(args.ena_tsv)

    # run accession for each pure single-source individual (host/donor source)
    pure_run: dict[str, str] = {}
    # mixtures grouped by patient; each entry: (sample_id, run, minor_pct, kind)
    two_person: dict[tuple[str, str], list[tuple]] = defaultdict(list)
    three_person: list[tuple] = []
    n_singles = 0

    for r in rows:
        run = r["run_accession"]
        name = r["sample_alias"].removesuffix(".bam")

        if m := PURE_RE.match(name):
            pure_run[m.group(1)] = run
        elif SINGLE_RE.match(name):
            n_singles += 1
        elif m := THREE_PERSON_RE.match(name):
            a, b = int(m.group(1)), int(m.group(2))
            c1, c2, c3 = m.group(3), m.group(4), m.group(5)
            three_person.append((name, run, (a, b), (c1, c2, c3)))
        elif (m := TWO_PERSON_RE.match(name)) or (m := DEGRADED_RE.match(name)):
            n = int(m.group(1))
            first, second = m.group(2), m.group(3)
            minor, major = (first, second) if MINOR_IS_FIRST else (second, first)
            key = (minor, major)
            two_person[key].append((name, run, minor_fraction(n) * 100.0))
        else:
            raise SystemExit(f"Unrecognised sample alias: {name!r}")

    os.makedirs(args.out_dir, exist_ok=True)
    manifest = [
        ("patient", "sample_id", "run_accession", "role", "individual", "known_minor_pct")
    ]

    def ref_rows(host: str, donors: list[str]) -> list[tuple[str, str, str]]:
        out = [(host, bam_path(args.bam_dir, pure_run[host], args.bam_suffix), "HOST")]
        for d in donors:
            out.append((d, bam_path(args.bam_dir, pure_run[d], args.bam_suffix), "DONOR"))
        return out

    def write_csv(patient: str, rows_out: list[tuple[str, str, str]]) -> None:
        path = os.path.join(args.out_dir, f"{patient}.csv")
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["sample_id", "bam_filename", "sample_type"])
            w.writerows(rows_out)

    n_patients = 0
    for (minor, major), mixes in sorted(two_person.items()):
        # Stable composition label "mix_<minor>_into_<major>" (spike-in into
        # background), independent of the HOST/DONOR role mapping below.
        patient = f"mix_{minor}_into_{major}"
        # Clinical mapping: minor (titrated) = HOST, major (background) = DONOR.
        rows_out = ref_rows(minor, [major])
        manifest.append((patient, minor, pure_run[minor], "HOST", minor, ""))
        manifest.append((patient, major, pure_run[major], "DONOR", major, ""))
        for sid, run, pct in sorted(mixes, key=lambda x: x[2]):
            rows_out.append((sid, bam_path(args.bam_dir, run, args.bam_suffix), "ADMIX"))
            manifest.append((patient, sid, run, "ADMIX", minor, f"{pct:g}"))
        write_csv(patient, rows_out)
        n_patients += 1

    # Three-person mixture: alias 1_3_5_F2-M1-M2 = ratio 1:3:5 of F2:M1:M2.
    # Clinical "host + up to 2 donors = 3 genomes" case: the smallest fraction
    # (the monitored minority = patient) = HOST, the two larger = DONORs. The
    # patient label keeps the stable composition form (smaller contributors
    # "into" the largest background) so it does not change with the role mapping.
    for name, run, (a, b), (c1, c2, c3) in three_person:
        parts = [(c1, 1), (c2, a), (c3, b)]  # tokens "1_3_5" -> 1:3:5 across c1,c2,c3
        total = sum(p for _, p in parts)
        parts_sorted = sorted(parts, key=lambda x: x[1])  # ascending fraction
        label = "_".join(p[0] for p in parts_sorted[:-1])
        patient = f"mix_{label}_into_{parts_sorted[-1][0]}"
        host = parts_sorted[0][0]  # smallest fraction = monitored host
        donors = [p[0] for p in parts_sorted[1:]]
        rows_out = ref_rows(host, donors)
        host_frac = next(p for ind, p in parts if ind == host) / total * 100.0
        manifest.append((patient, host, pure_run[host], "HOST", host, f"{host_frac:g}"))
        for d in donors:
            manifest.append((patient, d, pure_run[d], "DONOR", d, ""))
        rows_out.append((name, bam_path(args.bam_dir, run, args.bam_suffix), "ADMIX"))
        breakdown = ", ".join(f"{ind}={p / total * 100:g}%" for ind, p in parts_sorted)
        manifest.append((patient, name, run, "ADMIX", host, breakdown))
        write_csv(patient, rows_out)
        n_patients += 1

    with open(args.manifest, "w", newline="") as f:
        csv.writer(f, delimiter="\t").writerows(manifest)

    print(f"Pure single-source references: {len(pure_run)} ({', '.join(sorted(pure_run))})")
    print(f"Unused single*.bam runs (not in any mixture): {n_singles}")
    print(f"Two-person mixture patients: {len(two_person)}")
    print(f"Three-person mixture patients: {len(three_person)}")
    print(f"Wrote {n_patients} CSVs to {args.out_dir}")
    print(f"Wrote manifest to {args.manifest}")


if __name__ == "__main__":
    main()
