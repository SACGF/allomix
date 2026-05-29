#!/usr/bin/env python3
"""Batch allomix monitor over per-patient CSVs from the joint-calling pipeline.

For each CSV in ``samples_csv_dir`` (default ``pipeline/sample_csvs``), this
script:

  1. Reads HOST / DONOR / ADMIX sample IDs from the CSV.
  2. Locates the matching panel VCF at ``<vcf_dir>/<patient>.vcf.gz`` and
     admix VCF at ``<vcf_dir>/<patient>.admix.vcf.gz``.
  3. Runs ``allomix monitor`` once per patient with all admix timepoints,
     writing ``<output_dir>/<patient>.tsv``.
  4. Concatenates per-patient outputs into ``<output_dir>/batch.tsv``.

Patients with no ADMIX rows are skipped (nothing to estimate).

Usage:
    python scripts/run_csv_batch.py \\
        --samples-csv-dir pipeline/sample_csvs \\
        --vcf-dir output/joint_call \\
        --output-dir output/batch \\
        [--bias-table bias.tsv] \\
        [--error-table errors.tsv] \\
        [--extra-arg --min-dp=200 --extra-arg --min-gq=30]
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import subprocess
import sys


def _load_patient_csv(path: str) -> tuple[list[str], list[str], list[str]]:
    """Return (hosts, donors, admix) lists of sample IDs from one CSV."""
    hosts: list[str] = []
    donors: list[str] = []
    admix: list[str] = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        if "sample_type" not in (reader.fieldnames or []):
            sys.exit(f"{path}: missing required 'sample_type' column")
        for row in reader:
            stype = row["sample_type"].strip().upper()
            sid = row["sample_id"].strip()
            if stype == "HOST":
                hosts.append(sid)
            elif stype == "DONOR":
                donors.append(sid)
            elif stype == "ADMIX":
                admix.append(sid)
            else:
                sys.exit(f"{path}: unknown sample_type {stype!r} for {sid}")
    return hosts, donors, admix


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run allomix monitor for every patient CSV in samples-csv-dir, "
            "using the panel + admix VCFs produced by pipeline/Snakefile."
        )
    )
    parser.add_argument(
        "--samples-csv-dir",
        default="pipeline/sample_csvs",
        help="Directory of per-patient CSVs (default: pipeline/sample_csvs)",
    )
    parser.add_argument(
        "--vcf-dir",
        default="output/joint_call",
        help=(
            "Directory containing <patient>.vcf.gz and <patient>.admix.vcf.gz "
            "(default: output/joint_call)"
        ),
    )
    parser.add_argument(
        "--output-dir",
        default="output/batch",
        help="Directory for per-patient TSVs and batch.tsv (default: output/batch)",
    )
    parser.add_argument(
        "--panel-suffix",
        default=".vcf.gz",
        help=(
            "Filename suffix for the panel VCF: "
            "<vcf-dir>/<patient><panel-suffix>. Default '.vcf.gz'. Use "
            "e.g. '.union_sid_haem_vendor_probes.vcf.gz' when the pipeline "
            "was run with an --intervals BED."
        ),
    )
    parser.add_argument(
        "--admix-suffix",
        default=".admix.vcf.gz",
        help="Filename suffix for the admix VCF (default '.admix.vcf.gz')",
    )
    parser.add_argument(
        "--bias-table",
        default=None,
        help="Per-marker bias table TSV (passed through to allomix monitor)",
    )
    parser.add_argument(
        "--error-table",
        default=None,
        help="Per-site error table TSV (passed through to allomix monitor)",
    )
    parser.add_argument(
        "--allomix",
        default="allomix",
        help="Path to the allomix executable (default: allomix on $PATH)",
    )
    parser.add_argument(
        "--extra-arg",
        action="append",
        default=[],
        metavar="ARG",
        help=(
            "Extra argument forwarded verbatim to each `allomix monitor` "
            "invocation. Repeat for multiple, e.g. "
            "--extra-arg --min-dp=200 --extra-arg --min-gq=30"
        ),
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    csv_paths = sorted(glob.glob(os.path.join(args.samples_csv_dir, "*.csv")))
    if not csv_paths:
        sys.exit(f"No CSVs found in {args.samples_csv_dir}")

    n_run = 0
    n_skip_no_admix = 0
    n_fail = 0
    completed: list[str] = []

    for csv_path in csv_paths:
        patient = os.path.splitext(os.path.basename(csv_path))[0]
        hosts, donors, admix = _load_patient_csv(csv_path)

        if not admix:
            print(f"[{patient}] no ADMIX rows — skipping")
            n_skip_no_admix += 1
            continue
        if len(hosts) != 1:
            print(
                f"[{patient}] expected exactly 1 HOST row, got {len(hosts)} — skipping",
                file=sys.stderr,
            )
            n_fail += 1
            continue
        if not donors:
            print(f"[{patient}] no DONOR rows — skipping", file=sys.stderr)
            n_fail += 1
            continue

        panel_vcf = os.path.join(args.vcf_dir, f"{patient}{args.panel_suffix}")
        admix_vcf = os.path.join(args.vcf_dir, f"{patient}{args.admix_suffix}")
        for v in (panel_vcf, admix_vcf):
            if not os.path.exists(v):
                print(f"[{patient}] missing VCF: {v}", file=sys.stderr)
                n_fail += 1
                break
        else:
            output_tsv = os.path.join(args.output_dir, f"{patient}.tsv")
            cmd = [
                args.allomix, "monitor",
                "--panel-vcf", panel_vcf,
                "--admix-vcf", admix_vcf,
                "--host-sample", hosts[0],
                "--output", output_tsv,
            ]
            for d in donors:
                cmd += ["--donor-sample", d]
            for s in admix:
                cmd += ["--sample", s]
            if args.bias_table:
                cmd += ["--bias-table", args.bias_table]
            if args.error_table:
                cmd += ["--error-table", args.error_table]
            cmd += args.extra_arg

            print(f"[{patient}] {' '.join(cmd)}")
            result = subprocess.run(cmd, check=False)
            if result.returncode != 0:
                print(
                    f"[{patient}] allomix monitor failed (exit {result.returncode})",
                    file=sys.stderr,
                )
                n_fail += 1
            else:
                completed.append(output_tsv)
                n_run += 1

    if completed:
        batch_tsv = os.path.join(args.output_dir, "batch.tsv")
        header_line: str | None = None
        with open(batch_tsv, "w", encoding="utf-8") as out:
            for path in completed:
                with open(path, encoding="utf-8") as f:
                    lines = f.readlines()
                if not lines:
                    continue
                # allomix monitor writes a header per --sample, so a per-patient
                # TSV with N timepoints has N header lines interleaved with data.
                # Keep the first header globally and drop every subsequent one.
                if header_line is None:
                    header_line = lines[0]
                    out.write(header_line)
                for line in lines:
                    if line == header_line:
                        continue
                    out.write(line)
        print(f"\nCombined output: {batch_tsv}")

    print(
        f"\nDone: {n_run} run, {n_skip_no_admix} skipped (no admix), "
        f"{n_fail} failed."
    )
    if n_fail:
        sys.exit(1)


if __name__ == "__main__":
    main()
