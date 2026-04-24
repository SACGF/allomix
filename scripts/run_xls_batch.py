#!/usr/bin/env python3
"""Batch validation: run allomix monitor for each row in an XLS/XLSX file.

Usage:
    python scripts/run_xls_batch.py samples.xlsx \
        --vcf joint_called.vcf.gz \
        --host-column Host \
        --donor-column Donor \
        --test-sample-column Sample \
        --bias-table-tsv bias.tsv \
        --output-dir output/batch
"""

import argparse
import os
import subprocess
import sys

try:
    import openpyxl
except ImportError:
    sys.exit(
        "openpyxl is required for this script.\n"
        "Install with: pip install 'allomix[xls]'  (or: pip install openpyxl)"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run allomix monitor for each row in an XLS file."
    )
    parser.add_argument("xls_file", help="Path to the XLS/XLSX file")
    parser.add_argument("--vcf", required=True, help="Path to the joint-called VCF")
    parser.add_argument("--host-column", required=True, help="Column name for host sample")
    parser.add_argument("--donor-column", required=True, help="Column name for donor sample")
    parser.add_argument(
        "--test-sample-column", required=True, help="Column name for test (admixture) sample"
    )
    parser.add_argument("--bias-table-tsv", default=None, help="Path to bias table TSV")
    parser.add_argument(
        "--output-dir", default="output/batch", help="Directory for output files (default: output/batch)"
    )
    parser.add_argument(
        "--copy-columns",
        default=None,
        help="Comma-separated XLS column names to append to each row in batch.tsv",
    )
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    wb = openpyxl.load_workbook(args.xls_file, read_only=True, data_only=True)
    ws = wb.active

    rows = iter(ws.rows)
    headers = [cell.value for cell in next(rows)]

    def col_idx(name: str) -> int:
        try:
            return headers.index(name)
        except ValueError:
            sys.exit(f"Column '{name}' not found. Available columns: {headers}")

    host_idx = col_idx(args.host_column)
    donor_idx = col_idx(args.donor_column)
    test_idx = col_idx(args.test_sample_column)

    copy_col_names: list[str] = (
        [c.strip() for c in args.copy_columns.split(",")] if args.copy_columns else []
    )
    copy_col_idxs = [col_idx(name) for name in copy_col_names]

    n_run = 0
    n_skip = 0
    completed_files: list[tuple[str, list[str]]] = []  # (path, extra_values)

    for row in rows:
        values = [cell.value for cell in row]
        host = values[host_idx]
        donor = values[donor_idx]
        test_sample = values[test_idx]

        row_vals = {"host": host, "donor": donor, "test": test_sample}
        if any(str(v).strip().upper() == "N/A" or v is None for v in row_vals.values()):
            print(f"Skipping: {row_vals}")
            n_skip += 1
            continue

        extra_values = [str(values[i]) if values[i] is not None else "" for i in copy_col_idxs]

        output_file = os.path.join(args.output_dir, f"{test_sample}.tsv")
        cmd = [
            "allomix", "monitor",
            "--vcf", args.vcf,
            "--host-sample", str(host),
            "--donor-sample", str(donor),
            "--sample", str(test_sample),
            "--output", output_file,
        ]
        if args.bias_table_tsv:
            cmd += ["--bias-table", args.bias_table_tsv]

        print(f"Running [{test_sample}]: {' '.join(cmd)}")
        result = subprocess.run(cmd, check=False)
        if result.returncode != 0:
            print(
                f"ERROR: allomix monitor failed for {test_sample} (exit {result.returncode})",
                file=sys.stderr,
            )
        else:
            completed_files.append((output_file, extra_values))
            n_run += 1

    batch_tsv = os.path.join(args.output_dir, "batch.tsv")
    if completed_files:
        with open(batch_tsv, "w", encoding="utf-8") as out:
            for i, (path, extra_values) in enumerate(completed_files):
                with open(path, encoding="utf-8") as f:
                    lines = f.readlines()
                if not lines:
                    continue
                if extra_values:
                    extra_tsv = "\t" + "\t".join(extra_values)
                    if i == 0:
                        header = lines[0].rstrip("\n") + "\t" + "\t".join(copy_col_names) + "\n"
                        out.write(header)
                        for line in lines[1:]:
                            out.write(line.rstrip("\n") + extra_tsv + "\n")
                    else:
                        for line in lines[1:]:
                            out.write(line.rstrip("\n") + extra_tsv + "\n")
                else:
                    if i == 0:
                        out.writelines(lines)
                    else:
                        out.writelines(lines[1:])
        print(f"Combined output written to {batch_tsv}")

    print(f"\nDone: {n_run} run, {n_skip} skipped.")


if __name__ == "__main__":
    main()
