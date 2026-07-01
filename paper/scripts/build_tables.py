#!/usr/bin/env python3
"""Assemble display-ready table CSVs for the paper from the single-row fact CSVs.

The paper used to build several multi-row tables cell by cell, one Jinja fact
reference per cell. That is verbose in the template and easy to get wrong (a
mislabelled row, a swapped column, a forgotten filter). vibepaper can import a
CSV directly with an ``<!-- include-csv: ... -->`` directive, so instead we emit
one CSV per table, shaped exactly like the rendered table, and let vibepaper
lay it out.

Each table cell is pre-formatted here as a string (including the combined
"value ± sd" and "min-max" cells that a single CSV column cannot express through
include-csv's per-column formatting). Numbers are read from the existing fact
CSVs, so the rendered values are identical to the previous per-cell template;
only the mechanism changes.

Run after the depth, relatedness, and multidonor validation scripts have written
their facts:

    python paper/scripts/build_tables.py

Outputs (into output/facts/):
    table_depth.csv         supplementary depth-accuracy table
    table_relatedness.csv   results relatedness table
    table_multidonor.csv    results multi-donor accuracy table
"""

import csv
from pathlib import Path

FACTS_DIR = Path("output/facts")

# en dash for numeric ranges, matching the previous template text.
NDASH = "–"


def read_fact(name: str) -> dict[str, str]:
    """Read a single-row fact CSV as a dict of raw string values."""
    path = FACTS_DIR / f"{name}.csv"
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if len(rows) != 1:
        raise ValueError(f"{path}: expected exactly 1 data row, got {len(rows)}")
    return rows[0]


def dp(value: str, decimals: int) -> str:
    """Fixed-decimal format, matching vibepaper's ``dp`` filter."""
    return f"{float(value):.{decimals}f}"


def write_table(name: str, header: list[str], rows: list[list[str]]) -> None:
    """Write a display-ready table CSV with pre-formatted string cells."""
    path = FACTS_DIR / f"{name}.csv"
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows)
    print(f"Wrote {path} ({len(rows)} rows)")


def build_depth_table() -> None:
    """Supplementary depth-accuracy table (one row per sequencing depth)."""
    depths = [(50, "50x"), (100, "100x"), (200, "200x"), (500, "500x"), (1000, "1,000x")]
    header = [
        "Depth",
        "MAE (%)",
        "RMSE (%)",
        "Max Error (%)",
        "CI Coverage (%)",
        "Mean CI Width (%)",
    ]
    rows = []
    for depth, label in depths:
        d = read_fact(f"depth_{depth}")
        rows.append(
            [
                label,
                f"{dp(d['mean_abs_error_pct'], 2)} ± {dp(d['mean_abs_error_sd_pct'], 2)}",
                f"{dp(d['rmse_pct'], 2)} ± {dp(d['rmse_sd_pct'], 2)}",
                f"{dp(d['max_abs_error_pct'], 2)} ± {dp(d['max_abs_error_sd_pct'], 2)}",
                f"{d['ci_coverage_pct']} ± {d['ci_coverage_sd_pct']}",
                f"{dp(d['mean_ci_width_pct'], 2)} ± {dp(d['mean_ci_width_sd_pct'], 2)}",
            ]
        )
    write_table("table_depth", header, rows)


def build_relatedness_table() -> None:
    """Results relatedness table (one row per relatedness level)."""
    levels = [
        ("unrelated", "Unrelated"),
        ("cousin", "1st cousin"),
        ("half_sibling", "Half-sibling"),
        ("sibling", "Full sibling"),
    ]
    header = ["Relatedness", "Mean Informative", "Range", "MAE (%)", "RMSE (%)"]
    rows = []
    for fact, label in levels:
        r = read_fact(f"rel_{fact}")
        rows.append(
            [
                label,
                r["mean_informative"],
                f"{r['min_informative']}{NDASH}{r['max_informative']}",
                r["mean_mae_pct"],
                r["mean_rmse_pct"],
            ]
        )
    write_table("table_relatedness", header, rows)


def build_multidonor_table() -> None:
    """Results multi-donor accuracy table (metrics as rows, donors as columns)."""
    m = read_fact("multidonor")
    header = ["Metric", "Donor 1", "Donor 2", "Total"]
    rows = [
        ["MAE (%)", m["mae_d1_pct"], m["mae_d2_pct"], m["mae_total_pct"]],
        ["RMSE (%)", m["rmse_d1_pct"], m["rmse_d2_pct"], m["rmse_total_pct"]],
        ["Max error (%)", m["max_error_d1_pct"], m["max_error_d2_pct"], m["max_error_total_pct"]],
        ["CI coverage (%)", m["ci_coverage_d1_pct"], m["ci_coverage_d2_pct"], "n/a"],
    ]
    write_table("table_multidonor", header, rows)


def main() -> None:
    FACTS_DIR.mkdir(parents=True, exist_ok=True)
    build_depth_table()
    build_relatedness_table()
    build_multidonor_table()


if __name__ == "__main__":
    main()
