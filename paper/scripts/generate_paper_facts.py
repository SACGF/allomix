#!/usr/bin/env python3
"""Generate the paper's fixed facts CSVs: simulation calibration, panel specs and
the empirical panel stats.

The bias-correction facts (``bias_no_bias``, ``bias_with_bias``) are written by
``compare_bias_correction.py``, which runs both arms on the same simulated data.

Usage:
    python paper/scripts/generate_paper_facts.py
"""

import csv
from pathlib import Path

from facts_io import write_fact_csv

FACTS_DIR = Path("output/facts")


def write_fact(name: str, data: dict) -> None:
    """Write a typed facts CSV (see facts_io)."""
    path = write_fact_csv(FACTS_DIR / f"{name}.csv", data)
    print(f"  Wrote {path}")


def main():
    FACTS_DIR.mkdir(parents=True, exist_ok=True)

    # Simulation calibration inputs shared across the validation scripts
    # (run_depth_validation.py, run_relatedness_validation.py, etc.), rounded
    # from empirical panel_empirical
    # values (depth CV 0.429, no-call 1.6%). Exposing them as a fact keeps the
    # paper prose in sync with what the simulator actually ran.
    write_fact("sim_calibration", {
        "depth_cv": 0.43,
        "locus_dropout_pct": 1.6,
        "seq_error_pct": 1,
    })

    # The commercial-tool figures (tool_landscape.*) are literature values, not
    # measured here: they live in paper/frozen_facts.csv with their citations.

    # Panel specs
    write_fact("panel_specs", {
        "n_markers_panel": 76,
        "marker_type": "SNP",
        "panel_name": "IDT rhAmpSeq Sample ID",
        "typical_depth": ">1000x",
        "min_informative_markers": 3,
    })

    # Empirical panel stats (measured on the internal data by
    # scripts/measure_panel_bias.py, committed as a one-row CSV) rewritten typed.
    empirical_src = Path("paper/empirical_results/panel_empirical.csv")
    # Passed as raw strings, not through vibepaper's read_facts: for a one-row
    # (horizontal) CSV that upcasts every integer to float (1732 -> 1732.0).
    with open(empirical_src, newline="", encoding="utf-8") as fh:
        (empirical,) = csv.DictReader(fh)
    write_fact_csv(FACTS_DIR / "panel_empirical.csv", empirical)
    print(f"  Wrote {empirical_src} -> {FACTS_DIR / 'panel_empirical.csv'} (typed)")

    print("\nDone! Facts written to output/facts/")


if __name__ == "__main__":
    main()
