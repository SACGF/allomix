"""Shared utilities for allomix scripts."""

from __future__ import annotations

import csv
from pathlib import Path


def write_truth_table(
    rows: list[dict],
    path: Path,
    fieldnames: list[str] | None = None,
) -> None:
    """Write a truth table as a tab-separated file.

    Args:
        rows: List of dicts, one per sample.
        path: Output TSV path.
        fieldnames: Column order. Defaults to keys of the first row.
    """
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
