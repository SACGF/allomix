"""Tests for the paper's facts convention (paper/scripts/facts_io.py).

The writing mechanics (line endings, numpy, range checks, sidecars) belong to
vibepaper's ``FactsWriter`` and are tested there. This checks only the paper's
own rules: field names map to the types the prose expects, and a missing value
is written blank.
"""

import math
import sys
from pathlib import Path

import pytest

# Make the paper/scripts module importable (this file lives in paper/tests/).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from facts_io import FACTS  # noqa: E402


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("mae_reliable_pct", "percent"),
        ("est_100pct", "percent"),
        ("n_timepoints", "count"),
        ("zero_host_presence_falsepos_n", "count"),
        ("detect_rate_0p1", "fraction"),
        ("pooled_detect_0p2", "fraction"),
        ("mean_nocall_rate", "fraction"),
        ("three_host_lowest_ratio", "ratio"),
        ("two_person_offset_reduction_fold", "ratio"),
        ("fold_real_rho_vs_binomial", "ratio"),
        ("depth", None),
        ("mixture_label", None),
        ("markers_used_median", None),
    ],
)
def test_type_for(name: str, expected: str | None) -> None:
    cls = FACTS.type_for(name)
    assert (cls.fact_type if cls else None) == expected


def test_missing_value_is_blank(tmp_path: Path) -> None:
    FACTS.write(tmp_path / "grid.csv", {"lod_pct": math.nan}, tools=None)
    assert (tmp_path / "grid.csv").read_text().splitlines()[1:] == ["lod_pct,,"]
