"""Tests for paper/scripts/facts_io.py (typed facts CSVs for vibepaper)."""

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest
from vibepaper.facts import FactTypeError
from vibepaper.render import load_facts

# Make the paper/scripts module importable (this file lives in paper/tests/).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

import facts_io  # noqa: E402


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
def test_fact_type_for(name: str, expected: str | None) -> None:
    cls = facts_io.fact_type_for(name)
    assert (cls.fact_type if cls else None) == expected


def test_round_trip_through_vibepaper(tmp_path: Path) -> None:
    facts = {
        "mae_pct": np.float64(0.42),
        "n_samples": np.int64(56),
        "n_reps": 5.0,
        "detect_rate_0p2": "0.66",
        "lowest_ratio": 2.8,
        "label": "F2 into F1",
        "depth": 1000,
    }
    path = facts_io.write_fact_csv(tmp_path / "demo.csv", facts)
    assert path.read_text().splitlines()[0] == "field,value,type"
    assert b"\r\n" not in path.read_bytes()

    loaded = load_facts(tmp_path)["demo"]
    assert loaded["mae_pct"] == pytest.approx(0.42)
    assert loaded["mae_pct"].fact_type == "percent"
    assert loaded["n_samples"] == 56 and loaded["n_samples"].fact_type == "count"
    assert str(loaded["n_reps"]) == "5"
    assert loaded["detect_rate_0p2"].fact_type == "fraction"
    assert loaded["lowest_ratio"].fact_type == "ratio"
    assert loaded["label"] == "F2 into F1"
    assert str(loaded["depth"]) == "1000"

    prov = json.loads((tmp_path / "demo.provenance.json").read_text())
    assert set(prov["tools"]) == {"allomix"}


def test_non_finite_and_text_values_stay_untyped(tmp_path: Path) -> None:
    facts = {"lod_pct": math.nan, "inf_pct": float("inf"), "ceiling_pct": ">20"}
    facts_io.write_fact_csv(tmp_path / "grid.csv", facts, provenance=False)
    rows = (tmp_path / "grid.csv").read_text().splitlines()[1:]
    assert rows == ["lod_pct,nan,", "inf_pct,inf,", "ceiling_pct,>20,"]
    assert not (tmp_path / "grid.provenance.json").exists()


def test_out_of_range_typed_value_raises(tmp_path: Path) -> None:
    with pytest.raises(FactTypeError):
        facts_io.write_fact_csv(tmp_path / "bad.csv", {"detect_rate_0p1": 1.5})
    with pytest.raises(FactTypeError):
        facts_io.write_fact_csv(tmp_path / "bad.csv", {"n_markers": 3.5})


def test_read_fact_csv_both_layouts(tmp_path: Path) -> None:
    horizontal = tmp_path / "h.csv"
    horizontal.write_text("a_pct,n_x,label\r\n0.5,3,hello\r\n")
    vertical = facts_io.write_fact_csv(
        tmp_path / "v.csv", {"a_pct": 0.5, "n_x": 3, "label": "hello"}, provenance=False
    )
    assert facts_io.read_fact_csv(horizontal) == {"a_pct": "0.5", "n_x": "3", "label": "hello"}
    assert facts_io.read_fact_csv(vertical) == {"a_pct": "0.5", "n_x": "3", "label": "hello"}
