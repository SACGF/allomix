"""Tests for the quick-build draft render's placeholder filling
(paper/scripts/build_quick_paper.py)."""

import sys
from pathlib import Path

# Make the paper/scripts module importable (this file lives in paper/tests/).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

from build_quick_paper import fill_blank_facts  # noqa: E402
from facts_io import FACTS  # noqa: E402
from vibepaper.facts import Count, Percent, Ratio, read_facts  # noqa: E402


def test_blank_facts_get_typed_placeholders(tmp_path):
    FACTS.write(
        tmp_path / "lod.csv",
        {
            "mle_lod_pct": None,
            "n_pairs": None,
            "gain_fold": None,
            "label": None,
            "kept_pct": 1.5,
        },
        tools=[],
    )
    # A horizontal data table is left alone, blanks included.
    (tmp_path / "raw.csv").write_text("a,b\n1,\n2,3\n")

    filled = fill_blank_facts(tmp_path)

    assert len(filled) == 4
    facts = read_facts(tmp_path / "lod.csv")
    assert isinstance(facts["mle_lod_pct"], Percent) and facts["mle_lod_pct"] == 0
    assert isinstance(facts["n_pairs"], Count) and facts["n_pairs"] == 0
    assert isinstance(facts["gain_fold"], Ratio) and facts["gain_fold"] == 1
    assert facts["label"] == 0
    assert facts["kept_pct"] == 1.5
    assert (tmp_path / "raw.csv").read_text() == "a,b\n1,\n2,3\n"
