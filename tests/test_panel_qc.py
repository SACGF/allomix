"""Tests for allomix.qc.panel_qc and the panel-qc CLI subcommand (issue #37)."""

import pytest

from allomix.cli import main
from allomix.qc.panel_qc import (
    PanelQCThresholds,
    assign_verdicts,
    load_marker_stats,
    panel_p95_abs_bias,
    summarize,
    write_sites_bed,
    write_verdict_tsv,
)

# Column order matches scripts/measure_panel_bias.py per-marker TSV output.
HEADER = [
    "marker",
    "marker_index",
    "total_obs",
    "n_called",
    "n_nocall",
    "call_rate",
    "n_hom_ref",
    "n_het",
    "n_hom_alt",
    "het_ratio_vs_hwe",
    "mean_depth",
    "depth_cv",
    "n_het_for_bias",
    "median_bias",
    "mean_bias",
    "sd_within",
]


def _row(marker, call_rate, het_ratio, median_bias, depth_cv, index=0):
    """Build a per-marker TSV row dict with the fields panel-qc reads."""
    return {
        "marker": marker,
        "marker_index": str(index),
        "total_obs": "100",
        "n_called": "100",
        "n_nocall": "0",
        "call_rate": f"{call_rate:.4f}",
        "n_hom_ref": "25",
        "n_het": "50",
        "n_hom_alt": "25",
        "het_ratio_vs_hwe": f"{het_ratio:.4f}",
        "mean_depth": "1500",
        "depth_cv": f"{depth_cv:.4f}",
        "n_het_for_bias": "50",
        "median_bias": f"{median_bias:.6f}",
        "mean_bias": f"{median_bias:.6f}",
        "sd_within": "0.020000",
    }


def _write_tsv(path, rows, header=HEADER):
    lines = ["\t".join(header)]
    for r in rows:
        lines.append("\t".join(r[c] for c in header))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


@pytest.fixture
def sample_rows():
    # A clean marker plus one hit for each primary criterion. The bias marker sits
    # far above the panel p95 |bias|; a low-call marker also has high depth_cv.
    return [
        _row(
            "chr1:100:A:G", call_rate=1.00, het_ratio=1.00, median_bias=0.010, depth_cv=0.3, index=0
        ),
        _row(
            "chr1:200:C:T", call_rate=0.80, het_ratio=0.99, median_bias=0.010, depth_cv=1.5, index=1
        ),
        _row(
            "chr2:300:G:A", call_rate=1.00, het_ratio=0.50, median_bias=0.010, depth_cv=0.3, index=2
        ),
        _row(
            "chr2:400:T:C", call_rate=1.00, het_ratio=1.00, median_bias=0.400, depth_cv=0.3, index=3
        ),
    ]


def test_clean_marker_kept(sample_rows):
    verdicts, _ = assign_verdicts(sample_rows, PanelQCThresholds())
    assert verdicts[0].keep
    assert verdicts[0].reasons == []


def test_low_call_rate_dropped(sample_rows):
    verdicts, _ = assign_verdicts(sample_rows, PanelQCThresholds())
    assert not verdicts[1].keep
    assert "low_call_rate" in verdicts[1].reasons


def test_het_deficit_dropped(sample_rows):
    verdicts, _ = assign_verdicts(sample_rows, PanelQCThresholds())
    assert not verdicts[2].keep
    assert "het_deficit" in verdicts[2].reasons


def test_extreme_bias_uses_absolute_cap(sample_rows):
    # With a tight absolute cap the high-bias marker (0.40) is dropped.
    verdicts, _ = assign_verdicts(sample_rows, PanelQCThresholds(max_abs_bias=0.2))
    assert not verdicts[3].keep
    assert "extreme_bias" in verdicts[3].reasons


def test_depth_cv_is_secondary(sample_rows):
    # Marker 1 (low call rate) also has high depth_cv, so unstable_depth is
    # appended as a secondary reason.
    verdicts, _ = assign_verdicts(sample_rows, PanelQCThresholds())
    assert "unstable_depth" in verdicts[1].reasons
    # A marker with ONLY high depth_cv and no primary hit stays kept.
    clean_high_cv = [
        _row("chr9:1:A:G", call_rate=1.0, het_ratio=1.0, median_bias=0.01, depth_cv=5.0)
    ]
    v2, _ = assign_verdicts(clean_high_cv, PanelQCThresholds())
    assert v2[0].keep
    assert v2[0].reasons == []


def test_missing_stat_does_not_drop():
    # An NA het_ratio must not trigger het_deficit.
    row = _row("chr1:1:A:G", call_rate=1.0, het_ratio=1.0, median_bias=0.01, depth_cv=0.3)
    row["het_ratio_vs_hwe"] = "NA"
    verdicts, _ = assign_verdicts([row], PanelQCThresholds())
    assert verdicts[0].keep


def test_panel_p95_abs_bias(sample_rows):
    p95 = panel_p95_abs_bias(sample_rows)
    assert p95 == pytest.approx(0.400)


def test_summary_counts(sample_rows):
    verdicts, _ = assign_verdicts(sample_rows, PanelQCThresholds())
    summary = summarize(verdicts)
    assert summary["n_markers"] == 4
    assert summary["n_keep"] == 2
    assert summary["n_drop"] == 2
    assert summary["low_call_rate"] == 1
    assert summary["het_deficit"] == 1
    assert summary["unstable_depth"] == 1


def test_verdict_tsv_roundtrip(sample_rows, tmp_path):
    verdicts, _ = assign_verdicts(sample_rows, PanelQCThresholds())
    out = tmp_path / "verdicts.tsv"
    write_verdict_tsv(verdicts, out)
    reloaded = load_marker_stats(out)
    assert reloaded[0]["verdict"] == "keep"
    assert reloaded[1]["verdict"] == "drop"
    assert "low_call_rate" in reloaded[1]["reasons"]


def test_write_exclude_bed(sample_rows, tmp_path):
    verdicts, _ = assign_verdicts(sample_rows, PanelQCThresholds())
    bed = tmp_path / "exclude.bed"
    n = write_sites_bed(verdicts, bed, which="drop")
    assert n == 2
    lines = [ln.split("\t") for ln in bed.read_text().strip().splitlines()]
    # chr1:200 -> 0-based 199 200; chr2:300 -> 299 300
    assert ["chr1", "199", "200"] in lines
    assert ["chr2", "299", "300"] in lines


def test_bed_requires_coordinates(sample_rows, tmp_path):
    # Strip the coordinate column: BED emission must refuse rather than silently
    # drop markers.
    for r in sample_rows:
        r["marker"] = ""
    verdicts, _ = assign_verdicts(sample_rows, PanelQCThresholds())
    with pytest.raises(ValueError, match="coordinate column"):
        write_sites_bed(verdicts, tmp_path / "x.bed", which="drop")


def test_cli_panel_qc(sample_rows, tmp_path):
    tsv = _write_tsv(tmp_path / "per_marker.tsv", sample_rows)
    out = tmp_path / "verdicts.tsv"
    bed = tmp_path / "exclude.bed"
    rc = main(
        [
            "panel-qc",
            str(tsv),
            "--output",
            str(out),
            "--exclude-bed",
            str(bed),
            "--max-abs-bias",
            "0.2",
        ]
    )
    assert rc == 0
    assert out.exists()
    # Clean + bias-only-dropped: with the absolute cap the high-bias marker drops.
    reloaded = load_marker_stats(out)
    verdict_by_marker = {r["marker"]: r["verdict"] for r in reloaded}
    assert verdict_by_marker["chr1:100:A:G"] == "keep"
    assert verdict_by_marker["chr2:400:T:C"] == "drop"
    assert bed.read_text().strip()
