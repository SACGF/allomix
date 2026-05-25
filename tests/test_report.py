"""Tests for allomix.report — TSV and JSON output formatting."""

from __future__ import annotations

import inspect
import io
import json

import pytest

from allomix.qc import ChimerismResult, MarkerResult, QCReport
from allomix.report import timeline_json, to_json, to_tsv


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_marker_result(
    chrom: str = "chr1",
    pos: int = 100,
    marker_type: int = 0,
    expected_vaf: float = 0.10,
    observed_vaf: float = 0.105,
    residual: float = 0.005,
    ad_ref: int = 895,
    ad_alt: int = 105,
    dp: int = 1000,
    included: bool = True,
) -> MarkerResult:
    return MarkerResult(
        chrom=chrom,
        pos=pos,
        marker_type=marker_type,
        expected_vaf=expected_vaf,
        observed_vaf=observed_vaf,
        residual=residual,
        ad_ref=ad_ref,
        ad_alt=ad_alt,
        dp=dp,
        included=included,
    )


def _make_chimerism_result(
    donor_fraction: float = 0.1234,
    ci: tuple[float, float] = (0.1102, 0.1371),
    n_informative: int = 10,
    per_marker: list[MarkerResult] | None = None,
    lob_fraction: float = 0.0021,
    lod_fraction: float = 0.0045,
) -> ChimerismResult:
    if per_marker is None:
        per_marker = [_make_marker_result(pos=i * 100) for i in range(n_informative)]
    n_used = sum(1 for m in per_marker if m.included)
    return ChimerismResult(
        donor_fraction=donor_fraction,
        donor_fraction_ci=ci,
        host_fraction=1.0 - donor_fraction,
        log_likelihood=-100.0,
        n_informative=n_informative,
        n_markers_used=n_used,
        per_marker=per_marker,
        error_rate=0.01,
        lob_fraction=lob_fraction,
        lod_fraction=lod_fraction,
    )


def _make_qc_report(
    n_used: int = 10,
    mean_depth: float = 1500.0,
    gof_pval: float | None = 0.45,
    pass_: bool = True,
    warnings: list[str] | None = None,
) -> QCReport:
    return QCReport(
        n_total_markers=76,
        n_shared_markers=50,
        n_informative=10,
        n_used=n_used,
        n_excluded_depth=3,
        n_excluded_quality=0,
        n_excluded_outlier=0,
        mean_depth=mean_depth,
        median_depth=1400.0,
        min_depth=800,
        goodness_of_fit_pval=gof_pval,
        warnings=warnings if warnings is not None else [],
        pass_=pass_,
    )


# ---------------------------------------------------------------------------
# TSV output
# ---------------------------------------------------------------------------


class TestTsvOutput:
    """Test TSV summary output."""

    def test_tsv_has_header_and_data(self):
        result = _make_chimerism_result()
        qc = _make_qc_report()
        buf = io.StringIO()
        to_tsv(result, qc, buf)
        content = buf.getvalue()
        lines = content.strip().split("\n")
        assert len(lines) == 2
        assert lines[0].startswith("sample\t")

    def test_tsv_values_parseable(self):
        result = _make_chimerism_result(
            donor_fraction=0.1234,
            ci=(0.1102, 0.1371),
        )
        qc = _make_qc_report(n_used=10, mean_depth=1500.0, gof_pval=0.45)
        buf = io.StringIO()
        to_tsv(result, qc, buf)
        lines = buf.getvalue().strip().split("\n")
        fields = lines[1].split("\t")
        assert float(fields[1]) == pytest.approx(12.34, abs=0.01)
        assert float(fields[2]) == pytest.approx(11.02, abs=0.01)
        assert float(fields[3]) == pytest.approx(13.71, abs=0.01)
        assert float(fields[4]) == pytest.approx(0.21, abs=0.01)  # lob_pct
        assert float(fields[5]) == pytest.approx(0.45, abs=0.01)  # lod_pct
        assert int(fields[6]) == 10
        assert int(fields[7]) == 10
        assert float(fields[8]) == pytest.approx(1500, abs=1)
        assert float(fields[9]) == pytest.approx(0.45, abs=0.01)
        assert fields[10] == "PASS"

    def test_tsv_fail_status(self):
        result = _make_chimerism_result()
        qc = _make_qc_report(pass_=False)
        buf = io.StringIO()
        to_tsv(result, qc, buf)
        lines = buf.getvalue().strip().split("\n")
        fields = lines[1].split("\t")
        assert fields[10] == "FAIL"

    def test_tsv_warnings_column(self):
        result = _make_chimerism_result()
        qc = _make_qc_report(warnings=["Insufficient informative markers: 1 < 3", "Low depth"])
        buf = io.StringIO()
        to_tsv(result, qc, buf)
        lines = buf.getvalue().splitlines()
        assert lines[0].split("\t")[11] == "qc_warnings"
        assert lines[1].split("\t")[11] == "Insufficient informative markers: 1 < 3; Low depth"

    def test_tsv_warnings_column_empty_when_clean(self):
        result = _make_chimerism_result()
        qc = _make_qc_report(warnings=[])
        buf = io.StringIO()
        to_tsv(result, qc, buf)
        # Trailing empty cell: line ends with a tab then nothing.
        assert buf.getvalue().splitlines()[1].split("\t")[11] == ""

    def test_tsv_gof_na(self):
        result = _make_chimerism_result()
        qc = _make_qc_report(gof_pval=None)
        buf = io.StringIO()
        to_tsv(result, qc, buf)
        lines = buf.getvalue().strip().split("\n")
        fields = lines[1].split("\t")
        assert fields[9] == "NA"


class TestTsvVerbose:
    """Test verbose TSV output with per-marker detail lines."""

    def test_verbose_has_detail_lines(self):
        markers = [
            _make_marker_result(pos=100),
            _make_marker_result(pos=200),
            _make_marker_result(pos=300),
        ]
        result = _make_chimerism_result(n_informative=3, per_marker=markers)
        qc = _make_qc_report(n_used=3)
        buf = io.StringIO()
        to_tsv(result, qc, buf, verbose=True)
        content = buf.getvalue()
        lines = content.split("\n")

        # Find blank separator line
        blank_idx = None
        for i, line in enumerate(lines):
            if line.strip() == "":
                blank_idx = i
                break
        assert blank_idx is not None, "Expected blank separator line in verbose output"

        # Detail header follows the blank line
        detail_header = lines[blank_idx + 1]
        assert "chrom" in detail_header
        assert "observed_vaf" in detail_header

        # Detail data lines follow
        detail_lines = [line for line in lines[blank_idx + 2 :] if line.strip()]
        assert len(detail_lines) == 3

    def test_non_verbose_has_no_detail(self):
        result = _make_chimerism_result()
        qc = _make_qc_report()
        buf = io.StringIO()
        to_tsv(result, qc, buf, verbose=False)
        content = buf.getvalue()
        assert "observed_vaf" not in content


# ---------------------------------------------------------------------------
# JSON output
# ---------------------------------------------------------------------------


class TestJsonOutput:
    """Test JSON dict output."""

    def test_json_has_required_keys(self):
        result = _make_chimerism_result()
        qc = _make_qc_report()
        d = to_json(result, qc, sample_name="day30")

        required_keys = {
            "sample",
            "donor_pct",
            "ci_lo",
            "ci_hi",
            "n_informative",
            "n_used",
            "mean_depth",
            "gof_pval",
            "qc_pass",
            "warnings",
            "markers",
        }
        assert required_keys <= set(d.keys())

    def test_json_types(self):
        result = _make_chimerism_result()
        qc = _make_qc_report()
        d = to_json(result, qc, sample_name="day30")

        assert isinstance(d["sample"], str)
        assert isinstance(d["donor_pct"], float)
        assert isinstance(d["ci_lo"], float)
        assert isinstance(d["ci_hi"], float)
        assert isinstance(d["n_informative"], int)
        assert isinstance(d["n_used"], int)
        assert isinstance(d["mean_depth"], float)
        assert isinstance(d["qc_pass"], bool)
        assert isinstance(d["warnings"], list)
        assert isinstance(d["markers"], list)

    def test_json_values(self):
        result = _make_chimerism_result(
            donor_fraction=0.1234,
            ci=(0.1102, 0.1371),
        )
        qc = _make_qc_report()
        d = to_json(result, qc, sample_name="day30")
        assert d["sample"] == "day30"
        assert d["donor_pct"] == pytest.approx(12.34, abs=0.01)
        assert d["ci_lo"] == pytest.approx(11.02, abs=0.01)
        assert d["ci_hi"] == pytest.approx(13.71, abs=0.01)

    def test_json_serialisable(self):
        result = _make_chimerism_result()
        qc = _make_qc_report()
        d = to_json(result, qc)
        # Should not raise
        json_str = json.dumps(d)
        assert isinstance(json_str, str)

    def test_json_markers_present(self):
        markers = [_make_marker_result(pos=100)]
        result = _make_chimerism_result(n_informative=1, per_marker=markers)
        qc = _make_qc_report(n_used=1)
        d = to_json(result, qc)
        assert len(d["markers"]) == 1
        m = d["markers"][0]
        assert "chrom" in m
        assert "pos" in m
        assert "observed_vaf" in m
        assert "included" in m

    def test_json_gof_null(self):
        result = _make_chimerism_result()
        qc = _make_qc_report(gof_pval=None)
        d = to_json(result, qc)
        assert d["gof_pval"] is None
        # Verify it serialises to JSON null
        json_str = json.dumps(d)
        parsed = json.loads(json_str)
        assert parsed["gof_pval"] is None


# ---------------------------------------------------------------------------
# Timeline JSON
# ---------------------------------------------------------------------------


class TestTimelineJson:
    """Test multi-timepoint timeline output."""

    def test_timeline_structure(self):
        entries = []
        for i, name in enumerate(["day30", "day60", "day90"]):
            frac = 0.10 - i * 0.02
            result = _make_chimerism_result(
                donor_fraction=frac,
                ci=(frac - 0.02, frac + 0.02),
            )
            qc = _make_qc_report()
            entries.append((name, result, qc))

        tl = timeline_json(entries)
        assert "timepoints" in tl
        assert len(tl["timepoints"]) == 3

    def test_timeline_sample_names(self):
        entries = []
        for name in ["day30", "day60"]:
            result = _make_chimerism_result()
            qc = _make_qc_report()
            entries.append((name, result, qc))

        tl = timeline_json(entries)
        names = [tp["sample"] for tp in tl["timepoints"]]
        assert names == ["day30", "day60"]

    def test_timeline_values(self):
        result = _make_chimerism_result(
            donor_fraction=0.05,
            ci=(0.03, 0.07),
        )
        qc = _make_qc_report()
        tl = timeline_json([("day30", result, qc)])
        tp = tl["timepoints"][0]
        assert tp["donor_pct"] == pytest.approx(5.0, abs=0.01)
        assert tp["ci_lo"] == pytest.approx(3.0, abs=0.01)
        assert tp["ci_hi"] == pytest.approx(7.0, abs=0.01)
        assert tp["qc_pass"] is True

    def test_timeline_serialisable(self):
        result = _make_chimerism_result()
        qc = _make_qc_report()
        tl = timeline_json([("day30", result, qc)])
        json_str = json.dumps(tl)
        assert isinstance(json_str, str)

    def test_timeline_empty(self):
        tl = timeline_json([])
        assert tl == {"timepoints": []}


# ---------------------------------------------------------------------------
# TSV sample_name support
# ---------------------------------------------------------------------------


class TestTsvSampleName:
    """to_tsv should accept a sample_name parameter like to_json."""

    def test_default_sample_name_not_hardcoded(self):
        """TSV default should not be the literal string 'sample'."""
        result = _make_chimerism_result()
        qc = _make_qc_report()
        buf = io.StringIO()
        to_tsv(result, qc, buf)
        lines = buf.getvalue().strip().split("\n")
        first_col = lines[1].split("\t")[0]
        assert first_col != "sample"

    def test_accepts_sample_name_parameter(self):
        """to_tsv should accept a sample_name keyword argument."""
        sig = inspect.signature(to_tsv)
        assert "sample_name" in sig.parameters

    def test_sample_name_appears_in_output(self):
        """Passing sample_name should put it in the TSV data line."""
        result = _make_chimerism_result()
        qc = _make_qc_report()
        buf = io.StringIO()
        to_tsv(result, qc, buf, sample_name="patient_001")
        lines = buf.getvalue().strip().split("\n")
        first_col = lines[1].split("\t")[0]
        assert first_col == "patient_001"


class TestOutputConsistency:
    """JSON and TSV should produce consistent output for the same input."""

    def test_json_and_tsv_donor_fractions_match(self):
        result = _make_chimerism_result(donor_fraction=0.1234, ci=(0.10, 0.15))
        qc = _make_qc_report()

        json_out = to_json(result, qc, sample_name="patient_001")
        buf = io.StringIO()
        to_tsv(result, qc, buf, sample_name="patient_001")
        tsv_cols = buf.getvalue().strip().split("\n")[1].split("\t")

        assert json_out["donor_pct"] == pytest.approx(float(tsv_cols[1]), abs=0.01)
        assert json_out["sample"] == "patient_001"
        assert tsv_cols[0] == "patient_001"
