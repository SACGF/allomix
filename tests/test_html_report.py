"""Tests for the HTML chimerism report (issue #27).

Determinism is enforced on the HTML text, not on chart pixels: the only source
of nondeterminism is the generation timestamp, which is passed in and pinned
here. These tests cover the single-sample report (header, headline, host-presence
callout, footer), the CLI output-flag wiring (``--json`` / ``--html``), and the
standalone ``report`` subcommand that renders HTML from a saved JSON.
"""

import io
import json
from html.parser import HTMLParser
from pathlib import Path

import pytest

from allomix.cli import main
from allomix.contamination import ContaminationResult
from allomix.detect import HostPresenceResult
from allomix.qc import QCReport
from allomix.relatedness import AdmixConsistencyResult, RelatednessResult
from allomix.report import DonorMeta, ReportMeta, to_html
from allomix.results import ChimerismResult, MarkerResult, MultiDonorResult
from allomix.runmeta import RunUnitInfo

TEST_DATA_DIR = Path(__file__).resolve().parent / "test_data"
JOINT_VCF = TEST_DATA_DIR / "joint_single_donor.vcf"

FIXED_TS = "2026-06-28 12:00:00"


# ---------------------------------------------------------------------------
# Fixtures / helpers
# ---------------------------------------------------------------------------


def _markers(n: int = 12) -> list[MarkerResult]:
    return [
        MarkerResult(
            chrom="chr1",
            pos=100 + i,
            marker_type=0,
            expected_vaf=0.0,
            observed_vaf=0.005,
            residual=0.005,
            ad_ref=990,
            ad_alt=10,
            dp=1000,
            included=True,
        )
        for i in range(n)
    ]


def _host_presence(detected: bool = True, n_markers: int = 30) -> HostPresenceResult:
    return HostPresenceResult(
        n_markers=n_markers,
        n_donor_absent_reads=120 if detected else 5,
        expected_background=40.0,
        poisson_pval=1e-9 if detected else 0.8,
        lrt_pval=2e-8 if detected else 0.6,
        f_host_mle=0.004 if detected else 0.0,
        f_host_ci=(0.002, 0.006) if detected else (0.0, 0.0),
        used_per_site_error=True,
        error_rate_source="per-site",
        n_artifact_filtered=2,
    )


def _result(
    *,
    lob: float = 0.002,
    lod: float = 0.004,
    host_presence: HostPresenceResult | None = None,
) -> ChimerismResult:
    return ChimerismResult(
        donor_fraction=0.995,
        donor_fraction_ci=(0.990, 0.998),
        host_fraction=0.005,
        log_likelihood=-10.0,
        n_informative=12,
        n_markers_used=12,
        per_marker=_markers(),
        error_rate=0.001,
        lob_fraction=lob,
        lod_fraction=lod,
        host_presence=host_presence,
    )


def _qc(status: str = "REVIEW", warnings: list[str] | None = None) -> QCReport:
    return QCReport(
        n_total_markers=76,
        n_shared_markers=70,
        n_informative=12,
        n_used=12,
        n_excluded_depth=2,
        n_excluded_quality=0,
        n_excluded_outlier=0,
        mean_depth=1200.0,
        median_depth=1180.0,
        min_depth=900,
        goodness_of_fit_pval=0.6,
        status=status,
        warnings=warnings if warnings is not None else [],
    )


def _params() -> dict:
    return {
        "genotype_vcf": "/data/genotype.vcf.gz",
        "admix_vcf": "/data/admix.vcf.gz",
        "error_table": "/data/err.tsv",
        "bias_table": None,
        "contamination_table": None,
        "min_dp": 50,
        "min_gq": 20,
        "error_rate": 0.001,
        "robust": "auto",
        "robust_k": 4.0,
        "marker_type_overdispersion": True,
        "bias_correction": True,
        "estimate_bias": False,
        "error_correction": True,
        "contamination_correction": False,
        "host_presence": True,
        "artifact_filter": True,
        "use_sex_chroms": False,
    }


def _render(result, qc, *, sample_name="S1", meta=None, params=None) -> str:
    buf = io.StringIO()
    to_html(
        result,
        qc,
        buf,
        sample_name=sample_name,
        meta=meta,
        timestamp=FIXED_TS,
        params=params,
    )
    return buf.getvalue()


class _H1Counter(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.h1 = 0

    def handle_starttag(self, tag, attrs):
        if tag == "h1":
            self.h1 += 1


# ---------------------------------------------------------------------------
# Smoke / structure
# ---------------------------------------------------------------------------


class TestStructure:
    def test_parses_and_single_h1(self):
        html = _render(_result(host_presence=_host_presence()), _qc(), params=_params())
        counter = _H1Counter()
        counter.feed(html)
        assert counter.h1 == 1
        assert html.lstrip().startswith("<!DOCTYPE html>")
        assert html.rstrip().endswith("</html>")

    def test_self_contained(self):
        """No external network references: CSS and JS are inlined."""
        html = _render(_result(), _qc(), params=_params())
        assert "<style>" in html and "<script>" in html
        assert "http://" not in html and "https://" not in html
        assert "src=" not in html  # no external image/script src in phase 2


# ---------------------------------------------------------------------------
# Determinism
# ---------------------------------------------------------------------------


class TestDeterminism:
    def test_byte_stable_with_fixed_timestamp(self):
        meta = ReportMeta(recipient_id="DEMO-1")
        a = _render(_result(host_presence=_host_presence()), _qc(), meta=meta, params=_params())
        b = _render(_result(host_presence=_host_presence()), _qc(), meta=meta, params=_params())
        assert a == b


# ---------------------------------------------------------------------------
# Field fidelity
# ---------------------------------------------------------------------------


class TestHeadline:
    def test_donor_and_host_fractions_and_ci(self):
        html = _render(_result(), _qc(), params=_params())
        assert "99.50%" in html  # donor fraction
        assert "(95% CI 99.0 to 99.8)" in html  # donor CI
        assert "0.50%" in html  # host fraction
        # host CI is the donor CI reflected about 1
        assert "(95% CI 0.2 to 1.0)" in html

    def test_sensitivity_line(self):
        html = _render(_result(lob=0.002, lod=0.004), _qc(), params=_params())
        assert "limit of blank 0.20%" in html
        assert "limit of detection 0.40%" in html

    def test_infinite_lod_renders_na_not_inf(self):
        html = _render(_result(lob=float("inf"), lod=float("inf")), _qc(), params=_params())
        assert "limit of blank —" in html
        assert "limit of detection —" in html
        # The undetectable limits must not leak the literal float repr anywhere a
        # reader sees a value (the JS sorter's "Infinity" token is unrelated).
        assert "inf%" not in html.lower()
        assert "limit of detection inf" not in html.lower()

    def test_verdict_badge(self):
        html = _render(_result(), _qc(status="REVIEW"), params=_params())
        assert "badge-review" in html
        assert "REVIEW" in html

    def test_multi_donor_headline(self):
        per_marker = _markers()
        result = MultiDonorResult(
            donor_fractions=[0.60, 0.30],
            donor_fraction_cis=[(0.58, 0.62), (0.28, 0.32)],
            host_fraction=0.10,
            log_likelihood=-10.0,
            n_informative=12,
            n_markers_used=12,
            per_marker=per_marker,
            error_rate=0.001,
        )
        html = _render(result, _qc(status="PASS"), params=_params())
        assert "60.00%" in html
        assert "30.00%" in html
        assert "Combined donor" in html
        assert "90.00%" in html  # combined total
        assert "10.00%" in html  # host


# ---------------------------------------------------------------------------
# Host-presence callout (no silent omission)
# ---------------------------------------------------------------------------


class TestHostPresence:
    def test_detected_is_amber(self):
        html = _render(
            _result(host_presence=_host_presence(detected=True)), _qc(), params=_params()
        )
        assert "Low-level host signal detected" in html
        assert "callout-amber" in html

    def test_not_detected_is_neutral(self):
        html = _render(
            _result(host_presence=_host_presence(detected=False)), _qc(), params=_params()
        )
        assert "No host signal above background" in html

    def test_disabled_states_so(self):
        html = _render(_result(host_presence=None), _qc(), params=_params())
        assert "Detection disabled" in html

    def test_no_markers_states_not_assessable(self):
        hp = _host_presence(n_markers=0)
        html = _render(_result(host_presence=hp), _qc(), params=_params())
        assert "Not assessable" in html


# ---------------------------------------------------------------------------
# Metadata / footer
# ---------------------------------------------------------------------------


class TestMetaAndFooter:
    def test_renders_without_meta(self):
        html = _render(_result(), _qc(), meta=None, params=None)
        counter = _H1Counter()
        counter.feed(html)
        assert counter.h1 == 1

    def test_header_shows_metadata(self):
        meta = ReportMeta(
            recipient_id="DEMO-1",
            donors=[DonorMeta(donor_id="M1", relationship="unrelated")],
            transplant_date="2026-01-01",
            sample_dates={"S1": "2026-03-01"},
        )
        html = _render(_result(), _qc(), meta=meta, params=_params())
        assert "DEMO-1" in html
        assert "M1 (unrelated)" in html
        assert "+59" in html  # days post-transplant: 2026-01-01 -> 2026-03-01

    def test_days_post_transplant_omitted_when_dates_missing(self):
        meta = ReportMeta(recipient_id="DEMO-1")  # no transplant/sample dates
        html = _render(_result(), _qc(), meta=meta, params=_params())
        assert "Days post-transplant" not in html

    def test_footer_basenames_only_and_citation(self):
        html = _render(_result(), _qc(), params=_params())
        assert "genotype.vcf.gz" in html
        assert "/data/" not in html  # no patient-identifying paths
        assert "Crysup" in html

    def test_html_escaping(self):
        meta = ReportMeta(recipient_id="<script>alert(1)</script>")
        html = _render(_result(), _qc(), meta=meta, params=_params())
        assert "<script>alert(1)</script>" not in html
        assert "&lt;script&gt;" in html


# ---------------------------------------------------------------------------
# QC panel (phase 3)
# ---------------------------------------------------------------------------


def _contamination() -> ContaminationResult:
    return ContaminationResult(
        n_markers=25,
        contamination_fraction=0.003,
        median_minor_frac=0.003,
        error_floor=1e-4,
        floor_empirical=True,
        pooled_minor_frac=0.003,
        n_minor_reads=75,
        total_depth=25000,
        p_value=1e-5,
        n_excluded_high=1,
        used_per_site_error=True,
        error_rate_source="per-site",
    )


def _admix_consistency() -> AdmixConsistencyResult:
    return AdmixConsistencyResult(
        n_consensus_hom=40,
        n_discordant=1,
        discordant_fraction=0.025,
        swap_pval=0.2,
    )


def _relatedness() -> list[RelatednessResult]:
    return [
        RelatednessResult(
            a_name="HOST",
            b_name="DONOR",
            coefficient=0.02,
            ci_low=0.0,
            ci_high=0.05,
            confidence="high",
            relationship="unrelated",
            degree=None,
            n_sites=70,
            het_a=30,
            het_b=28,
            shared_hets=8,
            ibs0=5,
        )
    ]


class TestQCPanel:
    def test_marker_counts_and_depth(self):
        html = _render(_result(), _qc(), params=_params())
        assert "Quality control" in html
        assert "Total markers (input)" in html
        assert "Mean depth" in html
        assert "1200x" in html  # mean_depth formatted with unit

    def test_warnings_listed_verbatim(self):
        warn = "Robust refit excluded 3 markers as residual outliers"
        html = _render(_result(), _qc(status="REVIEW", warnings=[warn]), params=_params())
        assert warn in html

    def test_no_warnings_states_none(self):
        html = _render(_result(), _qc(status="PASS", warnings=[]), params=_params())
        assert "No QC warnings raised." in html

    def test_overall_verdict_badge_present(self):
        html = _render(_result(), _qc(status="FAIL"), params=_params())
        assert "Overall verdict:" in html
        assert "badge-fail" in html

    def test_optional_subresults_render_when_present(self):
        qc = _qc(status="REVIEW")
        qc.contamination = _contamination()
        qc.admix_consistency = _admix_consistency()
        qc.relatedness = _relatedness()
        qc.run_unit = RunUnitInfo(run_unit="FC1.L2", source="header", shares_run_with_host=False)
        html = _render(_result(), qc, params=_params())
        assert "Contamination fraction" in html
        assert "Swap discordant fraction" in html
        assert "Reference-sample relatedness" in html
        assert "HOST vs DONOR" in html
        assert "FC1.L2" in html

    def test_optional_subresults_state_not_assessed_when_absent(self):
        html = _render(_result(), _qc(), params=_params())
        assert "not assessed" in html  # contamination / swap default to absent

    def test_clinician_breakdown_summary(self):
        qc = _qc(status="REVIEW", warnings=["one flag", "two flag"])
        html = _render(_result(), qc, params=_params())
        # Marker accounting line for the clinician (informative / used).
        assert "12 of 76 input markers were informative" in html
        assert "12 used in the fit" in html
        # QC flag count summary.
        assert "2 QC checks flagged" in html

    def test_no_per_marker_table_in_html(self):
        """The per-marker detail is the CSV's job; HTML must not carry the table.

        ``<details>`` is allowed (the report uses it for the collapsible method
        help); what must be absent is the per-marker table itself.
        """
        html = _render(_result(), _qc(), params=_params())
        assert "Marker detail" not in html
        assert "Observed VAF" not in html  # the per-marker table header is gone

    def test_expandable_method_help_present(self):
        """The report carries collapsible 'how this works' help (feedback item 2)."""
        html = _render(_result(host_presence=_host_presence()), _qc(), params=_params())
        assert "<details class=\"help\">" in html
        assert "How host-presence detection works" in html


# ---------------------------------------------------------------------------
# Per-marker CSV (bioinformatician-facing detail)
# ---------------------------------------------------------------------------


class TestMarkerCSV:
    def _read_csv(self, text: str) -> tuple[list[str], list[list[str]]]:
        import csv

        rows = list(csv.reader(io.StringIO(text)))
        return rows[0], rows[1:]

    def test_single_sample_csv_columns_and_rows(self):
        from allomix.report import to_marker_csv

        buf = io.StringIO()
        to_marker_csv([("S1", _result())], buf)
        header, rows = self._read_csv(buf.getvalue())
        assert header[0] == "sample"
        assert "marker_type_label" in header
        assert len(rows) == 12  # one per marker
        assert all(r[0] == "S1" for r in rows)
        # Type-0 markers map to the readable host/donor label.
        label_idx = header.index("marker_type_label")
        assert rows[0][label_idx] == "host 0/0, donor 1/1"

    def test_excluded_flag_round_trips(self):
        from allomix.report import to_marker_csv

        markers = _markers(5)
        markers[0].included = False
        result = _result()
        result.per_marker = markers
        buf = io.StringIO()
        to_marker_csv([("S1", result)], buf)
        header, rows = self._read_csv(buf.getvalue())
        inc_idx = header.index("included")
        assert rows[0][inc_idx] == "False"
        assert rows[1][inc_idx] == "True"

    def test_multi_sample_csv_has_sample_column(self):
        from allomix.report import to_marker_csv

        buf = io.StringIO()
        to_marker_csv([("A", _result()), ("B", _result())], buf)
        _, rows = self._read_csv(buf.getvalue())
        samples = {r[0] for r in rows}
        assert samples == {"A", "B"}

    def test_cli_marker_csv_written(self, tmp_path):
        out = tmp_path / "report.html"
        csv_path = tmp_path / "markers.csv"
        rc = main(
            [
                "monitor",
                "--genotype-vcf",
                str(JOINT_VCF),
                "--admix-vcf",
                str(JOINT_VCF),
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR",
                "--sample",
                "ADMIX_F0.10",
                "--html",
                str(out),
                "--marker-csv",
                str(csv_path),
                "--min-dp",
                "0",
                "--min-gq",
                "0",
            ]
        )
        assert rc == 0
        assert csv_path.exists()
        text = csv_path.read_text()
        assert text.splitlines()[0].startswith("sample,chrom,pos,marker_type")
        # The HTML still carries no per-marker table.
        assert "Marker detail" not in out.read_text()


# ---------------------------------------------------------------------------
# Timeline report (phase 4)
# ---------------------------------------------------------------------------


def _tp(name: str, donor_frac: float, status: str = "PASS"):
    """One (name, result, qc) timeline tuple at a given donor fraction."""
    result = ChimerismResult(
        donor_fraction=donor_frac,
        donor_fraction_ci=(max(0.0, donor_frac - 0.01), min(1.0, donor_frac + 0.01)),
        host_fraction=1.0 - donor_frac,
        log_likelihood=-10.0,
        n_informative=12,
        n_markers_used=12,
        per_marker=_markers(),
        error_rate=0.001,
        host_presence=_host_presence(detected=False),
    )
    return (name, result, _qc(status=status))


def _render_timeline(results, *, meta=None, params=None, log_scale=False) -> str:
    from allomix.html.timeline import timeline_html

    buf = io.StringIO()
    timeline_html(
        results,
        buf,
        meta=meta,
        timestamp=FIXED_TS,
        params=params,
        log_scale=log_scale,
    )
    return buf.getvalue()


class TestTimeline:
    def test_renders_with_trend_chart(self):
        results = [_tp("T1", 0.10), _tp("T2", 0.30), _tp("T3", 0.55)]
        html = _render_timeline(results, params=_params())
        assert "<h2>Trend</h2>" in html
        assert "<h2>Timepoints</h2>" in html
        assert "data:image/png;base64," in html  # chart embedded, no external src
        assert "https://" not in html

    def test_headline_shows_latest_and_delta(self):
        results = [_tp("T1", 0.10), _tp("T2", 0.30)]
        html = _render_timeline(results, params=_params())
        assert "Donor fraction (latest)" in html
        assert "30.00%" in html  # latest
        assert "Previous" in html
        assert "20.00 pp" in html  # change 30 - 10

    def test_single_timepoint_note_not_chart(self):
        html = _render_timeline([_tp("T1", 0.10)], params=_params())
        assert "at least two timepoints" in html
        assert "data:image/png;base64," not in html

    def test_deterministic_with_fixed_timestamp(self):
        results = [_tp("T1", 0.10), _tp("T2", 0.30), _tp("T3", 0.55)]
        a = _render_timeline(results, params=_params())
        b = _render_timeline(results, params=_params())
        assert a == b  # includes the PNG bytes: same data, same process

    def test_empty_results_raises(self):
        with pytest.raises(ValueError):
            _render_timeline([])

    def test_cli_timeline_format_html(self, tmp_path):
        out = tmp_path / "timeline.html"
        rc = main(
            [
                "timeline",
                "--genotype-vcf",
                str(JOINT_VCF),
                "--admix-vcf",
                str(JOINT_VCF),
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR",
                "--sample",
                "ADMIX_F0.10",
                "--sample",
                "ADMIX_F0.50",
                "--sample",
                "ADMIX_F1.00",
                "--html",
                str(out),
                "--min-dp",
                "0",
                "--min-gq",
                "0",
            ]
        )
        assert rc == 0
        html = out.read_text()
        assert html.lstrip().startswith("<!DOCTYPE html>")
        assert "data:image/png;base64," in html
        counter = _H1Counter()
        counter.feed(html)
        assert counter.h1 == 1


# ---------------------------------------------------------------------------
# CLI wiring
# ---------------------------------------------------------------------------


class TestCLI:
    def test_monitor_format_html(self, tmp_path):
        out = tmp_path / "report.html"
        rc = main(
            [
                "monitor",
                "--genotype-vcf",
                str(JOINT_VCF),
                "--admix-vcf",
                str(JOINT_VCF),
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR",
                "--sample",
                "ADMIX_F0.10",
                "--html",
                str(out),
                "--min-dp",
                "0",
                "--min-gq",
                "0",
                "--recipient-id",
                "PATIENT-X",
                "--donor-relationship",
                "unrelated",
            ]
        )
        assert rc == 0
        html = out.read_text()
        assert html.lstrip().startswith("<!DOCTYPE html>")
        assert "Chimerism report" in html
        assert "PATIENT-X" in html
        counter = _H1Counter()
        counter.feed(html)
        assert counter.h1 == 1

    def test_report_timestamp_override_is_deterministic(self, tmp_path):
        """--report-timestamp pins the only wall-clock field, so output is stable."""
        argv = [
            "monitor",
            "--genotype-vcf",
            str(JOINT_VCF),
            "--admix-vcf",
            str(JOINT_VCF),
            "--host-sample",
            "HOST",
            "--donor-sample",
            "DONOR",
            "--sample",
            "ADMIX_F0.10",
            "--min-dp",
            "0",
            "--min-gq",
            "0",
            "--report-timestamp",
            "FIXED-TS-123",
        ]
        a = tmp_path / "a.html"
        b = tmp_path / "b.html"
        assert main(argv + ["--html", str(a)]) == 0
        assert main(argv + ["--html", str(b)]) == 0
        assert a.read_text() == b.read_text()
        assert "FIXED-TS-123" in a.read_text()

    def test_monitor_html_rejects_multiple_samples(self, tmp_path):
        out = tmp_path / "report.html"
        with pytest.raises(SystemExit):
            main(
                [
                    "monitor",
                    "--genotype-vcf",
                    str(JOINT_VCF),
                    "--admix-vcf",
                    str(JOINT_VCF),
                    "--host-sample",
                    "HOST",
                    "--donor-sample",
                    "DONOR",
                    "--sample",
                    "ADMIX_F0.10",
                    "--sample",
                    "ADMIX_F0.50",
                    "--html",
                    str(out),
                    "--min-dp",
                    "0",
                    "--min-gq",
                    "0",
                ]
            )


# ---------------------------------------------------------------------------
# Structured-data first: JSON is the artifact the report is rendered from
# ---------------------------------------------------------------------------


def _monitor_argv(*, sample="ADMIX_F0.10"):
    return [
        "monitor",
        "--genotype-vcf",
        str(JOINT_VCF),
        "--admix-vcf",
        str(JOINT_VCF),
        "--host-sample",
        "HOST",
        "--donor-sample",
        "DONOR",
        "--sample",
        sample,
        "--min-dp",
        "0",
        "--min-gq",
        "0",
        "--recipient-id",
        "PATIENT-X",
        "--report-timestamp",
        FIXED_TS,
    ]


class TestStructuredOutputs:
    """monitor emits structured data; the report is rendered from it."""

    def test_json_and_html_in_one_invocation(self, tmp_path):
        js = tmp_path / "r.json"
        html = tmp_path / "r.html"
        rc = main(_monitor_argv() + ["--json", str(js), "--html", str(html)])
        assert rc == 0
        # Structured artifact: the envelope, with analysis nested and provenance.
        data = json.loads(js.read_text())
        assert data["kind"] == "single"
        assert "donor_pct" in data["analysis"]
        assert data["meta"]["recipient_id"] == "PATIENT-X"
        # The HTML is a self-contained report.
        assert html.read_text().lstrip().startswith("<!DOCTYPE html>")
        assert "PATIENT-X" in html.read_text()

    def test_report_subcommand_matches_inline_html(self, tmp_path):
        """`report file.json` reproduces the HTML the monitor step would write,
        byte for byte, because both render the same envelope."""
        js = tmp_path / "r.json"
        inline = tmp_path / "inline.html"
        assert main(_monitor_argv() + ["--json", str(js), "--html", str(inline)]) == 0

        from_json = tmp_path / "from_json.html"
        assert main(["report", str(js), "--output", str(from_json)]) == 0
        assert from_json.read_text() == inline.read_text()

    def test_default_output_is_tsv_on_stdout(self, capsys):
        assert main(_monitor_argv()) == 0
        out = capsys.readouterr().out
        assert out.startswith("sample\t")


class TestReportSubcommandTimeline:
    def test_report_renders_timeline_from_json(self, tmp_path):
        js = tmp_path / "tl.json"
        rc = main(
            [
                "timeline",
                "--genotype-vcf",
                str(JOINT_VCF),
                "--admix-vcf",
                str(JOINT_VCF),
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR",
                "--sample",
                "ADMIX_F0.10",
                "--sample",
                "ADMIX_F0.50",
                "--sample",
                "ADMIX_F1.00",
                "--json",
                str(js),
                "--min-dp",
                "0",
                "--min-gq",
                "0",
            ]
        )
        assert rc == 0
        data = json.loads(js.read_text())
        assert data["kind"] == "timeline"
        assert len(data["timepoints"]) == 3

        html = tmp_path / "tl.html"
        assert main(["report", str(js), "--output", str(html)]) == 0
        text = html.read_text()
        assert "data:image/png;base64," in text  # trend chart drawn from the JSON
        assert "<h2>Trend</h2>" in text
