"""Tests for pipeline/scripts/extract_run_units.py — run-unit extraction (issue #12)."""

import pipeline.scripts.extract_run_units as eru
from pipeline.scripts.extract_run_units import (
    _normalise_unit,
    build_metadata,
    to_header_lines,
)


class TestNormaliseUnit:
    """PU tags reduce to flowcell:lane, dropping the per-sample barcode."""

    def test_picard_flowcell_lane_barcode(self):
        # flowcell.lane.barcode -> flowcell:lane (barcode dropped so co-loaded
        # samples on the same lane match).
        assert _normalise_unit("HJW='FC.3.ACGTACGT") == "HJW='FC:3"

    def test_flowcell_lane_only(self):
        assert _normalise_unit("FLOWCELL1.2") == "FLOWCELL1:2"

    def test_colon_separated(self):
        assert _normalise_unit("FC1:4:BARCODE") == "FC1:4"

    def test_flowcell_only(self):
        assert _normalise_unit("FLOWCELLX") == "FLOWCELLX"

    def test_empty(self):
        assert _normalise_unit("") == ""


def _csv(tmp_path, rows):
    p = tmp_path / "patient.csv"
    lines = ["sample_id,bam_filename,sample_type"]
    for sid, bam, stype in rows:
        lines.append(f"{sid},{bam},{stype}")
    p.write_text("\n".join(lines) + "\n")
    return p


class TestBuildMetadata:
    """The host-share flag and columns, with extraction stubbed per BAM path."""

    def _patch(self, monkeypatch, mapping):
        # mapping: bam path -> (set(run_units), source)
        monkeypatch.setattr(eru, "run_units_for_bam", lambda bam, samtools: mapping[str(bam)])

    def test_admix_sharing_host_run_unit_flagged(self, tmp_path, monkeypatch):
        self._patch(
            monkeypatch,
            {
                "h.bam": ({"FC1:1"}, "RG:PU"),
                "d.bam": ({"FC9:9"}, "RG:PU"),
                "a1.bam": ({"FC1:1"}, "RG:PU"),  # shares with host
                "a2.bam": ({"FC2:2"}, "readname"),  # does not
            },
        )
        csv = _csv(
            tmp_path,
            [
                ("h", "h.bam", "HOST"),
                ("d", "d.bam", "DONOR"),
                ("a1", "a1.bam", "ADMIX"),
                ("a2", "a2.bam", "ADMIX"),
            ],
        )
        rows = {r["sample_id"]: r for r in build_metadata(csv, "samtools")}
        assert rows["h"]["shares_run_with_host"] == "NA"
        assert rows["d"]["shares_run_with_host"] == "false"
        assert rows["a1"]["shares_run_with_host"] == "true"
        assert rows["a2"]["shares_run_with_host"] == "false"
        assert rows["a1"]["run_units"] == "FC1:1"
        assert rows["a2"]["run_unit_source"] == "readname"

    def test_unknown_host_unit_is_undetermined(self, tmp_path, monkeypatch):
        # Host run unit unrecoverable (SRA-like): the comparison cannot be made,
        # so the flag is "unknown", not a false "no".
        self._patch(
            monkeypatch,
            {"h.bam": ({"unknown"}, "unknown"), "a.bam": ({"FC1:1"}, "RG:PU")},
        )
        csv = _csv(tmp_path, [("h", "h.bam", "HOST"), ("a", "a.bam", "ADMIX")])
        rows = {r["sample_id"]: r for r in build_metadata(csv, "samtools")}
        assert rows["a"]["shares_run_with_host"] == "unknown"

    def test_header_lines_admix_only_known_only(self, tmp_path, monkeypatch):
        # Header fragment: only ADMIX samples with a recoverable run unit.
        self._patch(
            monkeypatch,
            {
                "h.bam": ({"FC1:1"}, "RG:PU"),
                "d.bam": ({"FC9:9"}, "RG:PU"),  # donor: excluded (not in admix VCF)
                "a1.bam": ({"FC1:1"}, "RG:PU"),  # admix, known, shares
                "a2.bam": ({"unknown"}, "unknown"),  # admix, unknown: excluded
            },
        )
        csv = _csv(
            tmp_path,
            [
                ("h", "h.bam", "HOST"),
                ("d", "d.bam", "DONOR"),
                ("a1", "a1.bam", "ADMIX"),
                ("a2", "a2.bam", "ADMIX"),
            ],
        )
        lines = to_header_lines(build_metadata(csv, "samtools"))
        assert lines == [
            "##allomixRunUnit=<ID=a1,RunUnit=FC1:1,Source=RG:PU,SharesRunWithHost=true>"
        ]

    def test_header_lines_empty_when_all_unknown(self, tmp_path, monkeypatch):
        # Nothing recoverable -> empty fragment -> annotate is a no-op (optional).
        self._patch(
            monkeypatch,
            {"h.bam": ({"unknown"}, "unknown"), "a.bam": ({"unknown"}, "unknown")},
        )
        csv = _csv(tmp_path, [("h", "h.bam", "HOST"), ("a", "a.bam", "ADMIX")])
        assert to_header_lines(build_metadata(csv, "samtools")) == []

    def test_multi_lane_bam_shares_on_intersection(self, tmp_path, monkeypatch):
        # A BAM merged across lanes carries several units; a share is any overlap.
        self._patch(
            monkeypatch,
            {
                "h.bam": ({"FC1:1", "FC1:2"}, "RG:PU"),
                "a.bam": ({"FC1:2", "FC3:3"}, "RG:PU"),  # overlaps on FC1:2
            },
        )
        csv = _csv(tmp_path, [("h", "h.bam", "HOST"), ("a", "a.bam", "ADMIX")])
        rows = {r["sample_id"]: r for r in build_metadata(csv, "samtools")}
        assert rows["a"]["shares_run_with_host"] == "true"
        # run_units column is sorted and ;-joined.
        assert rows["a"]["run_units"] == "FC1:2;FC3:3"
