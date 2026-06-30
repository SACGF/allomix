"""Output formatting for chimerism results.

Provides TSV, JSON, and clinician-facing HTML output for single-sample results
and multi-timepoint timelines, plus a bioinformatician-facing per-marker CSV
(the detail the HTML report omits). Supports both single-donor and multi-donor
results.
"""

import csv
import math
from pathlib import Path
from typing import TextIO

from allomix import __version__
from allomix.genotype import MarkerType
from allomix.qc.host_presence import HostPresenceResult
from allomix.qc.qc import QCReport
from allomix.qc.relatedness import AdmixConsistencyResult, RelatednessResult
from allomix.qc.runmeta import RunUnitInfo
from allomix.qc.sample_contamination import ContaminationResult
from allomix.report.html.meta import DonorMeta, ReportMeta
from allomix.report.html.render import render_single
from allomix.results import ChimerismResult, MultiDonorResult

# Re-exported so the public surface is ``from allomix.report.report import ReportMeta``,
# matching where the other output types live. The dataclasses themselves sit in
# ``allomix.report.html.meta`` to avoid a report.py <-> html.render import cycle.
__all__ = [
    "DonorMeta",
    "ReportMeta",
    "report_data",
    "timeline_json",
    "timeline_report_data",
    "to_html",
    "to_json",
    "to_marker_csv",
    "to_tsv",
]

# Version of the report-data JSON schema (the envelope ``report_data`` /
# ``timeline_report_data`` produce and the HTML renderer consumes). Bump when the
# envelope shape changes so downstream consumers can branch on it.
REPORT_SCHEMA_VERSION = "1"

# Column order for the per-marker CSV.
_MARKER_CSV_COLS = [
    "sample",
    "chrom",
    "pos",
    "marker_type",
    "marker_type_label",
    "ad_ref",
    "ad_alt",
    "dp",
    "observed_vaf",
    "expected_vaf",
    "residual",
    "included",
]


def to_marker_csv(
    results: list[tuple[str, ChimerismResult | MultiDonorResult]],
    output: Path | TextIO,
) -> None:
    """Write per-marker detail for one or more samples as a CSV.

    The bioinformatician-facing artifact: the HTML report omits the per-marker
    table on purpose, so the full marker accounting goes here. One row per marker
    per sample, with a leading ``sample`` column. Marker-type codes are expanded
    to a readable host/donor genotype label (see ``allomix.genotype.MarkerType``).

    Args:
        results: ``(sample_name, result)`` pairs, one-element for a single-sample
            report, one entry per timepoint for a timeline.
    """
    if isinstance(output, (str, Path)):
        with open(output, "w", encoding="utf-8", newline="") as fh:
            _write_marker_csv(results, fh)
    else:
        _write_marker_csv(results, output)


def _write_marker_csv(
    results: list[tuple[str, ChimerismResult | MultiDonorResult]],
    fh: TextIO,
) -> None:
    """Write the per-marker CSV header and one row per marker per sample."""
    writer = csv.writer(fh)
    writer.writerow(_MARKER_CSV_COLS)
    for sample_name, result in results:
        for m in result.per_marker:
            writer.writerow(
                [
                    sample_name,
                    m.chrom,
                    str(m.pos),
                    str(int(m.marker_type)),
                    MarkerType.label_for(m.marker_type),
                    str(m.ad_ref),
                    str(m.ad_alt),
                    str(m.dp),
                    f"{m.observed_vaf:.6f}",
                    f"{m.expected_vaf:.6f}",
                    f"{m.residual:.6f}",
                    str(m.included),
                ]
            )


# Sentinel for the host-presence cells when the detector did not run (e.g.
# --no-host-presence) or produced no usable markers. Keeps the TSV
# rectangular and parseable by downstream scripts.
_NA = "NA"


def _host_presence_tsv_cells(hp: HostPresenceResult | None) -> list[str]:
    """Format the host-presence columns (order matches ``_HOST_PRESENCE_TSV_COLS``).

    ``NA`` for every column when the detector did not run; per-column ``NA`` when
    there are no usable markers (the source flag still carries information).
    """
    if hp is None:
        return [_NA] * 7
    if hp.n_markers == 0:
        # Carry the source flag but leave numeric cells NA: the test did not run.
        # The artifact count stays informative (markers may have been filtered to
        # zero).
        return [
            _NA,
            _NA,
            _NA,
            _NA,
            str(hp.n_markers),
            hp.error_rate_source,
            str(hp.n_artifact_filtered),
        ]
    f_lo, f_hi = hp.f_host_ci
    return [
        f"{hp.lrt_pval:.4g}",
        f"{hp.f_host_mle:.6f}",
        f"{f_lo:.6f}",
        f"{f_hi:.6f}",
        str(hp.n_markers),
        hp.error_rate_source,
        str(hp.n_artifact_filtered),
    ]


def _host_presence_json(hp: HostPresenceResult | None) -> dict | None:
    """JSON view of a HostPresenceResult; mirrors the TSV columns."""
    if hp is None:
        return None
    f_lo, f_hi = hp.f_host_ci
    return {
        "lrt_pval": hp.lrt_pval,
        "poisson_pval": hp.poisson_pval,
        "f_host_mle": hp.f_host_mle,
        "f_host_ci_lo": f_lo,
        "f_host_ci_hi": f_hi,
        "n_markers": hp.n_markers,
        "n_donor_absent_reads": hp.n_donor_absent_reads,
        "expected_background": hp.expected_background,
        "used_per_site_error": hp.used_per_site_error,
        "error_rate_source": hp.error_rate_source,
        "n_artifact_filtered": hp.n_artifact_filtered,
    }


# Appended to the existing TSV summary. Downstream scripts parse by header, so
# renaming or reordering existing columns breaks them; new columns go at the end.
_HOST_PRESENCE_TSV_COLS = [
    "host_present_p",
    "host_f_est",
    "host_f_ci_lo",
    "host_f_ci_hi",
    "host_detect_markers",
    "host_err_source",
    "host_artifact_filtered",
]

# Relatedness columns, appended after the host-presence block. Multi-donor
# host-vs-donor pairs are joined with ";" within each cell to keep the row
# rectangular. The pass/fail verdict lives in qc_status / qc_warnings.
_RELATEDNESS_TSV_COLS = [
    "relatedness",
    "relatedness_ci",
    "relatedness_confidence",
    "relationship",
]


def _relatedness_tsv_cells(relatedness: list[RelatednessResult] | None) -> list[str]:
    """Format the relatedness columns (order matches ``_RELATEDNESS_TSV_COLS``).

    Reports the host-vs-donor pairs in donor order, joined with ";". ``NA`` cells
    when no relatedness was computed.
    """
    if not relatedness:
        return [_NA] * len(_RELATEDNESS_TSV_COLS)
    host_pairs = [r for r in relatedness if r.a_name == "host"]
    if not host_pairs:
        return [_NA] * len(_RELATEDNESS_TSV_COLS)

    coefs, cis, confs, rels = [], [], [], []
    for r in host_pairs:
        if r.coefficient is None:
            coefs.append(_NA)
            cis.append(_NA)
        else:
            coefs.append(f"{r.coefficient:.4f}")
            cis.append(f"{r.ci_low:.4f},{r.ci_high:.4f}")
        confs.append(r.confidence)
        rels.append(r.relationship)
    return [";".join(coefs), ";".join(cis), ";".join(confs), ";".join(rels)]


def _relatedness_json(relatedness: list[RelatednessResult] | None) -> list[dict] | None:
    """JSON view of all relatedness pairs.

    Carries per-pair names and split CI bounds (not just the TSV's joined ``pair``
    / ``ci``) so the HTML QC table renders from this dict alone, without the source
    ``RelatednessResult`` objects.
    """
    if not relatedness:
        return None
    out = []
    for r in relatedness:
        out.append(
            {
                "pair": r.pair,
                "a_name": r.a_name,
                "b_name": r.b_name,
                "coefficient": round(r.coefficient, 6) if r.coefficient is not None else None,
                "ci": [round(r.ci_low, 6), round(r.ci_high, 6)] if r.ci_low is not None else None,
                "ci_low": round(r.ci_low, 6) if r.ci_low is not None else None,
                "ci_high": round(r.ci_high, 6) if r.ci_high is not None else None,
                "confidence": r.confidence,
                "relationship": r.relationship,
                "n_sites": r.n_sites,
            }
        )
    return out


# Contamination columns, appended after the relatedness block. ``contamination_frac``
# is the estimated third-party floor (fraction above sequencing error), the in-data
# half of issue #12.
_CONTAMINATION_TSV_COLS = [
    "contamination_frac",
    "contamination_p",
    "contamination_markers",
]


def _contamination_tsv_cells(contamination: ContaminationResult | None) -> list[str]:
    """Format the three contamination columns. Order matches the cols above.

    ``NA`` for every column when the estimate did not run; per-column ``NA`` for
    the numeric cells when there were no usable consensus-hom markers.
    """
    if contamination is None:
        return [_NA] * len(_CONTAMINATION_TSV_COLS)
    if contamination.n_markers == 0:
        return [_NA, _NA, str(contamination.n_markers)]
    return [
        f"{contamination.contamination_fraction:.6f}",
        f"{contamination.p_value:.4g}",
        str(contamination.n_markers),
    ]


# Run-unit / index-hopping columns, appended after the contamination block. The
# metadata is optional (read from the admix VCF header), so these are ``NA`` when
# it is absent. ``index_hop_risk`` is the host-share flag.
_RUNMETA_TSV_COLS = [
    "run_unit",
    "index_hop_risk",
]


def _runmeta_tsv_cells(run_unit: RunUnitInfo | None) -> list[str]:
    """Format the run-unit columns. Order matches ``_RUNMETA_TSV_COLS``.

    ``NA`` throughout when no run metadata was present; ``index_hop_risk`` is
    ``true`` / ``false`` / ``NA`` (undetermined).
    """
    if run_unit is None:
        return [_NA] * len(_RUNMETA_TSV_COLS)
    unit = run_unit.run_unit if run_unit.run_unit else _NA
    shares = run_unit.shares_run_with_host
    risk = _NA if shares is None else ("true" if shares else "false")
    return [unit, risk]


def _runmeta_json(run_unit: RunUnitInfo | None) -> dict | None:
    """JSON view of the run-unit metadata; mirrors the TSV columns plus detail."""
    if run_unit is None:
        return None
    return {
        "run_unit": run_unit.run_unit,
        "source": run_unit.source,
        "shares_run_with_host": run_unit.shares_run_with_host,
    }


def _contamination_json(contamination: ContaminationResult | None) -> dict | None:
    """JSON view of a ContaminationResult; mirrors the TSV columns plus detail."""
    if contamination is None:
        return None
    return {
        "contamination_fraction": contamination.contamination_fraction,
        "median_minor_frac": contamination.median_minor_frac,
        "error_floor": contamination.error_floor,
        "floor_empirical": contamination.floor_empirical,
        "pooled_minor_frac": contamination.pooled_minor_frac,
        "p_value": contamination.p_value,
        "n_markers": contamination.n_markers,
        "n_minor_reads": contamination.n_minor_reads,
        "total_depth": contamination.total_depth,
        "n_excluded_high": contamination.n_excluded_high,
        "used_per_site_error": contamination.used_per_site_error,
        "error_rate_source": contamination.error_rate_source,
    }


def _admix_consistency_json(ac: AdmixConsistencyResult | None) -> dict | None:
    """JSON view of the consensus-homozygote swap check."""
    if ac is None:
        return None
    return {
        "n_consensus_hom": ac.n_consensus_hom,
        "n_discordant": ac.n_discordant,
        "discordant_fraction": round(ac.discordant_fraction, 6),
        "swap_pval": ac.swap_pval,
    }


def _is_multi_donor(result: object) -> bool:
    """Check whether result is a MultiDonorResult."""
    return hasattr(result, "donor_fractions")


def _warnings_cell(qc: QCReport) -> str:
    """Join QC warnings into a single tab-safe TSV cell (empty if none)."""
    return "; ".join(w.replace("\t", " ").replace("\n", " ") for w in qc.warnings)


def to_tsv(
    result: ChimerismResult | MultiDonorResult,
    qc: QCReport,
    output: Path | TextIO,
    verbose: bool = False,
    sample_name: str = "",
) -> None:
    """Write chimerism result and QC report as a TSV file.

    Writes a header and summary line; ``verbose`` appends a per-marker detail
    section after a blank separator. Handles single- and multi-donor results.
    """
    if isinstance(output, (str, Path)):
        with open(output, "w", encoding="utf-8") as fh:
            if _is_multi_donor(result):
                _write_tsv_multi(result, qc, fh, verbose, sample_name)
            else:
                _write_tsv(result, qc, fh, verbose, sample_name)
    else:
        if _is_multi_donor(result):
            _write_tsv_multi(result, qc, output, verbose, sample_name)
        else:
            _write_tsv(result, qc, output, verbose, sample_name)


def _write_tsv(
    result: ChimerismResult,
    qc: QCReport,
    fh: TextIO,
    verbose: bool,
    sample_name: str = "sample",
) -> None:
    """Write single-donor TSV content to an open file handle."""
    # Host-presence and later columns are appended at the end so existing parsers
    # keep finding ``donor_pct`` etc. by header.
    summary_cols = [
        "sample",
        "donor_pct",
        "ci_lo",
        "ci_hi",
        "lob_pct",
        "lod_pct",
        "n_total_markers",
        "n_informative",
        "n_used",
        "mean_depth",
        "gof_pval",
        "qc_status",
        "qc_warnings",
        *_HOST_PRESENCE_TSV_COLS,
        *_RELATEDNESS_TSV_COLS,
        *_CONTAMINATION_TSV_COLS,
        *_RUNMETA_TSV_COLS,
    ]
    fh.write("\t".join(summary_cols) + "\n")

    ci_lo, ci_hi = result.donor_fraction_ci
    gof_str = f"{qc.goodness_of_fit_pval:.4f}" if qc.goodness_of_fit_pval is not None else "NA"
    lob_str = f"{result.lob_fraction * 100:.3f}" if math.isfinite(result.lob_fraction) else "NA"
    lod_str = f"{result.lod_fraction * 100:.3f}" if math.isfinite(result.lod_fraction) else "NA"

    summary_vals = [
        sample_name,
        f"{result.donor_fraction * 100:.2f}",
        f"{ci_lo * 100:.2f}",
        f"{ci_hi * 100:.2f}",
        lob_str,
        lod_str,
        str(qc.n_total_markers),
        str(result.n_informative),
        str(qc.n_used),
        f"{qc.mean_depth:.0f}",
        gof_str,
        qc.status,
        _warnings_cell(qc),
        *_host_presence_tsv_cells(getattr(result, "host_presence", None)),
        *_relatedness_tsv_cells(getattr(result, "relatedness", None)),
        *_contamination_tsv_cells(getattr(result, "contamination", None)),
        *_runmeta_tsv_cells(getattr(result, "run_unit", None)),
    ]
    fh.write("\t".join(summary_vals) + "\n")

    if verbose and qc.warnings:
        fh.write("\n# warnings\n")
        for w in qc.warnings:
            fh.write(f"# {w}\n")

    if verbose and result.per_marker:
        fh.write("\n")
        detail_cols = [
            "chrom",
            "pos",
            "marker_type",
            "host_gt",
            "donor_gt",
            "ad_ref",
            "ad_alt",
            "observed_vaf",
            "expected_vaf",
            "residual",
            "included",
        ]
        fh.write("\t".join(detail_cols) + "\n")
        for m in result.per_marker:
            row = [
                m.chrom,
                str(m.pos),
                str(m.marker_type),
                ".",
                ".",
                str(m.ad_ref),
                str(m.ad_alt),
                f"{m.observed_vaf:.4f}",
                f"{m.expected_vaf:.4f}",
                f"{m.residual:.4f}",
                str(m.included),
            ]
            fh.write("\t".join(row) + "\n")


def _write_tsv_multi(
    result: MultiDonorResult,
    qc: QCReport,
    fh: TextIO,
    verbose: bool,
    sample_name: str = "sample",
) -> None:
    """Write multi-donor TSV content to an open file handle."""
    n_donors = len(result.donor_fractions)

    cols = ["sample"]
    for i in range(n_donors):
        d = i + 1
        cols.extend([f"donor{d}_pct", f"donor{d}_ci_lo", f"donor{d}_ci_hi"])
    cols.extend(
        [
            "host_pct",
            "n_total_markers",
            "n_informative",
            "n_used",
            "mean_depth",
            "gof_pval",
            "qc_status",
            "qc_warnings",
            *_HOST_PRESENCE_TSV_COLS,
            *_RELATEDNESS_TSV_COLS,
            *_CONTAMINATION_TSV_COLS,
            *_RUNMETA_TSV_COLS,
        ]
    )
    fh.write("\t".join(cols) + "\n")

    gof_str = f"{qc.goodness_of_fit_pval:.4f}" if qc.goodness_of_fit_pval is not None else "NA"

    vals = [sample_name]
    for i in range(n_donors):
        ci_lo, ci_hi = result.donor_fraction_cis[i]
        vals.extend(
            [
                f"{result.donor_fractions[i] * 100:.2f}",
                f"{ci_lo * 100:.2f}",
                f"{ci_hi * 100:.2f}",
            ]
        )
    vals.extend(
        [
            f"{result.host_fraction * 100:.2f}",
            str(qc.n_total_markers),
            str(result.n_informative),
            str(qc.n_used),
            f"{qc.mean_depth:.0f}",
            gof_str,
            qc.status,
            _warnings_cell(qc),
            *_host_presence_tsv_cells(getattr(result, "host_presence", None)),
            *_relatedness_tsv_cells(getattr(result, "relatedness", None)),
            *_contamination_tsv_cells(getattr(result, "contamination", None)),
            *_runmeta_tsv_cells(getattr(result, "run_unit", None)),
        ]
    )
    fh.write("\t".join(vals) + "\n")

    if verbose and qc.warnings:
        fh.write("\n# warnings\n")
        for w in qc.warnings:
            fh.write(f"# {w}\n")

    if verbose and result.per_marker:
        fh.write("\n")
        detail_cols = [
            "chrom",
            "pos",
            "marker_type",
            "ad_ref",
            "ad_alt",
            "observed_vaf",
            "expected_vaf",
            "residual",
            "included",
        ]
        fh.write("\t".join(detail_cols) + "\n")
        for m in result.per_marker:
            row = [
                m.chrom,
                str(m.pos),
                str(m.marker_type),
                str(m.ad_ref),
                str(m.ad_alt),
                f"{m.observed_vaf:.4f}",
                f"{m.expected_vaf:.4f}",
                f"{m.residual:.4f}",
                str(m.included),
            ]
            fh.write("\t".join(row) + "\n")


def _markers_json(result: ChimerismResult | MultiDonorResult) -> list[dict]:
    """Per-marker rows for the structured output (the HTML report omits these)."""
    return [
        {
            "chrom": m.chrom,
            "pos": m.pos,
            "marker_type": int(m.marker_type),
            "ad_ref": m.ad_ref,
            "ad_alt": m.ad_alt,
            "dp": m.dp,
            "observed_vaf": round(m.observed_vaf, 6),
            "expected_vaf": round(m.expected_vaf, 6),
            "residual": round(m.residual, 6),
            "included": m.included,
        }
        for m in result.per_marker
    ]


def _qc_common_json(result: ChimerismResult | MultiDonorResult, qc: QCReport) -> dict:
    """QC and sub-analysis fields shared by the single- and multi-donor views.

    The complete per-sample accounting the HTML report renders, so the report can
    be produced from the JSON alone (no result/QC objects needed).
    """
    return {
        "n_total_markers": qc.n_total_markers,
        "n_shared_markers": qc.n_shared_markers,
        "n_informative": result.n_informative,
        "n_used": qc.n_used,
        "per_donor_n_informative": qc.per_donor_n_informative,
        "n_excluded_depth": qc.n_excluded_depth,
        "n_excluded_quality": qc.n_excluded_quality,
        "n_excluded_outlier": qc.n_excluded_outlier,
        "n_robust_excluded": getattr(result, "n_robust_excluded", 0),
        "robust_drop_fraction": getattr(result, "robust_drop_fraction", 0.0),
        "mean_depth": round(qc.mean_depth, 1),
        "median_depth": round(qc.median_depth, 1),
        "min_depth": qc.min_depth,
        "gof_pval": (
            round(qc.goodness_of_fit_pval, 4) if qc.goodness_of_fit_pval is not None else None
        ),
        "gof_pval_pretrim": (
            round(qc.goodness_of_fit_pval_pretrim, 4)
            if qc.goodness_of_fit_pval_pretrim is not None
            else None
        ),
        "qc_pass": qc.pass_,
        "qc_status": qc.status,
        "warnings": list(qc.warnings),
        # host_presence is carried on the result (the QC report has no such field);
        # the other sub-analyses are owned by the QC report.
        "host_presence": _host_presence_json(getattr(result, "host_presence", None)),
        "relatedness": _relatedness_json(qc.relatedness),
        "admix_consistency": _admix_consistency_json(qc.admix_consistency),
        "contamination": _contamination_json(qc.contamination),
        "run_unit": _runmeta_json(qc.run_unit),
    }


def to_json(
    result: ChimerismResult | MultiDonorResult,
    qc: QCReport,
    sample_name: str = "",
) -> dict:
    """Convert a chimerism result and QC report to a JSON-serialisable dict.

    The per-sample payload the HTML report renders from: headline estimate, full
    QC accounting, every sub-analysis (host presence, relatedness, contamination,
    swap check, run unit), and the per-marker rows. Handles single-donor
    (``ChimerismResult``) and multi-donor (``MultiDonorResult``) results.

    Percent-valued headline fields use the ``*_pct`` suffix; nested sub-analyses
    keep their native fractions.
    """
    if _is_multi_donor(result):
        donors = []
        for i, frac in enumerate(result.donor_fractions):
            ci_lo, ci_hi = result.donor_fraction_cis[i]
            d: dict = {
                "label": f"donor{i + 1}",
                "donor_pct": round(frac * 100, 4),
                "ci_lo": round(ci_lo * 100, 4),
                "ci_hi": round(ci_hi * 100, 4),
            }
            if result.per_donor_n_informative:
                d["n_informative"] = result.per_donor_n_informative[i]
            donors.append(d)

        out: dict = {
            "sample": sample_name,
            "host_pct": round(result.host_fraction * 100, 4),
            "donors": donors,
            "total_donor_pct": round(sum(result.donor_fractions) * 100, 4),
        }
    else:
        ci_lo, ci_hi = result.donor_fraction_ci
        out = {
            "sample": sample_name,
            "donor_pct": round(result.donor_fraction * 100, 4),
            "ci_lo": round(ci_lo * 100, 4),
            "ci_hi": round(ci_hi * 100, 4),
            "host_pct": round(result.host_fraction * 100, 4),
            "lob_pct": round(result.lob_fraction * 100, 4)
            if math.isfinite(result.lob_fraction)
            else None,
            "lod_pct": round(result.lod_fraction * 100, 4)
            if math.isfinite(result.lod_fraction)
            else None,
        }

    out.update(_qc_common_json(result, qc))
    out["markers"] = _markers_json(result)
    return out


def report_data(
    result: ChimerismResult | MultiDonorResult,
    qc: QCReport,
    *,
    sample_name: str = "",
    meta: ReportMeta | None = None,
    params: dict | None = None,
    timestamp: str | None = None,
    version: str = __version__,
) -> dict:
    """Build the single-sample report envelope (``kind == "single"``).

    Wraps the per-sample analysis (``to_json``) with report provenance (metadata,
    params, version, generation time). Written by ``allomix detect --json`` and
    consumed by the HTML renderer, so a report can be produced in one step or later
    from the saved JSON.
    """
    return {
        "schema_version": REPORT_SCHEMA_VERSION,
        "kind": "single",
        "allomix_version": version,
        "generated": timestamp,
        "meta": (meta or ReportMeta()).to_dict(),
        "params": params or {},
        "analysis": to_json(result, qc, sample_name=sample_name),
    }


def to_html(
    result: ChimerismResult | MultiDonorResult,
    qc: QCReport,
    output: Path | TextIO,
    *,
    sample_name: str = "",
    meta: ReportMeta | None = None,
    timestamp: str | None = None,
    params: dict | None = None,
    template_dir: str | Path | None = None,
) -> None:
    """Write a single-sample chimerism report as a self-contained HTML file.

    Builds the envelope (``report_data``) and renders it; the HTML layer adds no
    analysis, only formatting. The document inlines all CSS and JavaScript so it
    opens from disk with no network access.

    Args:
        timestamp: Preformatted report generation time. Passed in (rather than
            read from the clock) so the output is deterministic and testable.
        template_dir: Optional directory of template overrides (the ``--template``
            flag), searched ahead of the built-in templates.
    """
    data = report_data(
        result,
        qc,
        sample_name=sample_name,
        meta=meta,
        params=params,
        timestamp=timestamp,
    )
    html = render_single(data, template_dir=template_dir)
    if isinstance(output, (str, Path)):
        with open(output, "w", encoding="utf-8") as fh:
            fh.write(html)
    else:
        output.write(html)


def timeline_json(
    results: list[tuple[str, ChimerismResult | MultiDonorResult, QCReport]],
) -> dict:
    """Build a timeline of chimerism results across multiple timepoints.

    Each timepoint is the full per-sample ``to_json`` dict, so the longitudinal
    report (trend chart, table, latest-timepoint detail) renders from this
    structure alone. ``results`` are in chronological order.
    """
    return {"timepoints": [to_json(result, qc, sample_name=name) for name, result, qc in results]}


def timeline_report_data(
    results: list[tuple[str, ChimerismResult | MultiDonorResult, QCReport]],
    *,
    meta: ReportMeta | None = None,
    params: dict | None = None,
    timestamp: str | None = None,
    version: str = __version__,
) -> dict:
    """Build the timeline report envelope (``kind == "timeline"``).

    The longitudinal counterpart of ``report_data``: wraps every timepoint's
    analysis with report provenance. Written by ``allomix timeline --json`` and
    consumed by the timeline HTML renderer. ``results`` are in chronological order.
    """
    return {
        "schema_version": REPORT_SCHEMA_VERSION,
        "kind": "timeline",
        "allomix_version": version,
        "generated": timestamp,
        "meta": (meta or ReportMeta()).to_dict(),
        "params": params or {},
        "timepoints": timeline_json(results)["timepoints"],
    }
