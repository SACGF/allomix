"""Output formatting for chimerism results.

Provides TSV and JSON output for single-sample results and
multi-timepoint timelines. Supports both single-donor and multi-donor results.
"""

import math
from pathlib import Path
from typing import TextIO

from allomix.contamination import ContaminationResult
from allomix.detect import HostPresenceResult
from allomix.qc import QCReport
from allomix.relatedness import AdmixConsistencyResult, RelatednessResult
from allomix.results import ChimerismResult, MultiDonorResult
from allomix.runmeta import RunUnitInfo

# Sentinel for the host-presence cells when the detector did not run (e.g.
# --no-host-presence) or produced no usable markers. Keeps the TSV
# rectangular and parseable by downstream scripts.
_NA = "NA"


def _host_presence_tsv_cells(hp: HostPresenceResult | None) -> list[str]:
    """Format the six host-presence columns for a single result.

    Order matches ``_HOST_PRESENCE_TSV_COLS`` below. Returns ``NA`` for every
    column when the detector did not run, and per-column ``NA`` when there
    are no usable markers (the source flag still carries information).
    """
    if hp is None:
        return [_NA] * 7
    if hp.n_markers == 0:
        # Carry the source flag (typically "none") but leave numeric cells
        # empty since the test was not actually exercised. The artifact count
        # is still informative (markers may have been filtered down to zero).
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


# Names APPENDED to the existing TSV summary. Order matters: downstream
# scripts parse by header so renaming or reordering existing columns would
# break them. New columns sit at the end.
_HOST_PRESENCE_TSV_COLS = [
    "host_present_p",
    "host_f_est",
    "host_f_ci_lo",
    "host_f_ci_hi",
    "host_detect_markers",
    "host_err_source",
    "host_artifact_filtered",
]

# Relatedness columns, appended after the host-presence block. For multi-donor
# samples the host-vs-donor pairs are joined with ";" within each cell, keeping
# the row rectangular. The pass/fail verdict itself lives in qc_status /
# qc_warnings, not a separate column.
_RELATEDNESS_TSV_COLS = [
    "relatedness",
    "relatedness_ci",
    "relatedness_confidence",
    "relationship",
]


def _relatedness_tsv_cells(relatedness: list[RelatednessResult] | None) -> list[str]:
    """Format the four relatedness columns from the host-vs-donor pairs.

    Order matches ``_RELATEDNESS_TSV_COLS``. Reports the host-vs-donor pairs (in
    donor order); multiple donors are joined with ";". Returns ``NA`` cells when
    no relatedness was computed.
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
    """JSON view of all relatedness pairs; mirrors the TSV columns."""
    if not relatedness:
        return None
    out = []
    for r in relatedness:
        out.append(
            {
                "pair": r.pair,
                "coefficient": round(r.coefficient, 6) if r.coefficient is not None else None,
                "ci": [round(r.ci_low, 6), round(r.ci_high, 6)] if r.ci_low is not None else None,
                "confidence": r.confidence,
                "relationship": r.relationship,
                "n_sites": r.n_sites,
            }
        )
    return out


# Contamination columns, appended after the relatedness block. ``contamination_frac``
# is the estimated third-party floor (fraction above sequencing error), the
# in-data half of issue #12.
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

    Writes a header and summary line. When verbose is True, also writes
    per-marker detail lines after a blank separator line.
    Handles both single-donor and multi-donor results.

    Args:
        result: Chimerism estimation result (single or multi-donor).
        qc: Quality control report.
        output: File path or writable text stream.
        verbose: If True, include per-marker detail section.
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
    """Write TSV content to an open file handle.

    Args:
        result: Chimerism estimation result.
        qc: Quality control report.
        fh: Open text file handle.
        verbose: If True, include per-marker detail section.
    """
    # Summary header and line. Host-presence columns are appended at the end
    # so existing parsers continue to find ``donor_pct`` etc. by header.
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
    qc_status_str = qc.status
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
        qc_status_str,
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

    # Build header dynamically
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

    # Build data line
    gof_str = f"{qc.goodness_of_fit_pval:.4f}" if qc.goodness_of_fit_pval is not None else "NA"
    qc_status_str = qc.status

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
            qc_status_str,
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


def to_json(
    result: ChimerismResult | MultiDonorResult,
    qc: QCReport,
    sample_name: str = "",
) -> dict:
    """Convert chimerism result and QC report to a JSON-serialisable dict.

    Handles both single-donor (ChimerismResult) and multi-donor
    (MultiDonorResult) results.

    Args:
        result: Chimerism estimation result.
        qc: Quality control report.
        sample_name: Sample identifier to include in output.

    Returns:
        Dictionary suitable for json.dumps() with summary and optional
        per-marker data.
    """
    markers_list = [
        {
            "chrom": m.chrom,
            "pos": m.pos,
            "marker_type": m.marker_type,
            "ad_ref": m.ad_ref,
            "ad_alt": m.ad_alt,
            "observed_vaf": round(m.observed_vaf, 6),
            "expected_vaf": round(m.expected_vaf, 6),
            "residual": round(m.residual, 6),
            "included": m.included,
        }
        for m in result.per_marker
    ]

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
            "n_total_markers": qc.n_total_markers,
            "n_informative": result.n_informative,
            "n_used": qc.n_used,
            "n_robust_excluded": getattr(result, "n_robust_excluded", 0),
            "mean_depth": round(qc.mean_depth, 1),
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
            "host_presence": _host_presence_json(getattr(result, "host_presence", None)),
            "relatedness": _relatedness_json(getattr(result, "relatedness", None)),
            "admix_consistency": _admix_consistency_json(
                getattr(result, "admix_consistency", None)
            ),
            "contamination": _contamination_json(getattr(result, "contamination", None)),
            "run_unit": _runmeta_json(getattr(result, "run_unit", None)),
            "markers": markers_list,
        }
    else:
        ci_lo, ci_hi = result.donor_fraction_ci
        out = {
            "sample": sample_name,
            "donor_pct": round(result.donor_fraction * 100, 4),
            "ci_lo": round(ci_lo * 100, 4),
            "ci_hi": round(ci_hi * 100, 4),
            "lob_pct": round(result.lob_fraction * 100, 4)
            if math.isfinite(result.lob_fraction)
            else None,
            "lod_pct": round(result.lod_fraction * 100, 4)
            if math.isfinite(result.lod_fraction)
            else None,
            "n_total_markers": qc.n_total_markers,
            "n_informative": result.n_informative,
            "n_used": qc.n_used,
            "n_robust_excluded": getattr(result, "n_robust_excluded", 0),
            "mean_depth": round(qc.mean_depth, 1),
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
            "host_presence": _host_presence_json(getattr(result, "host_presence", None)),
            "relatedness": _relatedness_json(getattr(result, "relatedness", None)),
            "admix_consistency": _admix_consistency_json(
                getattr(result, "admix_consistency", None)
            ),
            "contamination": _contamination_json(getattr(result, "contamination", None)),
            "run_unit": _runmeta_json(getattr(result, "run_unit", None)),
            "markers": markers_list,
        }
    return out


def timeline_json(
    results: list[tuple[str, ChimerismResult | MultiDonorResult, QCReport]],
) -> dict:
    """Build a timeline of chimerism results across multiple timepoints.

    Handles both single-donor and multi-donor results.

    Args:
        results: List of (sample_name, result, QCReport) tuples,
            one per timepoint.

    Returns:
        Dictionary with a 'timepoints' key containing a list of per-sample
        summary dicts.
    """
    timepoints = []
    for sample_name, result, qc in results:
        if _is_multi_donor(result):
            tp: dict = {
                "sample": sample_name,
                "host_pct": round(result.host_fraction * 100, 4),
                "donors": [],
                "n_total_markers": qc.n_total_markers,
                "n_informative": result.n_informative,
                "n_used": qc.n_used,
                "mean_depth": round(qc.mean_depth, 1),
                "gof_pval": (
                    round(qc.goodness_of_fit_pval, 4)
                    if qc.goodness_of_fit_pval is not None
                    else None
                ),
                "gof_pval_pretrim": (
                    round(qc.goodness_of_fit_pval_pretrim, 4)
                    if qc.goodness_of_fit_pval_pretrim is not None
                    else None
                ),
                "qc_pass": qc.pass_,
                "qc_status": qc.status,
                "host_presence": _host_presence_json(getattr(result, "host_presence", None)),
                "contamination": _contamination_json(getattr(result, "contamination", None)),
                "run_unit": _runmeta_json(getattr(result, "run_unit", None)),
            }
            for i, frac in enumerate(result.donor_fractions):
                ci_lo, ci_hi = result.donor_fraction_cis[i]
                tp["donors"].append(
                    {
                        "label": f"donor{i + 1}",
                        "donor_pct": round(frac * 100, 4),
                        "ci_lo": round(ci_lo * 100, 4),
                        "ci_hi": round(ci_hi * 100, 4),
                    }
                )
            timepoints.append(tp)
        else:
            ci_lo, ci_hi = result.donor_fraction_ci
            timepoints.append(
                {
                    "sample": sample_name,
                    "donor_pct": round(result.donor_fraction * 100, 4),
                    "ci_lo": round(ci_lo * 100, 4),
                    "ci_hi": round(ci_hi * 100, 4),
                    "n_total_markers": qc.n_total_markers,
                    "n_informative": result.n_informative,
                    "n_used": qc.n_used,
                    "mean_depth": round(qc.mean_depth, 1),
                    "gof_pval": (
                        round(qc.goodness_of_fit_pval, 4)
                        if qc.goodness_of_fit_pval is not None
                        else None
                    ),
                    "qc_pass": qc.pass_,
                    "qc_status": qc.status,
                    "host_presence": _host_presence_json(getattr(result, "host_presence", None)),
                    "contamination": _contamination_json(getattr(result, "contamination", None)),
                    "run_unit": _runmeta_json(getattr(result, "run_unit", None)),
                }
            )
    return {"timepoints": timepoints}
