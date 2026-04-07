"""Output formatting for chimerism results.

Provides TSV and JSON output for single-sample results and
multi-timepoint timelines.
"""

from __future__ import annotations

from pathlib import Path
from typing import TextIO

try:
    from allomix.chimerism import ChimerismResult
except ImportError:
    from allomix.qc import ChimerismResult  # type: ignore[assignment]

from allomix.qc import QCReport


def to_tsv(
    result: ChimerismResult,
    qc: QCReport,
    output: Path | TextIO,
    verbose: bool = False,
) -> None:
    """Write chimerism result and QC report as a TSV file.

    Writes a header and summary line. When verbose is True, also writes
    per-marker detail lines after a blank separator line.

    Args:
        result: Chimerism estimation result.
        qc: Quality control report.
        output: File path or writable text stream.
        verbose: If True, include per-marker detail section.
    """
    if isinstance(output, (str, Path)):
        with open(output, "w") as fh:
            _write_tsv(result, qc, fh, verbose)
    else:
        _write_tsv(result, qc, output, verbose)


def _write_tsv(
    result: ChimerismResult,
    qc: QCReport,
    fh: TextIO,
    verbose: bool,
) -> None:
    """Write TSV content to an open file handle.

    Args:
        result: Chimerism estimation result.
        qc: Quality control report.
        fh: Open text file handle.
        verbose: If True, include per-marker detail section.
    """
    # Summary header and line
    summary_header = (
        "sample\tdonor_pct\tci_lo\tci_hi\tn_informative\t"
        "n_used\tmean_depth\tgof_pval\tqc_pass"
    )
    fh.write(summary_header + "\n")

    ci_lo, ci_hi = result.donor_fraction_ci
    gof_str = f"{qc.goodness_of_fit_pval:.4f}" if qc.goodness_of_fit_pval is not None else "NA"
    qc_pass_str = "PASS" if qc.pass_ else "FAIL"

    summary_line = (
        f"sample\t"
        f"{result.donor_fraction * 100:.2f}\t"
        f"{ci_lo * 100:.2f}\t"
        f"{ci_hi * 100:.2f}\t"
        f"{result.n_informative}\t"
        f"{qc.n_used}\t"
        f"{qc.mean_depth:.0f}\t"
        f"{gof_str}\t"
        f"{qc_pass_str}"
    )
    fh.write(summary_line + "\n")

    if verbose and result.per_marker:
        fh.write("\n")
        detail_header = (
            "chrom\tpos\tmarker_type\thost_gt\tdonor_gt\t"
            "ad_ref\tad_alt\tobserved_vaf\texpected_vaf\tresidual\tincluded"
        )
        fh.write(detail_header + "\n")
        for m in result.per_marker:
            fh.write(
                f"{m.chrom}\t{m.pos}\t{m.marker_type}\t"
                f".\t.\t"
                f"{m.ad_ref}\t{m.ad_alt}\t"
                f"{m.observed_vaf:.4f}\t{m.expected_vaf:.4f}\t"
                f"{m.residual:.4f}\t{m.included}\n"
            )


def to_json(
    result: ChimerismResult,
    qc: QCReport,
    sample_name: str = "",
) -> dict:
    """Convert chimerism result and QC report to a JSON-serialisable dict.

    Args:
        result: Chimerism estimation result.
        qc: Quality control report.
        sample_name: Sample identifier to include in output.

    Returns:
        Dictionary suitable for json.dumps() with summary and optional
        per-marker data.
    """
    ci_lo, ci_hi = result.donor_fraction_ci
    out: dict = {
        "sample": sample_name,
        "donor_pct": round(result.donor_fraction * 100, 4),
        "ci_lo": round(ci_lo * 100, 4),
        "ci_hi": round(ci_hi * 100, 4),
        "n_informative": result.n_informative,
        "n_used": qc.n_used,
        "mean_depth": round(qc.mean_depth, 1),
        "gof_pval": (
            round(qc.goodness_of_fit_pval, 4)
            if qc.goodness_of_fit_pval is not None
            else None
        ),
        "qc_pass": qc.pass_,
        "warnings": list(qc.warnings),
        "markers": [
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
        ],
    }
    return out


def timeline_json(
    results: list[tuple[str, ChimerismResult, QCReport]],
) -> dict:
    """Build a timeline of chimerism results across multiple timepoints.

    Args:
        results: List of (sample_name, ChimerismResult, QCReport) tuples,
            one per timepoint.

    Returns:
        Dictionary with a 'timepoints' key containing a list of per-sample
        summary dicts.
    """
    timepoints = []
    for sample_name, result, qc in results:
        ci_lo, ci_hi = result.donor_fraction_ci
        timepoints.append({
            "sample": sample_name,
            "donor_pct": round(result.donor_fraction * 100, 4),
            "ci_lo": round(ci_lo * 100, 4),
            "ci_hi": round(ci_hi * 100, 4),
            "n_informative": result.n_informative,
            "n_used": qc.n_used,
            "mean_depth": round(qc.mean_depth, 1),
            "gof_pval": (
                round(qc.goodness_of_fit_pval, 4)
                if qc.goodness_of_fit_pval is not None
                else None
            ),
            "qc_pass": qc.pass_,
        })
    return {"timepoints": timepoints}
