"""Render context for the HTML report.

Turns the structured report envelope (``allomix.report.report_data`` /
``timeline_report_data``) into a flat, render-ready context dict for the Jinja
templates. This is where the small amount of presentational derivation lives
(days post-transplant, the host CI reflected about 100%, the methods-footer
parameter summary, the clinician marker-accounting sentence) so the templates
stay declarative: they place values, they do not compute them.

No analysis happens here either. Every number shown still comes from the
estimators and QC via the envelope; this layer only arranges and labels it. The
context is built from the envelope alone, so a report renders identically whether
produced in one step or later from a saved JSON file.
"""

import datetime
from pathlib import Path

from allomix.html.meta import ReportMeta
from allomix.qc import (
    CONTAMINATION_REVIEW_FRACTION,
    GOF_REVIEW_P,
    HOST_PRESENCE_REVIEW_P,
    LOW_MEAN_DEPTH_WARN,
    ROBUST_REVIEW_FRACTION,
    SWAP_REVIEW_P,
)

# Method citation shown in the footer (CLAUDE.md: cite Crysup & Woerner 2022).
CITATION = (
    "Maximum-likelihood mixture estimation following Crysup and Woerner (2022), "
    "Forensic Science International: Genetics."
)

# Reference-column threshold strings, formatted once for the QC tables. Plain
# text (">" / "<"); the templates autoescape them.
REVIEW_REFS = {
    "robust": f"review if > {ROBUST_REVIEW_FRACTION * 100:.0f}%",
    "depth": f"warn if < {LOW_MEAN_DEPTH_WARN}x",
    "gof": f"review if < {GOF_REVIEW_P}",
    "contamination": f"review if > {CONTAMINATION_REVIEW_FRACTION * 100:.0f}%",
    "swap": f"review if < {SWAP_REVIEW_P:g}",
}

# P-value below which the host-presence callout is treated as a detection.
HOST_PRESENCE_DETECT_P = HOST_PRESENCE_REVIEW_P


def _days_post_transplant(transplant_date: str | None, sample_date: str | None) -> int | None:
    """Days between transplant and sample collection, or None.

    Both dates are preformatted strings; this parses them as ISO dates only to
    derive the interval and returns None if either is absent or not ISO-parseable
    (the report never blocks on a date it cannot read).
    """
    if not transplant_date or not sample_date:
        return None
    try:
        t = datetime.date.fromisoformat(transplant_date)
        s = datetime.date.fromisoformat(sample_date)
    except ValueError:
        return None
    return (s - t).days


def header_rows(
    meta: ReportMeta, sample_name: str, *, version: str, timestamp: str | None
) -> list[tuple[str, str]]:
    """Build the (label, value) rows for the header identification band.

    Rows are included only when their value is present, so the band never shows
    empty fields. Days-post-transplant is derived only when both a transplant
    date and this sample's collection date are present.
    """
    rows: list[tuple[str, object | None]] = [
        ("Recipient ID", meta.recipient_id),
        ("Name", meta.recipient_name),
        ("Sex", meta.sex),
        ("Date of birth", meta.dob),
        ("Sample", sample_name),
        ("Transplant type", meta.transplant_type),
        ("Transplant date", meta.transplant_date),
    ]

    sample_date = meta.sample_dates.get(sample_name)
    rows.append(("Sample collected", sample_date))
    days = _days_post_transplant(meta.transplant_date, sample_date)
    if days is not None:
        rows.append(("Days post-transplant", f"+{days}"))

    for i, donor in enumerate(meta.donors, start=1):
        ident = donor.donor_id or f"donor {i}"
        rel = f" ({donor.relationship})" if donor.relationship else ""
        rows.append((f"Donor {i}", f"{ident}{rel}"))

    rows.append(("allomix version", version))
    rows.append(("Report generated", timestamp))

    return [(label, str(value)) for label, value in rows if value is not None and str(value) != ""]


def _marker_accounting(a: dict) -> str:
    """One-sentence informative-marker accounting line for the clinician."""
    excl: list[str] = []
    if a.get("n_excluded_outlier"):
        excl.append(f"{a['n_excluded_outlier']} outliers")
    if a.get("n_excluded_depth"):
        excl.append(f"{a['n_excluded_depth']} low-depth")
    if a.get("n_excluded_quality"):
        excl.append(f"{a['n_excluded_quality']} low-quality")
    excl_str = f" (excluded: {', '.join(excl)})" if excl else ""
    n_inf = a.get("n_informative")
    n_total = a.get("n_total_markers")
    n_used = a.get("n_used")
    return (
        f"{n_inf} of {n_total} input markers were informative; {n_used} used in the fit{excl_str}."
    )


def _flag_summary(a: dict) -> str:
    """One-line count of flagged QC checks for the clinician."""
    n = len(a.get("warnings") or [])
    if n:
        return f"{n} QC check{'s' if n != 1 else ''} flagged (listed below)."
    return "No QC checks flagged."


def _params_view(params: dict) -> dict:
    """Methods-footer view of the analysis parameters.

    Collapses the raw CLI parameter dict into the labelled lines the footer
    shows: input file basenames only (no patient-identifying paths), the
    error-rate and bias-correction descriptions, and the on/off toggles.
    """
    file_keys = [
        ("Panel VCF", "panel_vcf"),
        ("Admixture VCF", "admix_vcf"),
        ("Error table", "error_table"),
        ("Bias table", "bias_table"),
        ("Contamination table", "contamination_table"),
    ]
    input_files = [(label, Path(params[key]).name) for label, key in file_keys if params.get(key)]

    if params.get("no_error_correction"):
        error_model = "disabled"
    elif params.get("error_table"):
        error_model = f"table ({Path(params['error_table']).name})"
    else:
        error_model = f"flat rate {params.get('error_rate')}"

    if params.get("no_bias_correction"):
        bias = "disabled"
    elif params.get("bias_table"):
        bias = f"table ({Path(params['bias_table']).name})"
    elif params.get("estimate_bias"):
        bias = "estimated inline"
    else:
        bias = "none"

    def on_off(key: str) -> str:
        return "on" if params.get(key) else "off"

    lines = [
        ("Minimum depth", str(params.get("min_dp")) if params.get("min_dp") is not None else NA_),
        ("Minimum GQ", str(params.get("min_gq")) if params.get("min_gq") is not None else NA_),
        ("Error-rate model", error_model),
        ("Bias correction", bias),
        ("Robust refit", f"{params.get('robust', NA_)} (k={params.get('robust_k', NA_)})"),
        (
            "Overdispersion",
            "per marker type" if params.get("marker_type_overdispersion") else "shared (legacy)",
        ),
        ("Host-presence detection", on_off("host_presence")),
        ("Artifact filter", on_off("artifact_filter")),
        ("Contamination correction", on_off("contamination_correction")),
        (
            "Sex chromosomes",
            "included" if params.get("use_sex_chroms") else "excluded (autosomes only)",
        ),
    ]
    return {"input_files": input_files, "param_lines": lines}


# Local NA sentinel for footer parameter text (matches ``format.NA``); kept here
# so this module does not depend on the render-time environment globals.
NA_ = "—"


def base_context(data: dict) -> dict:
    """Provenance / footer context shared by the single and timeline reports."""
    params = data.get("params") or {}
    return {
        "version": data.get("allomix_version", ""),
        "timestamp": data.get("generated"),
        "params": _params_view(params),
        "citation": CITATION,
        "refs": REVIEW_REFS,
    }


def qc_context(analysis: dict) -> dict:
    """QC-panel context: the analysis dict plus the clinician summary lines."""
    return {
        "analysis": analysis,
        "marker_accounting": _marker_accounting(analysis),
        "flag_summary": _flag_summary(analysis),
        "refs": REVIEW_REFS,
    }


def host_presence_context(analysis: dict) -> dict:
    """Host-presence callout context: the detector sub-object and the verdict."""
    hp = analysis.get("host_presence")
    detected = (
        hp is not None and hp.get("n_markers", 0) > 0 and hp["lrt_pval"] < HOST_PRESENCE_DETECT_P
    )
    return {"hp": hp, "detected": detected}


def single_context(data: dict) -> dict:
    """Full render context for the single-sample report."""
    analysis = data["analysis"]
    meta = ReportMeta.from_dict(data.get("meta"))
    version = data.get("allomix_version", "")
    timestamp = data.get("generated")
    sample_name = analysis.get("sample", "")

    # Host CI for a single donor is the donor CI reflected about 100%
    # (host = 100 - donor); purely presentational, no new estimate.
    host_ci_lo = host_ci_hi = None
    if "donors" not in analysis:
        ci_lo, ci_hi = analysis.get("ci_lo"), analysis.get("ci_hi")
        if ci_lo is not None and ci_hi is not None:
            host_ci_lo, host_ci_hi = 100.0 - ci_hi, 100.0 - ci_lo

    ctx = {
        "analysis": analysis,
        "header_rows": header_rows(meta, sample_name, version=version, timestamp=timestamp),
        "host_ci_lo": host_ci_lo,
        "host_ci_hi": host_ci_hi,
        "qc": qc_context(analysis),
        "host_presence": host_presence_context(analysis),
    }
    ctx.update(base_context(data))
    ctx["title"] = _title("allomix chimerism report", sample_name, meta.recipient_id)
    return ctx


def _title(base: str, sample_name: str, recipient_id: str | None) -> str:
    """Document title: ``base: <label>`` when a label is available."""
    label = sample_name or (recipient_id or "")
    return f"{base}: {label}" if label else base
