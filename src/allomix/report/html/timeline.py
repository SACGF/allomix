"""Longitudinal (timeline) HTML report.

This module owns the timeline report and its trend chart. It imports
``allomix.report.html.charts`` (matplotlib) at the top level, so it deliberately sits
outside the import graph of the always-loaded ``report`` / ``cli`` modules: the
single-sample report and the tsv/json paths must keep working on the base
runtime deps alone. The CLI imports this module only after confirming matplotlib
is installed (the ``report`` optional extra), the same way the standalone
plotting scripts keep matplotlib out of the core package.

Like the single-sample renderer, this consumes the structured timeline envelope
(``allomix.report.report.timeline_report_data``), not the in-memory result/QC objects,
so a timeline report can be regenerated later from a saved JSON file.

The public entry point is ``timeline_html``, mirroring ``report.to_html`` for the
single-sample case.
"""

from pathlib import Path
from typing import TextIO

from allomix.qc.qc import QCReport
from allomix.report.html import charts, context
from allomix.report.html.engine import load_asset, make_environment
from allomix.report.html.meta import ReportMeta
from allomix.report.report import timeline_report_data
from allomix.results import ChimerismResult, MultiDonorResult

_TITLE_BASE = "allomix chimerism timeline"


def _tp_donor_pct(tp: dict) -> float | None:
    """Donor percentage for a timepoint (combined across donors if multi-donor)."""
    if "donors" in tp:
        return sum(d["donor_pct"] for d in tp["donors"])
    return tp.get("donor_pct")


def _tp_donor_ci(tp: dict) -> tuple[float | None, float | None]:
    """Donor CI bounds (percent) for a single-donor timepoint, else (None, None)."""
    if "donors" in tp:
        return (None, None)  # no combined CI is defined for multi-donor
    return (tp.get("ci_lo"), tp.get("ci_hi"))


def _headline(timepoints: list[dict]) -> dict:
    """Latest-value / previous-value / delta context for the timeline headline."""
    latest = timepoints[-1]
    prev = timepoints[-2] if len(timepoints) >= 2 else None
    cur = _tp_donor_pct(latest)
    ci_lo, ci_hi = _tp_donor_ci(latest)

    h = {
        "status": latest.get("qc_status", ""),
        "latest_pct": cur,
        "latest_ci_lo": ci_lo,
        "latest_ci_hi": ci_hi,
        "has_prev": prev is not None,
        "prev_pct": None,
        "delta": None,
    }
    if prev is not None:
        prev_val = _tp_donor_pct(prev)
        h["prev_pct"] = prev_val
        if cur is not None and prev_val is not None and cur == cur:  # cur finite
            delta = cur - prev_val
            arrow = "▲" if delta > 0 else ("▼" if delta < 0 else "→")
            h["delta"] = f"{arrow} {abs(delta):.2f} pp"
    return h


def _table_rows(timepoints: list[dict], sample_dates: dict[str, str]) -> list[dict]:
    """Per-timepoint summary-table rows (label + the values the columns show)."""
    rows = []
    for tp in timepoints:
        sample = tp.get("sample", "")
        ci_lo, ci_hi = _tp_donor_ci(tp)
        rows.append(
            {
                "label": sample_dates.get(sample, sample),
                "donor_pct": _tp_donor_pct(tp),
                "ci_lo": ci_lo,
                "ci_hi": ci_hi,
                "depth": tp.get("mean_depth"),
                "gof": tp.get("gof_pval"),
                "status": tp.get("qc_status", ""),
            }
        )
    return rows


def _build_chart(
    timepoints: list[dict], sample_dates: dict[str, str], log_scale: bool
) -> tuple[str | None, str]:
    """Build the trend-chart data URI, or a note explaining why there is none.

    Returns:
        ``(chart_uri, note)``. ``chart_uri`` is None and ``note`` is set when a
        chart cannot be drawn (fewer than two timepoints), so the section is
        never silently blank.
    """
    if len(timepoints) < 2:
        return None, "A trend chart needs at least two timepoints; one was supplied."

    x_labels = [sample_dates.get(tp.get("sample", ""), tp.get("sample", "")) for tp in timepoints]

    if "donors" in timepoints[0]:
        n_donors = len(timepoints[0]["donors"])
        series = [
            {
                "name": f"Donor {i + 1}",
                "y": [tp["donors"][i]["donor_pct"] for tp in timepoints],
                "ci_lo": [tp["donors"][i]["ci_lo"] for tp in timepoints],
                "ci_hi": [tp["donors"][i]["ci_hi"] for tp in timepoints],
                "qc_pass": [tp.get("qc_pass", True) for tp in timepoints],
            }
            for i in range(n_donors)
        ]
    else:
        series = [
            {
                "name": "Donor",
                "y": [tp["donor_pct"] for tp in timepoints],
                "ci_lo": [tp["ci_lo"] for tp in timepoints],
                "ci_hi": [tp["ci_hi"] for tp in timepoints],
                "qc_pass": [tp.get("qc_pass", True) for tp in timepoints],
            }
        ]

    return charts.trend_png(x_labels, series, log_scale=log_scale), ""


def _timeline_context(data: dict, log_scale: bool) -> dict:
    """Build the full render context for the timeline report."""
    timepoints = data["timepoints"]
    meta = ReportMeta.from_dict(data.get("meta"))
    version = data.get("allomix_version", "")
    timestamp = data.get("generated")
    latest = timepoints[-1]
    latest_name = latest.get("sample", "")

    chart_uri, note = _build_chart(timepoints, meta.sample_dates, log_scale)

    label = meta.recipient_id or latest_name or ""
    title = f"{_TITLE_BASE}: {label}" if label else _TITLE_BASE

    ctx = {
        "title": title,
        "header_rows": context.header_rows(meta, latest_name, version=version, timestamp=timestamp),
        "headline": _headline(timepoints),
        "chart_uri": chart_uri,
        "chart_note": note,
        "table_rows": _table_rows(timepoints, meta.sample_dates),
        "host_presence": context.host_presence_context(latest),
        "qc": context.qc_context(latest),
    }
    ctx.update(context.base_context(data))
    return ctx


def render_timeline(
    data: dict, *, log_scale: bool = False, template_dir: str | Path | None = None
) -> str:
    """Render a longitudinal report envelope to a complete HTML string.

    The headline is the most recent timepoint (value, previous value, and
    change), followed by the trend chart and a per-timepoint table. The full
    host-presence callout and QC panel are shown for the most recent sample,
    reusing the single-sample macros. Everything is read from the envelope, so no
    result/QC objects are needed.

    Args:
        data: A timeline report envelope (``kind == "timeline"``) from
            ``allomix.report.report.timeline_report_data``.
        log_scale: Use a logarithmic y axis on the trend chart.
        template_dir: Optional directory of template overrides (the ``--template``
            flag), searched ahead of the built-in templates.

    Returns:
        A complete, self-contained HTML document.

    Raises:
        ValueError: If the envelope carries no timepoints.
    """
    if not data.get("timepoints"):
        raise ValueError("render_timeline requires at least one timepoint")

    env = make_environment(template_dir)
    ctx = _timeline_context(data, log_scale)
    template = env.get_template("timeline.html")
    return template.render(
        css=load_asset(env, "styles.css"),
        js=load_asset(env, "report.js"),
        **ctx,
    )


def timeline_html(
    results: list[tuple[str, ChimerismResult | MultiDonorResult, QCReport]],
    output: Path | TextIO,
    *,
    meta: ReportMeta | None = None,
    timestamp: str | None = None,
    params: dict | None = None,
    log_scale: bool = False,
    template_dir: str | Path | None = None,
) -> None:
    """Write a longitudinal chimerism report as a self-contained HTML file.

    Builds the timeline envelope (``timeline_report_data``) and renders it. The
    most recent timepoint supplies the headline and the detailed QC / marker
    sections; all timepoints feed the trend chart and table.

    Args:
        results: ``(sample_name, result, QCReport)`` per timepoint, in
            chronological order. Must be non-empty.
        output: Output file path or a writable text stream.
        meta: Optional recipient / transplant metadata for the header band.
        timestamp: Preformatted generation time. Passed in (rather than read from
            the clock) so the output is deterministic and testable.
        params: Analysis parameters for the methods footer.
        log_scale: Use a logarithmic y axis on the trend chart.
        template_dir: Optional directory of template overrides (``--template``).

    Raises:
        ValueError: If ``results`` is empty.
    """
    if not results:
        raise ValueError("timeline_html requires at least one timepoint")
    data = timeline_report_data(results, meta=meta, params=params, timestamp=timestamp)
    html = render_timeline(data, log_scale=log_scale, template_dir=template_dir)
    if isinstance(output, (str, Path)):
        with open(output, "w", encoding="utf-8") as fh:
            fh.write(html)
    else:
        output.write(html)
