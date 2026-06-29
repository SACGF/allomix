"""Single-sample report rendering.

Fills the ``report.html`` Jinja template from the structured report envelope
produced by ``allomix.report.report_data`` (the same JSON written by
``monitor --json``), not the in-memory result/QC objects. Rendering therefore
adds no analysis: every value shown comes straight from the envelope, so the HTML
and the JSON always agree and a report can be regenerated later from a saved JSON
file.

This module is matplotlib-free so it stays importable on the base runtime; only
the timeline path (``allomix.html.timeline``) pulls in the charting extra.
"""

from pathlib import Path

from allomix.html import context
from allomix.html.engine import load_asset, make_environment


def render_single(data: dict, *, template_dir: str | Path | None = None) -> str:
    """Render a single-sample report envelope to a complete HTML string.

    Args:
        data: A single-sample report envelope (``kind == "single"``) from
            ``allomix.report.report_data``: the per-sample ``analysis`` dict plus
            the ``meta`` / ``params`` / ``allomix_version`` / ``generated``
            provenance.
        template_dir: Optional directory of template overrides (the ``--template``
            flag), searched ahead of the built-in templates. A lab can drop in a
            ``styles.css`` to restyle or a ``report.html`` to restructure.

    Returns:
        A complete, self-contained HTML document.
    """
    env = make_environment(template_dir)
    ctx = context.single_context(data)
    template = env.get_template("report.html")
    return template.render(
        css=load_asset(env, "styles.css"),
        js=load_asset(env, "report.js"),
        **ctx,
    )
