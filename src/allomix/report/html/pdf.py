"""PDF rendering of the chimerism report (issue #35).

Renders the same report envelope as the HTML path to PDF via WeasyPrint
(pure-Python, no headless browser). The layout reuses the shared section macros
(``macros.html``) through a PDF-specific template set (``pdf/report.html``,
``pdf/timeline.html``, ``pdf/styles.css``): the on-screen "how it works"
disclosures move to a "Methods and notes" appendix so nothing is hidden in a
printed document, and ``pdf/styles.css`` adds the paged-media running header and
page numbers on top of the shared look.

This is a demo-grade converter, not a LIMS replacement. A lab that wants a
branded or signed report overrides the ``pdf/`` template set with ``--template``,
exactly as for HTML.

WeasyPrint is an optional dependency (the ``[pdf]`` extra, which also pulls in
the timeline chart's matplotlib via ``[report]``). This module imports both at
the top level; the CLI imports this module only after confirming WeasyPrint is
installed (a ``find_spec`` check), the same deferred-after-capability-check
pattern the timeline HTML path uses for matplotlib, so the base install stays
light and never imports WeasyPrint.
"""

from pathlib import Path
from typing import BinaryIO

from markupsafe import Markup
from weasyprint import HTML

from allomix.report.html import context
from allomix.report.html.engine import load_asset, make_environment
from allomix.report.html.timeline import _timeline_context


def _render_pdf(
    template_name: str,
    ctx: dict,
    output: str | Path | BinaryIO,
    *,
    template_dir: str | Path | None,
) -> None:
    """Render a PDF template with the shared + PDF stylesheets and write it out.

    The CSS is the base ``styles.css`` and the additive ``pdf/styles.css``
    concatenated (base first), so the PDF cannot drift from the shared look. No
    JavaScript is inlined (``js=""``): the PDF is static.
    """
    env = make_environment(template_dir)
    css = Markup(f"{load_asset(env, 'styles.css')}\n{load_asset(env, 'pdf/styles.css')}")
    template = env.get_template(template_name)
    html = template.render(css=css, js="", **ctx)
    target = str(output) if isinstance(output, (str, Path)) else output
    HTML(string=html).write_pdf(target)


def to_pdf_single(
    data: dict, output: str | Path | BinaryIO, *, template_dir: str | Path | None = None
) -> None:
    """Render a single-sample report envelope to PDF.

    Args:
        data: A single-sample report envelope (``kind == "single"``) from
            ``allomix.report.report.report_data``.
        output: Output file path or a writable binary stream.
        template_dir: Optional directory of template overrides (``--template``),
            searched ahead of the built-in ``pdf/`` templates.
    """
    _render_pdf("pdf/report.html", context.single_context(data), output, template_dir=template_dir)


def to_pdf_timeline(
    data: dict,
    output: str | Path | BinaryIO,
    *,
    log_scale: bool = False,
    template_dir: str | Path | None = None,
) -> None:
    """Render a timeline report envelope to PDF.

    Args:
        data: A timeline report envelope (``kind == "timeline"``) from
            ``allomix.report.report.timeline_report_data``.
        output: Output file path or a writable binary stream.
        log_scale: Use a logarithmic y axis on the trend chart.
        template_dir: Optional directory of template overrides (``--template``),
            searched ahead of the built-in ``pdf/`` templates.

    Raises:
        ValueError: If the envelope carries no timepoints.
    """
    if not data.get("timepoints"):
        raise ValueError("to_pdf_timeline requires at least one timepoint")
    ctx = _timeline_context(data, log_scale)
    _render_pdf("pdf/timeline.html", ctx, output, template_dir=template_dir)
