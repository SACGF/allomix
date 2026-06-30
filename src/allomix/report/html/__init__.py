"""HTML chimerism report rendering (issue #27).

The public entry points live in ``allomix.report.report`` (``to_html`` and
``timeline_html``), mirroring the existing ``to_json`` / ``to_tsv`` surface. This
subpackage holds the rendering machinery:

- ``meta``: optional report metadata dataclasses (``ReportMeta``, ``DonorMeta``).
- ``format``: number / percentage / p-value / CI formatters and the status badge
  helper. Single source of formatting consistency, registered as Jinja
  filters/globals by ``engine``.
- ``engine``: builds the Jinja2 environment and registers the formatting
  filters/globals. A user template directory (the ``--template`` flag) is
  searched ahead of the built-in templates.
- ``context``: turns the structured report envelope into a render-ready context
  dict for the templates (the small presentational derivations live here).
- ``templates/``: the Jinja templates (``base.html``, ``report.html``,
  ``timeline.html``, ``macros.html``) and the inlined ``styles.css`` / ``report.js``.
- ``charts``: matplotlib chart builders that return base64 PNG data URIs.
- ``render``: single-sample page rendering.
- ``timeline``: longitudinal page rendering and its trend chart (matplotlib).

No analysis logic lives here; the layer only formats fields the estimators and
QC already produced.
"""
