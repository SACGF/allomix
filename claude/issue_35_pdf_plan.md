# Issue #35: PDF report output (WeasyPrint, demo-grade)

Plan for adding a `--pdf` output option to `detect`, `timeline`, and `report`,
rendering the existing HTML report to PDF via WeasyPrint (pure-Python, no headless
browser). Demo-grade converter, not a LIMS replacement; labs that want a fully
branded/signed report use `--template` overrides.

## What's already in place (nothing to build)

- HTML renders from the canonical JSON envelope via Jinja2, fully self-contained
  (CSS/JS inlined, timeline chart is a base64 PNG). `render_single` /
  `render_timeline` are the render entry points.
- `styles.css` already has an `@media print` block and
  `@page { size: A4; margin: 16mm }` (`styles.css:194-200`).
- The `--template` override is a `ChoiceLoader` (user `FileSystemLoader` first,
  then `PackageLoader` at `report/html/templates`), resolving by filename
  (`engine.py:42-45`).
- The matplotlib optional-extra pattern (guard with `find_spec` in the CLI, then
  deferred import of a module that imports the heavy dep at top level) is
  established in `cmd_timeline`/`cmd_report` and `timeline.py`. WeasyPrint follows
  the same pattern, so it stays within CLAUDE.md's "no lazy imports / no
  try-except import guards" rule (top-level import lives in a new module; the CLI
  does the `find_spec` check then imports it, with an explanatory comment like the
  matplotlib one).

## Design decisions

**Template namespacing -> `templates/pdf/` subdir.** New `pdf/report.html`,
`pdf/timeline.html`, `pdf/styles.css`. `macros.html` stays shared and
single-sourced (imported as `"macros.html"`; Jinja import paths are
loader-root-relative, so `pdf/report.html` imports it unchanged). A lab overrides
the PDF set by dropping `pdf/report.html` etc. into their `--template` dir, exactly
like HTML. Cleaner than distinct filenames (`report_pdf.html`) and gives labs a
directory to copy as the worked example.

**Separate `[pdf]` extra** (decided with the user). `weasyprint` goes in a new
`[pdf]` extra, added to `[dev]` so tests run. Keeps the core install light; PDF is
opt-in like the `report` extra.

**PDF CSS is additive, not standalone.** `render_pdf` inlines `styles.css` +
`pdf/styles.css` concatenated (base first), so `pdf/styles.css` holds only the
print-specific additions and can't drift from the base look. `base.html` unchanged;
the PDF shells extend it and pass `js=""` (no table-sort JS needed in PDF).

**`--pdf PATH` takes a real path, not `-`.** PDF is binary; stdout streaming isn't
worth it for a demo-grade converter. Error clearly if `-` is passed.

## Refactor: single-source the "how it works" prose

Three `{% call m.help(...) %}` blocks exist today, all with inline bodies:

- `report.html:32-41` - "How the donor and host fractions are estimated"
- `macros.html:63-73` (inside `host_presence`) - "How host-presence detection works"
- `macros.html:207-228` (inside `qc_panel`) - "How the quality-control checks work"

Extract each body into its own named macro in `macros.html`: `help_fractions()`,
`help_host_presence()`, `help_qc()` (each wrapping the `help(summary)` disclosure).
Then:

- Parameterize `host_presence(hp, detected, show_help=True)` and
  `qc_panel(ctx, show_help=True)` - when `show_help` is true (HTML default) they
  emit `{{ help_host_presence() }}` / `{{ help_qc() }}` inline as now; PDF passes
  `show_help=False`.
- `report.html` calls `{{ m.help_fractions() }}` inline where the block is now.
- `footer(..., show_run_command=True)` - PDF passes `False` and renders the
  run-command in the appendix instead.

This keeps HTML output byte-identical (existing tests such as
`test_expandable_method_help_present` and `test_byte_stable_with_fixed_timestamp`
guard this) while letting PDF relocate the same prose.

## PDF shell templates

`pdf/report.html` and `pdf/timeline.html` extend `base.html` and:

1. Render header, result, host-presence (`show_help=False`), QC (`show_help=False`),
   footer (`show_run_command=False`) up top.
2. Add a final `<section class="section methods-appendix">` "Methods and notes"
   that calls `help_fractions()` (single only), `help_host_presence()`,
   `help_qc()`, and re-emits the footer's run-command / invocation block.

## PDF stylesheet (`pdf/styles.css`)

- `@page` margin boxes: `@top-*` = recipient/sample ID (via `string-set` on a
  running-header element + `content: string(...)`), `@bottom-right` =
  `"page " counter(page) " of " counter(pages)`.
- Flatten `<details>` so the appendix prints: force `details.help .help-body`,
  `.run-command pre` to `display: block` and keep summaries visible (WeasyPrint
  collapses `<details>` otherwise).

## New render module + CLI wiring

- New `src/allomix/report/html/pdf.py`: top-level `from weasyprint import HTML`;
  `to_pdf(data, output, *, template_dir=None)` and a timeline variant, building the
  HTML string from the `pdf/` templates (reusing the context builders in
  `render.py` / `timeline.py`, refactored so context-building is shared and only
  the template name + css differ) then calling `HTML(string=html).write_pdf(output)`.
- `_add_output_args` (`cli.py:644`): add `--pdf PATH` to `detect` and `timeline`;
  add `--pdf` to the `report` subparser (`cli.py:1217`).
- In `cmd_detect`/`cmd_timeline`/`cmd_report`: guard with
  `find_spec("weasyprint")` -> `SystemExit("... needs WeasyPrint: pip install
  'allomix[pdf]'")`; timeline also keeps the matplotlib guard (chart PNG still
  needs it). Then deferred-import `pdf.py` and write.

## pyproject

Add `pdf = ["weasyprint>=60"]` extra (`pyproject.toml:32`), add it to `dev` so
tests can run, and note the system Pango/Cairo/HarfBuzz dependency.

## Docs

- `docs/reports.md`: document `--pdf` and the `[pdf]` extra.
- `docs/custom_report_template.md`: note WeasyPrint's partial flexbox/grid support
  (the headline uses a grid; PDF may differ from browser) and the system-library
  requirement (Pango/Cairo/HarfBuzz, usually present on Linux, occasionally a
  hurdle in locked-down clinical IT).

## Tests (`tests/test_html_report.py`, new PDF section)

- Round-trip JSON -> PDF for single-sample and timeline: assert a non-empty file
  starting with the `%PDF-` magic bytes.
- Extra guard: monkeypatch `find_spec` to return None and assert `cmd_*` raises
  `SystemExit` with the `allomix[pdf]` install hint (mirrors the matplotlib guard
  test).
- HTML path unchanged: existing byte-stability test guards that the macro refactor
  doesn't alter HTML output.

## Out of scope (per issue)

No sign-out / authorization / e-signature block. allomix stays an analysis tool;
the LIMS owns the legal report. A lab that wants a fileable signed report adds it
via a `--template` override.
