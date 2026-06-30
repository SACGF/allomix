# Customising the HTML report template

The HTML report (`allomix monitor --html`, `allomix timeline --html`, and
`allomix report`) is rendered from a small set of [Jinja2](https://jinja.palletsprojects.com/)
templates that ship inside the package. A laboratory can override any of them to
restyle or restructure the report (its own colours, logo, header wording, or
section layout) without touching the allomix source, by pointing `--template` at
a directory of overrides.

## How overrides work

The renderer loads templates from two places, in order:

1. your override directory (the `--template DIR` argument), if given;
2. the built-in templates that ship with allomix.

Any file you place in your directory replaces the built-in of the same name;
anything you leave out falls back to the built-in. So you can override just the
stylesheet, or just one section, and keep everything else as-is. The output is
still a single self-contained HTML file (your CSS and JS are inlined), so it
needs no network access and survives being emailed or filed.

## The template files

These are the files you can override. Copy the ones you want to change out of the
installed package (they live in `src/allomix/html/templates/` in the source
tree) into your own directory and edit the copies.

| File | What it controls |
|---|---|
| `styles.css` | All styling: colours, fonts, spacing, print layout. Inlined verbatim, so it may contain any characters. **Start here for a lab restyle.** |
| `report.js` | The small table-sort script (timeline report only). Inlined verbatim. |
| `base.html` | The page shell: `<head>`, the `<style>`/`<script>` wrappers, the overall page frame. |
| `report.html` | The single-sample report: which sections appear and in what order. |
| `timeline.html` | The longitudinal (timeline) report layout. |
| `macros.html` | The reusable section fragments (header band, host-presence callout, QC panel, methods footer, the collapsible help blocks) shared by `report.html` and `timeline.html`. |

The templates do no analysis and no number formatting: they only place
already-formatted, HTML-safe values. Formatting helpers (`pct`, `pval`, `ci`,
`count`, `badge`, ...) are available as Jinja filters and globals, and the values
themselves come from the structured report data (the same JSON `--json` writes).

## Example: restyle for your lab

Restyling is usually just CSS. Override only `styles.css`:

```bash
mkdir mylab-template
# copy the shipped stylesheet as a starting point, then edit it
python -c "import allomix.html, pathlib, shutil; \
  src = pathlib.Path(allomix.html.__file__).parent / 'templates' / 'styles.css'; \
  shutil.copy(src, 'mylab-template/styles.css')"
```

Edit `mylab-template/styles.css` (for example, change the accent colour):

```css
:root { --accent: #00558c; }   /* your lab's brand colour */
```

Then run with `--template`:

```bash
allomix monitor \
  --panel-vcf panel.vcf.gz --admix-vcf admix.vcf.gz \
  --host-sample HOST --donor-sample DONOR --sample S1 \
  --html report.html \
  --template mylab-template
```

The report renders with your stylesheet and the built-in layout. To add a logo or
change the header wording, copy `macros.html` as well and edit the `header`
macro; to change which sections appear or their order, copy `report.html` (and
`timeline.html` for the longitudinal report).

## Notes

- `--template` works the same on `monitor`, `timeline`, and `report`.
- Because `report` regenerates the HTML from a saved `--json` file, you can apply
  or change a template after the analysis has run, without recomputing anything.
- Keep your override directory under version control alongside your pipeline so
  the report style is reproducible.
