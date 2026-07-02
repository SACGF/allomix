# Customising the report templates

The HTML report (`allomix detect --html`, `allomix timeline --html`, and
`allomix report`) and the PDF report (`--pdf`) are rendered from a small set of
[Jinja2](https://jinja.palletsprojects.com/) templates that ship inside the
package. A laboratory can override any of them to restyle or restructure the
report (its own colours, logo, header wording, or section layout) without
touching the allomix source, by pointing `--template` at a directory of
overrides. The PDF reuses the same section content as the HTML (see [Customising
the PDF](#customising-the-pdf) below).

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
| `macros.html` | The reusable section fragments (header band, host-presence callout, QC panel, methods footer, the collapsible help blocks) shared by `report.html`, `timeline.html`, and the PDF templates. |
| `pdf/report.html`, `pdf/timeline.html`, `pdf/styles.css` | The PDF layout (used by `--pdf`). They reuse `macros.html`, so the section content cannot drift from the HTML; `pdf/styles.css` adds only the print running header, page numbers, and the "Methods and notes" appendix styling on top of `styles.css`. |

The templates do no analysis and no number formatting: they only place
already-formatted, HTML-safe values. Formatting helpers (`pct`, `pval`, `ci`,
`count`, `badge`, ...) are available as Jinja filters and globals, and the values
themselves come from the structured report data (the same JSON `--json` writes).

## Example: restyle for your lab

Restyling is usually just CSS. Override only `styles.css`:

```bash
mkdir mylab-template
# copy the shipped stylesheet as a starting point, then edit it
python -c "import allomix.report.html, pathlib, shutil; \
  src = pathlib.Path(allomix.report.html.__file__).parent / 'templates' / 'styles.css'; \
  shutil.copy(src, 'mylab-template/styles.css')"
```

Edit `mylab-template/styles.css` (for example, change the accent colour):

```css
:root { --accent: #00558c; }   /* your lab's brand colour */
```

Then run with `--template`:

```bash
allomix detect \
  --genotype-vcf panel.vcf.gz --admix-vcf admix.vcf.gz \
  --host-sample HOST --donor-sample DONOR --sample S1 \
  --html report.html \
  --template mylab-template
```

The report renders with your stylesheet and the built-in layout. To add a logo or
change the header wording, copy `macros.html` as well and edit the `header`
macro; to change which sections appear or their order, copy `report.html` (and
`timeline.html` for the longitudinal report).

## Customising the PDF

The PDF (`--pdf`) is its own template set under a `pdf/` subdirectory:
`pdf/report.html`, `pdf/timeline.html`, and `pdf/styles.css`. Override them the
same way, by placing a `pdf/` subdirectory in your `--template` directory. They
import the shared `macros.html`, so a change to a section (say the header macro)
flows to both the HTML and the PDF; `pdf/styles.css` is layered on top of the
base `styles.css`, so it holds only print-specific rules (the `@page` running
header and page numbers, and flattening the collapsible blocks for the appendix).

A lab that wants a fully branded or signed PDF (letterhead, a sign-out block)
copies the `pdf/` templates and edits them. allomix itself stays an analysis
tool and does not add a sign-out or e-signature section.

Two WeasyPrint caveats when editing PDF templates:

- WeasyPrint's support for CSS flexbox and grid is partial. The headline uses a
  flex row and the header band uses a grid; these may lay out differently in the
  PDF than in a browser. If a custom layout looks wrong in the PDF, prefer simple
  block/table CSS over flex/grid for the print path.
- WeasyPrint needs the system Pango, Cairo, and HarfBuzz libraries (see the
  [reports guide](reports.md#pdf-report)); they are usually present on Linux but
  can be a hurdle on a locked-down clinical workstation.

## Notes

- `--template` works the same on `detect`, `timeline`, and `report`, for both
  HTML and PDF.
- Because `report` regenerates the report from a saved `--json` file, you can
  apply or change a template after the analysis has run, without recomputing
  anything.
- Keep your override directory under version control alongside your pipeline so
  the report style is reproducible.
