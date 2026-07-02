# Reports and structured output

`allomix` emits several output formats: TSV (machine-readable), JSON (the
structured artifact), a self-contained HTML report, and an optional PDF. The JSON
is canonical: the HTML and PDF reports are rendered from it and nothing else, so
the report and the data always agree.

## Structured JSON

The `--json` output is a self-describing envelope carrying the chimerism
estimate, the full QC accounting, every sub-analysis (host presence,
relatedness, contamination, swap check, run unit), the per-marker rows, and the
report provenance (recipient metadata, analysis parameters, version, timestamp).
Bioinformaticians can process it, store it, or upload it elsewhere.

The provenance also records the exact analysis invocation under
`params.command`, for example `allomix detect --genotype-vcf ... --no-host-presence`.
Output-destination and presentation flags (`--json`, `--html`, `--tsv`,
`--marker-csv`, `--output`, `--template`, `--report-timestamp`, `--log-scale`,
`--verbose`) are stripped, so the recorded command captures only what determined
the result and stays byte-stable across runs to different output paths.

Because the report is rendered from exactly this structure,
`allomix report saved.json --output report.html` regenerates the report from
saved data without re-running the analysis (on any machine).

## HTML report

`--html` writes a single self-contained HTML file (all CSS and JavaScript
inlined, no network access needed) suitable for review or attaching to a record.
It is written for the clinician: the headline chimerism fractions with CIs, the
host-presence callout, a plain-language QC breakdown (informative-marker
accounting and the number of QC flags), the QC panel, and a methods/provenance
footer. The footer lists the analysis settings (thresholds, error/bias model,
and the on/off toggles) and ends with a collapsed "Run command" section holding
the recorded invocation; that section is the only place full file paths appear,
so it is hidden by default.

It does not include the per-marker table; that detail is for bioinformaticians
and is written separately with `--marker-csv PATH` (one row per marker per
sample: allele depths, observed vs expected VAF, residual, and the include
flag).

`allomix timeline --html` adds a trend chart across timepoints; it needs
matplotlib, installed with the `report` extra:

```bash
pip install 'allomix[report]'
```

## PDF report

`--pdf PATH` writes the same report as a PDF, for attaching to a record or a
LIMS. It is available on `detect`, `timeline`, and `report`, and can be produced
alongside the other outputs in one run:

```bash
allomix detect ... --json result.json --pdf report.pdf
allomix report result.json --pdf report.pdf     # from a saved JSON
```

The PDF carries the same sections as the HTML, laid out for print: a running
header, page numbers, and a "Methods and notes" appendix that gathers the "how it
works" explanations and the run command (in the HTML these are collapsible, which
a printed page cannot show). PDF is binary, so `--pdf -` (stdout) is not
supported; give a file path. When `allomix report` is given only `--pdf` (no
`--output`), it writes just the PDF and does not also print HTML to stdout.

PDF output needs WeasyPrint, installed with the `pdf` extra (which also pulls in
the `report` extra for the timeline chart):

```bash
pip install 'allomix[pdf]'
```

WeasyPrint depends on the system Pango, Cairo, and HarfBuzz libraries. These are
usually already present on Linux; on a locked-down clinical workstation they may
need installing (for example `apt install libpango-1.0-0 libpangocairo-1.0-0`).

This is a demo-grade converter, not a LIMS replacement. A lab that wants a fully
branded or signed report overrides the PDF templates with `--template` (see
below).

## Customising the report

The report layout and styling can be customised for your own laboratory (logo,
colours, headers, wording) by supplying your own templates with `--template`;
see the [custom report template guide](custom_report_template.md).

## Worked examples

These are built from the public SRP434573 mixtures (a public-data
demonstration, not patient data). The titrated minor contributor (F2) is
assigned the host role, so the monitored quantity is the host fraction:

The reports below are hosted on [GitHub Pages](https://sacgf.github.io/allomix/)
so they render in the browser (the raw `.html` files in `docs/examples/` show as
source if opened directly on GitHub).

- [Single-sample report, 1%](https://sacgf.github.io/allomix/srp434573_single_sample_1pct.html)
  (the headline example): the 1% titration `1_99_F2-M1`, with its
  [per-marker CSV](https://sacgf.github.io/allomix/srp434573_single_sample_1pct.markers.csv).
  The estimate recovers the 1% host fraction cleanly.
- [Single-sample report, 0.5%](https://sacgf.github.io/allomix/srp434573_single_sample.html)
  (near the panel's contamination floor): the 0.5% titration `1_199_F2-M1`,
  with its
  [per-marker CSV](https://sacgf.github.io/allomix/srp434573_single_sample.markers.csv). The
  MLE reads slightly low because of a donor-homozygous contamination background
  in this public dataset; the host-presence test still detects the residual host
  signal.
- [Dilution series timeline](https://sacgf.github.io/allomix/srp434573_dilution_series.html):
  the whole F2-into-M1 titration ladder fed to the timeline mode to show the
  trend chart. These are a titration series, not serial timepoints from one
  patient.

All regenerate deterministically with `scripts/gen_example_report.sh` (which
shows the exact `allomix detect` / `allomix timeline` command lines).
