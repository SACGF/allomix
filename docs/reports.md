# Reports and structured output

`allomix` emits three output formats: TSV (machine-readable), JSON (the
structured artifact), and a self-contained HTML report. The JSON is canonical:
the HTML report is rendered from it and nothing else, so the report and the data
always agree.

## Structured JSON

The `--json` output is a self-describing envelope carrying the chimerism
estimate, the full QC accounting, every sub-analysis (host presence,
relatedness, contamination, swap check, run unit), the per-marker rows, and the
report provenance (recipient metadata, analysis parameters, version, timestamp).
Bioinformaticians can process it, store it, or upload it elsewhere.

The provenance also records the exact analysis invocation under
`params.command`, for example `allomix monitor --genotype-vcf ... --no-host-presence`.
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

## Customising the report

The report layout and styling can be customised for your own laboratory (logo,
colours, headers, wording) by supplying your own templates with `--template`;
see the [custom report template guide](custom_report_template.md).

## Worked examples

These are built from the public SRP434573 mixtures (a public-data
demonstration, not patient data). The titrated minor contributor (F2) is
assigned the host role, so the monitored quantity is the host fraction:

- [Single-sample report, 1%](examples/srp434573_single_sample_1pct.html)
  (the headline example): the 1% titration `1_99_F2-M1`, with its
  [per-marker CSV](examples/srp434573_single_sample_1pct.markers.csv).
  The estimate recovers the 1% host fraction cleanly.
- [Single-sample report, 0.5%](examples/srp434573_single_sample.html)
  (near the panel's contamination floor): the 0.5% titration `1_199_F2-M1`,
  with its
  [per-marker CSV](examples/srp434573_single_sample.markers.csv). The
  MLE reads slightly low because of a donor-homozygous contamination background
  in this public dataset; the host-presence test still detects the residual host
  signal.
- [Dilution series timeline](examples/srp434573_dilution_series.html):
  the whole F2-into-M1 titration ladder fed to the timeline mode to show the
  trend chart. These are a titration series, not serial timepoints from one
  patient.

All regenerate deterministically with `scripts/gen_example_report.sh` (which
shows the exact `allomix monitor` / `allomix timeline` command lines).
