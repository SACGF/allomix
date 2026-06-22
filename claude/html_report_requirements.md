# HTML Report Requirements

Requirements for a standalone HTML chimerism report produced by allomix (issue #27).

## Motivation

allomix currently emits a CSV/TSV summary and optional plot files. For clinical
adoption the output needs to be a single, self-explanatory document that a
laboratory scientist or clinician can read, sign off, and file without needing
to know how the tool works internally. A well presented report is one of the
main things that distinguishes a research script from a tool people trust for
patient reporting. The report should make the result, its uncertainty, and its
quality status legible at a glance, and put the supporting detail one scroll
away for anyone who wants to check it.

## Design Principles

- **Single file, no dependencies.** The report is one `.html` file with all CSS,
  JavaScript, and images (charts) inlined or embedded as data URIs. It must
  render correctly when opened directly from disk with no network access, and
  survive being emailed or attached to a record as a single artifact.
- **Print and PDF friendly.** The layout must paginate sensibly through the
  browser "Print to PDF" path. A print stylesheet should set page breaks so
  sections do not split awkwardly, hide interactive-only controls, and keep
  colour-coded status legible in greyscale (use shape/text, not colour alone).
- **Clinical-first hierarchy.** The headline result and QC verdict come first.
  Methodological detail (per-marker tables, parameters) comes last. A reader
  should be able to answer "what is the chimerism, can I trust it, what changed
  since last time" from the first screen.
- **No silent omissions.** If a section has no data (for example, only one
  timepoint so no trend can be drawn, or host-presence detection was disabled),
  say so explicitly rather than dropping the section. A blank where the reader
  expects information is worse than a stated "not available".
- **Faithful to the computation.** Every number shown must trace to a field the
  tool actually produced. Do not invent, round away meaningful precision, or
  imply a precision the estimate does not have. Confidence intervals,
  detection limits, and QC thresholds are shown alongside point estimates, not
  hidden.

## Report Variants

Two report shapes are needed, sharing the same components and styling:

1. **Single-sample report.** One admixture sample analysed against its host and
   donor genotype(s). Used for a one-off measurement.
2. **Timeline (longitudinal) report.** Multiple serial timepoints for one
   recipient. This is the primary clinical view: it adds the trend chart and a
   per-timepoint table, and is the format most relevant to post-transplant
   monitoring where direction of change matters more than any single value.

A batch-level wrapper (an index linking to per-recipient reports) is a possible
later addition but is out of scope for the first version.

## Report Structure

### 1. Header / identification band

A compact band at the top of every report:

- Recipient identifier and optional display name.
- Optional demographics if supplied: sex, date of birth.
- Transplant context if supplied: transplant type (default HSCT but allow
  others), transplant date, and derived days-post-transplant for each sample.
- Donor(s): identifier per donor, and the declared relationship to the recipient
  (related/unrelated, and degree if known).
- Report metadata: report generation date/time, allomix version, and a clear
  marker that this is the analysis output (not a substitute for clinical
  judgement). A short disclaimer line is appropriate.

Recipient/donor/transplant metadata is optional input. The report must render
cleanly when it is absent and must never block on missing demographics.

### 2. Headline result

The single most important block. For each admixture sample:

- **Donor fraction** and **host (recipient) fraction** as percentages, shown
  prominently. For multi-donor cases, show each donor fraction separately plus
  the combined donor total and the host fraction.
- The **95% confidence interval** for each fraction, shown next to the point
  estimate (for example "94.8% (95% CI 94.0 to 95.6)").
- **Sample-specific sensitivity**: the limit of blank and limit of detection for
  this sample, expressed as percentages. This is a genuine strength to surface,
  because it tells the clinician the smallest fraction this particular sample
  could have detected rather than quoting a fixed assay-wide number.
- The overall **QC verdict** (PASS / REVIEW / FAIL) as a clear status badge,
  with REVIEW and FAIL drawing the eye.

For a timeline report, the headline shows the most recent timepoint's result,
with the previous value and the change since the prior timepoint.

### 3. Host-presence / detection callout

allomix runs a dedicated host-presence detector that can flag a small recurring
host signal before it shows up as a meaningful fraction. When host presence is
detected this should be a distinct, prominent callout (this is an early-warning
signal clinically and is worth more than a buried table cell):

- Detection p-value and the estimated host fraction with its confidence
  interval.
- Number of markers the test used and the error-rate source.
- Plain-language interpretation (for example "low-level host signal detected" vs
  "no host signal above background").

If detection was disabled, state that explicitly.

### 4. Trend chart (timeline reports)

A chimerism-over-time chart, embedded as an inline image or self-contained
vector/canvas:

- X axis: time (calendar date and/or days post-transplant).
- Y axis: percentage. Support a log scale option, since clinically the
  interesting movement is often at the low end (sub-1% host signal returning).
- Plot donor and host fractions across timepoints, with confidence intervals
  shown as error bars or a shaded band.
- Mark points that did not pass QC distinctly so a reviewer does not read trend
  into an unreliable value.
- If only one timepoint exists, replace the chart with a clear note rather than
  drawing a single dot.

### 5. Quality control panel

A structured QC section that mirrors how a lab would audit the run. Present each
QC check as a row with: check name, the value the tool measured, the
pass/review/fail thresholds, a status indicator, and a one-line plain-language
description of what the check means and why it matters. Checks to surface, drawn
from what allomix computes:

- **Marker counts**: total markers, shared across all inputs, informative
  (differ between host and donor), and used in the final fit. For multi-donor
  runs, informative-marker count per donor.
- **Depth / coverage**: mean, median, and minimum depth across used markers,
  with a low-coverage warning threshold.
- **Markers excluded**: counts excluded for low depth, for quality, and as
  residual outliers, plus the robust-trim fraction.
- **Goodness of fit**: the model fit p-value (and the pre-trim value), with a
  review threshold, so a reader can see whether the data actually match the
  fitted mixture.
- **Contamination**: the estimated third-party contamination fraction, its
  p-value, and the number of markers used, against warn/review thresholds.
- **Sample-swap / consistency check**: the admixture-vs-(host+donor)
  consistency result, including the discordant-marker fraction and its p-value.
- **Relatedness / identity**: the measured kinship between each reference pair,
  with confidence and the inferred relationship, compared against the declared
  relationship. Flag a mismatch (a possible sample swap or mislabelling).
- **Run provenance / index-hop risk**: the sequencing run unit and whether the
  admixture sample shares a run with the host (an index-hopping contamination
  risk), when this can be determined.

Each check shows its own PASS/REVIEW/FAIL state and the overall verdict is the
worst of them. All free-text QC warnings the tool produced should be listed.

### 6. Marker-level detail

A collapsible / appendix table for reviewers who want to inspect the evidence.
This can be hidden by default to keep the top of the report clean. Per
informative marker:

- Chromosome and position.
- Marker type (the host/donor genotype configuration).
- Reference and alternate allele depths, total depth, and observed VAF.
- Expected VAF at the fitted fraction and the residual.
- Whether the marker was included or excluded as an outlier.

For multi-donor and genotype review, also support showing the host and donor
genotypes per marker. Given panels can be small (tens of markers) up to larger
sets, the table should be sortable and should make excluded/outlier markers
visually distinct.

### 7. Methods / provenance footer

A footer block that makes the report reproducible and auditable:

- allomix version.
- Input files referenced (names only; no patient-identifying paths).
- All analysis parameters used: minimum depth, minimum GQ, error rate, robust
  setting, whether bias/error tables were applied, host-presence and
  artifact-filter settings, sex-chromosome handling.
- The MLE method citation.
- Generation timestamp.

This supports the audit-trail expectation for clinical use: a reader should be
able to see exactly how the numbers were produced.

## Visual / Layout Requirements

- Clean, neutral, clinical aesthetic. Generous whitespace, a readable serif or
  sans body font, clear section headings, no decorative chrome.
- A consistent status colour scheme: pass (green), review/warn (amber), fail
  (red), reinforced with an icon or text label so it survives greyscale and is
  accessible to colour-blind readers.
- Responsive enough to read on a laptop screen and to print to A4/Letter.
- Numbers right-aligned in tables; consistent significant figures (percentages
  to a sensible decimal place, p-values in a consistent format).
- The headline result and QC badge must be visible without scrolling on a
  typical screen.

## Technical Requirements

- Generated by a new path in the existing report module, reusing the same result
  objects that already feed the TSV/JSON output (the chimerism result, QC
  report, host-presence, relatedness, contamination, run-unit, and per-marker
  data). No new analysis logic in the report layer; it only formats existing
  fields.
- Exposed through the CLI (a `--format html` option on the monitor/timeline
  output, or an equivalent), writing to a user-specified path.
- Charts produced without requiring a running browser or heavy JS framework:
  either render static images server-side and embed them, or emit a small
  self-contained inline chart. Keep the dependency footprint minimal and
  consistent with the project's existing optional-extras model.
- Templating kept simple and dependency-light; prefer a single template module
  over pulling in a large framework.
- Deterministic output for a given input, so reports can be regression-tested.

## Out of Scope (first version)

- Interactive editing of patient/sample metadata inside the report.
- A batch dashboard or multi-recipient index page.
- Live filtering of marker tables beyond simple sort/collapse.
- Any write-back to a LIMS or external system.
