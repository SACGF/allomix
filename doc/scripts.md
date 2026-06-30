# Scripts

Utility and validation scripts live in `scripts/`. They are **not** part of the
installed `allomix` package; run them from the repo root with `python
scripts/<name>.py`. They depend on `allomix` being installed (`pip install -e
".[dev]"`); the plotting script also needs `matplotlib`.

> The example sample codes below (`AAAA`, `BBBB`, ...) are placeholders. Use your
> own sample names, and keep real sample sheets and `/tau` paths out of version
> control.

## Overview

| Script | Purpose |
|--------|---------|
| `run_csv_batch.py` | Run `allomix monitor` for every per-patient CSV in `pipeline/sample_csvs/` against the two-VCF pipeline outputs, and combine the per-patient results into one `batch.tsv`. |
| `plot_chimerism_comparison.py` | Plot whole-blood NGS chimerism (with CIs) against flow-sorted lineage values, optionally overlaying one or more other runs. |
| `diagnose_sample.py` | Per-marker residuals and noise model for one admixture sample: localise a goodness-of-fit failure (CNV/LOH) by chromosome and show why a per-sample LOD is what it is (fitted overdispersion `rho`, SE, LoB, LoD). |
| `plot_host_presence_per_marker.py` | Per-marker host-presence structure for selected admixture samples: the host fraction each donor-homozygous marker implies, with binomial CIs and the pooled MLE line, to see whether the host signal is spread evenly or carried by a few markers (a CNV/LOH clue). |
| `host_presence_manhattan.py` | Manhattan plot of the same per-marker host signal along the genome, so a host CNV/LOH region shows up as a contiguous run of markers lifted above the pooled host fraction (UPREG markers ringed) rather than a lone artifact. |
| `host_presence_markers_vcf.py` | Write a VCF of the donor-homozygous markers carrying host-presence signal, flagging the "upregulated" ones, with INFO fields so it can be annotated by VEP or intersected with driver-gene panels for CNV/LOH follow-up. |
| `plot_informative_karyogram.py` | Karyogram of informative (IBS0/IBS1) markers along the genome for host vs donor(s); IBD blocks read as gaps in the IBS0 ticks, a coarse relatedness sanity check. |
| `run_validation.py` | Run allomix on synthetic test data against a truth table and produce a validation report. |
| `generate_test_data.py` | Generate synthetic joint-called VCFs (host + donor + admixture) for tests. |
| `generate_multidonor_test_data.py` | Generate a multi-donor synthetic dataset (host + 2 related donors). |
| `generate_timeline_data.py` | Generate a synthetic post-HSCT timeline (serial chimeric VCFs). |
| `make_synthetic_genotypes.py` | Create synthetic host/donor genotype VCFs (100 biallelic SNPs). |
| `measure_panel_bias.py` | Measure per-marker bias/characteristics from joint-called genotyping VCFs. |
| `make_midpoint_bed.py` | Derive a thin background-site BED (one position per amplicon, the interval midpoint) from a capture BED, for the host-presence error-table background. Equivalent to the pipeline's `midpoint_bed` rule; use standalone to inspect or commit the positions. |
| `qc_bias_samples.py` | Sample-level QC for bias-estimation training samples. |
| `mix_bams.sh` | Mix N BAMs at explicit target fractions into one synthetic admixture sample, depth-normalized on on-target reads (host + 1 or 2 donors). |

The two most commonly used scripts are documented in detail below.

## `run_csv_batch.py` — batch the joint-calling pipeline outputs

Drives `allomix monitor` across every per-patient CSV in
`pipeline/sample_csvs/`. For each CSV it locates the matching
`<patient>.vcf.gz` (host/donor genotypes) and `<patient>.admix.vcf.gz`
(raw pileup AD) in `--vcf-dir`, runs `allomix monitor --panel-vcf ...
--admix-vcf ...` once per patient with all ADMIX timepoints, and
concatenates the per-patient TSVs into `batch.tsv`. Patients with no
ADMIX rows are skipped automatically.

### Example

```bash
python scripts/run_csv_batch.py \
    --samples-csv-dir pipeline/sample_csvs \
    --vcf-dir output/genotypes \
    --output-dir output/validation_run2 \
    --bias-table output/bias_training/bias_table.tsv \
    --error-table output/error_training/error_table.tsv
```

This produces `output/validation_run2/<patient>.tsv` for each patient
plus `output/validation_run2/batch.tsv`.

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `--samples-csv-dir` | `pipeline/sample_csvs` | Directory of per-patient CSVs (one per file). |
| `--vcf-dir` | `output/genotypes` | Directory containing `<patient>.vcf.gz` and `<patient>.admix.vcf.gz` from `pipeline/Snakefile`. |
| `--output-dir` | `output/batch` | Where per-patient TSVs and the combined `batch.tsv` are written. |
| `--bias-table` | none | Per-marker bias table passed through to `allomix monitor`. |
| `--error-table` | none | Per-site error table passed through to `allomix monitor`. |
| `--allomix` | `allomix` | Path to the `allomix` executable. |
| `--extra-arg ARG` | none | Forward an extra argument to every `allomix monitor` call. Repeat for multiple, e.g. `--extra-arg --min-dp=200 --extra-arg --min-gq=30`. |

### `batch.tsv` columns

The combined file has one row per successfully run admix sample, with the
columns `allomix monitor` emits:

```
sample  donor_pct  ci_lo  ci_hi  lob_pct  lod_pct  n_informative  n_used
        mean_depth  gof_pval  qc_status  qc_warnings
```

- `lob_pct` / `lod_pct` — per-sample limit of blank / limit of detection.
- `qc_status` — `PASS`, `REVIEW`, or `FAIL`. `FAIL` means the result is unusable
  (e.g. too few informative markers). `REVIEW` means it was computed but a
  reliability check failed (poor model fit or wide CI), so it needs manual
  interpretation rather than being trusted or discarded automatically.
- `qc_warnings` — `; `-joined QC warnings (empty when clean). For a `FAIL` or
  `REVIEW` sample this names the cause, e.g. which input starved the informative
  markers, or the poor-fit reason.

> Note: `plot_chimerism_comparison.py` reads the flow lineage column
> (and donor-type column) directly from `batch.tsv`. The old
> `run_xls_batch.py` driver merged those columns in via
> `--copy-columns`. `run_csv_batch.py` does not do this yet — see the
> follow-up in the plotting section below if you need flow overlay.

## `plot_chimerism_comparison.py` — NGS vs flow, and run-to-run

Draws a per-sample plot of the NGS donor estimate (with confidence interval)
against the flow-sorted lineage values parsed from a `--copy-columns` column.
The y axis is donor %, log-spaced by distance from 100 % and inverted so 100 %
sits at the top, which keeps the low-level signal near full chimerism readable.
QC-FAIL samples are skipped (and listed on the console).

### Example 1: single run vs flow

```bash
python scripts/plot_chimerism_comparison.py output/validation_run2/batch.tsv \
    --flow-column "Chimerism result TP2" \
    --label-code \
    --output output/chimerism_comparison.png
```

### Example 2: two or more runs compared (and flag explicit-donor samples)

```bash
python scripts/plot_chimerism_comparison.py output/validation_run3/batch.tsv \
    --compare-tsv output/validation_run1/batch.tsv output/validation_run2/batch.tsv \
    --flow-column "Chimerism result TP2" \
    --labels run1 run2 run3 \
    --label-code \
    --explicit-donor REDACTED \
    --title "run1 vs run2 vs run3 (explicit donor ★)" \
    --output output/run1_vs_run2_vs_run3.png
```

The primary file (the positional argument) is plotted on the right of each
sample and drawn filled; the `--compare-tsv` files are drawn hollow to its left
in the order given, so list them oldest first and time reads left to right.
Drop the extra `--compare-tsv` path (and a label) for a two-run plot, or omit
`--compare-tsv` entirely for a single run. Flow lineage markers are drawn once
(from the primary file). Each run gets its own x-axis label row, coloured to
match, showing that run's marker count and mean depth (e.g. `M:49 D:1356x`); the
patient code (with any `★`) sits on the primary row.

### Options

| Option | Default | Description |
|--------|---------|-------------|
| `batch_tsv` (positional) | — | Primary `batch.tsv`. Drawn filled and rightmost; flow, donor type and the `★` are read from it. |
| `--compare-tsv` | none | One or more `batch.tsv` files to overlay (hollow), oldest first. |
| `--flow-column` | `Chimerism result TP2` | Column holding flow lineage strings like `CD45 100%; CD3 98%; CD13 88%`. |
| `--labels LABEL ...` | `run1 run2 ...` | One legend label per run, left to right (compare runs first, then primary). Count must match the number of runs. |
| `--label-code` | off | Shorten x labels to the patient code (last all-uppercase token; robust to the code being at name field 3 or 4). |
| `--label-field N` | none | Alternative: shorten by splitting the sample name on `_` and taking field `N` (0-based). |
| `--explicit-donor TOK1,TOK2` | none | Sample tokens (matched as substrings) that had an explicit donor genotype; their primary-run label gets a `★`. Pair with a `★` in `--title`. |
| `--sort {name,tsv,chimerism}` | `name` | X-axis sample order. `name` (alphabetical) and `tsv` (file order) are stable across runs so two plots line up; `chimerism` orders by measured donor fraction (reshuffles between runs). |
| `--hide-lod` | off | Suppress the per-sample LOD band. The band is drawn by default when the primary `batch.tsv` has an `lod_pct` column. |
| `--anonymize` | off | Replace sample names with `S1, S2, ...`. |
| `--floor` | `0.02` | Host-% floor for the log axis. |
| `--title` | "Whole-blood NGS chimerism vs flow lineages" | Plot title. |
| `--output` | `output/chimerism_comparison.png` | Output PNG path. |

### Reading the plot

- **Blue (filled)** = primary run; compare runs are open circles in a colour
  cycle (orange, green, red, ...) matching the legend. The run name is carried
  by colour, not repeated in every label.
- Flow `CD45` is the whole-blood comparator; the `CD3`/`CD13` spread brackets
  where the true whole-blood value can sit. A correct NGS estimate falls inside
  that spread.
- The grey band at the top of each sample is the primary run's per-sample LOD:
  a point inside it is below the limit of detection, i.e. not a reportable
  detection (statistically consistent with full donor) even if its CI excludes
  100%. Use `--hide-lod` to suppress it.
- A red ring around a primary-run point marks a QC-REVIEW sample (e.g. poor
  model fit), so a confident-looking estimate is not mistaken for a clean one.
  QC-FAIL samples are skipped entirely (and listed on the console). On older
  `batch.tsv` files without a `qc_status` column, REVIEW is inferred from a
  failing `gof_pval`.
- A `★` on a sample's primary-run label marks an explicit-donor sample.
- The donor type (e.g. "Matched sibling donor") is printed in grey beneath each
  sample, read once from the primary file.
