# CLI Usage

`allomix` has four subcommands: `monitor` (single timepoint), `timeline`
(serial timepoints), `report` (render HTML from a saved JSON), and
`estimate-bias` (build a per-marker bias table).

## Two-VCF input

`allomix` takes two VCFs:

- a **panel VCF** with host/donor genotypes, typically from GATK joint calling
  of the reference samples, and
- a separate **admix VCF** with per-timepoint AD counts, typically from forced
  `bcftools mpileup` at the panel sites.

Joint calling of HOST + DONOR ensures ALT alleles discovered in the donor are
propagated to the panel even when one sample is hom-ref. Pileup of the ADMIX
samples preserves raw per-allele counts at the panel sites, which is essential
for detecting host fractions below ~5% (GATK's GVCF mode strips minority ALT
reads at hom-ref blocks).

A ready-to-use Snakemake pipeline that produces both files is included in
`pipeline/`. See the [Joint Calling Guide](joint_calling.md) for the two-phase
rationale and how to run it. When a new timepoint arrives, re-run the admix-only
pileup for it (the panel does not need rebuilding), then re-run allomix on the
updated admix VCF.

## monitor

```bash
# Calculate chimerism for a single timepoint (TSV to stdout by default)
allomix monitor \
    --genotype-vcf patient001_panel.vcf.gz \
    --admix-vcf patient001_admix.vcf.gz \
    --host-sample HOST_001 \
    --donor-sample DONOR_001 \
    --sample TP1_20240101 \
    --tsv results.tsv

# Multi-donor (2 donors)
allomix monitor \
    --genotype-vcf patient001_panel.vcf.gz \
    --admix-vcf patient001_admix.vcf.gz \
    --host-sample HOST_001 \
    --donor-sample DONOR1_001 \
    --donor-sample DONOR2_001 \
    --sample TP1_20240101 \
    --tsv results.tsv

# Structured JSON (the artifact the HTML report is rendered from)
allomix monitor \
    --genotype-vcf patient001_panel.vcf.gz \
    --admix-vcf patient001_admix.vcf.gz \
    --host-sample HOST_001 \
    --donor-sample DONOR_001 \
    --sample TP1_20240101 \
    --json results.json

# Structured JSON and the HTML report in one run, plus the per-marker CSV
# (bioinformatician-facing detail the report omits). Any output flags combine.
allomix monitor \
    --genotype-vcf patient001_panel.vcf.gz \
    --admix-vcf patient001_admix.vcf.gz \
    --host-sample HOST_001 \
    --donor-sample DONOR_001 \
    --sample TP1_20240101 \
    --json report.json \
    --html report.html \
    --marker-csv report.markers.csv
```

## timeline

```bash
# Timeline across multiple timepoints (JSON by default, --html for a trend chart)
allomix timeline \
    --genotype-vcf patient001_panel.vcf.gz \
    --admix-vcf patient001_admix.vcf.gz \
    --host-sample HOST_001 \
    --donor-sample DONOR_001 \
    --sample TP1_20240101 \
    --sample TP2_20240201 \
    --sample TP3_20240301 \
    --json timeline.json
```

## report

```bash
# Render the HTML report later from a saved JSON
allomix report report.json --output report.html
```

## estimate-bias

```bash
# Estimate bias from per-sample VCFs
allomix estimate-bias \
    --genotype-vcfs sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz \
    --output bias_table.tsv

# Estimate bias from named samples within a joint-called VCF
allomix estimate-bias \
    --genotype-vcf joint_called.vcf.gz \
    --samples DONOR_001 DONOR_002 DONOR_003 \
    --output bias_table.tsv

# Use bias correction during monitoring
allomix monitor \
    --genotype-vcf patient001_panel.vcf.gz \
    --admix-vcf patient001_admix.vcf.gz \
    --host-sample HOST_001 \
    --donor-sample DONOR_001 \
    --sample TP1_20240101 \
    --bias-table bias_table.tsv \
    --tsv results.tsv
```

If you do not yet have enough donor VCFs to train a bias table, `estimate-bias`
can also be driven from archived BAMs on the same panel via a joint-calling
pipeline plus sample-level QC. See
[Building a training cohort from BAMs](estimate_bias.md#building-a-training-cohort-from-bams)
in the bias guide.

## Common options

Both `monitor` and `timeline` accept these additional options:

| Option | Default | Description |
|---|---|---|
| `--min-dp` | 100 | Minimum read depth to use a marker |
| `--min-gq` | 20 | Minimum genotype quality for host/donor genotyping |
| `--error-rate` | 0.01 | Sequencing error rate for the likelihood model |
| `--bias-table` | none | Per-marker bias table TSV (from `estimate-bias`; see [Bias Estimation Guide](estimate_bias.md)) |
| `--no-bias-correction` | off | Disable bias correction even when a bias table is provided |
| `--verbose` | off | Include per-marker detail in output |

Output is selected by per-artifact flags that can be combined in one run:
`monitor` accepts `--tsv PATH`, `--json PATH`, and `--html PATH` (plus
`--marker-csv PATH`); `timeline` accepts `--json PATH` and `--html PATH`. Each
accepts `-` for stdout. With no output flag, `monitor` writes TSV and `timeline`
writes JSON, both to stdout.

## Inputs and outputs

### Inputs

The tool works with VCFs from any variant calling pipeline that supports joint
calling (GATK GenomicsDBImport + GenotypeGVCFs) as long as GT and AD fields are
present. Higher depth improves sensitivity; panels with >1000x coverage give the
best results at low chimerism fractions. Sample names are specified on the
command line via `--host-sample`, `--donor-sample`, and `--sample`.

### Outputs

| Output | Description |
|---|---|
| % chimerism | Estimated fraction of donor cells (per donor if multi-donor) |
| Confidence interval | 95% CI on the chimerism estimate |
| QC metrics | Number of informative markers used, mean depth, markers excluded and why, goodness-of-fit |
| Per-marker detail | Allele depths, expected vs observed VAF, residual, and the include flag for each marker (per-marker CSV, or the verbose TSV / JSON) |
| Timeline report | Chimerism trend across serial timepoints for a patient |

Output formats are TSV (machine-readable), JSON (the structured artifact, for
programmatic consumption and as the report source), and a self-contained HTML
report. See [Reports and structured output](reports.md) for the JSON envelope,
the HTML report, and worked examples.
