# CLI Usage

`allomix` has these subcommands: `detect` (single timepoint), `timeline`
(serial timepoints), `report` (render HTML/PDF from a saved JSON), `estimate-bias`
(build a per-marker bias table), and `panel-qc` (per-marker inclusion QC). The
calibration commands `estimate-errors` and `build-contamination-table` are
covered in the [Panel Guide](panel_guide.md).

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

`detect` and `timeline` sniff each VCF header for the caller that produced it and
print a stderr warning in two cases (soft warnings only, they do not stop the
run):

- the admix VCF looks GATK-called rather than `bcftools mpileup` (GATK's
  low-level filters drop the minority ALT reads a low host fraction needs), or
- bias correction is on and the bias table was estimated from a different caller
  than the admix data. Per-marker amplification bias is caller-specific, so a
  mismatched table can make the estimate worse. `estimate-bias` records its
  source caller in the table header (`# allomixCaller=`) so this check can run;
  building the table with `estimate-bias --both-het` from the same mpileup admix
  data avoids the mismatch entirely.

## detect

```bash
# Calculate chimerism for a single timepoint (TSV to stdout by default)
allomix detect \
    --genotype-vcf patient001_panel.vcf.gz \
    --admix-vcf patient001_admix.vcf.gz \
    --host-sample HOST_001 \
    --donor-sample DONOR_001 \
    --sample TP1_20240101 \
    --tsv results.tsv

# Multi-donor (2 donors)
allomix detect \
    --genotype-vcf patient001_panel.vcf.gz \
    --admix-vcf patient001_admix.vcf.gz \
    --host-sample HOST_001 \
    --donor-sample DONOR1_001 \
    --donor-sample DONOR2_001 \
    --sample TP1_20240101 \
    --tsv results.tsv

# Structured JSON (the artifact the HTML/PDF reports are rendered from)
allomix detect \
    --genotype-vcf patient001_panel.vcf.gz \
    --admix-vcf patient001_admix.vcf.gz \
    --host-sample HOST_001 \
    --donor-sample DONOR_001 \
    --sample TP1_20240101 \
    --json results.json

# Structured JSON, the HTML and PDF reports in one run, plus the per-marker CSV
# (bioinformatician-facing detail the report omits). Any output flags combine.
# --pdf needs the pdf extra: pip install 'allomix[pdf]'.
allomix detect \
    --genotype-vcf patient001_panel.vcf.gz \
    --admix-vcf patient001_admix.vcf.gz \
    --host-sample HOST_001 \
    --donor-sample DONOR_001 \
    --sample TP1_20240101 \
    --json report.json \
    --html report.html \
    --pdf report.pdf \
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

# Or render a PDF (needs the pdf extra); with only --pdf, no HTML is written to stdout
allomix report report.json --pdf report.pdf
```

## estimate-bias

```bash
# Estimate bias from per-sample VCFs (positional, globbable)
allomix estimate-bias \
    sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz \
    --output bias_table.tsv

# Estimate bias from named samples within a joint-called VCF
allomix estimate-bias \
    joint_called.vcf.gz \
    --sample DONOR_001 --sample DONOR_002 --sample DONOR_003 \
    --output bias_table.tsv

# Use bias correction during monitoring
allomix detect \
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

## panel-qc

Apply per-marker inclusion cutoffs to the panel characterization from
`scripts/measure_panel_bias.py`, emitting an auditable keep/drop verdict per
marker. See [Panel Guide, step 4](panel_guide.md#4-set-marker-inclusion-thresholds)
for the criteria and defaults.

```bash
allomix panel-qc output/panel_stats_per_marker.tsv \
    --output output/panel_qc_verdicts.tsv \
    --exclude-bed output/panel_exclude.bed

# Feed the dropped-marker BED straight into analysis
allomix detect ... --exclude-sites output/panel_exclude.bed
```

## Common options

Both `detect` and `timeline` accept these additional options:

| Option | Default | Description |
|---|---|---|
| `--min-dp` | 100 | Minimum read depth to use a marker |
| `--min-gq` | 20 | Minimum genotype quality for host/donor genotyping |
| `--error-rate` | 0.01 | Sequencing error rate for the likelihood model |
| `--bias-table` | none | Per-marker bias table TSV (from `estimate-bias`; see [Bias Estimation Guide](estimate_bias.md)) |
| `--no-bias-correction` | off | Disable bias correction even when a bias table is provided |
| `--exclude-sites` | none | BED of marker sites to drop before analysis (e.g. `panel-qc --exclude-bed`) |
| `--include-sites` | none | BED of marker sites to restrict analysis to (mutually exclusive with `--exclude-sites`) |
| `--recipient-sex` | none | Declared recipient sex. `F`/`female` or `M`/`male` (any case) is checked against the sex inferred from chrX; other text is shown in the report header only |
| `--donor-sex` | none | Declared donor sex, one per `--donor-sample` in order (repeat to match; `NA` for none). Parsed and checked like `--recipient-sex` |
| `--contig-policy` | `sex_aware` | Which contigs may enter the estimate. `sex_aware` routes non-PAR chrX on the inferred/declared sex pair (used for a sex-matched pair, homozygous sites only for a male pair; excluded otherwise); `autosomes_only` never uses chrX. PAR, chrY, MT and non-primary contigs are excluded under both |
| `--admix-depth-vcf` | none | Forced-pileup VCF with per-site `FORMAT/DP` for the admixture samples at the panel's interval midpoints (the pipeline's `<patient>.admix.midpoints.vcf.gz`). Enables the experimental chrY depth readout of the sex-mismatch cross-check; needs `--ref-depth-vcf`. Samples missing from the file are skipped with a warning |
| `--ref-depth-vcf` | none | The same forced-pileup VCF for the host and donor reference samples (the pipeline's `<patient>/refs/midpoints.vcf.gz`). Also gives sex inference its chrY-depth secondary signal (see below); usable without `--admix-depth-vcf` |
| `--verbose` | off | Include per-marker detail in output |

**Sex chromosomes.** allomix infers the sex of the recipient and of each donor
from the heterozygosity of their non-PAR chrX genotypes (see
`allomix.qc.sex`), compares it with any declared sex, and reports both. A
confident inference that contradicts the declaration fails QC as a likely
sample mix-up; a declaration resolves an ambiguous or unavailable inference
(an autosome-only panel always gives `unavailable`). With `--ref-depth-vcf`,
chrY relative depth (`max(non-PAR chrY DP) / median(autosomal DP)`, graded
with the lab's amplicon-pipeline thresholds M > 0.5, M* > 0.15, F* > 0.05,
else F) is a secondary signal: a full `M` grade turns an ambiguous or
unavailable chrX call into male (shown as "inferred (chrY depth)"). It never
changes a confident chrX call, and low chrY depth is not used to call a female,
because a panel that lists chrY targets its capture does not pull down would
read every sample as female. The resulting pair status
routes the non-PAR chrX markers: a female/female pair uses them like autosomal
markers, a male/male pair uses the sites where both are homozygous and drops
any site with a het call (a genotyping error on a hemizygous chromosome,
counted as `n_chrx_male_het_dropped`), and a mismatched or unresolved pair
excludes them (counted as `n_informative_sex_chrom_excluded`, with the reason
on stderr). Pseudoautosomal, chrY, mitochondrial and non-primary-contig
markers are never used. `--contig-policy autosomes_only` turns chrX off for a
run. The old `--use-sex-chroms` flag is retired: giving it exits with an
explanation, and it will be removed in a later release.

**Sex-mismatch cross-check.** For a sex-mismatched single-donor pair the chrX
markers the estimate excludes are fitted on their own, under the
copy-number-weighted expectation (a female contributes two X copies, a male
one), as an independent estimate of the donor fraction with its own CI
(`allomix.qc.sex_mismatch`). When `--admix-depth-vcf` and `--ref-depth-vcf`
are given and the panel has non-PAR chrY regions, the chrY depth ratio
(admixture over the male reference sample) takes priority as the basis;
this readout is experimental and has not been exercised on real data. The
result is reported beside the MLE (TSV `sexchrom_*` columns, JSON
`sex_mismatch`, an HTML panel, a stderr line) and compared with it: a
discordant pair (no CI overlap and a gap above 2 percentage points) is a QC
warning. It does not change the PASS/REVIEW/FAIL status while the readout's
real-sample behaviour is being mapped. The headline `donor_pct` is never
blended with it. Multi-donor runs skip the check (`NA`).

Output is selected by per-artifact flags that can be combined in one run:
`detect` accepts `--tsv PATH`, `--json PATH`, `--html PATH`, and `--pdf PATH`
(plus `--marker-csv PATH`); `timeline` accepts `--json PATH`, `--html PATH`, and
`--pdf PATH`. Each text artifact accepts `-` for stdout; `--pdf` is binary and
needs a file path. With no output flag, `detect` writes TSV and `timeline` writes
JSON, both to stdout. `--pdf` needs the `pdf` extra (`pip install
'allomix[pdf]'`); see the [reports guide](reports.md#pdf-report).

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
programmatic consumption and as the report source), a self-contained HTML report,
and an optional PDF (the `pdf` extra). See [Reports and structured
output](reports.md) for the JSON envelope, the HTML/PDF report, and worked
examples.
