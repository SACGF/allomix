# allomix

NGS-based donor chimerism monitoring for hematopoietic stem cell transplantation (HSCT).

## Overview

**allomix** calculates donor chimerism percentages from NGS data, replacing traditional STR-based analysis with a higher-sensitivity SNP-based approach. It works with any panel of bi-allelic markers (SNPs or indels) and is designed for clinical laboratories monitoring engraftment and relapse after HSCT.

The tool is panel-agnostic: it operates on whatever markers are present in the input VCFs. Bring your own panel: whether that's 24 indels, 76 SNPs, 202 SNPs, or any other set of bi-allelic loci with sufficient depth.

Chimerism MLE methodology based on [Crysup & Woerner](https://pubmed.ncbi.nlm.nih.gov/36152508/) (2022)

## Clinical Context

After HSCT, patients carry a mixture of their own (host) and transplanted (donor) cells. Monitoring the ratio of donor to host cells over time (**chimerism monitoring**) is critical for detecting graft rejection or disease relapse. Early detection allows early intervention.

Current STR-based methods have limited sensitivity (~3-5% LOD) and require separate workflows. allomix aims to:

- Achieve **<1% sensitivity** for detecting minority cell populations
- Support **up to 3 genomes** (host + 2 donors) for patients with multiple transplants
- Provide **timeline tracking** of chimerism across serial post-HSCT timepoints

## Installation

```bash
pip install allomix
```

For development:

```bash
git clone https://github.com/dlawrence/allomix.git
cd allomix
uv pip install -e ".[dev]"
```

## Workflow

```
1. Genotyping            Sequence host and each donor individually
   (upstream)              → per-sample GVCFs at marker loci

2. allomix estimate-bias (optional) Estimate per-marker amplification
                         bias from genotyping VCFs
                           → bias table TSV

3. Sequencing            Sequence post-HSCT admixture samples at
   (upstream)            serial timepoints (>=3 per patient)
                           → per-sample GVCFs at marker loci

4. Joint calling         Combine all GVCFs (host + donor + all
   (upstream)            timepoints) with GenomicsDBImport +
                         GenotypeGVCFs → one joint-called VCF

5. allomix monitor       Calculate chimerism for each sample
                           → per-sample TSV or JSON

   allomix timeline      Track chimerism across timepoints
                           → multi-timepoint JSON
```

**Two-VCF input.** allomix takes a panel VCF (host/donor genotypes, typically from GATK joint calling of the reference samples) and a separate admix VCF (per-timepoint AD counts, typically from forced `bcftools mpileup` at the panel sites). Joint calling of HOST + DONOR ensures ALT alleles discovered in the donor are propagated to the panel even when one sample is hom-ref. Pileup of the ADMIX samples preserves raw per-allele counts at the panel sites, which is essential for detecting host fractions below ~5% (GATK's GVCF mode strips minority ALT reads at hom-ref blocks).

A ready-to-use Snakemake pipeline that produces both files is included in `pipeline/`. See [Joint Calling Guide](doc/joint_calling.md) for the two-phase rationale and how to run it.

When a new timepoint arrives, re-run the admix-only pileup for it (the panel does not need rebuilding), then re-run allomix on the updated admix VCF.

## Usage

```bash
# Calculate chimerism for a single timepoint
allomix monitor \
    --panel-vcf patient001_panel.vcf.gz \
    --admix-vcf patient001_admix.vcf.gz \
    --host-sample HOST_001 \
    --donor-sample DONOR_001 \
    --sample TP1_20240101 \
    --output results.tsv

# Multi-donor (2 donors)
allomix monitor \
    --panel-vcf patient001_panel.vcf.gz \
    --admix-vcf patient001_admix.vcf.gz \
    --host-sample HOST_001 \
    --donor-sample DONOR1_001 \
    --donor-sample DONOR2_001 \
    --sample TP1_20240101 \
    --output results.tsv

# JSON output with per-marker detail
allomix monitor \
    --panel-vcf patient001_panel.vcf.gz \
    --admix-vcf patient001_admix.vcf.gz \
    --host-sample HOST_001 \
    --donor-sample DONOR_001 \
    --sample TP1_20240101 \
    --format json --verbose \
    --output results.json

# Timeline across multiple timepoints (always JSON)
allomix timeline \
    --panel-vcf patient001_panel.vcf.gz \
    --admix-vcf patient001_admix.vcf.gz \
    --host-sample HOST_001 \
    --donor-sample DONOR_001 \
    --sample TP1_20240101 \
    --sample TP2_20240201 \
    --sample TP3_20240301 \
    --output timeline.json

# Estimate bias from per-sample VCFs
allomix estimate-bias \
    --vcfs sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz \
    --output bias_table.tsv

# Estimate bias from named samples within a joint-called VCF
allomix estimate-bias \
    --vcf joint_called.vcf.gz \
    --samples DONOR_001 DONOR_002 DONOR_003 \
    --output bias_table.tsv

# Use bias correction during monitoring
allomix monitor \
    --panel-vcf patient001_panel.vcf.gz \
    --admix-vcf patient001_admix.vcf.gz \
    --host-sample HOST_001 \
    --donor-sample DONOR_001 \
    --sample TP1_20240101 \
    --bias-table bias_table.tsv \
    --output results.tsv
```

If you do not yet have enough donor VCFs to train a bias table, `estimate-bias` can also be driven from archived BAMs on the same panel via a joint-calling pipeline plus sample-level QC. See [Building a training cohort from BAMs](doc/estimate_bias.md#building-a-training-cohort-from-bams) in the bias guide.

### Common Options

Both `monitor` and `timeline` accept these additional options:

| Option | Default | Description |
|---|---|---|
| `--min-dp` | 100 | Minimum read depth to use a marker |
| `--min-gq` | 20 | Minimum genotype quality for host/donor genotyping |
| `--error-rate` | 0.01 | Sequencing error rate for the likelihood model |
| `--bias-table` | none | Per-marker bias table TSV (from `estimate-bias`; see [Bias Estimation Guide](doc/estimate_bias.md)) |
| `--no-bias-correction` | off | Disable bias correction even when a bias table is provided |
| `--verbose` | off | Include per-marker detail in output |

`monitor` also accepts `--format tsv|json` (default: tsv). `timeline` always outputs JSON.

## Input / Output

### Inputs

| Input | Format | Description |
|---|---|---|
| Joint-called VCF | VCF (.vcf/.vcf.gz) | Multi-sample VCF from GATK joint calling containing host, donor(s), and admixture samples. Must contain GT and AD fields. |

All samples for a patient (host, donor(s), and all post-HSCT timepoints) must be joint-called together in a single VCF. Sample names are specified on the command line via `--host-sample`, `--donor-sample`, and `--sample`.

The tool works with VCFs from any variant calling pipeline that supports joint calling (GATK GenomicsDBImport + GenotypeGVCFs) as long as GT and AD fields are present. Higher depth improves sensitivity. Panels with >1000x coverage will give the best results at low chimerism fractions.

### Outputs

| Output | Description |
|---|---|
| % chimerism | Estimated fraction of donor cells (per donor if multi-donor) |
| Confidence interval | 95% CI on the chimerism estimate |
| Per-marker details | Allele depths, expected vs observed VAF, and informativeness flag for each marker |
| QC metrics | Number of informative markers used, mean depth, markers excluded and why, goodness-of-fit |
| Timeline report | Chimerism trend across serial timepoints for a patient |

Output formats: TSV (machine-readable), JSON (for programmatic consumption), and optionally a summary plot.

## Comparison with Commercial Products

| Feature | allomix | AlloSeq HCT (CareDx) | Devyser Chimerism (Thermo Fisher) |
|---|---|---|---|
| Markers | Any bi-allelic panel | 202 SNPs | 24 indels |
| Max genomes | 3 (host + 2 donors) | 3 | 3 |
| Sensitivity | Depends on panel/depth | 0.22% LOD | 0.05% LOD |
| Additional wet-lab | None (uses existing data) | Dedicated kit | Dedicated kit |
| Software | Open-source CLI | Web-based (HCT Software) | Desktop (Advyser) |

## Validation Strategy

Validation follows a two-phase approach:

1. **In silico validation** (current): Synthetic chimeric VCFs with realistic noise models (per-marker bias, depth coefficient of variation, and locus dropout) calibrated from empirical panel data. All experiments use multiple independent replicates (N>=5) with different random seeds to capture sampling variability.
2. **Wet-lab validation** (planned): Real patient samples and controlled dilution series.

## Project Status

This project is under active development. The repository is currently private; an empty package has been published to [PyPI](https://pypi.org/project/allomix/) to reserve the name.

Single-donor and multi-donor (up to 2 donors) chimerism estimation is implemented and validated (in silico):

- MLE-based estimation using Crysup & Woerner (2023) likelihood framework with known genotypes
- Grid search + Brent refinement, profile likelihood 95% confidence intervals
- Multi-donor support with triangular grid search under the f1 + f2 <= 1 constraint
- QC assessment (marker counts, depth, goodness-of-fit, outlier detection)
- TSV and JSON output, including multi-timepoint timeline
- 274 automated tests, in-silico validation at 0-100% donor fractions (RMSE ~0.3%)

**Not yet implemented:** VariantGrid integration.

## Project Structure

```
src/allomix/          # Installable library and CLI, the shipped product
scripts/              # Development and validation utilities
paper/scripts/        # Publication-specific analysis and figures
tests/                # pytest tests
data/                 # De-identified example VCFs
```

**`src/allomix/`** contains everything a user gets when they `pip install allomix`: the core library modules (genotyping, chimerism estimation, simulation, QC, reporting) and the CLI entry point. See the [Architecture Guide](doc/architecture.md) for a module-by-module code map and the data flow through the package.

**`scripts/`** contains developer-facing tools that support building and testing allomix: generating synthetic test data, measuring panel bias from empirical data, and running validation suites. These are not part of the installed package. See the [Scripts Guide](doc/scripts.md) for what each script does and how to run it.

**`paper/scripts/`** contains scripts that produce the specific figures, validation experiments, and statistics for the publication. They use allomix as a library and are intended to make the paper's results fully reproducible.

### Building the Paper

The paper build is orchestrated by Snakemake. All validation and figure scripts run in parallel, then vibepaper renders the final document from the facts they produce.

The paper dependencies (snakemake >=8) require **Python 3.11 or newer**, even though the core tool runs on 3.10+. Create a 3.11 virtual environment for the build:

```bash
uv venv --python '>=3.11'                       # create a venv on Python 3.11 or newer (.venv)
source .venv/bin/activate
uv pip install -e ".[paper]"                   # install paper dependencies (matplotlib, snakemake, vibepaper)
snakemake -s paper/Snakefile --cores $(nproc)              # run all scripts in parallel, then build paper
snakemake -s paper/Snakefile --cores $(nproc) --forceall   # force rerun everything from scratch
snakemake -s paper/Snakefile --cores 1 paper               # just render the paper (assumes facts already exist)
snakemake -s paper/Snakefile clean             # remove all generated output
```

Snakemake tracks file timestamps, so editing a script or its input data reruns only the affected rule and the downstream paper build.

#### Output formats and system dependencies

The build produces a Word document (`output/allomix_paper_<date>.docx`) and rendered Markdown (`.md`). The DOCX step needs **pandoc >= 2.11** (for `--citeproc`). The apt pandoc on some distros is too old (Ubuntu 22.04 ships 2.9.2.1); install a newer one with conda-forge (`conda install -c conda-forge "pandoc>=3"`), the official `.deb` from the [pandoc releases](https://github.com/jgm/pandoc/releases), or `pip install pypandoc-binary` (bundles a modern pandoc, no system install).

The Snakemake build does **not** produce a PDF. The PDF render needs weasyprint's system libraries (pango, cairo, gdk-pixbuf). To make one, install those and run vibepaper yourself:

```bash
# Debian/Ubuntu (other distros: install the equivalent pango/cairo/gdk-pixbuf packages)
sudo apt-get install -y libpango-1.0-0 libpangocairo-1.0-0 libgdk-pixbuf-2.0-0 libcairo2

vibepaper build --md --pdf       # DOCX + Markdown + PDF
```

## License

MIT
