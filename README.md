# allomix

[![PyPi version](https://img.shields.io/pypi/v/allomix.svg)](https://pypi.org/project/allomix/) [![Python versions](https://img.shields.io/pypi/pyversions/allomix.svg)](https://pypi.org/project/allomix/)
[![tests](https://github.com/SACGF/allomix/actions/workflows/tests.yml/badge.svg)](https://github.com/SACGF/allomix/actions/workflows/tests.yml)

NGS-based donor chimerism monitoring for hematopoietic stem cell transplantation (HSCT).

**allomix** calculates donor chimerism percentages from NGS data, replacing traditional STR-based analysis with a higher-sensitivity SNP-based approach. It is panel-agnostic: it operates on whatever bi-allelic markers (SNPs or indels) are present in the input VCFs, whether that is 24 indels, 76 SNPs, 202 SNPs, or any other set of loci with sufficient depth.

> **Results are highly panel specific: do your own validation.** Sensitivity and limit of detection depend on your marker set, sequencing depth, and noise profile. Qualify the tool on your own panel before clinical use (see the [Panel guide](https://github.com/SACGF/allomix/blob/main/docs/panel_guide.md)).

Chimerism MLE methodology is based on [Crysup & Woerner](https://pubmed.ncbi.nlm.nih.gov/36152508/) (2022).

## Clinical context

After HSCT, patients carry a mixture of their own (host) and transplanted (donor) cells. Monitoring the donor-to-host ratio over time detects graft rejection or disease relapse early enough to intervene. Current STR-based methods have limited sensitivity (~3-5% LOD) and require separate workflows. allomix aims to:

- achieve **<1% sensitivity** for detecting minority cell populations,
- support **up to 3 genomes** (host + 2 donors) for patients with multiple transplants, and
- provide **timeline tracking** of chimerism across serial post-HSCT timepoints.

## Installation

```bash
pip install allomix
```

For development:

```bash
git clone https://github.com/SACGF/allomix.git
cd allomix
uv pip install -e ".[dev]"
```

## Quickstart

allomix takes two VCFs: a **panel VCF** with host/donor genotypes (from GATK joint calling of the reference samples) and an **admix VCF** with per-timepoint AD counts (from forced `bcftools mpileup` at the panel sites). A ready-to-use Snakemake pipeline that produces both is in `pipeline/`.

```bash
allomix detect \
    --genotype-vcf patient001_panel.vcf.gz \
    --admix-vcf patient001_admix.vcf.gz \
    --host-sample HOST_001 \
    --donor-sample DONOR_001 \
    --sample TP1_20240101 \
    --html report.html
```

See the [CLI usage guide](https://github.com/SACGF/allomix/blob/main/docs/cli.md) for multi-donor runs, timelines, bias correction, output options, and input/output reference.

## Workflow

```
1. Genotyping            Sequence host and each donor individually
   (upstream)              → per-sample GVCFs at marker loci

2. estimate-bias         (optional) Estimate per-marker amplification
                         bias from genotyping VCFs → bias table TSV

3. Sequencing            Sequence post-HSCT admixture samples at serial
   (upstream)            timepoints (>=3 per patient) → per-sample GVCFs

4. Joint calling         Combine HOST + DONOR GVCFs (GenomicsDBImport +
   (upstream)            GenotypeGVCFs) → panel VCF; pileup admix samples
                         at the panel sites → admix VCF

5. allomix detect       Chimerism for one timepoint → TSV / JSON / HTML / PDF
   allomix timeline      Track chimerism across timepoints → JSON / HTML / PDF
   allomix report        Render HTML/PDF from a saved detect/timeline JSON
```

Joint calling of HOST + DONOR propagates donor ALT alleles to the panel even when one sample is hom-ref; pileup of the admix samples preserves raw per-allele counts needed for host fractions below ~5%. See the [Joint Calling Guide](https://github.com/SACGF/allomix/blob/main/docs/joint_calling.md) for the full rationale.

## Documentation

- [Panel guide](https://github.com/SACGF/allomix/blob/main/docs/panel_guide.md) — qualifying your own panel for chimerism use (start here for a new panel)
- [Marker types](https://github.com/SACGF/allomix/blob/main/docs/marker_types.md) — how allomix classifies markers and what each class is used for
- [CLI usage](https://github.com/SACGF/allomix/blob/main/docs/cli.md) — all subcommands, options, and input/output reference
- [Reports and structured output](https://github.com/SACGF/allomix/blob/main/docs/reports.md) — the JSON envelope, HTML/PDF report, and worked examples ([live example reports](https://sacgf.github.io/allomix/))
- [Joint Calling Guide](https://github.com/SACGF/allomix/blob/main/docs/joint_calling.md) — two-phase upstream pipeline and rationale
- [Bias Estimation Guide](https://github.com/SACGF/allomix/blob/main/docs/estimate_bias.md) — per-marker bias tables (and building a training cohort from BAMs)
- [Custom report templates](https://github.com/SACGF/allomix/blob/main/docs/custom_report_template.md) — branding the HTML/PDF report for your lab
- [Architecture](https://github.com/SACGF/allomix/blob/main/docs/architecture.md) — module-by-module code map and data flow
- [Scripts](https://github.com/SACGF/allomix/blob/main/docs/scripts.md) — developer and validation utilities
- [Building the paper](https://github.com/SACGF/allomix/blob/main/docs/paper.md) — Snakemake validation and figure build

## Comparison with commercial products

| Feature | allomix | AlloSeq HCT (CareDx) | Devyser Chimerism (Thermo Fisher) |
|---|---|---|---|
| Markers | Any bi-allelic panel | 202 SNPs | 24 indels |
| Max genomes | 3 (host + 2 donors) | 3 | 3 |
| Sensitivity | Depends on panel/depth | 0.22% LOD | 0.05% LOD |
| Additional wet-lab | None (uses existing data) | Dedicated kit | Dedicated kit |
| Software | Open-source CLI | Web-based (HCT Software) | Desktop (Advyser) |

## Validation and status

allomix has been validated in silico (synthetic chimeric VCFs with realistic noise models: per-marker bias, depth CV, locus dropout) and on real reads from a public dataset of titrated DNA mixtures (SRA study [SRP434573](https://www.ebi.ac.uk/ena/browser/view/PRJNA960854)). On the real mixtures it recovered known host fractions from 10% down to 1%, resolved a three-person mixture, and called residual host with no false positives on the pure-donor controls. Full validation, including the real-data limit of detection, is in the [paper build guide](https://github.com/SACGF/allomix/blob/main/docs/paper.md).

These are analytical bounds, not wet-lab limits. Wet-lab validation against STR chimerism on real patient samples is planned, and is a per-laboratory step for any new panel.

This project is under active development.

## Project structure

```
src/allomix/          # Installable library and CLI, the shipped product
scripts/              # Development and validation utilities
paper/scripts/        # Publication-specific analysis and figures
tests/                # pytest tests
```

`src/allomix/` is everything a user gets from `pip install allomix`: the core library (genotyping, chimerism estimation, simulation, QC, reporting) and the CLI entry point. `scripts/` and `paper/scripts/` are developer-facing and not part of the installed package. See the [Architecture Guide](https://github.com/SACGF/allomix/blob/main/docs/architecture.md), [Scripts Guide](https://github.com/SACGF/allomix/blob/main/docs/scripts.md), and [paper build guide](https://github.com/SACGF/allomix/blob/main/docs/paper.md).

## License

allomix is distributed under the MIT licence, and comes with no warranty. Results are highly panel specific: do your own validation.
