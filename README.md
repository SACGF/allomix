# allomix

NGS-based donor chimerism monitoring for hematopoietic stem cell transplantation (HSCT).

## Overview

**allomix** calculates donor chimerism percentages from NGS data, replacing traditional STR-based analysis with a higher-sensitivity SNP-based approach. It is designed for clinical laboratories monitoring engraftment and relapse after HSCT.

The tool uses the 76 bi-allelic SNPs from the IDT xGen Human ID Hybridization Capture Panel (rhAmpSeq Sample ID), which are already sequenced as part of the standard haematology capture panel at >1000x depth. No additional wet-lab work is required.

## Clinical Context

After HSCT, patients carry a mixture of their own (host) and transplanted (donor) cells. Monitoring the ratio of donor to host cells over time — **chimerism monitoring** — is critical for detecting graft rejection or disease relapse. Early detection allows early intervention.

Current STR-based methods have limited sensitivity (~3-5% LOD) and require a separate workflow outside GMP. allomix aims to:

- Achieve **<1% sensitivity** for detecting minority cell populations
- Support **up to 3 genomes** (host + 2 donors) for patients with multiple transplants
- Run entirely **within GMP**, eliminating the risk of donor/host sample mix-ups during data analysis
- Provide **timeline tracking** of chimerism across serial post-HSCT timepoints

## Workflow

```
1. GENOTYPE    Sequence host and each donor individually
                 → VCF with genotypes at 76 SNP markers

2. COMPARE     Identify informative markers (sites where
               donor and host genotypes differ)

3. MONITOR     Sequence post-HSCT admixture samples at
               serial timepoints (≥3 per patient)
                 → VCF with allele depths at 76 markers

4. CALCULATE   At each informative marker, use observed
               allele frequencies to estimate the fraction
               of donor vs host DNA

5. REPORT      Output % chimerism with confidence intervals,
               per-marker details, and QC metrics
```

## Input / Output

### Inputs

| Input | Format | Description |
|---|---|---|
| Host genotype | VCF (.vcf.gz) | Per-sample GATK-called VCF at 76 SNP loci. Must contain GT and AD fields. |
| Donor genotype(s) | VCF (.vcf.gz) | One VCF per donor (up to 2 donors). Same format as host. |
| Admixture sample(s) | VCF (.vcf.gz) | Post-HSCT monitoring samples. One or more timepoints. Must contain AD (allele depth) fields. |

VCF files are expected to be produced by the existing GATK per-sample calling pipeline from rhAmpSeq amplicon data, with depth >>1000x. The AF tag (VAF) should be present.

### Outputs

| Output | Description |
|---|---|
| % chimerism | Estimated fraction of donor cells (per donor if multi-donor) |
| Confidence interval | 95% CI on the chimerism estimate |
| Per-marker details | Allele depths, expected vs observed VAF, and informativeness flag for each of the 76 markers |
| QC metrics | Number of informative markers used, mean depth, markers excluded and why, goodness-of-fit |
| Timeline report | Chimerism trend across serial timepoints for a patient |

Output formats: TSV (machine-readable), JSON (for programmatic consumption / VariantGrid integration), and optionally a summary plot.

## SNP Panel

The 76 SNPs from the IDT xGen Human ID Hybridization Capture Panel provide:

- **Discrimination power**: >1 in 5 million individuals
- **Bi-allelic markers**: Simplifies genotype comparison and mixture modelling
- **High depth**: >1000x coverage from the existing haematology capture panel
- **Population independence**: Markers selected for low population bias (low FST)

Panel details: [IDT xGen Human ID Hyb Cap Panel](https://sg.idtdna.com/pages/products/next-generation-sequencing/workflow/xgen-ngs-hybridization-capture/pre-designed-hyb-cap-panels/human-id-hyb-panel)

## Installation

```bash
pip install allomix
```

For development:

```bash
git clone https://github.com/dlawrence/allomix.git
cd allomix
pip install -e ".[dev]"
```

## Usage

```bash
# Calculate chimerism for a single timepoint
allomix monitor \
    --host host_genotype.vcf.gz \
    --donor donor_genotype.vcf.gz \
    --sample post_hsct_day30.vcf.gz \
    --output results.tsv

# Multi-donor (2 donors)
allomix monitor \
    --host host_genotype.vcf.gz \
    --donor donor1_genotype.vcf.gz \
    --donor donor2_genotype.vcf.gz \
    --sample post_hsct_day30.vcf.gz \
    --output results.tsv

# Timeline across multiple timepoints
allomix timeline \
    --host host_genotype.vcf.gz \
    --donor donor_genotype.vcf.gz \
    --sample day30.vcf.gz \
    --sample day60.vcf.gz \
    --sample day90.vcf.gz \
    --output timeline.json
```

## Comparison with Commercial Products

| Feature | allomix | AlloSeq HCT (CareDx) | Devyser Chimerism (Thermo Fisher) |
|---|---|---|---|
| Markers | 76 SNPs | 202 SNPs | 24 indels |
| Marker type | Bi-allelic SNPs | Bi-allelic SNPs | Indels |
| Max genomes | 3 (host + 2 donors) | 3 | 3 |
| Sensitivity target | <1% | 0.22% LOD | 0.05% LOD |
| Additional wet-lab | None (uses existing panel) | Dedicated kit | Dedicated kit |
| Software | Open-source CLI + VariantGrid | Web-based (HCT Software) | Desktop (Advyser) |
| GMP compatible | Yes (by design) | Yes | Yes |

## Project Status

This project is in early development (v0.0.1). The core chimerism calculation algorithm is not yet implemented.

## License

MIT
