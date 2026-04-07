# allomix

NGS-based donor chimerism monitoring for hematopoietic stem cell transplantation (HSCT).

## Overview

**allomix** calculates donor chimerism percentages from NGS data, replacing traditional STR-based analysis with a higher-sensitivity SNP-based approach. It works with any panel of bi-allelic markers (SNPs or indels) and is designed for clinical laboratories monitoring engraftment and relapse after HSCT.

The tool is panel-agnostic — it operates on whatever markers are present in the input VCFs. Bring your own panel: whether that's 24 indels, 76 SNPs, 202 SNPs, or any other set of bi-allelic loci with sufficient depth.

## Clinical Context

After HSCT, patients carry a mixture of their own (host) and transplanted (donor) cells. Monitoring the ratio of donor to host cells over time — **chimerism monitoring** — is critical for detecting graft rejection or disease relapse. Early detection allows early intervention.

Current STR-based methods have limited sensitivity (~3-5% LOD) and require separate workflows. allomix aims to:

- Achieve **<1% sensitivity** for detecting minority cell populations
- Support **up to 3 genomes** (host + 2 donors) for patients with multiple transplants
- Provide **timeline tracking** of chimerism across serial post-HSCT timepoints

## Workflow

```
1. GENOTYPE    Sequence host and each donor individually
                 → VCF with genotypes at marker loci

2. COMPARE     Identify informative markers (sites where
               donor and host genotypes differ)

3. MONITOR     Sequence post-HSCT admixture samples at
               serial timepoints (≥3 per patient)
                 → VCF with allele depths at marker loci

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
| Host genotype | VCF (.vcf.gz) | Per-sample VCF at marker loci. Must contain GT and AD fields. |
| Donor genotype(s) | VCF (.vcf.gz) | One VCF per donor (up to 2 donors). Same format as host. |
| Admixture sample(s) | VCF (.vcf.gz) | Post-HSCT monitoring samples. One or more timepoints. Must contain AD (allele depth) fields. |

The tool works with VCFs from any variant calling pipeline (GATK, DeepVariant, etc.) as long as GT and AD fields are present. Higher depth improves sensitivity — panels with >1000x coverage will give the best results at low chimerism fractions.

### Outputs

| Output | Description |
|---|---|
| % chimerism | Estimated fraction of donor cells (per donor if multi-donor) |
| Confidence interval | 95% CI on the chimerism estimate |
| Per-marker details | Allele depths, expected vs observed VAF, and informativeness flag for each marker |
| QC metrics | Number of informative markers used, mean depth, markers excluded and why, goodness-of-fit |
| Timeline report | Chimerism trend across serial timepoints for a patient |

Output formats: TSV (machine-readable), JSON (for programmatic consumption), and optionally a summary plot.

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
| Markers | Any bi-allelic panel | 202 SNPs | 24 indels |
| Max genomes | 3 (host + 2 donors) | 3 | 3 |
| Sensitivity | Depends on panel/depth | 0.22% LOD | 0.05% LOD |
| Additional wet-lab | None (uses existing data) | Dedicated kit | Dedicated kit |
| Software | Open-source CLI | Web-based (HCT Software) | Desktop (Advyser) |

## Project Status

Single-donor chimerism estimation is implemented and validated:

- MLE-based estimation using Crysup & Woerner (2022) likelihood framework with known genotypes
- Grid search + Brent refinement, profile likelihood 95% confidence intervals
- QC assessment (marker counts, depth, goodness-of-fit, outlier detection)
- TSV and JSON output, including multi-timepoint timeline
- 192 automated tests, in-silico validation at 0–100% donor fractions (RMSE ~0.3%)

**Not yet implemented:** multi-donor support (host + 2 donors), per-marker bias correction (Vynck et al.), VariantGrid integration.

## License

MIT

Chimerism MLE methodology based on Crysup & Woerner (2022) — independently reimplemented, not derived from Demixtify source code.
