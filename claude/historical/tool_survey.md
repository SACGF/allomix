# Tool Survey: In-House NGS Chimerism Monitoring

## Executive Summary

We surveyed >30 tools across 6 categories (dedicated chimerism, contamination detection, forensic mixture, tumour purity, fingerprinting, and broader methods) to determine whether an existing open-source tool can quantify % donor chimerism from our rhAmpSeq panel data (76 SNPs, >1000x depth, VCFs with AF tags).

**Key finding: No existing open-source tool directly solves our problem.** All commercial chimerism solutions bundle proprietary software with proprietary assay kits. However, several tools from adjacent domains are strong candidates for adaptation, and the underlying math is simple enough that a custom solution is also viable.

**Top 3 candidates:**
1. **Demixtify** — forensic MLE mixture deconvolution, designed for biallelic SNPs, works with targeted panels
2. **EuroForMix** — already validated for 2-person SNP mixtures from NGS read-depth data
3. **Custom script** — the math is straightforward; contamination/tumour purity tools all converge on the same VAF-to-fraction calculation

---

## The Core Math

The chimerism problem is simpler than most tools are designed for:

1. At informative markers (where donor genotype ≠ host genotype), the observed VAF directly encodes the mixture fraction
2. Example: Host=AA, Donor=BB → post-transplant VAF for B allele ≈ donor fraction / 2 (for het donor) or ≈ donor fraction (for hom-alt donor)
3. Average across all informative markers → % chimerism
4. With >1000x depth and 76 SNPs, expect ~15-25 informative sites per pair (literature confirms ≥20 is sufficient)

Both genomes are diploid (no copy number complications), and pre-transplant genotypes are known — this makes it far simpler than tumour purity or forensic blind deconvolution.

---

## Detailed Tool Evaluation

### Category 1: Dedicated NGS Chimerism Tools

All are **commercial/proprietary** — included as benchmarks.

| Tool | Vendor | Markers | Type | LOD | Multi-donor | Software |
|------|--------|---------|------|-----|-------------|----------|
| AlloSeq HCT | CareDx | 202 SNPs, 22 autosomes | SNP | 0.3% | Up to 3 genomes | HCT Software (web) |
| Devyser Chimerism NGS | Thermo Fisher / One Lambda | 24 indels, 17 chromosomes | Indel | 0.05% | Yes | Advyser (desktop) |
| NGStrack + TRKengine | GenDx | 34 indels + X/Y, 18 chromosomes | Indel | 0.1% | Yes | TRKengine |
| ScisGo Chimerism MD | Scisco Genetics | SNPs + indels | Mixed | 0.2% (single), 0.5% (multi) | Yes | GeMS-UI |

**Notable:** Devyser and GenDx use **indels** rather than SNPs, arguing they are less prone to sequencing noise at very low fractions. Our panel uses SNPs, which may limit sensitivity at the lowest chimerism levels.

**Key publications:**
- Kakodkar et al. (2023) Front Genet 14:1282947 — AlloSeq HCT validation
- Blouin et al. (2024) J Mol Diagn 26:995-1006 — ScisGo/AlloSeq comparison
- Pedini et al. (2021) Transplant Cell Ther 27:89.e1 — Devyser validation

---

### Category 2: Open-Source Chimerism-Adjacent Tools

| Tool | What it does | Quantifies mixture? | 76-SNP compatible? | Multi-donor? | Language | GitHub |
|------|-------------|--------------------|--------------------|--------------|----------|--------|
| **Rchimerism** | Automates STR chimerism analysis | Yes (STR only) | No — STR input | Yes (double donor) | R/Shiny | [BioHPC/Rchimerism](https://github.com/BioHPC/Rchimerism) |
| **chimerism_smmip** | Ultrasensitive chimerism via smMIPs | Yes (0.01% LOD) | No — requires smMIP capture | Yes (up to 4 donors) | Perl | [salipante/chimerism_smmip](https://github.com/salipante/chimerism_smmip) |
| **Chimerism-Bias** | Bias correction for NGS chimerism assays | Improves quantification | **Yes** — directly applicable | N/A | R | [mvynck/Chimerism-Bias](https://github.com/mvynck/Chimerism-Bias) |
| **Chimerism-FABCASE** | Assay design: how many markers needed | No — design tool | **Yes** — evaluate our 76-SNP panel | Accounts for relatedness | R/Shiny | [mvynck/Chimerism-FABCASE](https://github.com/mvynck/Chimerism-FABCASE) |
| **Chimerism-nMarkers** | Calculate marker informativity | No — design tool | **Yes** | N/A | R | [mvynck/Chimerism-nMarkers](https://github.com/mvynck/Chimerism-nMarkers) |

**The Vynck et al. tools are essential.** Their bias correction methodology addresses preferential allele observation — a systematic error that directly impacts chimerism accuracy at low fractions. The FABCASE tool can evaluate whether our 76 SNPs provide sufficient informativity. Key papers:
- Vynck et al. (2022) Clin Chim Acta 532:123-129 — How many markers?
- Vynck et al. (2023) Clin Chim Acta — Bias reduction
- Vynck et al. (2025) Int J Lab Hematol — FABCASE

---

### Category 3: Forensic Mixture Deconvolution Tools

| Tool | Mixture quantification? | SNP/NGS support? | 76-SNP compatible? | Open source? | Chimerism suitability |
|------|------------------------|------------------|--------------------|--------------|-----------------------|
| **EuroForMix** | Yes — MLE mixture weights | **Yes — validated for 2-3 person SNP mixtures from NGS** | Yes — custom kit definitions | Yes (LGPL-3.0, R/C++) | **HIGH** |
| **MixDeR** | Yes (via EuroForMix) | Yes — ForenSeq SNP data | Needs panel adapter | Yes | MODERATE (workflow template) |
| **Demixtify** | Yes — MLE mixture fraction | Yes — biallelic SNPs, targeted panels | **Yes** | Yes (AGPL-3.0) | **HIGH** |
| **DNAmixtures/Lite** | Yes — peak height model | Primarily STR | Needs adaptation | Partially (HUGIN dependency) | LOW |
| **STRmix NGS** | Yes | STR + some SNP | No | No (commercial) | LOW |
| **K-means FMAR** | Yes — clustering approach | Yes — 90 autosomal SNPs | Conceptually yes | Unclear (paper supplementary?) | MODERATE |

#### EuroForMix — Top Candidate

Already validated for exactly our scenario: 2-person SNP mixtures from NGS read-depth data (PMID 28942111, using 134 autosomal SNPs from HID-Ion AmpliSeq). The quantitative read-depth model massively outperformed qualitative methods.

To adapt for chimerism:
1. Format our VCF allele counts into EuroForMix's evidence/reference CSV format
2. Define a custom "kit" for our 76 SNP markers
3. Use the estimated mixture weight as the chimerism percentage
4. MixDeR demonstrates how to build a preprocessing pipeline around EuroForMix

GitHub: [oyvble/euroformix](https://github.com/oyvble/euroformix) | Web: euroformix.com

#### Demixtify — Top Candidate

MLE of mixture fraction from biallelic SNP read counts — mathematically identical to chimerism quantification. Explicitly supports targeted panels (ships with WES panel files). Handles extreme ratios (≤1:99). Includes empirical sequence error estimation.

Key advantage over EuroForMix: simpler, more focused — does mixture fraction estimation without the full probabilistic genotyping framework.

GitHub: [Ahhgust/Demixtify](https://github.com/Ahhgust/Demixtify)
Publication: Vohr et al. 2023, Forensic Sci Int Genet (PMID 38016331)

---

### Category 4: Contamination Detection Tools

| Tool | Quantifies fraction? | Min markers | 76-SNP feasible? | Sensitivity | Open source? | Chimerism suitability |
|------|---------------------|-------------|------------------|-------------|--------------|----------------------|
| **VerifyBamID2** | Yes (FREEMIX) | 5,000-10,000 | **No** (MSE=0.69 at 1000 markers, 1% contam.) | ~1% | Yes (BSD) | LOW — assumes random contaminant |
| **Conpair** | Yes | ~7,000 (ships with) | Uncertain | 0.1% | Yes (CC-BY-4.0) | MODERATE — paired model fits |
| **ART-DeCo** | Yes | **≥30 SNPs** (MAF 30-70%) | **Yes** | ~1% | Yes (SourceForge) | MODERATE — panel-compatible |
| **ContEst/GATK4** | Yes | Flexible | Moderate | <0.1% (ContEst) | Yes (BSD-3) | MODERATE |
| **CHARR** | Yes | ≥500 hom-alt | **No** | 0.5% | Yes | LOW |
| **NGSCheckMate** | No (binary) | N/A | N/A | N/A | Yes | NONE |
| **NGSTroubleFinder** | Yes (ML regression) | 164,767 | **No** | ~1% | Yes (MIT) | LOW — fixed large panel |

#### ART-DeCo — Compatible with Our Panel

Explicitly designed for targeted gene panels. Requires ≥30 SNPs with MAF 30-70% at ≥200x depth — our 76 SNPs at >1000x exceed all requirements. Provides quantitative contamination percentage per sample. Two-stage approach: screening then contaminant identification.

Limitation: sensitivity floor of ~1%, and assumes contaminant is another sample in the same sequencing run.

GitHub: [SourceForge](https://sourceforge.net/projects/ngs-art-deco/) | PMID: 30683922

#### Conpair — Best Head-to-Head Performance

Ranked #1 in a 2024 comparison of 9 contamination tools (PMID: 38479675). Its paired-sample model (tumour=chimeric sample, normal=pre-transplant host) maps well to chimerism. Detects contamination as low as 0.1%.

Concern: ships with ~7,000 markers; unclear how well it performs with only 76.

GitHub: [nygenome/Conpair](https://github.com/nygenome/Conpair)

---

### Category 5: Tumour Purity Estimation Tools

| Tool | VAF-based? | Targeted panel support? | Open source? | Chimerism suitability |
|------|-----------|------------------------|--------------|----------------------|
| **All-FIT** | Yes — weighted least squares | **Yes** — designed for high-depth targeted | Yes (Python) | **MODERATE-HIGH** |
| **PureCN** | Yes | Yes (Bioconductor) | Yes (R) | MODERATE — deeply coupled with CNV |
| **TPES** | Yes — TP = observed_VAF / expected_VAF | No — needs WES/WGS | Yes (R, CRAN) | MODERATE (formula transferable) |
| **ABSOLUTE** | Yes | No — needs genome-wide CNV | Yes (R) | LOW |
| **PurityEst** | Yes | No — needs WGS | Perl | LOW |
| **AbsCN-seq** | Partially | No | R | LOW |

#### All-FIT — Simplest Adaptable Codebase

Small Python tool: takes VAF + depth + ploidy per variant, estimates mixture fraction using iterative weighted least squares. For chimerism, all loci are diploid (ploidy=2), simplifying the model. The somatic mutation model would need to be replaced with a germline informative-SNP model, but the codebase is small and modifiable.

GitHub: [KhiabanianLab/All-FIT](https://github.com/KhiabanianLab/All-FIT) | PMC7141867

---

### Category 6: Fingerprinting / Relatedness Tools

| Tool | Quantifies mixture? | Relevance |
|------|---------------------|-----------|
| **Somalier** | Experimental (commented-out contamination estimation) | LOW-MOD — genotype extraction useful, mixture quantification not supported |
| **Peddy** | No (idr_baf proxy only) | LOW |
| **CrosscheckFingerprints** (Picard) | No (LOD score, binary) | LOW |
| **bcftools gtcheck** | No (concordance only) | VERY LOW |
| **read_haps** | Detection only, needs dense SNP clusters | LOW — 76 spread-out SNPs incompatible |

None of these are suitable for chimerism quantification, but somalier's genotype extraction at polymorphic sites could feed into a custom pipeline.

---

## Additional Key Publications

| Reference | Year | Key Contribution |
|-----------|------|-----------------|
| Lee et al., J Clin Med 8:2077 | 2019 | 121-SNP chimerism algorithm piggybacking on existing MDS mutation panel |
| Aloisio et al., Mol Med Rep 14:2967 | 2016 | Custom 44-amplicon SNP panel, 1% LOD, R²=0.999 vs STR |
| Wu et al., Clin Chem 64:938 | 2018 | smMIP chimerism, 0.01% LOD, 4-donor deconvolution |
| Wu et al., J Mol Diagn 24:167 | 2022 | Improved smMIP with copy-number neutral control loci |
| Caulier et al., Front Immunol (2023) | 2023 | Comprehensive review comparing all NGS chimerism methods |
| PMID 35802296 | 2022 | SNP-NGS chimerism with sensitivity 0.01-0.05%, only 8-200ng DNA |

---

## Recommendation

### Path 1: Adapt Demixtify (Recommended — Fastest to Production)

**Demixtify** is purpose-built for quantitative 2-person mixture estimation from biallelic SNPs using MLE. It already:
- Supports targeted panels
- Handles extreme ratios (≤1:99)
- Includes empirical sequence error estimation
- Takes BAM + BCF panel input

**Work needed:**
- Test with our 76-SNP panel VCFs
- Validate against known mixture samples
- Add multi-donor support (currently 2-person)
- Build a wrapper for serial timepoint tracking and reporting

**Concern:** AGPL-3.0 license requires consideration for clinical deployment.

### Path 2: Adapt EuroForMix (Best Statistical Framework)

**EuroForMix** is the most statistically rigorous option, already validated for 2-3 person SNP mixtures from NGS. It provides:
- Full probabilistic framework with confidence intervals
- Validated for 134 autosomal SNPs from NGS (close to our 76)
- Mixture weight estimation = chimerism %

**Work needed:**
- Format VCF data into EuroForMix evidence/reference CSVs
- Define custom kit for our 76 markers
- Build preprocessing pipeline (MixDeR is a template)
- Validate against known mixtures

### Path 3: Custom Script (Simplest, Most Maintainable)

Given the mathematical simplicity (average VAF at informative markers → % chimerism), a custom pipeline may be the most practical:

1. **Genotype comparison** — identify informative markers where donor/host differ
2. **VAF extraction** — pull AF from post-transplant VCFs at informative sites
3. **Chimerism calculation** — convert VAF to donor % based on genotype configuration
4. **Bias correction** — apply Vynck et al. methods
5. **Reporting** — track across timepoints

**Advantages:** Full control, no license concerns, integrates with existing VG/TAU infrastructure, easy to audit for GMP.

**Use FABCASE** to evaluate whether 76 SNPs provide sufficient informativity, and **Chimerism-Bias** to correct systematic errors.

### Path 4: Increase Panel Size

If 76 SNPs proves insufficient for <1% sensitivity (especially with related donors), consider:
- Adding more common SNPs to the capture panel (the IDT xGen panel targets 76, but custom panels can include hundreds)
- Commercial panels use 24-202 markers; literature suggests ≥20 informative markers is sufficient, but more is better for related donors
- The FABCASE tool can model exactly how many markers are needed for a given donor-recipient relatedness level

---

## Assessment Matrix

| Criterion | Demixtify | EuroForMix | Custom Script | ART-DeCo | All-FIT |
|-----------|-----------|------------|---------------|----------|---------|
| Works with 76-SNP VCFs | Yes (with BAM) | Yes | Yes | Yes | Yes (with adaptation) |
| Quantifies mixture fraction | Yes (MLE) | Yes (MLE) | Yes | Yes | Yes (IWLS) |
| Multi-donor (3 genomes) | No (2-person) | Yes (up to 6) | Yes (custom) | No | No |
| Sensitivity <1% | Yes (≤1:99) | Yes (validated) | Depends on implementation | ~1% floor | Depends |
| Serial monitoring / trending | No | No | Yes (custom) | No | No |
| Genotype database | No | No | Yes (via VG) | No | No |
| Integration effort | Medium | Medium | Medium-High | Low | Medium |
| Actively maintained | Yes (2023) | Yes (active) | N/A | 2019 | 2020 |
| License | AGPL-3.0 | LGPL-3.0 | N/A | Unclear | Open |
| Statistical rigour | High | Highest | Depends | Moderate | Moderate |

---

## Conclusion

**No off-the-shelf open-source tool exists for general-purpose NGS chimerism from VCF/BAM data.** This is a genuine gap in the field — all commercial solutions are proprietary, and the only open-source NGS chimerism pipeline (smMIP) requires a specialized capture method.

The **recommended approach** is a hybrid:
1. Use **FABCASE** to validate that our 76-SNP panel has sufficient informativity
2. Use **Chimerism-Bias** methods to correct systematic allele observation biases
3. Build the chimerism calculation using either **Demixtify's MLE framework** (adapted for known donor/host genotypes) or a **custom implementation** of the same math
4. Integrate with **VariantGrid** for genotype storage and matching
5. Build timepoint tracking and reporting on top

The math is well-understood and simple for our use case. The harder engineering problems are: (a) robust genotype database management (the GMP constraint), (b) handling related donors (fewer informative markers), and (c) achieving <1% sensitivity with only 76 SNPs (may require panel expansion).
