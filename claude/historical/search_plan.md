# Search Plan: In-House NGS Chimerism Monitoring Tools

## Goal Recap

We want to calculate **% donor chimerism** from post-HSCT admixture samples using our existing rhAmpSeq/capture panel data (76 SNPs, >1000x depth, VCFs with AF/VAF tags). The core math is:

1. Genotype donor and host individually at informative SNP markers (where their genotypes differ)
2. In the post-HSCT mixed sample, measure VAF at those informative markers
3. Convert VAF to % chimerism (e.g., if host is AA and donor is BB, a post-transplant VAF of 0.03 for B implies ~3% donor)
4. Track across serial timepoints

This is essentially **quantifying a DNA mixture at known distinguishing sites** — the same problem addressed by contamination detection tools, forensic mixture analysis, and tumour purity estimation.

## Search Strategy

### Phase 1: Dedicated NGS Chimerism Tools

Search for tools explicitly designed for chimerism from NGS data.

**Web/GitHub searches:**
- `NGS chimerism monitoring tool software`
- `chimerism SNP NGS bioinformatics`
- `hematopoietic stem cell transplant chimerism NGS analysis software`
- `donor chimerism calculation from VCF`
- `chimerism from amplicon sequencing`

**PubMed searches:**
- `chimerism NGS software` 
- `chimerism next-generation sequencing bioinformatics`
- `"chimerism monitoring" AND ("next-generation sequencing" OR NGS) AND (software OR tool OR pipeline OR method)`
- `"donor chimerism" AND SNP AND (algorithm OR calculation OR quantification)`
- `"mixed chimerism" AND "variant allele frequency"`

### Phase 2: Sample Contamination / Mixture Quantification Tools

Contamination detection tools solve a closely related problem: quantifying the fraction of a second genome mixed into a sample. The math is nearly identical to chimerism.

**Tools to investigate:**
- **VerifyBamID / VerifyBamID2** — estimates contamination fraction from BAM files using population allele frequencies; could potentially be repurposed
- **ContEst** (Broad/GATK) — contamination estimation in tumour sequencing
- **Conpair** — concordance and contamination estimator for tumour-normal pairs
- **ART-DeCo** — contamination detection tool
- **CHARR** (gnomAD) — contamination estimate from sequencing data
- **NGSCheckMate** — sample identity verification via correlation of VAFs at polymorphic sites

**Web/GitHub searches:**
- `sample contamination quantification NGS tool`
- `DNA mixture fraction estimation bioinformatics`
- `cross-sample contamination detection sequencing`
- `verify sample identity NGS VAF`

**PubMed searches:**
- `"sample contamination" AND "next-generation sequencing" AND (quantification OR estimation OR fraction)`
- `"DNA mixture" AND "allele frequency" AND (deconvolution OR quantification)`
- `contamination estimation sequencing genotype`

### Phase 3: Forensic / Mixture Deconvolution Tools

Forensic genetics has extensive literature on DNA mixture interpretation, which is mathematically the same problem.

**Tools to investigate:**
- **EuroForMix** — continuous DNA mixture interpretation
- **STRmix** — probabilistic genotyping (STR-based but methodology applicable)
- **DNAmixtures** — R package for DNA mixture analysis

**PubMed searches:**
- `"DNA mixture" AND SNP AND deconvolution AND "next-generation sequencing"`
- `forensic mixture interpretation NGS SNP`
- `"mixture proportion" AND SNP AND sequencing`

### Phase 4: Tumour Purity / Clonality Tools (Analogous Math)

Tumour purity estimation from VAF data solves a very similar problem (what fraction of cells are tumour vs normal?).

**Tools to investigate:**
- **PureCN** — tumour purity estimation from targeted sequencing
- **ABSOLUTE** — tumour purity and ploidy
- **THetA** — tumour heterogeneity analysis
- **PyClone** — clonal population structure

**PubMed searches:**
- `"tumour purity" AND "variant allele frequency" AND "targeted sequencing"`
- `"allele frequency" AND "mixture proportion" AND estimation AND sequencing`

### Phase 5: Fingerprinting / Relatedness Tools (for the Genotype Matching Component)

Even if these don't do chimerism calculation, they solve the genotype storage/matching problem.

**Tools to investigate:**
- **Somalier** — already identified; fingerprinting and relatedness
- **Peddy** — pedigree/sex/ancestry checks
- **CrosscheckFingerprints** (Picard) — sample identity via genotype concordance
- **GTcheck** (bcftools) — genotype concordance checking

**Web/GitHub searches:**
- `sample fingerprinting genotype database tool`
- `donor recipient genotype matching transplant`

### Phase 6: Broader / Creative Searches

Catch things that don't fit neatly into the above categories.

**PubMed searches:**
- `"engraftment monitoring" AND NGS AND (tool OR software OR pipeline)`
- `"microchimerism" AND "next-generation sequencing" AND detection`
- `"post-transplant monitoring" AND SNP AND "allele frequency"`
- `chimerism AND (amplicon OR "targeted sequencing" OR "capture sequencing") AND quantification`
- `"minimal residual" AND chimerism AND NGS` (MRD and chimerism are often discussed together)
- `"lineage-specific chimerism" AND NGS`

**Web searches:**
- `open source chimerism analysis software NGS 2023 2024`
- `bioconda chimerism tool`
- `chimerism NGS pipeline github`

## Evaluation Criteria

For each tool found, assess:

1. **Input compatibility** — Can it work with our existing VCFs (76 SNPs, AF tag, >1000x depth)? Or does it need BAM/CRAM?
2. **Mixture quantification** — Does it calculate a mixture fraction / % chimerism, or only binary identity?
3. **Multi-donor support** — Can it handle host + 2 donors (3-way mixture)?
4. **Sensitivity** — Can it detect <1% minority fraction? (our requirement)
5. **Serial monitoring** — Does it support tracking over timepoints, or is it single-sample?
6. **Genotype database** — Does it store/retrieve donor-host genotype profiles?
7. **Ease of integration** — How much custom development would be needed to plug it into our workflow?
8. **Actively maintained** — Is it still supported?

## Expected Deliverable

A write-up (`claude/tool_survey.md`) with:
- Table of all tools found, scored against the evaluation criteria
- Recommendation on which tool(s) to adopt or adapt
- Assessment of whether a simple custom script (VAF at informative markers → % chimerism) might be simpler than adapting an existing tool
- Estimate of what additional markers (beyond the 76) might be needed for robust chimerism at <1% sensitivity
