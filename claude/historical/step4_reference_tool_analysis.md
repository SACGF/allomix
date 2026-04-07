# Step 4: Reference Tool Analysis for allomix

Deep analysis of open-source tools relevant to NGS-based chimerism monitoring.
Research conducted 2026-04-07 by examining GitHub repositories (READMEs, source code, algorithms).

---

## 1. mvynck/Chimerism-Bias

**Repository**: https://github.com/mvynck/Chimerism-Bias
**Language**: R | **Stars**: 0 | **License**: No explicit license file (academic research code)

### Core Algorithm

The tool implements **per-marker allelic bias correction** for NGS chimerism measurement. It does NOT estimate mixture fractions via MLE -- instead it computes chimerism as a simple per-marker algebraic calculation from VAF, then corrects for systematic amplification bias.

The key insight is classifying markers into 6 types based on donor/recipient genotype constellations:

| Type | Donor | Recipient | Formula (raw host %) |
|------|-------|-----------|---------------------|
| 0 | AA | aa | 1 - VAF |
| 1 | aa | AA | VAF |
| 10 | AA | Aa | 2*(1 - VAF) |
| 11 | aa | Aa | 2*VAF |
| 20 | Aa | AA | 2*(VAF - 0.5) |
| 21 | Aa | aa | 2*(0.5 - VAF) |

Types 0/1 are "fully informative" (both homozygous, different alleles). Types 10/11 use markers where the donor is homozygous and recipient is heterozygous. Types 20/21 are "potentially informative" -- donor het, recipient hom.

**Bias correction**: For each marker, bias `b` is estimated as the median deviation of heterozygous VAF from 0.5 across a training set of samples. Then the raw host percentage is corrected analytically:
- For type 0: `HC_corrected = (HC - 2*HC*b) / (-4*b*HC + 2*b + 1)`
- Similar closed-form expressions for each type

This is a **first-order correction** for systematic reference/alt allele amplification imbalance at each locus.

### Input Format

Tab-delimited text files with columns: marker_id, total_reads, ref_count, alt_count. One file each for recipient screening, donor screening, and follow-up samples.

### Error Model

- **Amplification bias**: Explicitly modeled and corrected per marker using median heterozygous deviation from 0.5
- **Sequencing error**: Not explicitly modeled
- **Stochastic sampling**: Analyzed as binomial variance `Var(VAF) = p(1-p)/n`; theoretical 95% quantile `= 1.96 * sqrt(4*0.5*0.5/depth)` matches observations well
- **Run-specific bias**: Investigated via ANOVA (marker + run effects); found marker effects dominate

### Edge Cases

- **No informative markers**: Returns NA
- **Extreme ratios**: At 0% or 100% host, the approach is algebraic so it works, but "potentially informative" markers (types 20/21) contribute more noise at extremes
- **Low depth**: Stochastic effects analyzed extensively; with 3+ markers at >40,000x depth, deviation typically <0.5%

### Multi-donor Support

Not supported -- single donor assumption throughout.

### Output

Per-sample summary: mean, median, SD, and count of host percentages, separately for informative/potentially-informative markers, and for bias-corrected/uncorrected. 24-element vector per follow-up sample.

### Code Quality/Reusability

Research code quality -- functional but not packaged. Functions are well-documented and self-contained. The `biasFuns.R` file is ~200 lines and could be easily ported. No tests.

### Key Insight

**The per-marker amplification bias correction is essential and directly applicable to allomix.** Their data shows that heterozygous markers consistently deviate from 0.5 VAF in a marker-specific way (not random, not run-specific). Bias correction reduces median per-marker absolute error from 0.27% to 0.19% and dramatically improves type-II (potentially informative) markers. This is the single most important practical correction for VAF-based chimerism.

---

## 2. mvynck/Chimerism-FABCASE

**Repository**: https://github.com/mvynck/Chimerism-FABCASE
**Language**: R/Shiny | **Stars**: 0 | **License**: MIT

### Core Algorithm

FABCASE = "Forensic/clinical Analysis of Biallelic marker panel Case Assessment and Sample Evaluation". It is a **panel sufficiency calculator** -- given a set of marker allele frequencies (MAFs), it computes the probability of having at least N informative markers for a given donor-recipient pair.

The math uses the Poisson-Binomial distribution. For each marker, the probability of being informative is computed from the MAF using HWE genotype frequencies:

For **unrelated** donor-recipient (type-I markers only):
```
P(informative) = 2 * (p^2*(1-p)^2 + p*(1-p)^3 + p^3*(1-p))
```

For **siblings** (reduced by IBD sharing):
```
P(informative) = 0.5 * (3*p^2*(1-p)^2 + 2*p*(1-p)^3 + 2*p^3*(1-p))
```

Then the probability of having >= N informative markers across the whole panel follows a Poisson-Binomial distribution (each marker has different success probability). Confidence intervals are obtained by bootstrapping the MAF estimates.

### Input Format

CSV with marker MAF values (can be estimated from observed allele counts via `MAFvec` function which resamples using binomial).

### Multi-donor Support

The `simDoubleTransplant.R` companion script handles dual transplants -- simulates three genomes and counts markers informative for each donor vs. host.

### Key Insight

**This provides the theoretical framework for evaluating whether our 76-SNP panel has enough informative markers.** With 76 markers at MAF ~0.5, even sibling pairs should have >99.5% probability of getting >=3 informative markers. This validates our panel design. The functions for computing informativity probabilities from MAFs can be directly reused.

---

## 3. mvynck/Chimerism-nMarkers

**Repository**: https://github.com/mvynck/Chimerism-nMarkers
**Language**: R | **Stars**: 0 | **License**: No explicit license

### Core Algorithm

Monte Carlo simulation study to answer: "How many biallelic markers do you need for chimerism analysis?" Simulates genotype draws under HWE for various MAF distributions, donor-recipient relatedness levels (sibling, unrelated, parent-child), and marker types (informative, homozygous informative, all informative).

Key findings from the code:
- At MAF 0.5: ~10 markers gives >95% chance of >=3 informative (unrelated); ~20 for siblings
- The Devyser 24-plex (real MAFs) achieves >99% for unrelated, ~97% for siblings
- With 76 markers, informativity is essentially guaranteed for all relatedness levels
- Hardy-Weinberg testing on 218 samples showed no significant departures (q-value adjusted)

### Key Insight

**Our 76-marker panel is massively oversized for informativity -- the real question is whether each individual marker performs well at the sequencing/VAF level.** This tool confirms we should focus engineering effort on per-marker quality rather than panel size. Also validates that we should use type-I and type-II markers together for maximum sensitivity.

---

## 4. salipante/chimerism_smmip

**Repository**: https://github.com/salipante/chimerism_smmip
**Language**: Perl | **Stars**: 0 | **License**: Academic/non-commercial (UW license -- NOT compatible with MIT)

### Core Algorithm

Uses single-molecule molecular inversion probes (smMIPs) targeting copy number deletion polymorphisms. The approach is fundamentally different from VAF-based methods:

1. **Normalize** read counts per MIP probe to each sample's geometric mean
2. **Model** the expected count for each probe as: `model[i] = donor_norm[i] * (1-spike/100) + recipient_norm[i] * (spike/100)`
3. **Binary search** on `spike` parameter: User manually guesses recipient %, runs the model, checks if output ratio is >1 or <1, then iterates (bracketing + weighted average)
4. Uses 14 control probes to calibrate the normalization factor between model and observed

The 2-donor variant (`copy_num2_2donor_calculate.pl`) extends the model to sum two donor contributions.

### Input Format

Custom MIP count files (tab-delimited, 51+ columns per row), one per sample.

### Error Model

- Geometric mean normalization across probes (assumes log-normal count distribution)
- Control probes for inter-sample calibration
- Minimum read count thresholds for both genotyping (5) and samples (5)
- No probabilistic error model -- purely deterministic ratio approach

### Edge Cases

- Manual binary search (!!) -- user must iterate manually to find the mixture fraction
- Hardcoded control probe indices (positions 85-95, 170-172)
- Hardcoded exclusion of two "bad control" probes by name

### Multi-donor Support

Yes -- `copy_num2_2donor_calculate.pl` handles host + 2 donors. The model becomes:
`model[i] = (donor1_norm[i] + donor2_norm[i]) * (1-spike/100) + recipient_norm[i] * (spike/100)`

### Code Quality

Low. Hardcoded array indices, manual iteration, no error handling. ~200 lines of functional Perl.

### Key Insight

**The geometric mean normalization for probe-level bias correction is conceptually sound** but the implementation is extremely rough. The control-probe calibration idea is interesting -- we could use non-informative markers as internal controls. However, the overall approach (deletion polymorphisms, manual iteration) is too crude for our needs. The 2-donor model structure is worth noting for reference.

---

## 5. Ahhgust/Demixtify

**Repository**: https://github.com/Ahhgust/Demixtify
**Language**: C (with htslib) | **Stars**: 1 | **License**: AGPL-3.0 (NOT compatible with MIT for linked code)

### Core Algorithm

This is the **closest mathematical match to our problem**. Demixtify performs MLE-based mixture fraction estimation from biallelic SNPs.

**Likelihood function** (from demix.h and demix.c, referencing Crysup & Woerner 2022):

For each SNP locus with `N_ref` ref-supporting reads and `N_alt` alt-supporting reads, and a proposed mixture fraction `mf`, there are 9 possible two-person genotype combinations (AA|AA, AA|AB, AA|BB, AB|AA, ..., BB|BB).

For each genotype combination, the expected "weight" of allele A is:
```
w = (g1[0]=='A')*mf/2 + (g1[1]=='A')*mf/2 + (g2[0]=='A')*(1-mf)/2 + (g2[1]=='A')*(1-mf)/2
```

**Simplified likelihood** (Formula 5 of Crysup & Woerner):
```
log_like = N_ref * log(w*(1-e) + (1-w)*e/3) + N_alt * log((1-w)*(1-e) + w*e/3)
```

**Full per-base likelihood** (Formula 3):
```
For each base b at the site:
  like_b = w * Pr(b|A,q) + (1-w) * Pr(b|B,q)
  where Pr(b|match,q) = 1 - 10^(-q/10) and Pr(b|mismatch,q) = 10^(-q/10) / 3
```

The per-locus likelihood is marginalized over all 9 genotype combinations using population allele frequencies (with optional FST correction):
```
L(mf | locus) = sum_g [ P(g|AF,FST) * L(data|g,mf,e) ]
```

The overall likelihood is the product across all loci. Optimization is via grid search over `mf in [0, 0.5]` (101 grid points by default).

**Detection**: Likelihood ratio test comparing mixture hypothesis (best mf) vs single-source (mf=0), evaluated against chi-square(1 df).

**Confidence interval**: Chi-square approximation -- mf values where `2*(logL_max - logL(mf)) < chi_sq_critical`.

### Input Format

BAM file + VCF of known polymorphic sites (with AF tag). Directly reads from BAM using htslib pileup. Also supports known contributor BCF.

### Error Model

- **Sequencing error**: Estimated empirically from "other" allele counts (neither ref nor alt), or from base quality scores. Recalibration option.
- **Population structure**: FST correction on genotype priors via `GENOPROBFST(c, a, f) = (c[0]=='A'&&c[1]=='A' ? (1-a)^2 + (1-a)*a*f : ...)`
- **Quality filtering**: Min mapping quality (20), min base quality (20), max base quality cap (30), read length filter, indel-adjacent filter, duplicate/secondary read exclusion

### Edge Cases

- **Balanced mixtures**: Acknowledged down-bias as MF approaches 0.5 (inherent limitation of 2-unknown model)
- **Ultra-low coverage**: Recommends disabling theta/FST correction at <0.2x
- **Single-source**: If LLR does not favor mixture, reports as single source

### Multi-donor Support

**Two contributors only.** README explicitly states this limitation. However, notes that 3-person mixtures will likely still be flagged.

### Output

- Grid of mf values vs log-likelihoods
- Point estimate (MLE of mixture fraction)
- Confidence interval (chi-square)
- Per-site allele counts and genotype likelihoods
- Full VCF with deconvolved genotypes (PL, GQ, AD fields)
- Chi-square test LR for mixture vs single-source

### Code Quality

Moderate. ~4000 lines of C. Well-commented, references the paper. Uses htslib directly. Thread support. However, hardcoded to 2 contributors, and the code mixes C and C++ idioms. Would need significant refactoring to generalize to 3+ contributors.

### Key Insight

**The Crysup & Woerner likelihood framework is exactly what allomix needs.** The key formula is the marginalization over unknown genotypes weighted by population allele frequencies. For our case (known genotypes from screening), this simplifies dramatically -- we don't need to marginalize, we just compute the expected allele weight directly from the known genotypes. The FST correction and empirical error estimation are also directly applicable. The confidence interval via chi-square on the likelihood surface is the right approach for clinical reporting.

---

## 6. oyvble/euroformix

**Repository**: https://github.com/oyvble/euroformix
**Language**: R (with C++ backend) | **Stars**: 17 | **License**: LGPL-3.0 (compatible with MIT for linking, but derivative works must be LGPL)

### Core Algorithm

EuroForMix is the gold standard for forensic DNA mixture interpretation using STR peak heights. It uses a **continuous probabilistic model** that jointly estimates:
- Mixture proportions (N contributors)
- Per-contributor degradation slope
- Mean peak height (mu) and coefficient of variation (sigma)
- Back-stutter and forward-stutter ratios
- Drop-in probability

The likelihood function models observed peak heights as gamma-distributed:
```
h_a ~ Gamma(shape = (mu_a/sigma)^2, scale = sigma^2/mu_a)
```
where `mu_a` is the expected peak height for allele `a`, computed as:
```
mu_a = mu * sum_c(n_{c,a} * mx_c * deg_c(bp))
```
with `n_{c,a}` being the allele count for contributor c at allele a, `mx_c` the mixture proportion, and `deg_c(bp)` a degradation function.

Optimization is via `nlm()` with a presearch strategy that restricts the genotype outcome space. Supports MCMC for Bayesian inference and numerical integration for Bayes Factors.

### Input Format

Forensic-style CSV tables with alleles and peak heights per marker per sample. Population allele frequency databases. Kit-specific metadata (fragment sizes, stutter ratios).

### Error Model

Extremely comprehensive:
- Degradation (exponential decline with fragment length)
- Back-stutter and forward-stutter
- Allele drop-in (exponential model)
- Allele drop-out (probability from peak height distribution)
- Population structure (FST/theta correction on allele frequencies)
- Related individuals (IBD coefficients)

### Multi-donor Support

Yes -- arbitrary number of contributors. Each adds one mixture proportion parameter. Uses nlm optimization with multiple random restarts and genotype restriction for efficiency.

### Code Quality

High. Well-maintained R package with C++ backend, comprehensive documentation, test data, vignettes. Active development since ~2014. The `calcMLE.R` function is ~350 lines including presearch, optimization, deconvolution, and model validation.

### Key Insight

**The mathematical framework is elegant but over-engineered for our problem.** EuroForMix handles STR artifacts (stutters, degradation, drop-in/out) that don't apply to our SNP data. However, **the multi-contributor mixture proportion estimation framework is exactly right**: each contributor has a mixture proportion, genotypes are either known or marginalized, and the likelihood is product over loci. The presearch/restriction strategy for optimization with many unknowns is clever. We should adopt the overall architecture but with a much simpler per-locus model.

---

## 7. KhiabanianLab/All-FIT

**Repository**: https://github.com/KhiabanianLab/All-FIT
**Language**: Python | **Stars**: 18 | **License**: MIT (fully compatible)

### Core Algorithm

All-FIT (Allele Frequency-based Imputation of Tumor Purity) estimates tumor purity from variant allele frequencies. The core approach:

1. For each candidate purity `p` (grid from 0.01 to 0.99, step 0.01):
   - For each variant, compute the Cancer Cell Fraction (CCF) under various mutation models (somatic/germline, with/without LOH, various copy numbers)
   - Weight models by AIC-based model selection
   - Compute `L(p) = sum_variants( sum_models( weight * (CCF - 1)^2 ) )`
   - The best purity minimizes the aggregate weighted squared deviation of CCF from 1.0

2. Iteratively remove:
   - Germline heterozygous mutations (high weight for non-LOH germline models)
   - Subclonal mutations (VAF significantly below expected for purity)

3. Confidence interval: Computed as purity values where a lower bound (using mutation-wise variance) falls below the minimum aggregate score

The CCF for a somatic mutation with `c` mutant copies is:
```
CCF = VAF / (c*p / (2*(1-p) + ploidy*p))
```

### Input Format

Tab-delimited file with columns: unique_ID, variant_allele_frequency (as percentage), sequencing depth, ploidy.

### Error Model

- **Binomial confidence intervals** on VAF using beta distribution
- **AIC-based model weighting** across mutation models
- **Iterative outlier removal** (germline, subclonal)
- No explicit sequencing error model

### Edge Cases

- **Discontinuous CI**: Detected and warned about
- **Low variant count**: Works but with wide CI
- **Ploidy variation**: Supported per-variant

### Multi-donor Support

No -- single tumor/normal pair assumed.

### Output

Point estimate of purity with confidence interval. Detailed likelihood plots. Per-variant CCF assignments.

### Code Quality

Moderate. Single 400-line Python file. Clear structure but no tests. Input parsing is basic. Good visualization output.

### Key Insight

**The CCF-based approach of "which purity makes all mutations look clonal" is conceptually analogous to "which mixture fraction makes all informative markers consistent".** The grid search + iterative outlier removal is practical and robust. We should consider a similar approach: for each candidate chimerism %, compute the residual at each marker, then minimize the aggregate weighted residual. The AIC-based model weighting for handling different marker types (hom-hom vs hom-het) is a clean pattern.

---

## 8. nygenome/Conpair

**Repository**: https://github.com/nygenome/Conpair
**Language**: Python | **Stars**: 59 | **License**: Non-commercial academic license (NOT compatible with MIT)

### Core Algorithm

Conpair estimates contamination in tumor-normal pairs using a **Bayesian likelihood model** over a set of ~8000 pre-selected biallelic SNPs (MAF ~0.4-0.5, low LD).

**Contamination model**: For each marker with known genotype `G` (AA, AB, or BB) in the "pure" sample, and a contamination fraction `x`, the expected allele frequency in the contaminated sample follows one of 9 patterns (AAAA, AABB, AABA, ABAB, ABAA, etc. -- the 9 possible true_genotype|contaminant_genotype combinations).

For example, if the true genotype is AA and contamination source is BB at fraction `x`:
```
P(ref_base | AABB, x, bq) = (1-x)*(1-e) + x*e
P(alt_base | AABB, x, bq) = x*(1-e) + (1-x)*e
```
where `e = 10^(-bq/10) / 3`.

The per-marker likelihood marginalizes over the 9 genotype pairs using HWE priors from population allele frequencies. The overall log-likelihood is summed across markers. Optimization uses grid search (0.01 steps) followed by Brent's method for refinement.

**Concordance**: Separately computed using homozygous markers only -- if >95% match, samples are from the same individual.

### Input Format

GATK pileup files (generated from BAM via included script). Pre-selected marker set in BED/TXT format.

### Error Model

- **Base quality scores**: Used per-read (phred-to-probability conversion)
- **Genotype priors**: HWE from population allele frequencies
- **Homozygosity threshold**: P(AA) > 0.999 for calling homozygous
- **Downsampling**: Random downsample to 450x to prevent float underflow (!)
- No amplification bias correction

### Edge Cases

- **Copy number changes**: Concordance analysis supports `-H` flag to use only normal-homozygous markers (robust to CNV)
- **Zero coverage**: Markers with no reads at both alleles are skipped
- **Optimization**: Brent's method with bounds checking

### Multi-donor Support

No -- single contamination source model.

### Output

Contamination percentage for tumor and normal. Concordance percentage.

### Code Quality

Moderate-good. Clean Python modules with clear separation of concerns. Good documentation. However, Python 2.7 requirement is dated. The `ContaminationModel.py` at ~80 lines is elegant and reusable.

### Key Insight

**The per-base quality-aware likelihood is the right way to handle sequencing error**, and the Brent's algorithm refinement after grid search is more efficient than pure grid search. The 9-genotype-pair marginalization pattern (identical to Demixtify) is the canonical approach. The downsampling-at-450x hack is a warning: at our >1000x depth, numerical underflow in product-of-probabilities is a real concern -- we must work in log-space.

---

## 9. brentp/somalier

**Repository**: https://github.com/brentp/somalier
**Language**: Nim | **Stars**: 307 | **License**: MIT (fully compatible)

### Core Algorithm

Somalier extracts genotype-like information at ~17,000 pre-selected polymorphic sites and computes relatedness between all sample pairs using bit-vector operations.

**Extract phase**: For each site, count ref and alt alleles. Classify as HOM_REF, HET, HOM_ALT, or UNKNOWN using allele balance thresholds (default: AB < 0.01 = hom_ref, AB > 0.3 and < 0.7 = het).

**Relate phase**: For each sample pair, compute:
- IBS0: sites where one is hom-ref and other is hom-alt
- IBS2: sites with identical genotypes
- shared-hets: both heterozygous
- shared-hom-alts: both homozygous alt

**Relatedness** = `2 * (shared_hets - 2*IBS0) / het_ab` (Pedersen et al. 2020)

Uses 3 bit-vectors per sample (64-bit integers) with popcount hardware instructions for O(1) per-pair computation.

### Input Format

BAM/CRAM/VCF/GVCF + sites VCF + reference FASTA. Extracts to compact binary `.somalier` files.

### Error Model

- **Allele balance threshold**: Configurable (default 0.3 for het calls, adjustable to 0.2 for RNA-seq)
- **Minimum depth**: 7 reads default
- **Quality filters**: PASS/RefCall filter only
- No bias correction -- relies on threshold-based genotype calling

### Multi-donor Support

Not applicable -- this is a relatedness/QC tool, not a mixture analysis tool.

### Output

Per-sample QC metrics (depth, het rate, ancestry). Per-pair relatedness metrics. Interactive HTML visualization. TSV output for programmatic use.

### Code Quality

Excellent. Well-engineered Nim code, extensive documentation, binary releases, bioconda package, Docker image, active maintenance. ~300 GitHub stars. Clean architecture with extract/relate/ancestry subcommands.

### Key Insight

**The site selection strategy and extract-then-analyze architecture are directly applicable to allomix.** Somalier's pre-selected sites have MAF ~0.5 for maximum discrimination power -- exactly what we want. The binary extract format (compact per-sample files) is an elegant pattern for building a genotype database. For allomix, we could use a similar architecture: extract allele counts at our 76 SNPs into compact files, then run chimerism analysis on these extracts. The bit-vector relatedness approach could serve as a fast pre-check for donor-host identity confirmation before chimerism analysis. The `find-sites` subcommand's site selection criteria (MAF, LD, quality) is also useful reference for panel evaluation.

---

## Synthesis and Recommendations

### Mathematical Approach for allomix

**Recommended: MLE with known-genotype simplification, drawing from Demixtify's framework**

The mathematical framework should be:

1. **Likelihood function**: For each informative marker `i` with known host genotype `G_h` and donor genotype `G_d`, at a proposed chimerism fraction `f` (fraction host):

```
w_i(f) = f * allele_dose(G_h) / 2 + (1-f) * allele_dose(G_d) / 2
```

where `allele_dose` is the count of the reference allele (0, 1, or 2).

The per-marker log-likelihood is:
```
LL_i(f) = n_ref_i * log(w_i * (1-e) + (1-w_i) * e/3) + n_alt_i * log((1-w_i) * (1-e) + w_i * e/3)
```

This is Demixtify's Formula 5 but with **known genotypes** (no marginalization needed), making it much simpler and faster.

2. **Multi-donor extension**: For host + 2 donors with proportions `f_h, f_d1, f_d2` (where `f_h + f_d1 + f_d2 = 1`):

```
w_i(f_h, f_d1) = f_h * dose(G_h)/2 + f_d1 * dose(G_d1)/2 + (1-f_h-f_d1) * dose(G_d2)/2
```

Optimize over 2D grid (f_h, f_d1) with constraint f_h + f_d1 <= 1.

3. **Optimization**: Grid search (0.1% steps from 0% to 100%) for point estimate, then Brent refinement. For 2-donor: 2D grid then Nelder-Mead.

4. **Confidence interval**: Profile likelihood -- chi-square(1df) for single donor, chi-square(2df) for two donors.

**Reasoning**: Known genotypes eliminate the most complex part of forensic mixture analysis (genotype marginalization). This makes the model simple, fast, and identifiable. The MLE approach gives proper statistical inference including CIs, unlike the mean/median approaches in Chimerism-Bias.

### Error/Bias Correction

Implement in this order of priority:

1. **Per-marker amplification bias** (from Chimerism-Bias): Estimate bias as `median(VAF_het - 0.5)` per marker across training samples. Apply analytic correction before or incorporate into likelihood as `w_corrected = w + bias_term`. This is the single biggest accuracy improvement.

2. **Sequencing error rate** (from Demixtify/Conpair): Either estimate from base qualities or empirically from "other" allele counts. Use in the likelihood function's `e` parameter.

3. **Depth-based weighting**: Markers with higher depth contribute more to the likelihood naturally (more reads = tighter binomial). This is handled automatically by the MLE framework -- no explicit weighting needed.

4. **Outlier detection** (from All-FIT): After initial MLE, compute per-marker residuals. Flag markers with residuals >3 SD for review (possible genotyping error, CNV, or sample issue).

### What We Can Learn from Each Tool

| Tool | Key Lesson for allomix |
|------|----------------------|
| **Chimerism-Bias** | Per-marker allelic bias correction is essential; estimate from het samples; ~0.3% improvement |
| **Chimerism-FABCASE** | Panel sufficiency framework; confirms 76 markers is more than enough |
| **Chimerism-nMarkers** | Theoretical informativity rates by relatedness; validates our design |
| **chimerism_smmip** | 2-donor model structure; geometric mean normalization; control probe calibration |
| **Demixtify** | MLE likelihood framework (Crysup & Woerner); grid search + chi-sq CI; FST correction; error estimation |
| **EuroForMix** | Multi-contributor architecture; presearch/restriction for optimization; proper Bayesian alternative |
| **All-FIT** | Grid search + iterative outlier removal; CCF residual minimization pattern |
| **Conpair** | Per-base-quality likelihood; Brent refinement; 9-genotype marginalization; log-space arithmetic |
| **somalier** | Extract-then-analyze architecture; site selection criteria; compact binary format; QC metrics |

### Ranked List: Most Useful References

1. **Demixtify** (CRITICAL) -- The MLE likelihood function is our starting point. Formulas 3 and 5 of Crysup & Woerner are the core math. Grid search + chi-square CI is our inference strategy. Cannot directly use code (AGPL), but the math is public domain.

2. **Chimerism-Bias** (HIGH) -- Per-marker bias correction is the single most impactful practical improvement. The marker type classification (6 types) and correction formulas should be implemented. Code is simple R, easy to port.

3. **Conpair** (HIGH) -- The contamination model implementation is a clean, compact reference for the 9-genotype likelihood computation. The Brent's algorithm refinement is more efficient than pure grid search. Code is non-commercial license, but the math is standard.

4. **somalier** (HIGH) -- Architecture model for extract/analyze workflow. Site selection criteria. MIT license. We should consider somalier's `.somalier` format or similar for our genotype database.

5. **All-FIT** (MEDIUM) -- The grid search + outlier removal pattern is practical. MIT license means we could reference the code directly. The CCF residual minimization is a useful conceptual model.

6. **EuroForMix** (MEDIUM) -- The multi-contributor optimization strategy (presearch + restriction) is relevant if we need to handle >2 donors or unknown genotypes. LGPL license is workable.

7. **Chimerism-FABCASE** (MEDIUM) -- Panel sufficiency calculator. Worth running once to validate our 76-marker panel. MIT license.

8. **Chimerism-nMarkers** (LOW) -- Confirms panel design. Reference for publication only.

9. **chimerism_smmip** (LOW) -- Different assay approach (deletion polymorphisms). The 2-donor model structure is a reference, but the implementation is too crude to learn from. Non-commercial license.

### License Compatibility Summary

| Tool | License | MIT Compatible? |
|------|---------|----------------|
| Chimerism-Bias | None stated | Cannot reuse code; math only |
| Chimerism-FABCASE | MIT | Yes |
| Chimerism-nMarkers | None stated | Cannot reuse code; math only |
| chimerism_smmip | UW Academic | No -- non-commercial only |
| Demixtify | AGPL-3.0 | No -- copyleft, linking triggers AGPL |
| EuroForMix | LGPL-3.0 | Partial -- can link but not incorporate |
| All-FIT | MIT | Yes |
| Conpair | NYGC Non-commercial | No -- academic use only |
| somalier | MIT | Yes |

**For allomix code**: Implement the math independently (it's published science), do not copy code from AGPL/non-commercial repos. Can reference/adapt code from MIT-licensed repos (All-FIT, FABCASE, somalier).
