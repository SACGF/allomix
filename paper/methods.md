## Materials and Methods

### Overview

allomix is implemented in Python (version 3.10 or later) and operates on standard Variant Call Format (VCF) files. The workflow comprises four stages: (1) VCF parsing and marker classification, (2) maximum likelihood chimerism estimation, (3) optional per-marker bias correction, and (4) quality control assessment. The tool is designed to work with any set of biallelic markers (SNPs or indels) and imposes no requirements on the specific panel used.

### Input Requirements and VCF Parsing

allomix requires three sets of VCF files: host (recipient) pre-transplant genotyping, donor pre-transplant genotyping, and post-transplant admixture sample(s). VCF parsing is performed using the cyvcf2 library. For each biallelic site, the tool extracts the genotype (GT), allele depth (AD), total depth (DP), and genotype quality (GQ) FORMAT fields.

For low-fraction donor detection (below approximately 5%), admixture samples should be joint-called alongside donor and host genotyping samples in the same variant calling pipeline (e.g., GATK GenomicsDBImport followed by GenotypeGVCFs). This ensures that alternative alleles discovered in donor samples produce two-element AD fields (reference, alternative) in admixture samples even when those sites are called homozygous-reference, preventing loss of informative read count data.

### Marker Classification

Markers are classified as informative when the host and donor differ in genotype at that locus. Following the classification system of De Vynck et al.,[@Vynck2023bias] informative markers are assigned to one of six types based on the host and donor alternative allele dosage:

- **Type 0**: Host homozygous-reference (0/0), donor homozygous-alternative (1/1)
- **Type 1**: Host homozygous-alternative (1/1), donor homozygous-reference (0/0)
- **Type 10**: Host heterozygous (0/1), donor homozygous-reference (0/0)
- **Type 11**: Host heterozygous (0/1), donor homozygous-alternative (1/1)
- **Type 20**: Host homozygous-reference (0/0), donor heterozygous (0/1)
- **Type 21**: Host homozygous-alternative (1/1), donor heterozygous (0/1)

Types 0 and 1 are fully informative (maximum allelic contrast between host and donor), while types 10, 11, 20, and 21 are partially informative (one contributor is heterozygous, providing half the allelic contrast). Markers where host and donor share the same genotype are non-informative and excluded from analysis. Filtering is applied based on minimum genotype quality for host and donor samples (default GQ >= 20) and minimum read depth for the admixture sample (default DP >= 100).

### Maximum Likelihood Estimation

The donor fraction is estimated by maximum likelihood. For a proposed donor fraction *f*, the expected reference allele weight at each informative marker is:

$$w_i(f) = (1 - f) \cdot \frac{g_{h,i}}{2} + f \cdot \frac{g_{d,i}}{2}$$

where $g_{h,i}$ and $g_{d,i}$ are the reference allele doses (0, 1, or 2) for the host and donor at marker *i*, respectively.

To account for sequencing errors (base substitutions, polymerase errors), the observed allele probabilities are modeled as:

$$p_{ref,i} = w_i(1 - \varepsilon) + (1 - w_i)\frac{\varepsilon}{3}$$

$$p_{alt,i} = (1 - w_i)(1 - \varepsilon) + w_i\frac{\varepsilon}{3}$$

where $\varepsilon$ is the per-base sequencing error rate (default 0.01). The factor of 3 distributes error probability among the three non-reference (or non-alternative) bases.

The per-marker log-likelihood is:

$$\ell_i(f) = n_{ref,i} \cdot \log(p_{ref,i}) + n_{alt,i} \cdot \log(p_{alt,i})$$

where $n_{ref,i}$ and $n_{alt,i}$ are the observed reference and alternative allele read counts at marker *i*. The total log-likelihood across all *M* informative markers is:

$$\mathcal{L}(f) = \sum_{i=1}^{M} \ell_i(f)$$

This formulation is adapted from the mixture deconvolution model of Crysup and Woerner,[@CrysupWoerner2022] simplified for the case of known contributor genotypes.

Optimization proceeds in two stages. First, a grid search evaluates the likelihood at 1,001 evenly spaced points across the interval [0, 1], identifying the approximate maximum. Second, bounded Brent optimization (via scipy.optimize.minimize_scalar) refines the estimate within a ±1% window around the grid maximum, yielding the MLE point estimate $\hat{f}$.

### Confidence Intervals

A 95% profile likelihood confidence interval is constructed by identifying the bounds where the log-likelihood drops by a threshold derived from the chi-squared distribution:

$$2[\mathcal{L}(\hat{f}) - \mathcal{L}(f)] = \chi^2_{1, 0.95} \approx 3.84$$

The lower and upper bounds are found by scanning outward from the MLE in steps of 0.001, following standard profile likelihood methodology.[@Wilks1938]

### Multi-Donor Extension

The single-donor model generalises naturally to multiple donors. For two donors with fractions $f_1$ and $f_2$ (subject to the constraint $f_1 + f_2 \leq 1$), the expected reference allele weight at marker $i$ becomes:

$$w_i(f_1, f_2) = (1 - f_1 - f_2) \cdot \frac{g_{h,i}}{2} + f_1 \cdot \frac{g_{d1,i}}{2} + f_2 \cdot \frac{g_{d2,i}}{2}$$

where $g_{d1,i}$ and $g_{d2,i}$ are the reference allele doses for donors 1 and 2. The per-marker and total log-likelihoods follow identically from the single-donor case, with $w_i$ now a function of two parameters. A marker is considered informative if the host genotype differs from that of any donor; per-donor informative counts are tracked separately.

Optimization proceeds by triangular grid search over the simplex $\{(f_1, f_2) : f_1, f_2 \geq 0,\; f_1 + f_2 \leq 1\}$ at 101 steps per dimension (~5,150 evaluations), followed by Nelder-Mead refinement from the grid maximum. Profile likelihood 95% confidence intervals are computed per donor using a $\chi^2$ threshold with 1 degree of freedom ($\chi^2_{1,0.95} \approx 3.84$), since each CI profiles one donor fraction while optimising the other.

### Per-Marker Bias Correction

Capture and amplicon-based sequencing panels exhibit systematic per-marker amplification biases that cause observed variant allele frequencies to deviate from their true values.[@Vynck2023bias] allomix supports optional per-marker bias correction following the approach of De Vynck et al.

Bias is estimated from a set of training samples (typically genotyping controls with known heterozygous genotypes). For each marker, the bias $b_i$ is computed as the median deviation of observed variant allele frequency from 0.5 across all heterozygous observations:

$$b_i = \text{median}(\text{VAF}_{het,i} - 0.5)$$

During chimerism estimation, the expected reference allele weight at each marker is adjusted:

$$w'_i = w_i - b_i$$

clamped to the interval [$10^{-6}$, $1 - 10^{-6}$] to prevent numerical instability. The corrected weight is then used in the likelihood calculation. Bias correction is optional and requires a pre-computed bias table; when no bias table is provided, correction is not applied.

### Quality Control

allomix performs several quality control assessments for each chimerism estimate:

1. **Marker sufficiency**: A minimum of 3 informative markers is required (configurable).
2. **Depth assessment**: Mean and median sequencing depth across informative markers are reported, with a warning if mean depth falls below 100-fold.
3. **Confidence interval width**: A warning is issued if the 95% CI exceeds 20 percentage points.
4. **Goodness-of-fit**: A chi-squared test is performed on the per-marker Pearson residuals (observed minus expected variant allele frequency). A significant result (p < 0.01) may indicate genotype errors, copy number alterations, or other systematic model violations.
5. **Outlier detection**: Markers with standardized residuals exceeding 3 standard deviations from the mean are flagged.

Each sample receives an overall pass/fail QC assessment based on these criteria.

### Simulation Framework

For validation, allomix includes a simulation module that generates synthetic chimeric VCFs by blending two genotype VCFs at a specified donor fraction. For each marker, the expected alternative allele frequency is calculated from the mixture model:

$$\text{VAF}_{expected} = \frac{(1-f) \cdot a_h + f \cdot a_d}{2}$$

where $a_h$ and $a_d$ are the alternative allele doses (0, 1, or 2) for host and donor, respectively. The simulation incorporates four sources of measurement noise calibrated from empirical data:

1. **Per-marker amplification bias**: Each marker receives a fixed bias modelling systematic allele capture efficiency differences. We measured bias from {{ panel_empirical.n_het_total | commas }} heterozygous observations across {{ panel_empirical.n_bias_markers | fmt('g') }} markers in {{ panel_empirical.n_vcfs | fmt('g') }} joint-called VCFs ({{ panel_empirical.n_samples | commas }} samples) from a 76-SNP rhAmpSeq sample identification panel, obtaining an overall $\sigma_{bias}$ = {{ panel_empirical.sd_bias }}. The empirical bias distribution is heavy-tailed (median |bias| = {{ panel_empirical.median_abs_bias }}, 95th percentile = {{ panel_empirical.p95_abs_bias }}, maximum = {{ panel_empirical.max_abs_bias }}), so biases are drawn from a Gaussian mixture: 95% of markers from $\mathcal{N}(0, 0.012)$ and 5% from $\mathcal{N}(0, 0.08)$, yielding an overall SD of ~0.018 matching the empirical measurement.

2. **Non-uniform depth across markers**: In real panels, sequencing depth varies substantially across markers due to differences in primer/probe efficiency. Empirical characterisation showed a per-sample depth coefficient of variation of {{ panel_empirical.mean_sample_depth_cv }} (mean depth {{ panel_empirical.mean_depth | commas }}x, range {{ panel_empirical.min_depth | fmt('g') }}–{{ panel_empirical.max_depth | fmt('g') }}x). Per-marker depths are drawn from a log-normal distribution parameterised to match the target mean depth and empirical CV.

3. **Sequencing errors**: Each read is mis-called with probability $\varepsilon$ = 0.01, matching the error rate used in the likelihood model.

4. **Locus dropout**: Each marker has a probability of producing zero reads, set to {{ panel_empirical.mean_nocall_pct }}% based on the empirical no-call rate.

Alternative allele counts are drawn from a binomial distribution with the biased, error-adjusted expected frequency and per-marker depth. Empirical characterisation also showed a mean observed-to-expected heterozygosity ratio of {{ panel_empirical.mean_het_ratio }}, indicating negligible allele dropout at these depths.

To evaluate longitudinal monitoring, we simulated a six-timepoint post-HSCT engraftment trajectory (day +14 to day +365) with true donor fractions ranging from 15% (early engraftment) to 97% (full donor chimerism), including a clinically relevant 3-percentage-point dip at day +180. Each timepoint was generated at 500x depth with the same noise model parameters as the depth validation, and five independent replicates were run with different random seeds.

{# TODO: Real sequencing data validation methods #}
{# Add subsection "### Clinical Sample Validation" describing: #}
{# - Sample cohort (retrospective post-HSCT patients from /tau) #}
{# - STR chimerism comparison methodology #}
{# - Concordance analysis approach #}
{# - LOD characterisation with dilution series if available #}

### Software Availability

allomix is implemented in Python with dependencies on cyvcf2, NumPy, and SciPy. It is available under the MIT license at https://github.com/SACGF/allomix. Installation is via pip (`pip install allomix`). The command-line interface provides three subcommands: `monitor` for single-sample or multi-timepoint analysis, `timeline` for consolidated multi-timepoint reporting, and `estimate-bias` for panel bias calibration.
