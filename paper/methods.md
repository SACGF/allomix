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

The donor fraction is estimated by maximum likelihood using a beta-binomial model that accounts for extra-binomial variance (overdispersion) arising from per-marker amplification bias and depth variability.

For a proposed donor fraction *f*, the expected reference allele weight at each informative marker is:

$$w_i(f) = (1 - f) \cdot \frac{g_{h,i}}{2} + f \cdot \frac{g_{d,i}}{2}$$

where $g_{h,i}$ and $g_{d,i}$ are the reference allele doses (0, 1, or 2) for the host and donor at marker *i*, respectively.

To account for sequencing errors (base substitutions, polymerase errors), the observed allele probabilities are modeled using a 4-state (trinucleotide) error model:

$$p_{alt,i} = (1 - w_i)(1 - \varepsilon) + w_i\frac{\varepsilon}{3}, \qquad p_{ref,i} = w_i(1 - \varepsilon) + (1 - w_i)\frac{\varepsilon}{3}$$

where $\varepsilon$ is the per-base sequencing error rate (default 0.01). The factor of 3 distributes error probability among the three non-observed bases. Since VCF allele-depth fields count only reference and alternative alleles, the conditional probability used in the likelihood is $\tilde{p}_i = p_{alt,i} / (p_{ref,i} + p_{alt,i})$.

This formulation uses the mixture genotype likelihood of Crysup and Woerner,[@CrysupWoerner2023] applied in the inverse direction. Crysup and Woerner derived Formula 5 for genotyping unknown contributors at a known mixture fraction; here we apply the same likelihood to estimate the mixture fraction given known contributor genotypes, a simplification afforded by the clinical chimerism setting where host and donor are genotyped independently before transplant.

#### Beta-binomial likelihood

A standard binomial model assumes all variance in allele counts comes from random sampling. In practice, per-marker amplification bias and depth variability produce overdispersion: the observed variance exceeds the binomial prediction. allomix uses a beta-binomial likelihood, the standard conjugate model for overdispersed count data,[@HindeDemetrio1998] to account for this. The beta-binomial models the true success probability at each marker as drawn from a Beta distribution centred on the expected probability, rather than being fixed.

The model is parameterised by the donor fraction *f* and a shared concentration parameter $\rho > 0$. For marker *i* with alternative allele count $k_i$ out of $n_i = n_{ref,i} + n_{alt,i}$ total reads, the per-marker log-likelihood (up to a constant) is:

$$\ell_i(f, \rho) = \log\Gamma(k_i + \alpha_i) + \log\Gamma(n_i - k_i + \beta_i) - \log\Gamma(n_i + \rho) - \log\Gamma(\alpha_i) - \log\Gamma(\beta_i) + \log\Gamma(\rho)$$

where $\alpha_i = \tilde{p}_i \cdot \rho$ and $\beta_i = (1 - \tilde{p}_i) \cdot \rho$. As $\rho \to \infty$, the beta-binomial converges to the binomial (no overdispersion); smaller values of $\rho$ produce flatter likelihoods and wider confidence intervals.

The total log-likelihood across all *M* informative markers is:

$$\mathcal{L}(f, \rho) = \sum_{i=1}^{M} \ell_i(f, \rho)$$

#### Optimization

Both *f* and $\rho$ are estimated jointly from the data. Optimization proceeds in two stages. First, a grid search evaluates the likelihood at 1,001 evenly spaced values of *f* across the interval [0, 1], with $\rho$ profiled out (optimised on the log-scale via bounded Brent) at each grid point. Second, Nelder-Mead refinement over both (*f*, log $\rho$) from the grid maximum yields the MLE point estimates $\hat{f}$ and $\hat{\rho}$.

### Confidence Intervals

A 95% profile likelihood confidence interval for *f* is constructed by profiling out $\rho$ at each candidate value of *f*. The profile log-likelihood is:

$$\mathcal{L}_P(f) = \max_{\rho} \mathcal{L}(f, \rho)$$

The CI bounds are the values of *f* where the profile log-likelihood ratio reaches the chi-squared threshold:

$$2[\mathcal{L}_P(\hat{f}) - \mathcal{L}_P(f)] = \chi^2_{1, 0.95} \approx 3.84$$

The lower and upper bounds are found using Brent's root-finding method (scipy.optimize.brentq), following standard profile likelihood methodology.[@Wilks1938] By profiling out $\rho$, the confidence interval automatically adapts to the level of overdispersion in the data: panels with high marker-to-marker variability produce wider intervals, while well-behaved panels produce intervals close to the binomial baseline.

### Multi-Donor Extension

The single-donor model generalises naturally to multiple donors. For two donors with fractions $f_1$ and $f_2$ (subject to the constraint $f_1 + f_2 \leq 1$), the expected reference allele weight at marker $i$ becomes:

$$w_i(f_1, f_2) = (1 - f_1 - f_2) \cdot \frac{g_{h,i}}{2} + f_1 \cdot \frac{g_{d1,i}}{2} + f_2 \cdot \frac{g_{d2,i}}{2}$$

where $g_{d1,i}$ and $g_{d2,i}$ are the reference allele doses for donors 1 and 2. The per-marker and total log-likelihoods follow identically from the single-donor case, with $w_i$ now a function of two parameters. A marker is considered informative if the host genotype differs from that of any donor; per-donor informative counts are tracked separately.

Optimization proceeds by triangular grid search over the simplex $\{(f_1, f_2) : f_1, f_2 \geq 0,\; f_1 + f_2 \leq 1\}$ at 101 steps per dimension (~5,150 evaluations), followed by Nelder-Mead refinement from the grid maximum. Profile likelihood 95% confidence intervals are computed per donor using a $\chi^2$ threshold with 1 degree of freedom ($\chi^2_{1,0.95} \approx 3.84$), since each CI profiles one donor fraction while optimising the other.

### Per-Marker Bias Correction

Capture and amplicon-based sequencing panels exhibit systematic per-marker amplification biases that cause observed variant allele frequencies to deviate from their true values.[@Vynck2023bias] allomix supports optional per-marker bias correction following the approach of De Vynck et al.

Bias is estimated from a set of training samples (typically genotyping controls with known heterozygous genotypes). For each marker, the bias $b_i$ is computed as the median deviation of observed variant allele frequency from 0.5 across all heterozygous observations:

$$b_i = \text{median}(\text{VAF}_{het,i} - 0.5)$$

The bias is estimated at heterozygous sites, where the expected reference allele weight is 0.5. Applying it as a flat additive shift $w'_i = w_i - b_i$ is only valid near 0.5; at informative markers whose expected weight is near 0 or 1 (the common case at low chimerism) a fixed additive shift overcorrects and degrades the fit. The correction is therefore applied multiplicatively, in logit space, where it remains valid at any expected weight:

$$w'_i = \text{expit}\!\left(\text{logit}(w_i) - \text{logit}(0.5 + b_i)\right)$$

clamped to the interval [$10^{-6}$, $1 - 10^{-6}$] to prevent numerical instability. At a heterozygous site ($w_i = 0.5$) this reduces to $w'_i = 0.5 - b_i$, matching the estimate; at an extreme expected weight it is a small proportional shift rather than a large additive jump. The corrected weight is then used in the likelihood calculation. Bias correction is optional and uses a pre-computed per-marker table. A marker is only correctable if its bias was measured where it was heterozygous; the markers that are informative for a given host/donor pair are homozygous in both contributors, so their bias cannot be measured from that pair and must come from a table built across other samples. The table can be estimated from a cohort of reference samples called the same way as the admixture (so that per-marker bias, which is caller-specific, transfers), or from the admixture samples of a patient cohort at markers where the host and every donor are heterozygous (true VAF 0.5 regardless of mixing). When no table is provided, correction is not applied.

### Quality Control

allomix performs several quality control assessments for each chimerism estimate:

1. **Marker sufficiency**: A minimum of 3 informative markers is required (configurable).
2. **Depth assessment**: Mean and median sequencing depth across informative markers are reported, with a warning if mean depth falls below 100-fold.
3. **Confidence interval width**: A warning is issued if the 95% CI exceeds 20 percentage points.
4. **Goodness-of-fit**: A chi-squared test is performed on the per-marker Pearson residuals (observed minus expected variant allele frequency). A significant result (p < 0.01) may indicate genotype errors, copy number alterations, or other systematic model violations.
5. **Outlier detection**: Markers with standardized residuals exceeding 3 standard deviations from the mean are flagged.

Each sample receives an overall pass/fail QC assessment based on these criteria.

### Simulation Framework

For validation, allomix includes a simulation module that generates synthetic chimeric VCFs by blending two genotype VCFs at a specified donor fraction. The simulated fraction corresponds to the donor proportion in the analysed DNA, regardless of whether that DNA was extracted from unfractionated whole blood or from a sorted cell subset; the same statistical framework applies in either case, only the clinical interpretation changes. For each marker, the expected alternative allele frequency is calculated from the mixture model:

$$\text{VAF}_{expected} = \frac{(1-f) \cdot a_h + f \cdot a_d}{2}$$

where $a_h$ and $a_d$ are the alternative allele doses (0, 1, or 2) for host and donor, respectively. The simulation incorporates four sources of measurement noise calibrated from empirical data:

1. **Per-marker amplification bias**: Each marker receives a fixed bias drawn from a heavy-tailed Gaussian mixture (95% from $\mathcal{N}(0, 0.012)$, 5% from $\mathcal{N}(0, 0.08)$; overall SD ~0.018), calibrated from {{ panel_empirical.n_het_total | commas }} heterozygous observations across {{ panel_empirical.n_bias_markers | fmt('g') }} markers in {{ panel_empirical.n_vcfs | fmt('g') }} joint-called VCFs from a 76-SNP rhAmpSeq sample identification panel ($\sigma_{bias}$ = {{ panel_empirical.sd_bias }}; Supplementary Table S1).

2. **Non-uniform depth**: Per-marker depths are drawn from a log-normal distribution matching the empirically observed CV of {{ panel_empirical.mean_sample_depth_cv }} (Supplementary Table S1).

3. **Sequencing errors**: The same 4-state (trinucleotide) error model used in the likelihood function is applied in simulation. Each read has probability $\varepsilon$ = 0.01 of being mis-called, with errors distributed uniformly among the three non-observed bases. Since only reference and alternative alleles are counted in VCF allele-depth fields, the observed alternative probability is the conditional $p_{alt} / (p_{ref} + p_{alt})$, where $p_{ref}$ and $p_{alt}$ are as defined in the likelihood model. This ensures the simulator and estimator use a consistent generative model.

4. **Locus dropout**: Each marker has a {{ panel_empirical.mean_nocall_pct }}% probability of producing zero reads, based on the empirical no-call rate.

Alternative allele counts are drawn from a binomial distribution with the biased, error-adjusted, and conditionally normalised expected frequency and per-marker depth. The simulator can optionally draw reads from a beta-binomial with a specified overdispersion concentration $\rho$ instead (binomial is the $\rho \to \infty$ limit); this is used to characterise the dependence of the limit of detection on overdispersion (Supplementary Figures S7, S8) but not in the main validation, which uses binomial sampling.

To evaluate longitudinal monitoring, we simulated a six-timepoint post-HSCT engraftment trajectory (day +14 to day +365) with true donor fractions ranging from 15% (early engraftment) to 97% (full donor chimerism), including a clinically relevant 3-percentage-point dip at day +180. Each timepoint was generated at 500x depth with the same noise model parameters as the depth validation, and five independent replicates were run with different random seeds.

### Limit of Detection

Limit of detection (LoD) and limit of blank (LoB) are defined following the Clinical and Laboratory Standards Institute guideline EP17-A2,[@CLSIEP17A2; @PiersonPerry2012] which is the framework used by published evaluations of comparable NGS chimerism assays.[@Blouin2024comparison; @Qama2026devyser] For a given combination of donor-host relatedness, sequencing depth, and panel size:

- **Limit of blank (LoB)** is the 95th percentile of the estimated donor fraction across replicates of a pure-host sample (true donor fraction = 0). Because the estimator is bounded at zero, LoB captures the upper tail of nonzero estimates produced by sampling noise, per-marker bias, and locus dropout on a blank input.
- **Limit of detection (LoD)** is the lowest true donor fraction at which at least 95% of sequencing replicates yield an estimate exceeding LoB. Empirical detection rates at each tested fraction are fitted with a 2-parameter logistic in $\log_{10}(f)$, $P(\text{detected} \mid f) = [1 + \exp(-(a + b\log_{10} f))]^{-1}$, and LoD is the fraction at which the fitted curve equals 0.95. LoB and LoD are computed separately for each donor/host pair (over that pair's sequencing replicates); we report the median LoD across pairs as the curve and the 10th-90th percentile across pairs as a band.

LoD was characterised across a sweep of 2 relatedness levels (unrelated, full sibling), 5 mean sequencing depths (100x, 250x, 500x, 1,000x, 2,000x), 6 panel sizes (25, 50, 75, 100, 200, 400 markers), and 7 true donor fractions (0, 0.1%, 0.2%, 0.5%, 1%, 2%, 5%). The design separates the two sources of variation that the previous pooled estimate conflated. For each relatedness level we draw multiple donor/host pairs (10 unrelated, 40 sibling); each pair's genotypes (from population allele frequencies, MAF 0.2--0.5) and per-marker amplification bias are fixed and reused across all depths and panel sizes, and markers are nested so a smaller panel is a strict prefix of a larger one. For each pair we then run 30 sequencing replicates per cell that vary only the read-sampling noise. Holding the pair fixed isolates sequencing noise within a pair and makes each pair's LoD curve monotone in panel size; the across-pair spread (driven by identity-by-descent sharing, which varies markedly between sibling pairs and little between unrelated pairs) is reported as the band rather than leaking into the central estimate as it did when a single pooled LoD mixed pairs of differing informative-marker counts. The estimator's profile-likelihood lower confidence bound is constrained to non-negative values, so at low true fractions a substantial proportion of CIs will touch zero and a CI-lower-bound detection rule would behave conservatively at small panel sizes; the EP17 LoB rule used here is independent of this constraint. Sibling-donor simulations use Mendelian segregation from independent parental haplotypes and therefore do not model non-random identity-by-descent around the HLA locus, which is irrelevant for sample-identification panels that avoid HLA but would shift effective informative-marker counts on panels overlapping the major histocompatibility complex.

### Recipient Copy-Number Aberrations

To assess sensitivity to somatic copy-number changes in the recipient clone, the simulator applies a per-marker recipient aberration before read sampling. The recipient cell population is modelled as a mixture of normal diploid cells and an aberrant clone present at a clonal fraction; the expected allele fraction is then a copy-number-weighted average over the normal recipient, the recipient clone, and the donor, rather than the diploid mean. Three aberration types are modelled by mutating one randomly chosen germline homolog of the clone: copy-neutral LoH (the retained homolog is duplicated, two copies, affecting only heterozygous markers), deletion (one homolog lost, one copy), and gain (one homolog duplicated, three copies). Deletion and gain change the locus DNA contribution and so affect homozygous as well as heterozygous markers. A burden parameter sets the fraction of eligible markers carrying the aberration. The recipient genotype reference is left as the clean germline, so the aberration affects only the admixture sample, matching the two-phase upstream workflow in which recipient genotypes come from a reference sample.

Using this generator we computed the EP17-A2 LoD (as above) in two directions at a pure clone, {{ cnv_loh_headline.ref_markers }} markers and {{ cnv_loh_headline.ref_depth }}x depth, sweeping burden over 0, 10%, 25%, and 50%. In the relapse direction the recipient clone is the minor component being detected against a pure-donor blank; in the donor direction the aberration-bearing recipient is the major background and the donor is the minor component. Depth and panel size can be swept as in the depth-by-markers LoD analysis; the figure fixes them at the reference operating point because the donor-LoD inflation is a systematic-bias floor rather than a sampling-noise limit, so it does not improve with depth (the same saturation behaviour as overdispersion). To limit the influence of aberrant markers, the estimator applies an iterative median/MAD outlier-resistant refit by default: per-marker residuals beyond a robust cutoff are dropped and the fraction re-estimated, gated so that it engages only when the number of outliers exceeds the chance expectation (leaving clean samples unchanged) and floored so it does not over-trim small panels. When the excluded fraction is large, the result is flagged for review rather than reported as a confident estimate.

{# TODO: Real sequencing data validation methods #}
{# Add subsection "### Clinical Sample Validation" describing: #}
{# - Sample cohort (retrospective post-HSCT patients from /tau) #}
{# - STR chimerism comparison methodology #}
{# - Concordance analysis approach #}
{# - LOD characterisation with dilution series if available #}

### Software Availability

allomix is implemented in Python with dependencies on cyvcf2, NumPy, and SciPy. It is available under the MIT license at https://github.com/SACGF/allomix. Installation is via pip (`pip install allomix`). The command-line interface provides three subcommands: `monitor` for single-sample or multi-timepoint analysis, `timeline` for consolidated multi-timepoint reporting, and `estimate-bias` for panel bias calibration.
