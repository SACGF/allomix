# Supplementary Data

## Supplementary Methods

This section gives the statistical detail summarised in plain language in the main Methods. The notation is shared across subsections: for informative marker *i*, $g_{h,i}$ and $g_{d,i}$ are the reference-allele doses (0, 1, or 2) of host and donor, $n_i$ is the total read count, and $k_i$ is the alternative-allele count.

### S1. Mixture model and marker classification

For a proposed donor fraction *f*, the expected reference-allele weight at marker *i* is a blend of the two known genotypes:

$$w_i(f) = (1 - f)\,\frac{g_{h,i}}{2} + f\,\frac{g_{d,i}}{2}$$

Markers are informative when host and donor genotypes differ. Following De Vynck et al.,[@Vynck2023bias] each informative marker is assigned one of six types by host/donor alternative-allele dose:

- **Type 0**: host 0/0, donor 1/1 (fully informative)
- **Type 1**: host 1/1, donor 0/0 (fully informative)
- **Type 10**: host 0/1, donor 0/0 (partially informative)
- **Type 11**: host 0/1, donor 1/1 (partially informative)
- **Type 20**: host 0/0, donor 0/1 (partially informative)
- **Type 21**: host 1/1, donor 0/1 (partially informative)

Types 0 and 1 give the maximum allelic contrast (the minority allele has a single possible source); the heterozygous types give half the contrast. Markers where host and donor share a genotype are non-informative. Default filters: host and donor GQ $\geq$ 20, admixture DP $\geq$ 100, and at least three informative markers.

### S2. Sequencing-error model

To account for base substitutions and polymerase errors, the observed allele probabilities use a 4-state (trinucleotide) model:

$$p_{alt,i} = (1 - w_i)(1 - \varepsilon) + w_i\frac{\varepsilon}{3}, \qquad p_{ref,i} = w_i(1 - \varepsilon) + (1 - w_i)\frac{\varepsilon}{3}$$

where $\varepsilon$ is the per-base error rate (default 0.01) and the factor of 3 distributes error among the three non-observed bases. Because VCF allele-depth fields count only reference and alternative, the likelihood uses the conditional $\tilde{p}_i = p_{alt,i} / (p_{ref,i} + p_{alt,i})$. This is the mixture genotype likelihood of Crysup and Woerner[@CrysupWoerner2023] applied in the inverse direction (estimating the fraction from known genotypes rather than genotyping at a known fraction).

### S3. Beta-binomial likelihood and optimization

A binomial model assumes all variance comes from read sampling; in practice per-marker amplification bias and depth variability produce overdispersion. allomix uses a beta-binomial,[@HindeDemetrio1998] parameterised by *f* and a shared concentration $\rho > 0$. The per-marker log-likelihood (up to a constant) is:

$$\ell_i(f, \rho) = \log\Gamma(k_i + \alpha_i) + \log\Gamma(n_i - k_i + \beta_i) - \log\Gamma(n_i + \rho) - \log\Gamma(\alpha_i) - \log\Gamma(\beta_i) + \log\Gamma(\rho)$$

with $\alpha_i = \tilde{p}_i\,\rho$ and $\beta_i = (1 - \tilde{p}_i)\,\rho$. As $\rho \to \infty$ this converges to the binomial; smaller $\rho$ flattens the likelihood and widens intervals. The total log-likelihood is $\mathcal{L}(f,\rho) = \sum_{i=1}^{M}\ell_i(f,\rho)$. Both parameters are fit jointly: a grid search over 1,001 evenly spaced values of *f* in [0, 1] with $\rho$ profiled out (bounded Brent on log-scale) at each point, followed by Nelder-Mead refinement over $(f, \log\rho)$ from the grid maximum.

### S4. Profile-likelihood confidence intervals

The 95% interval for *f* inverts the profile log-likelihood $\mathcal{L}_P(f) = \max_\rho \mathcal{L}(f,\rho)$, with bounds where

$$2\left[\mathcal{L}_P(\hat{f}) - \mathcal{L}_P(f)\right] = \chi^2_{1,\,0.95} \approx 3.84$$

found by Brent root-finding.[@Wilks1938] The reference maximum is re-derived from the same profiled optimizer (the profile value at $\hat{f}$, not the joint optimum) so the root-finder brackets a sign change. Bounds are pinned at 0 and 1, so an estimate near a boundary does not produce an interval running past it.

### S5. Multi-donor extension

For two donors with fractions $f_1, f_2$ (host the remainder), the weight at marker *i* is

$$w_i(f_1, f_2) = (1 - f_1 - f_2)\,\frac{g_{h,i}}{2} + f_1\,\frac{g_{d1,i}}{2} + f_2\,\frac{g_{d2,i}}{2}$$

A marker is informative if the host differs from any donor; per-donor informative counts are tracked separately. Optimization is a triangular grid over the simplex $\{(f_1, f_2): f_1, f_2 \geq 0,\; f_1 + f_2 \leq 1\}$ at 101 steps per dimension (~5,150 evaluations), then Nelder-Mead refinement. Per-donor 95% profile intervals use $\chi^2_{1,0.95}$ (each interval profiles one donor while optimising the other).

### S6. Per-marker bias correction

Per-marker bias is estimated at heterozygous training observations as $b_i = \text{median}(\text{VAF}_{het,i} - 0.5)$. A flat additive shift $w_i - b_i$ is valid only near 0.5 and overcorrects at the extreme expected weights that dominate low-fraction samples, so the correction is multiplicative in logit space:

$$w'_i = \text{expit}\!\left(\text{logit}(w_i) - \text{logit}(0.5 + b_i)\right)$$

clamped to $[10^{-6}, 1 - 10^{-6}]$. At a heterozygous site ($w_i = 0.5$) this reduces to $0.5 - b_i$; at an extreme weight it is a small proportional shift. A marker is correctable only where its bias was measured (where it was heterozygous), so the table is built across other samples, either a reference cohort called the same way as the admixture, or admixture samples at markers where host and every donor are heterozygous (true VAF 0.5 regardless of mixing).

### S7. Residual-host presence test

The presence test uses only markers where the donor is homozygous and the host carries the donor-absent allele (Vynck types 0, 1, 10, 11). Let $y_i$ be the donor-absent allele count out of $n_i$ reads, $e_i$ the per-marker error background in that direction, and $h_i$ the host dose of the donor-absent allele. Under a host fraction $f_h$ the expected donor-absent allele probability is

$$q_i(f_h) = e_i + \frac{h_i}{2}\,f_h$$

Two statistics are reported. A pooled one-sided Poisson test uses $Y = \sum_i y_i$ against $\Lambda = \sum_i n_i e_i$, with $p = P(\text{Poisson}(\Lambda) \geq Y)$. A bounded-MLE likelihood-ratio test maximises a per-marker binomial likelihood in $q_i(f_h)$ over $f_h \geq 0$. Because the null sits on the boundary $f_h = 0$, the LRT p-value uses a chi-bar-square reference (a 50:50 mixture of a point mass at 0 and $\chi^2_1$) rather than naive Wilks, and the reported confidence interval is the profile-likelihood interval clipped at 0. The test thus returns a p-value, a host-fraction estimate $\hat{f}_h$, and a CI, and is calibrated against the per-marker error background ($e_i$); in the present work a symmetric global error rate sets $e_i = \varepsilon/3$ per marker, pending a per-site, per-direction empirical error table.

### S8. In-data contamination estimation

Contamination is measured at consensus-homozygous markers, where host and every donor are homozygous for the same allele, so the minor allele can only be sequencing error or foreign DNA. The headline estimate is the background-subtracted median per-site minor-allele fraction: the median is used rather than a pooled mean so a few gross miscall sites do not dominate, sites above 10% minor fraction are capped as miscalls, and the error floor is the 10th percentile of per-site minor fractions (the no-carrier/error sites), so contamination is reported as the heterogeneous excess over a uniform error floor (a uniform error elevation lifts the floor too and is correctly not called contamination). Contamination is distinguished from real low-level chimerism by marker geometry rather than magnitude: a dose-response in which the minor fraction rises with the number of co-pooled panel individuals carrying that allele indicates foreign reads, whereas a flat elevation indicates error. A separate sample-swap / third-genome test runs at the same consensus sites: a per-site binomial tail at the error rate flags sites where the minor allele is individually significant, combined into a swap p-value over discordant sites, catching a wrong-patient VCF that the informative-marker goodness-of-fit never sees. An optional `##allomixRunUnit` VCF header (flowcell:lane) supports a pure-metadata index-hopping flag, kept separate from the in-data estimate.

The optional per-marker contamination correction (off by default) acts on the magnitude estimate rather than the contamination report. On a co-pooled run a donor-homozygous informative marker carries extra reads on the host (donor-absent) allele from co-pooled genomes that happen to carry it, scaling with the number of those co-pooled carriers; the host signal is the same at every such marker while this contamination scales with carrier dose. The correction subtracts a dose term, `slope * n_carriers * depth`, from each donor-homozygous host-allele count before the fit, leaving the flat error floor to the per-site error model so it is not double-counted. Two quantities are measured per run, not assumed. The gate is the per-flowcell consensus-homozygous dose-response: the minor-allele fraction at consensus-homozygous sites is regressed on the co-pooled carrier count (weighted by depth, pooled across the run's serial timepoints with a per-timepoint intercept), and the correction is applied only when that slope is significantly positive; a clean run has a flat slope and the correction is a no-op. The magnitude is calibrated separately on the informative donor-homozygous markers themselves (the same weighted dose regression), because the consensus-homozygous slope predicts the informative-marker slope well enough to gate on but not to transfer one-for-one. The carrier counts come from the cohort's joint-called genotypes, the same input as the per-site error table. Only donor-homozygous markers (genotype-contrast types where the donor is homozygous) are corrected; the magnitude slope is clamped non-negative, so the correction can only lower the estimate, which is why it is gated rather than applied unconditionally.

### S9. Simulation, limit of detection, and copy-number model

The simulator draws each marker's expected alternative-allele frequency from the same mixture and error model used by the estimator, $\text{VAF}_{expected} = [(1-f)a_h + f a_d]/2$ (with $a_h, a_d$ the alternative doses), then applies per-marker bias, log-normal depth, sequencing error, and locus dropout (main Methods), and draws counts from a binomial (or, for the overdispersion characterisation, a beta-binomial at concentration $\rho$). The LoD sweep uses a nested design: for each relatedness level, multiple donor/host pairs ({{ lod_headline.n_pairs_unrelated }} unrelated, {{ lod_headline.n_pairs_sibling }} sibling) have fixed genotypes (population allele frequencies, MAF 0.2--0.5) and per-marker bias reused across all depths and panel sizes, with panels strictly nested (a smaller panel is a bit-identical prefix of a larger one); each pair then has {{ lod_headline.n_seq_reps }} sequencing replicates per cell that vary only read-sampling noise. Holding the pair fixed isolates sequencing noise within a pair and makes each pair's LoD curve monotone in panel size, while the across-pair identity-by-descent spread is reported as a band rather than leaking into the central estimate. Seeds are SHA-256-derived for process-stable reproducibility. For recipient copy-number aberrations, the recipient is modelled as a mixture of normal diploid cells and an aberrant clone at a clonal fraction, with the expected allele fraction a copy-number-weighted average over normal recipient, recipient clone, and donor; three aberration types are produced by mutating one germline homolog of the clone (copy-neutral LoH: retained homolog duplicated, two copies, heterozygous sites only; deletion: one copy; gain: three copies), with deletion and gain also changing the locus DNA contribution at homozygous sites. The aberration is applied only to the admixture sample, matching the two-phase workflow.

## Supplementary Table S1. Empirical Panel Characterisation

Per-marker amplification bias, depth distribution, locus dropout, and allele dropout were measured from {{ panel_empirical.n_vcfs | fmt('g') }} joint-called VCFs ({{ panel_empirical.n_samples | commas }} samples) generated from the {{ panel_specs.n_markers_panel }}-SNP IDT rhAmpSeq Sample ID panel as part of routine clinical sequencing. All {{ panel_empirical.n_bias_markers | fmt('g') }} biallelic markers with heterozygous observations were included. Simulation parameters used throughout this study were calibrated from these measurements.

| Parameter | Empirical Value | Simulation Default | Notes |
|:---|:---:|:---:|:---|
| **Amplification bias** | | | |
| Per-marker bias SD ($\sigma_{bias}$) | {{ panel_empirical.sd_bias }} | 0.02 | Close match |
| Mean \|bias\| | {{ panel_empirical.mean_abs_bias }} | | |
| Median \|bias\| | {{ panel_empirical.median_abs_bias }} | | |
| 95th percentile \|bias\| | {{ panel_empirical.p95_abs_bias }} | | Heavy tail |
| Max \|bias\| | {{ panel_empirical.max_abs_bias }} | | Single outlier marker |
| **Sequencing depth** | | | |
| Mean depth | {{ panel_empirical.mean_depth | commas }}x | 50–1,000x | Sims test lower depths |
| Median depth | {{ panel_empirical.median_depth | commas }}x | | |
| Min mean depth | {{ panel_empirical.min_depth | fmt('g') }}x | | Weakest marker |
| Max mean depth | {{ panel_empirical.max_depth | fmt('g') }}x | | |
| Per-sample depth CV | {{ panel_empirical.mean_sample_depth_cv }} | 0 (uniform) | Not yet modelled |
| **Locus dropout** | | | |
| Mean no-call rate | {{ panel_empirical.mean_nocall_pct }}% | 0% | Not yet modelled |
| Markers with >5% no-call | {{ panel_empirical.markers_gt5pct_nocall | fmt('g') }}/{{ panel_empirical.n_bias_markers | fmt('g') }} | | Single problematic marker |
| **Allele dropout** | | | |
| Mean het/HWE ratio | {{ panel_empirical.mean_het_ratio }} | | 1.0 = no ADO |
| Markers with ratio < 0.8 | {{ panel_empirical.markers_low_het | fmt('g') }}/{{ panel_empirical.n_bias_markers | fmt('g') }} | | Negligible at high depth |
| Estimated ADO rate | {{ panel_empirical.ado_estimate }} | 0 | Negligible |

## Supplementary Table S2. Per-Marker Detail

Per-marker statistics are available in the allomix repository at `paper/empirical_results/panel_per_marker.tsv`. Fields include: number of observations, call rate, genotype counts (hom-ref, het, hom-alt), observed-to-expected heterozygosity ratio (HWE), mean depth, depth CV, and median amplification bias. Marker identities are anonymised (sequential index only).

## Supplementary Table S3. Per-Sample Validation Results by Depth

Detailed per-sample validation results for each sequencing depth (50x, 100x, 200x, 500x, 1,000x) are available in the allomix repository at `output/depth_validation/`. For each depth, the true donor fraction, estimated fraction, error, and 95% confidence interval bounds are reported for all simulated mixture levels.

## Supplementary Figures: Simulation Model Validation

### S1. Amplification Bias Distribution

![Figure S1]({{ facts_dir }}/figS1_bias_distributions.png)

**Figure S1.** Per-marker amplification bias distribution. (A) Histogram of empirical per-marker bias (median het VAF deviation from 0.5) measured across {{ supp_synthetic.n_empirical_markers | dp(0) }} markers, with kernel density estimates from a simple Gaussian model and the heavy-tailed mixture model used in allomix simulations. The mixture model (95% N(0, 0.012), 5% N(0, 0.08)) captures the heavy tails observed in the empirical data, where the 95th percentile of |bias| reaches {{ supp_synthetic.empirical_p95_abs_bias }}. (B) Cumulative distribution of |bias| showing the mixture model tracks the empirical tail while the simple Gaussian underestimates extreme values.

### S2. Depth Distribution

![Figure S2]({{ facts_dir }}/figS2_depth_distributions.png)

**Figure S2.** Per-marker sequencing depth. (A) Empirical mean depth per marker vs log-normal model draws at the same mean and CV. (B) Within-marker depth CV across samples for each marker, showing the range of per-marker depth variability.

### S3. Heterozygous VAF Comparison

![Figure S3]({{ facts_dir }}/figS3_het_vaf.png)

**Figure S3.** Violin plots of median heterozygous VAF per marker: empirical measurements ({{ supp_synthetic.n_empirical_markers | dp(0) }} markers from {{ panel_specs.n_markers_panel }}-SNP rhAmpSeq panel) vs simulated values drawn from the heavy-tailed mixture bias model. Both distributions are centred on 0.5 with comparable spread, confirming that the simulation reproduces the per-marker VAF displacement observed in real sequencing data.

### S4. Noise Component Ablation

![Figure S4]({{ facts_dir }}/figS4_ablation.png)

**Figure S4.** Effect of individual noise components on estimation accuracy (500x depth, 10 replicates per condition, 7 conditions). (A) Overall RMSE by noise condition. Under ideal conditions RMSE is {{ supp_synthetic.ablation_rmse_ideal_pct }}%; amplification bias alone raises it to {{ supp_synthetic.ablation_rmse_bias_only_pct }}%, and bias correction leaves it essentially unchanged at this depth ({{ supp_synthetic.ablation_rmse_bias_corrected_pct }}%), since the injected biases average near zero and their spread is largely absorbed by the overdispersion term. The full realistic model with binomial read sampling (all noise sources, bias corrected) produces {{ supp_synthetic.ablation_rmse_full_pct }}% RMSE; this is the no-overdispersion baseline. Adding per-marker overdispersion to that full model (beta-binomial read sampling at the fitted concentration rho = 100, applied at intermediate-VAF markers where amplification jitter is physical) raises RMSE to {{ supp_synthetic.ablation_rmse_overdispersion_pct }}%, a larger effect than any single bias, depth, or sequencing-error component. This is consistent with overdispersion, rather than depth, being the dominant control on accuracy and on the limit of detection at clinical coverage (Figures S7, S8). (B) Mean absolute error by true donor fraction for each condition. Dashed lines indicate conditions with bias correction applied.

### S5. Confidence Interval Calibration

![Figure S5]({{ facts_dir }}/figS5_ci_calibration.png)

**Figure S5.** CI calibration under the full noise model (100 replicates per fraction, run as 10 parallel batches of 10 via Snakemake). (A) Observed 95% CI coverage rate by true donor fraction; overall coverage is {{ supp_synthetic.cal_coverage_pct }}%. (B) Mean CI width by true donor fraction, with standard deviation bars.

### S6. Per-Marker Residuals

![Figure S6]({{ facts_dir }}/figS6_residuals.png)

**Figure S6.** Per-marker residuals (observed minus expected VAF) from a simulated 30% donor mixture at 500x. (A) Residual histogram with normal fit. (B) Residuals plotted against expected VAF, showing no systematic trend across the VAF range.

### S7. Limit of Detection vs Depth (Saturation)

![Figure S7]({{ facts_dir }}/fig_lod_saturation.png)

**Figure S7.** Limit of detection as a function of mean depth for the simulated unrelated panels in the LoD sweep. Points are the in silico LoD per panel size; lines are the LoD model $\mathrm{LoD} = (A/\sqrt{M})\sqrt{(n + \rho)/(n(\rho + 1))}$ fitted across panels ($M$ = informative markers, $n$ = depth, $\rho$ = beta-binomial overdispersion concentration). The simulator draws reads from a binomial (Methods), so the fit returns a near-infinite $\rho$ and the LoD falls close to $1/\sqrt{n}$ with no floor (dashed binomial reference). Under real, finite $\rho$ the per-marker variance instead approaches $p(1-p)/(\rho + 1)$, so the LoD saturates at a floor and depth beyond an effective cap of $\rho + 1$ reads yields diminishing returns. The in silico LoD reported elsewhere is therefore an analytical best case under near-binomial sampling.

### S8. Effect of Overdispersion on the Limit of Detection

![Figure S8]({{ facts_dir }}/fig_overdispersion_lod.png)

**Figure S8.** In silico LoD as a function of the beta-binomial overdispersion concentration $\rho$, at {{ overdispersion_lod_headline.depth }}x depth with {{ overdispersion_lod_headline.n_markers }} informative markers (unrelated donor). Reads were simulated beta-binomial across a grid of $\rho$ and the donor fraction estimated with the standard pipeline; the analytic and simulated (tool) LoD agree closely. The LoD rises from {{ overdispersion_lod_headline.lod_binomial_pct }}% under pure-binomial sampling ($\rho \to \infty$) to {{ overdispersion_lod_headline.lod_rho100_pct }}% at $\rho = 100$ (a {{ overdispersion_lod_headline.fold_rho100_vs_binomial }}-fold increase) and {{ overdispersion_lod_headline.lod_rho30_pct }}% at $\rho = 30$. At clinical coverage the overdispersion, not the depth, is the dominant control on the achievable LoD, which is why a simulated $\rho$ calibrated from real per-sample fits is needed before the in silico LoD is taken as a performance figure.

### S9. Fixed-Bias-Per-Marker Stability

![Figure S9]({{ facts_dir }}/fig_bias_stability.png)

**Figure S9.** Validation of the fixed-bias-per-marker assumption used by the simulator and the bias-correction model. Each point is one of the {{ supp_synthetic.n_empirical_markers | dp(0) }} panel markers: the x-axis is its absolute median amplification bias (the systematic, marker-specific component), and the y-axis is its within-marker standard deviation of heterozygous VAF across samples (the random, sample-to-sample component). The two are only weakly correlated (r = {{ supp_synthetic.bias_stability_r | dp(2) }}), so a marker's systematic bias does not predict its sample-to-sample scatter. Bias behaves as a stable per-marker offset rather than a quantity that grows with marker noise, which supports modelling it as a fixed offset (Supplementary Methods S6) and absorbing the residual scatter separately through the overdispersion term.

### S10. Absolute Error by Depth (Boxplots)

![Figure S10]({{ facts_dir }}/fig2_depth_boxplots.png)

**Figure S10.** Distribution of absolute estimation error by sequencing depth (N={{ depth_50.n_replicates | dp(0) }} replicates per depth). Boxes show median and interquartile range for interior fractions (excluding 0% and 100% donor). Whiskers extend to 1.5× IQR. This is the per-fraction distribution behind the summary metrics in main-text Table 2.

### S11. Depth-Performance Summary

![Figure S11]({{ facts_dir }}/fig3_depth_summary.png)

**Figure S11.** allomix performance as a function of sequencing depth (mean ± SD, N={{ depth_50.n_replicates | dp(0) }} replicates). Left: accuracy metrics (MAE, RMSE, maximum error). Centre: 95% profile-likelihood CI coverage versus the nominal 95% level (dashed). Right: mean CI width.

![Figure S12]({{ facts_dir }}/figS12_srp434573_logy.png)

**Figure S12.** Confidence-interval view of the SRP434573 two-person dilution series (main-text Figure 4A shows the same data as a log-log scatter). Each timepoint is plotted on a log host-fraction axis grouped by mixture, with the maximum-likelihood estimate (filled circle, 100 minus donor%, with the per-marker contamination correction applied) and the residual-host presence-test estimate (open square) each shown with its 95% confidence interval, against the known fraction (grey diamond). The dashed horizontal line in each mixture is that mixture's independent in-data contamination level, measured at consensus-homozygous markers (a marker class the magnitude estimate never reads), so it is a floor estimated from different sites than the ones being corrected: an estimate below its mixture's line is not separable from contamination. The pure-donor (true-0%-host) endpoints, at the right of each group, sit at or near the 0 row after correction, at or below their contamination line. At the higher titration levels the intervals are tight and bracket the known value; at the 0.5% level they widen toward the contamination line, where the residual floor competes with the true host signal (Results, Figure 4 caption).

![Figure S13]({{ facts_dir }}/figS13_srp_contam.png)

**Figure S13.** Co-pooled contamination floor in SRP434573 as a dose-response, the figure behind the median values in Results. Consensus-homozygous sites are those where host and donor are both homozygous for the same allele, so the minor allele at that site cannot come from either contributor. For each such site, reads were pooled across the dilution samples (sites with pooled depth below 500 were dropped for a stable per-site fraction) and the per-site minor-allele fraction was computed. Sites are binned on the x-axis by the number of the five other co-pooled panel individuals that carry the minor allele, het or homozygous counting equally as one carrier: {{ srp_contam.n_nocarrier_sites | dp(0) }} sites have no carrier (the no-carrier bin) and {{ srp_contam.n_carrier_sites | dp(0) }} have at least one. Boxes show the per-site distribution on a log y-axis (median, interquartile range, 1.5x IQR whiskers, individual outliers; exact-zero sites are drawn at a 0.001% floor so they render on the log axis); the red line connects the per-bin medians and n below each box is the number of sites. The median rises monotonically from the no-carrier floor (sequencing error, {{ srp_contam.nocarrier_floor_pct | dp(3) }}%) through the carrier bins, the signature of real reads from co-pooled material (most plausibly index hopping) rather than flat sequencing error, which would not scale with co-pooled dose.
