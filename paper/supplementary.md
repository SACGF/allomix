# Supplementary Data

## Supplementary Table S1. Empirical Panel Characterisation

Per-marker amplification bias, depth distribution, locus dropout, and allele dropout were measured from {{ panel_empirical.n_vcfs | fmt('g') }} joint-called VCFs ({{ panel_empirical.n_samples | commas }} samples) generated from the 76-SNP IDT rhAmpSeq Sample ID panel as part of routine clinical sequencing. All {{ panel_empirical.n_bias_markers | fmt('g') }} biallelic markers with heterozygous observations were included. Simulation parameters used throughout this study were calibrated from these measurements.

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

![Figure S1](output/facts/figS1_bias_distributions.png)

**Figure S1.** Per-marker amplification bias distribution. (A) Histogram of empirical per-marker bias (median het VAF deviation from 0.5) measured across {{ supp_synthetic.n_empirical_markers | dp(0) }} markers, with kernel density estimates from a simple Gaussian model and the heavy-tailed mixture model used in allomix simulations. The mixture model (95% N(0, 0.012), 5% N(0, 0.08)) captures the heavy tails observed in the empirical data, where the 95th percentile of |bias| reaches {{ supp_synthetic.empirical_p95_abs_bias }}. (B) Cumulative distribution of |bias| showing the mixture model tracks the empirical tail while the simple Gaussian underestimates extreme values.

### S2. Depth Distribution

![Figure S2](output/facts/figS2_depth_distributions.png)

**Figure S2.** Per-marker sequencing depth. (A) Empirical mean depth per marker vs log-normal model draws at the same mean and CV. (B) Within-marker depth CV across samples for each marker, showing the range of per-marker depth variability.

### S3. Heterozygous VAF Comparison

![Figure S3](output/facts/figS3_het_vaf.png)

**Figure S3.** Violin plots of median heterozygous VAF per marker: empirical measurements ({{ supp_synthetic.n_empirical_markers | dp(0) }} markers from 76-SNP rhAmpSeq panel) vs simulated values drawn from the heavy-tailed mixture bias model. Both distributions are centred on 0.5 with comparable spread, confirming that the simulation reproduces the per-marker VAF displacement observed in real sequencing data.

### S4. Noise Component Ablation

![Figure S4](output/facts/figS4_ablation.png)

**Figure S4.** Effect of individual noise components on estimation accuracy (500x depth, 10 replicates per condition, 7 conditions). (A) Overall RMSE by noise condition. Under ideal conditions RMSE is {{ supp_synthetic.ablation_rmse_ideal_pct }}%; amplification bias alone raises it to {{ supp_synthetic.ablation_rmse_bias_only_pct }}%, and bias correction leaves it essentially unchanged at this depth ({{ supp_synthetic.ablation_rmse_bias_corrected_pct }}%), since the injected biases average near zero and their spread is largely absorbed by the overdispersion term. The full realistic model with binomial read sampling (all noise sources, bias corrected) produces {{ supp_synthetic.ablation_rmse_full_pct }}% RMSE; this is the no-overdispersion baseline. Adding per-marker overdispersion to that full model (beta-binomial read sampling at the fitted concentration rho = 100, applied at intermediate-VAF markers where amplification jitter is physical) raises RMSE to {{ supp_synthetic.ablation_rmse_overdispersion_pct }}%, a larger effect than any single bias, depth, or sequencing-error component. This is consistent with overdispersion, rather than depth, being the dominant control on accuracy and on the limit of detection at clinical coverage (Figures S7, S8). (B) Mean absolute error by true donor fraction for each condition. Dashed lines indicate conditions with bias correction applied.

### S5. Confidence Interval Calibration

![Figure S5](output/facts/figS5_ci_calibration.png)

**Figure S5.** CI calibration under the full noise model (100 replicates per fraction, run as 10 parallel batches of 10 via Snakemake). (A) Observed 95% CI coverage rate by true donor fraction; overall coverage is {{ supp_synthetic.cal_coverage_pct }}%. (B) Mean CI width by true donor fraction, with standard deviation bars.

### S6. Per-Marker Residuals

![Figure S6](output/facts/figS6_residuals.png)

**Figure S6.** Per-marker residuals (observed minus expected VAF) from a simulated 30% donor mixture at 500x. (A) Residual histogram with normal fit. (B) Residuals plotted against expected VAF, showing no systematic trend across the VAF range.

### S7. Limit of Detection vs Depth (Saturation)

![Figure S7](output/facts/fig_lod_saturation.png)

**Figure S7.** Limit of detection as a function of mean depth for the simulated unrelated panels in the LoD sweep. Points are the in silico LoD per panel size; lines are the LoD model $\mathrm{LoD} = (A/\sqrt{M})\sqrt{(n + \rho)/(n(\rho + 1))}$ fitted across panels ($M$ = informative markers, $n$ = depth, $\rho$ = beta-binomial overdispersion concentration). The simulator draws reads from a binomial (Methods), so the fit returns a near-infinite $\rho$ and the LoD falls close to $1/\sqrt{n}$ with no floor (dashed binomial reference). Under real, finite $\rho$ the per-marker variance instead approaches $p(1-p)/(\rho + 1)$, so the LoD saturates at a floor and depth beyond an effective cap of $\rho + 1$ reads yields diminishing returns. The in silico LoD reported elsewhere is therefore an analytical best case under near-binomial sampling.

### S8. Effect of Overdispersion on the Limit of Detection

![Figure S8](output/facts/fig_overdispersion_lod.png)

**Figure S8.** In silico LoD as a function of the beta-binomial overdispersion concentration $\rho$, at {{ overdispersion_lod_headline.depth }}x depth with {{ overdispersion_lod_headline.n_markers }} informative markers (unrelated donor). Reads were simulated beta-binomial across a grid of $\rho$ and the donor fraction estimated with the standard pipeline; the analytic and simulated (tool) LoD agree closely. The LoD rises from {{ overdispersion_lod_headline.lod_binomial_pct }}% under pure-binomial sampling ($\rho \to \infty$) to {{ overdispersion_lod_headline.lod_rho100_pct }}% at $\rho = 100$ (a {{ overdispersion_lod_headline.fold_rho100_vs_binomial }}-fold increase) and {{ overdispersion_lod_headline.lod_rho30_pct }}% at $\rho = 30$. At clinical coverage the overdispersion, not the depth, is the dominant control on the achievable LoD, which is why a simulated $\rho$ calibrated from real per-sample fits is needed before the in silico LoD is taken as a performance figure.

### S9. Fixed-Bias-Per-Marker Stability

![Figure S9](output/facts/fig_bias_stability.png)

**Figure S9.** Validation of the fixed-bias-per-marker assumption used by the simulator and the bias-correction model. Each point is one of the {{ supp_synthetic.n_empirical_markers | dp(0) }} panel markers: the x-axis is its absolute median amplification bias (the systematic, marker-specific component), and the y-axis is its within-marker standard deviation of heterozygous VAF across samples (the random, sample-to-sample component). The two are only weakly correlated (r = {{ supp_synthetic.bias_stability_r | dp(2) }}), so a marker's systematic bias does not predict its sample-to-sample scatter. Bias behaves as a stable per-marker offset rather than a quantity that grows with marker noise, which supports modelling it as a fixed offset (Methods) and absorbing the residual scatter separately through the overdispersion term.
