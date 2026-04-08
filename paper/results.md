## Results

### In Silico Validation Across Sequencing Depths

To assess allomix accuracy under realistic conditions, we generated synthetic chimeric VCFs across donor fractions from 0% to 100% using 100 biallelic SNP markers (80 informative) with per-marker capture bias drawn from the empirically calibrated heavy-tailed mixture model (see Methods), non-uniform depth (CV = 0.43), and 1.6% locus dropout. We repeated the validation at five sequencing depths (50x, 100x, 200x, 500x, and 1,000x) spanning the range encountered in clinical panels, from low-depth whole-exome capture to high-depth targeted amplicon sequencing. At each depth, {{ depth_50.n_replicates | dp(0) }} independent replicates were generated with different random seeds for per-marker biases and sampling noise.

allomix accurately estimated donor fraction across all depths and mixture levels (Figure 1). Estimation precision improved with increasing sequencing depth. At 50x depth, the mean absolute error (MAE) was {{ depth_50.mean_abs_error_pct | dp(2) }} ± {{ depth_50.mean_abs_error_sd_pct | dp(2) }}% (mean ± SD across replicates) with RMSE {{ depth_50.rmse_pct | dp(2) }} ± {{ depth_50.rmse_sd_pct | dp(2) }}%, while at 1,000x depth these improved to {{ depth_1000.mean_abs_error_pct | dp(2) }} ± {{ depth_1000.mean_abs_error_sd_pct | dp(2) }}% MAE and {{ depth_1000.rmse_pct | dp(2) }} ± {{ depth_1000.rmse_sd_pct | dp(2) }}% RMSE (Table 1, Figure 2). MAE remained below 1% at all depths tested (Figure 3).

![**Figure 1.** In silico validation of allomix across sequencing depths. Each panel shows true donor fraction (x-axis) versus estimated donor fraction (y-axis) at the indicated depth. Points represent all replicates (N={{ depth_50.n_replicates | dp(0) }}). Synthetic chimeric samples were generated with 100 markers (80 informative), empirically calibrated per-marker bias, non-uniform depth (CV = 0.43), and 1.6% locus dropout. Dashed line indicates perfect agreement.]({{ facts_dir }}/fig1_depth_scatter.png)

![**Figure 2.** Distribution of absolute estimation error by sequencing depth (N={{ depth_50.n_replicates | dp(0) }} replicates per depth). Boxes show median and interquartile range for interior fractions (excluding 0% and 100% donor). Whiskers extend to 1.5× IQR.]({{ facts_dir }}/fig2_depth_boxplots.png)

| Depth | MAE (%) | RMSE (%) | Max Error (%) | CI Coverage (%) | Mean CI Width (%) |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 50x | {{ depth_50.mean_abs_error_pct | dp(2) }} ± {{ depth_50.mean_abs_error_sd_pct | dp(2) }} | {{ depth_50.rmse_pct | dp(2) }} ± {{ depth_50.rmse_sd_pct | dp(2) }} | {{ depth_50.max_abs_error_pct | dp(2) }} ± {{ depth_50.max_abs_error_sd_pct | dp(2) }} | {{ depth_50.ci_coverage_pct }} ± {{ depth_50.ci_coverage_sd_pct }} | {{ depth_50.mean_ci_width_pct | dp(2) }} ± {{ depth_50.mean_ci_width_sd_pct | dp(2) }} |
| 100x | {{ depth_100.mean_abs_error_pct | dp(2) }} ± {{ depth_100.mean_abs_error_sd_pct | dp(2) }} | {{ depth_100.rmse_pct | dp(2) }} ± {{ depth_100.rmse_sd_pct | dp(2) }} | {{ depth_100.max_abs_error_pct | dp(2) }} ± {{ depth_100.max_abs_error_sd_pct | dp(2) }} | {{ depth_100.ci_coverage_pct }} ± {{ depth_100.ci_coverage_sd_pct }} | {{ depth_100.mean_ci_width_pct | dp(2) }} ± {{ depth_100.mean_ci_width_sd_pct | dp(2) }} |
| 200x | {{ depth_200.mean_abs_error_pct | dp(2) }} ± {{ depth_200.mean_abs_error_sd_pct | dp(2) }} | {{ depth_200.rmse_pct | dp(2) }} ± {{ depth_200.rmse_sd_pct | dp(2) }} | {{ depth_200.max_abs_error_pct | dp(2) }} ± {{ depth_200.max_abs_error_sd_pct | dp(2) }} | {{ depth_200.ci_coverage_pct }} ± {{ depth_200.ci_coverage_sd_pct }} | {{ depth_200.mean_ci_width_pct | dp(2) }} ± {{ depth_200.mean_ci_width_sd_pct | dp(2) }} |
| 500x | {{ depth_500.mean_abs_error_pct | dp(2) }} ± {{ depth_500.mean_abs_error_sd_pct | dp(2) }} | {{ depth_500.rmse_pct | dp(2) }} ± {{ depth_500.rmse_sd_pct | dp(2) }} | {{ depth_500.max_abs_error_pct | dp(2) }} ± {{ depth_500.max_abs_error_sd_pct | dp(2) }} | {{ depth_500.ci_coverage_pct }} ± {{ depth_500.ci_coverage_sd_pct }} | {{ depth_500.mean_ci_width_pct | dp(2) }} ± {{ depth_500.mean_ci_width_sd_pct | dp(2) }} |
| 1,000x | {{ depth_1000.mean_abs_error_pct | dp(2) }} ± {{ depth_1000.mean_abs_error_sd_pct | dp(2) }} | {{ depth_1000.rmse_pct | dp(2) }} ± {{ depth_1000.rmse_sd_pct | dp(2) }} | {{ depth_1000.max_abs_error_pct | dp(2) }} ± {{ depth_1000.max_abs_error_sd_pct | dp(2) }} | {{ depth_1000.ci_coverage_pct }} ± {{ depth_1000.ci_coverage_sd_pct }} | {{ depth_1000.mean_ci_width_pct | dp(2) }} ± {{ depth_1000.mean_ci_width_sd_pct | dp(2) }} |

**Table 1.** allomix accuracy and confidence interval performance across sequencing depths (mean ± SD, N={{ depth_50.n_replicates | dp(0) }} replicates). MAE = mean absolute error; RMSE = root mean square error. Error metrics are computed on interior fractions (excluding 0% and 100%). CI coverage and width are computed across all fractions.

![**Figure 3.** allomix performance as a function of sequencing depth (mean ± SD, N={{ depth_50.n_replicates | dp(0) }} replicates). Left: accuracy metrics (MAE, RMSE, maximum error). Centre: 95% profile likelihood CI coverage versus nominal 95% level (dashed). Right: mean CI width.]({{ facts_dir }}/fig3_depth_summary.png)

### Effect of Per-Marker Bias Correction

To evaluate the impact of per-marker amplification bias correction, we generated synthetic chimeric samples at 2,000x depth with known per-marker biases ($\sigma_{bias}$ = 0.02) and ran allomix both with and without bias correction using the true bias values.

Without bias correction, the mean absolute error on interior fractions was {{ bias_no_bias.mean_abs_error_pct | dp(2) }}% with RMSE {{ bias_no_bias.rmse_pct | dp(2) }}%. With bias correction, the mean absolute error improved to {{ bias_with_bias.mean_abs_error_pct | dp(2) }}% with RMSE {{ bias_with_bias.rmse_pct | dp(2) }}%. Bias correction also improved boundary fraction estimates: the 0% donor sample was estimated at {{ bias_no_bias.est_0pct }}% without correction versus {{ bias_with_bias.est_0pct }}% with correction, and the 100% donor sample improved from {{ bias_no_bias.est_100pct }}% to {{ bias_with_bias.est_100pct }}%.

### Confidence Interval Calibration

The 95% profile likelihood confidence intervals showed coverage well below the nominal level, ranging from {{ depth_1000.ci_coverage_pct }}% at 1,000x depth to {{ depth_50.ci_coverage_pct }}% at 50x depth (Table 1, Figure 3). This undercoverage reflects the impact of unmodeled systematic noise sources in the simulation (non-uniform depth, heavy-tailed capture biases, locus dropout) that the binomial likelihood model does not account for. Coverage was higher at lower depths, where stochastic sampling noise is larger relative to systematic biases, widening the CIs. At higher depths, the narrower CIs become sensitive to these unmodeled effects, reducing coverage. Approaches to improving CI calibration include incorporating bias correction, overdispersion modelling via beta-binomial likelihoods, and empirical CI recalibration from training data.

### Comparison with Existing Tools

Table 2 summarises available NGS-based chimerism tools. allomix is the only open-source option that works with arbitrary marker panels from standard VCF files.

| Tool | Markers | LOD | Open Source | Panel Agnostic | Input |
|:---|:---:|:---:|:---:|:---:|:---:|
| AlloSeq HCT | 202 SNPs | 0.3% | No | No | Proprietary |
| Devyser Chimerism | 24 indels | 0.05% | No | No | Proprietary |
| NGStrack | 34 indels | 0.1% | No | No | Proprietary |
| ScisGo Chimerism MD | SNPs + indels | 0.2–0.5% | No | No | Proprietary |
| **allomix** | **Any biallelic** | **~0.6% MAE (in silico)** | **Yes (MIT)** | **Yes** | **VCF** |

**Table 2.** Comparison of NGS-based chimerism monitoring tools. LOD = limit of detection.

### Effect of Donor-Host Relatedness

In clinical HSCT, donors may be unrelated, or may be relatives of the host, with siblings being the most common related donor type. Increased relatedness reduces the number of informative markers (where donor and host genotypes differ), potentially degrading chimerism estimation. To evaluate this, we generated synthetic donor-host pairs at four relatedness levels (unrelated, first cousin, half-sibling, full sibling) using population allele frequencies (MAF 0.2–0.5), with {{ rel_unrelated.n_replicates }} replicate pairs per level and {{ rel_unrelated.n_markers }} markers at 500x depth (Figure 4, Table 3).

The mean number of informative markers decreased with increasing relatedness, from {{ rel_unrelated.mean_informative }} (unrelated) to {{ rel_sibling.mean_informative }} (full sibling). Despite this reduction, allomix maintained sub-2% mean absolute error across all relatedness levels: {{ rel_unrelated.mean_mae_pct }}% (unrelated), {{ rel_cousin.mean_mae_pct }}% (cousin), {{ rel_half_sibling.mean_mae_pct }}% (half-sibling), and {{ rel_sibling.mean_mae_pct }}% (sibling). Even in the worst case (sibling donors), the minimum number of informative markers observed was {{ rel_sibling.min_informative }}, well above the minimum of 3 required for estimation.

| Relatedness | Mean Informative | Range | MAE (%) | RMSE (%) |
|:---|:---:|:---:|:---:|:---:|
| Unrelated | {{ rel_unrelated.mean_informative }} | {{ rel_unrelated.min_informative }}–{{ rel_unrelated.max_informative }} | {{ rel_unrelated.mean_mae_pct }} | {{ rel_unrelated.mean_rmse_pct }} |
| 1st cousin | {{ rel_cousin.mean_informative }} | {{ rel_cousin.min_informative }}–{{ rel_cousin.max_informative }} | {{ rel_cousin.mean_mae_pct }} | {{ rel_cousin.mean_rmse_pct }} |
| Half-sibling | {{ rel_half_sibling.mean_informative }} | {{ rel_half_sibling.min_informative }}–{{ rel_half_sibling.max_informative }} | {{ rel_half_sibling.mean_mae_pct }} | {{ rel_half_sibling.mean_rmse_pct }} |
| Full sibling | {{ rel_sibling.mean_informative }} | {{ rel_sibling.min_informative }}–{{ rel_sibling.max_informative }} | {{ rel_sibling.mean_mae_pct }} | {{ rel_sibling.mean_rmse_pct }} |

**Table 3.** Effect of donor-host relatedness on marker informativity and chimerism accuracy. Each relatedness level was tested with {{ rel_unrelated.n_replicates }} replicate donor-host pairs, {{ rel_unrelated.n_markers }} markers, 500x depth.

![**Figure 4.** Effect of donor-host relatedness on allomix performance. Left: number of informative markers by relatedness level (dots = individual replicates, bars = means). Centre: mean absolute error. Right: truth versus estimated donor fraction across all replicates. Simulated with {{ rel_unrelated.n_markers }} markers, 500x mean depth (CV = 0.43), 1% sequencing error rate, empirically calibrated per-marker bias, and 1.6% locus dropout.]({{ facts_dir }}/fig4_relatedness.png)

### Multi-Donor Estimation with Sibling Donors

To evaluate the multi-donor extension, we generated a three-sibling scenario: host and two donors sharing both parents, with genotypes produced by Mendelian segregation from parental haplotypes across {{ multidonor.n_markers | dp(0) }} biallelic markers at {{ multidonor.depth | commas }}x depth. Of {{ multidonor.n_markers | dp(0) }} total markers, {{ multidonor.n_informative_any | dp(0) }} were informative for at least one donor ({{ multidonor.n_informative_d1 | dp(0) }} for donor 1, {{ multidonor.n_informative_d2 | dp(0) }} for donor 2), reflecting the expected reduction due to full-sibling relatedness.

allomix was run on {{ multidonor.n_samples | dp(0) }} chimeric samples spanning the simplex of (donor 1, donor 2) fractions, including pure host, pure single-donor, balanced mixes, and asymmetric mixes (Table 4, Figure 5). Per-donor mean absolute error was {{ multidonor.mae_d1_pct }}% for donor 1 (RMSE {{ multidonor.rmse_d1_pct }}%) and {{ multidonor.mae_d2_pct }}% for donor 2 (RMSE {{ multidonor.rmse_d2_pct }}%), both well below 2%. The total donor fraction (f1 + f2) was estimated with {{ multidonor.mae_total_pct }}% MAE. All {{ multidonor.n_asymmetric | dp(0) }} asymmetric mixes were correctly ranked (the donor with the higher true fraction was estimated as the larger contributor in every case). Profile likelihood CI coverage was {{ multidonor.ci_coverage_d1_pct }}% for donor 1 and {{ multidonor.ci_coverage_d2_pct }}% for donor 2, showing the same undercoverage pattern observed in single-donor estimation.

| Metric | Donor 1 | Donor 2 | Total |
|:---|:---:|:---:|:---:|
| MAE (%) | {{ multidonor.mae_d1_pct }} | {{ multidonor.mae_d2_pct }} | {{ multidonor.mae_total_pct }} |
| RMSE (%) | {{ multidonor.rmse_d1_pct }} | {{ multidonor.rmse_d2_pct }} | {{ multidonor.rmse_total_pct }} |
| Max error (%) | {{ multidonor.max_error_d1_pct }} | {{ multidonor.max_error_d2_pct }} | {{ multidonor.max_error_total_pct }} |
| CI coverage (%) | {{ multidonor.ci_coverage_d1_pct }} | {{ multidonor.ci_coverage_d2_pct }} | — |

**Table 4.** Multi-donor chimerism accuracy with sibling donors. {{ multidonor.n_markers | dp(0) }} markers ({{ multidonor.n_informative_any | dp(0) }} informative), {{ multidonor.depth | commas }}x depth. Error metrics computed on interior fractions (excluding 0% and 100%).

![**Figure 5.** Multi-donor chimerism estimation. (A) Per-donor accuracy: true versus estimated donor fraction for donor 1 (circles, blue) and donor 2 (triangles, orange), with 95% profile likelihood CIs. (B) Two-dimensional log-likelihood surface for a representative mixture (60% host, 30% donor 1, 10% donor 2). Coloured contours show delta log-likelihood from maximum; dashed red line marks the 95% joint CI (chi-squared, df=2). Grey region is infeasible (f1 + f2 > 1). Star = true value, circle = MLE.]({{ facts_dir }}/fig_multidonor.png)

### Longitudinal Timeline Monitoring

To demonstrate allomix's support for serial chimerism monitoring, the primary clinical use case, we simulated a post-HSCT engraftment trajectory across {{ timeline.n_timepoints | dp(0) }} timepoints from day +14 to day +365 (Figure 6). The scenario includes early engraftment (15% donor at day +14), progressive engraftment through day +100 (95% donor), a clinically relevant dip at day +180 (92% donor, representing possible mixed chimerism that would trigger clinical attention), and recovery to 97% donor at day +365. Synthetic chimeric VCFs were generated at {{ timeline.depth | commas }}x depth with empirically calibrated noise and {{ timeline.n_replicates | dp(0) }} independent replicates. allomix tracked the trajectory accurately, with mean absolute error {{ timeline.mae_pct | dp(2) }} ± {{ timeline.mae_sd_pct | dp(2) }}% across all timepoints (maximum {{ timeline.max_error_pct | dp(2) }} ± {{ timeline.max_error_sd_pct | dp(2) }}%). The day +180 dip was correctly detected in all replicates, with a mean absolute error of {{ timeline.dip_abs_error_pct | dp(2) }}% at that timepoint.

![**Figure 6.** Simulated post-HSCT engraftment monitoring. True donor fraction trajectory (grey squares, dashed) versus allomix estimates (blue circles, solid) across six timepoints. Thin blue lines show individual replicates (N={{ timeline.n_replicates | dp(0) }}); thick blue line is the mean estimate. Shaded band shows the mean 95% profile likelihood CI. The trajectory includes a clinically relevant dip at day +180. Simulated with 100 markers, {{ timeline.depth | commas }}x depth, empirically calibrated per-marker bias, and 1.6% locus dropout.]({{ facts_dir }}/fig_timeline.png)

{# TODO: Add subsection "### Validation with Clinical Sequencing Data" #}
{# - Concordance with STR-based chimerism results (scatter, Bland-Altman) #}
{# - LOD characterisation from dilution series or low-fraction samples #}
{# - Multi-donor validation with real cord blood / dual-transplant cases if available #}
{# - Lineage-specific chimerism (CD3, CD33, etc.) if sorted fractions available #}
