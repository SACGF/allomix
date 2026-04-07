## Results

### In Silico Validation

To assess the accuracy of allomix, we generated synthetic chimeric VCFs spanning the full range of donor fractions (0% to 100% in 10% increments) using {{ validation_summary.n_markers | fmt('g') }} biallelic SNP markers ({{ validation_summary.n_informative | fmt('g') }} informative) at {{ validation_summary.depth | fmt(',.0f') }}-fold sequencing depth with per-marker capture bias ($\sigma_{bias}$ = {{ validation_summary.bias_sd }}).

allomix accurately estimated donor fraction across all tested levels (Figure 1, Table 1). The mean signed error was {{ validation_summary.mean_signed_error_pct | dp(4) }}%, indicating minimal systematic bias. The mean absolute error was {{ validation_summary.mean_abs_error_pct | dp(2) }}% and the root mean square error (RMSE) was {{ validation_summary.rmse_pct | dp(2) }}%. The maximum absolute error observed was {{ validation_summary.max_abs_error_pct | dp(2) }}%, occurring at the boundary fractions (0% and 100% donor) where per-marker capture bias has the largest relative effect.

{# TODO: Table 1. Per-sample validation results. For each simulated donor fraction, #}
{# the estimated fraction, error, 95% CI, and number of informative markers are shown. #}
{# Generate from output/validation/validation_results.tsv #}

| True Donor % | Estimated % | Error (%) | 95% CI | CI Width (%) | Informative Markers |
|:---:|:---:|:---:|:---:|:---:|:---:|
| {{ val_host_100_donor_0.true_pct }} | {{ val_host_100_donor_0.est_pct }} | {{ val_host_100_donor_0.error_pct | dp(2) }} | {{ val_host_100_donor_0.ci_lo_pct }}–{{ val_host_100_donor_0.ci_hi_pct }} | {{ val_host_100_donor_0.ci_width_pct }} | {{ val_host_100_donor_0.n_informative }} |
| {{ val_host_90_donor_10.true_pct }} | {{ val_host_90_donor_10.est_pct }} | {{ val_host_90_donor_10.error_pct | dp(2) }} | {{ val_host_90_donor_10.ci_lo_pct }}–{{ val_host_90_donor_10.ci_hi_pct }} | {{ val_host_90_donor_10.ci_width_pct }} | {{ val_host_90_donor_10.n_informative }} |
| {{ val_host_80_donor_20.true_pct }} | {{ val_host_80_donor_20.est_pct }} | {{ val_host_80_donor_20.error_pct | dp(2) }} | {{ val_host_80_donor_20.ci_lo_pct }}–{{ val_host_80_donor_20.ci_hi_pct }} | {{ val_host_80_donor_20.ci_width_pct }} | {{ val_host_80_donor_20.n_informative }} |
| {{ val_host_70_donor_30.true_pct }} | {{ val_host_70_donor_30.est_pct }} | {{ val_host_70_donor_30.error_pct | dp(2) }} | {{ val_host_70_donor_30.ci_lo_pct }}–{{ val_host_70_donor_30.ci_hi_pct }} | {{ val_host_70_donor_30.ci_width_pct }} | {{ val_host_70_donor_30.n_informative }} |
| {{ val_host_60_donor_40.true_pct }} | {{ val_host_60_donor_40.est_pct }} | {{ val_host_60_donor_40.error_pct | dp(2) }} | {{ val_host_60_donor_40.ci_lo_pct }}–{{ val_host_60_donor_40.ci_hi_pct }} | {{ val_host_60_donor_40.ci_width_pct }} | {{ val_host_60_donor_40.n_informative }} |
| {{ val_host_50_donor_50.true_pct }} | {{ val_host_50_donor_50.est_pct }} | {{ val_host_50_donor_50.error_pct | dp(2) }} | {{ val_host_50_donor_50.ci_lo_pct }}–{{ val_host_50_donor_50.ci_hi_pct }} | {{ val_host_50_donor_50.ci_width_pct }} | {{ val_host_50_donor_50.n_informative }} |
| {{ val_host_40_donor_60.true_pct }} | {{ val_host_40_donor_60.est_pct }} | {{ val_host_40_donor_60.error_pct | dp(2) }} | {{ val_host_40_donor_60.ci_lo_pct }}–{{ val_host_40_donor_60.ci_hi_pct }} | {{ val_host_40_donor_60.ci_width_pct }} | {{ val_host_40_donor_60.n_informative }} |
| {{ val_host_30_donor_70.true_pct }} | {{ val_host_30_donor_70.est_pct }} | {{ val_host_30_donor_70.error_pct | dp(2) }} | {{ val_host_30_donor_70.ci_lo_pct }}–{{ val_host_30_donor_70.ci_hi_pct }} | {{ val_host_30_donor_70.ci_width_pct }} | {{ val_host_30_donor_70.n_informative }} |
| {{ val_host_20_donor_80.true_pct }} | {{ val_host_20_donor_80.est_pct }} | {{ val_host_20_donor_80.error_pct | dp(2) }} | {{ val_host_20_donor_80.ci_lo_pct }}–{{ val_host_20_donor_80.ci_hi_pct }} | {{ val_host_20_donor_80.ci_width_pct }} | {{ val_host_20_donor_80.n_informative }} |
| {{ val_host_10_donor_90.true_pct }} | {{ val_host_10_donor_90.est_pct }} | {{ val_host_10_donor_90.error_pct | dp(2) }} | {{ val_host_10_donor_90.ci_lo_pct }}–{{ val_host_10_donor_90.ci_hi_pct }} | {{ val_host_10_donor_90.ci_width_pct }} | {{ val_host_10_donor_90.n_informative }} |
| {{ val_host_0_donor_100.true_pct }} | {{ val_host_0_donor_100.est_pct }} | {{ val_host_0_donor_100.error_pct | dp(2) }} | {{ val_host_0_donor_100.ci_lo_pct }}–{{ val_host_0_donor_100.ci_hi_pct }} | {{ val_host_0_donor_100.ci_width_pct }} | {{ val_host_0_donor_100.n_informative }} |

**Table 1.** In silico validation results across the full range of donor fractions. Synthetic chimeric samples were generated with {{ validation_summary.n_markers | fmt('g') }} markers ({{ validation_summary.n_informative | fmt('g') }} informative) at {{ validation_summary.depth | fmt(',.0f') }}x depth with per-marker bias $\sigma$ = {{ validation_summary.bias_sd }}.

### Effect of Per-Marker Bias Correction

To evaluate the impact of per-marker amplification bias correction, we generated synthetic chimeric samples with known per-marker biases ($\sigma_{bias}$ = 0.02) and ran allomix both with and without bias correction using the true bias values (Table 2).

Without bias correction, the mean absolute error on interior fractions (excluding 0% and 100%) was {{ bias_no_bias.mean_abs_error_pct | dp(2) }}% with RMSE {{ bias_no_bias.rmse_pct | dp(2) }}%. With bias correction, the mean absolute error was {{ bias_with_bias.mean_abs_error_pct | dp(2) }}% with RMSE {{ bias_with_bias.rmse_pct | dp(2) }}%. The most notable improvement from bias correction was at the boundary fractions: without correction, the 0% donor sample was estimated at 0.44%, while with correction the estimate was 0.00%. Similarly, the 100% donor sample improved from 99.37% to 100.00% (Figure 2).

![**Figure 1.** In silico validation of allomix. True donor fraction (x-axis) versus estimated donor fraction (y-axis) for synthetic chimeric samples with {{ validation_summary.n_markers | fmt('g') }} markers at {{ validation_summary.depth | fmt(',.0f') }}x depth. Dashed line indicates perfect agreement.](paper/figures/fig1_scatter.png)

![**Figure 2.** Effect of per-marker bias correction. Left: truth versus estimated donor fraction without (red) and with (blue) bias correction. Centre: residuals. Right: 95% confidence intervals. Per-marker biases were simulated with $\sigma_{bias}$ = 0.02.](paper/figures/fig2_bias_correction.png)

### Confidence Interval Calibration

The 95% profile likelihood confidence intervals showed nominal coverage of {{ validation_summary.ci_coverage_pct }}% in the standard validation. This under-coverage relative to the 95% nominal level is attributable to systematic per-marker capture biases that are not accounted for in the standard error model. When bias correction is applied, CI calibration is expected to improve, as the residual error after correction is dominated by stochastic sampling noise that is well-modeled by the binomial likelihood.

![**Figure 3.** Confidence interval coverage. Vertical bars show 95% profile likelihood CIs for each simulated donor fraction. Blue bars contain the true value; red bars miss it. Black dots are point estimates.](paper/figures/fig3_ci_coverage.png)

### Comparison with Existing Tools

Table 3 summarizes the landscape of available NGS-based chimerism tools. allomix is the only open-source option that is panel-agnostic and works from standard VCF files.

| Tool | Markers | LOD | Open Source | Panel Agnostic | Input |
|:---|:---:|:---:|:---:|:---:|:---:|
| AlloSeq HCT | 202 SNPs | 0.3% | No | No | Proprietary |
| Devyser Chimerism | 24 indels | 0.05% | No | No | Proprietary |
| NGStrack | 34 indels | 0.1% | No | No | Proprietary |
| ScisGo Chimerism MD | SNPs + indels | 0.2–0.5% | No | No | Proprietary |
| **allomix** | **Any biallelic** | **<1%** | **Yes (MIT)** | **Yes** | **VCF** |

**Table 3.** Comparison of NGS-based chimerism monitoring tools. LOD = limit of detection.

{# TODO: Validation with real sequencing data #}
{# Once real patient/control data is processed through the pipeline, add: #}
{# - Concordance with STR-based chimerism results #}
{# - Dilution series from cell line mixing experiments #}
{# - Multi-donor detection capability #}
{# - Lineage-specific chimerism (CD3, CD33, etc.) if sorted fractions available #}
