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
