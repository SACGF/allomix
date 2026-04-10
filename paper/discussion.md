## Discussion

allomix enables laboratories to add post-HSCT chimerism monitoring to existing NGS workflows by repurposing polymorphic markers they are already sequencing, without a dedicated assay or proprietary software.

### Accuracy and Performance

The <1% MAE across all depths tested is competitive with published performance data for commercial tools: Kakodkar et al. reported 0.3--1.5% MAE for AlloSeq HCT,[@Kakodkar2023alloseq] Pedini et al. demonstrated comparable precision for the Devyser system,[@Pedini2021devyser] and Blouin et al. reported R^2 = 0.9987 across the full chimerism range for the ScisGo assay.[@Blouin2024comparison] Qama et al. validated the Devyser assay with a limit of detection of 0.06% and high concordance with STR (R^2 = 0.998), and showed that NGS detected residual host DNA (>0.1%) in 85% of samples classified as full donor chimerism (>95%) by STR, compared with only 5% detection by STR.[@Qama2026devyser] This finding suggests that current definitions of full donor chimerism based on STR sensitivity thresholds may underestimate residual host haematopoiesis, and highlights the clinical value of sensitive chimerism monitoring at all engraftment levels. The largest errors in our simulations occurred at boundary fractions (0% and 100% donor), where per-marker amplification biases have the greatest relative impact, consistent with De Vynck et al.[@Vynck2023bias] Bias correction reduced these boundary errors (e.g., 0% donor: {{ bias_no_bias.est_0pct }}% uncorrected vs {{ bias_with_bias.est_0pct }}% corrected).

### High-Depth Regime and Overdispersion

The likelihood model of Crysup and Woerner was evaluated at read depths of 2--100, where sampling noise dominates. Clinical targeted panels operate at 500--2,000x or higher, a regime where per-marker systematic biases (amplification efficiency, GC content, capture probe affinity) become the dominant source of variance. This motivates the per-marker bias correction and overdispersion modelling in allomix, which are not addressed in the original derivation.

### Confidence Interval Calibration

The observed CI undercoverage ({{ depth_1000.ci_coverage_pct }}--{{ depth_50.ci_coverage_pct }}% versus a nominal 95%) reflects systematic noise sources that the binomial likelihood model does not capture, a known limitation of MLE-based approaches.[@Vynck2023bias] Bias correction (demonstrated here) partially addresses this; beta-binomial likelihoods and empirical recalibration are planned improvements.

### Repurposing Existing Panels

The central advantage of allomix is that it works with markers laboratories are already sequencing. Sample identification SNPs, pharmacogenomic markers, and other polymorphic loci included in clinical NGS panels for quality control or diagnostic purposes can serve double duty for chimerism monitoring. This eliminates the cost and logistical overhead of running a separate dedicated chimerism assay. Lee et al. demonstrated this principle by extracting chimerism from 121 SNPs embedded in a myeloid neoplasm panel,[@Lee2019snp] but the analysis required custom scripting with no reusable tool. allomix generalises this approach into a tool any laboratory can deploy.

De Vynck et al. have shown that as few as 3 informative markers are sufficient for chimerism quantification, though accuracy improves with additional markers.[@Vynck2022markers] Their simulation study found that panels of approximately 20 markers with MAFs near 0.5 provide a >95% probability of yielding at least 3 informative markers even for sibling donor-host pairs, and their FABCASE tool can be used to prospectively assess panel sufficiency for a specific donor-host pair.[@Vynck2025fabcase] Panels with tens of polymorphic markers, typical of sample ID marker sets, are therefore expected to provide adequate informativity for most clinical scenarios.

### Multi-Donor Chimerism

Multi-donor transplants (cord blood, sequential transplants) are increasingly common and require simultaneous quantification of multiple donor fractions. Blouin et al. validated multi-donor chimerism with the ScisGo assay using clinical samples, reporting analytical sensitivity of 0.5% for double-donor detection and successful quantification of triple-donor transplants, though with reduced informative marker counts (average 21 for double donors, 8 for triple donors from >200 markers).[@Blouin2024comparison] Our in silico multi-donor validation with sibling donor pairs demonstrated <2% per-donor MAE at 500x depth, with correct ranking of asymmetric donor fractions in all cases. The reduced number of informative markers when donors are related to each other (as siblings sharing both parents) represents the most challenging scenario for multi-donor estimation; unrelated multi-donor settings such as cord blood transplants would yield more informative markers and correspondingly better precision.

### Clinical Workflow Considerations

By accepting standard VCF files, allomix decouples chimerism analysis from upstream alignment and variant calling, allowing laboratories to use their existing bioinformatics infrastructure. One implementation requirement is that admixture samples should be joint-called alongside donor and host samples (see Methods); this is straightforward where GATK pipelines are already in place but does require pipeline configuration.

### Limitations and Future Directions

Several limitations of the current work should be noted. First, the validation presented here is entirely in silico; clinical validation against STR-based chimerism results and with controlled cell-line dilution series is required before clinical deployment. Blouin et al. describe a practical framework for clinical validation of NGS chimerism assays, including run-level quality metrics (cluster density, base call quality, coverage uniformity) and sample-level acceptance criteria (minimum coverage, LOD thresholds, confidence interval width), which provides a useful model for future allomix validation studies.[@Blouin2024comparison] Second, the limit of detection has not been formally characterized following AMP guidelines, though the mathematical framework supports detection below 1% with sufficient depth and informative markers.

{# TODO: Add results from clinical validation once available #}
{# Planned validation studies: #}
{# - Concordance with STR chimerism on retrospective patient samples #}
{# - Controlled dilution series using cell lines or DNA mixtures #}
{# - Multi-donor detection with cord blood transplant samples #}
{# - Comparison of LOD with commercial tools on matched samples #}
{# - VariantGrid integration for production deployment #}

Future development priorities include formal analytical validation following AMP guidelines, integration with the VariantGrid clinical genomics platform for production deployment, and longitudinal monitoring features including trend analysis and alerting for clinically significant chimerism changes.

## Conclusions

allomix enables laboratories to repurpose polymorphic markers already present in their clinical NGS panels for donor chimerism monitoring after HSCT, achieving <1% mean absolute error in silico without requiring a dedicated assay, additional reagents, or proprietary software. Clinical validation studies are underway.
