## Discussion

allomix enables laboratories to add post-HSCT chimerism monitoring to existing NGS workflows by repurposing polymorphic markers they are already sequencing, without a dedicated assay or proprietary software.

### Accuracy and Performance

The <1% MAE across all depths tested is competitive with published performance data for commercial tools: Kakodkar et al. reported 0.3--1.5% MAE for AlloSeq HCT,[@Kakodkar2023alloseq] and Pedini et al. demonstrated comparable precision for the Devyser system.[@Pedini2021devyser] The largest errors occurred at boundary fractions (0% and 100% donor), where per-marker amplification biases have the greatest relative impact, consistent with De Vynck et al.[@Vynck2023bias] Bias correction reduced these boundary errors (e.g., 0% donor: {{ bias_no_bias.est_0pct }}% uncorrected vs {{ bias_with_bias.est_0pct }}% corrected).

### Confidence Interval Calibration

The observed CI undercoverage ({{ depth_1000.ci_coverage_pct }}--{{ depth_50.ci_coverage_pct }}% versus a nominal 95%) reflects systematic noise sources that the binomial likelihood model does not capture, a known limitation of MLE-based approaches.[@Vynck2023bias] Bias correction (demonstrated here) partially addresses this; beta-binomial likelihoods and empirical recalibration are planned improvements.

### Repurposing Existing Panels

The central advantage of allomix is that it works with markers laboratories are already sequencing. Sample identification SNPs, pharmacogenomic markers, and other polymorphic loci included in clinical NGS panels for quality control or diagnostic purposes can serve double duty for chimerism monitoring. This eliminates the cost and logistical overhead of running a separate dedicated chimerism assay. Lee et al. demonstrated this principle by extracting chimerism from 121 SNPs embedded in a myeloid neoplasm panel,[@Lee2019snp] but the analysis required custom scripting with no reusable tool. allomix generalises this approach into a tool any laboratory can deploy.

De Vynck et al. have shown that as few as 3 informative markers are sufficient for chimerism quantification, though accuracy improves with additional markers.[@Vynck2022markers] Their FABCASE tool can be used to prospectively assess whether a given panel contains enough informative markers for a specific donor-host pair.[@Vynck2025fabcase] For panels with tens of polymorphic markers, typical of sample ID marker sets, the probability of having sufficient informative markers exceeds 99% even for sibling donor-host pairs.

### Clinical Workflow Considerations

By accepting standard VCF files, allomix decouples chimerism analysis from upstream alignment and variant calling, allowing laboratories to use their existing bioinformatics infrastructure. One implementation requirement is that admixture samples should be joint-called alongside donor and host samples (see Methods); this is straightforward where GATK pipelines are already in place but does require pipeline configuration.

### Limitations and Future Directions

Several limitations of the current work should be noted. First, the validation presented here is entirely in silico; clinical validation against STR-based chimerism results and with controlled cell-line dilution series is required before clinical deployment. Second, the limit of detection has not been formally characterized following AMP guidelines, though the mathematical framework supports detection below 1% with sufficient depth and informative markers.

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
