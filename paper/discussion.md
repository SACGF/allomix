## Discussion

We have presented allomix, an open-source tool that enables donor chimerism monitoring by repurposing biallelic polymorphic markers already present in clinical NGS panels. By extracting chimerism information from existing sequencing data, allomix allows laboratories to add post-HSCT monitoring to their workflows without a dedicated assay, additional reagents, or proprietary software.

### Accuracy and Performance

In silico validation demonstrated sub-percentage-point mean absolute error across all sequencing depths tested (50x to 1,000x), with MAE ranging from {{ depth_50.mean_abs_error_pct | dp(2) }}% at 50x to {{ depth_1000.mean_abs_error_pct | dp(2) }}% at 1,000x depth. Maximum single-sample errors ranged from {{ depth_50.max_abs_error_pct | dp(2) }}% at 50x to {{ depth_1000.max_abs_error_pct | dp(2) }}% at 1,000x. This level of accuracy is competitive with published performance metrics from commercial tools. Kakodkar et al. reported mean absolute errors of 0.3–1.5% for AlloSeq HCT across various donor fractions,[@Kakodkar2023alloseq] while Pedini et al. demonstrated comparable precision for the Devyser system at high sequencing depths.[@Pedini2021devyser]

The largest estimation errors were observed at boundary fractions (0% and 100% donor), where per-marker amplification biases have the greatest relative impact. This is consistent with the observations of De Vynck et al., who demonstrated that systematic capture biases cause the largest quantitative distortion when the true minor component is near zero.[@Vynck2023bias] When bias correction was applied with known per-marker biases, boundary estimates improved (e.g., 0% donor: {{ bias_no_bias.est_0pct }}% uncorrected vs {{ bias_with_bias.est_0pct }}% corrected).

### Statistical Framework

The MLE framework employed by allomix, adapted from Crysup and Woerner's mixture deconvolution model,[@CrysupWoerner2022] offers several advantages for the clinical chimerism setting. The use of known donor and host genotypes — available from pre-transplant samples — dramatically simplifies the estimation problem compared to forensic mixture analysis, where contributor profiles must be jointly inferred. This simplification yields faster computation, more interpretable results, and tighter confidence intervals.

The two-stage optimization strategy (coarse grid search followed by Brent refinement) provides robust global optimization without the risk of converging to local optima, which can occur with gradient-based methods alone. Profile likelihood confidence intervals provide a principled measure of estimation uncertainty that naturally accounts for the information content of the marker panel and sequencing depth.

### Confidence Interval Calibration

The observed CI coverage ({{ depth_1000.ci_coverage_pct }}–{{ depth_50.ci_coverage_pct }}% across depths, versus a nominal 95%) reflects the impact of systematic noise sources — per-marker capture biases, non-uniform depth, and locus dropout — that the binomial likelihood model does not account for. The current model assumes that allele counts are drawn from a binomial distribution with a single error rate parameter, but does not model marker-specific systematic shifts in allele efficiency. This is a known limitation shared with other MLE-based approaches.[@Vynck2023bias] Approaches to improving CI calibration include incorporating bias correction (demonstrated here), overdispersion modeling via beta-binomial likelihoods, and empirical CI recalibration from training data.

### Repurposing Existing Panels

The central advantage of allomix is that it works with markers laboratories are already sequencing. Sample identification SNPs, pharmacogenomic markers, and other polymorphic loci included in clinical NGS panels for quality control or diagnostic purposes can serve double duty for chimerism monitoring. This eliminates the cost and logistical overhead of running a separate dedicated chimerism assay. Lee et al. demonstrated this principle by extracting chimerism from 121 SNPs embedded in a myeloid neoplasm panel,[@Lee2019snp] but the analysis required custom scripting with no reusable tool. allomix generalises this approach into a tool any laboratory can deploy.

De Vynck et al. have shown that as few as 3 informative markers are sufficient for chimerism quantification, though accuracy improves with additional markers.[@Vynck2022markers] Their FABCASE tool can be used to prospectively assess whether a given panel contains enough informative markers for a specific donor-host pair.[@Vynck2025fabcase] For panels with tens of polymorphic markers — typical of sample ID marker sets — the probability of having sufficient informative markers exceeds 99% even for sibling donor-host pairs.

### Clinical Workflow Considerations

allomix is designed for integration into standard clinical genomics pipelines. By accepting VCF files as input rather than raw sequencing data, it decouples chimerism analysis from upstream alignment and variant calling, allowing laboratories to use their established bioinformatics infrastructure. The tool's output formats (TSV and JSON) facilitate integration with laboratory information management systems and clinical reporting workflows.

A critical implementation consideration is the requirement for joint variant calling of admixture samples alongside donor and host genotyping samples. When admixture samples are called independently, sites where the minor component is below the variant calling threshold may be reported as homozygous-reference with only a single-element AD field, losing the informative read count data. Joint calling ensures that all markers genotyped in the donor are represented with full AD information in the admixture sample, even at very low donor fractions.

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

allomix enables laboratories to repurpose polymorphic markers already present in their clinical NGS panels for donor chimerism monitoring after HSCT, achieving sub-percentage-point accuracy in silico without requiring a dedicated assay, additional reagents, or proprietary software. Clinical validation studies are underway.
