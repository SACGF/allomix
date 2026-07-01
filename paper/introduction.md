## Introduction

Allogeneic hematopoietic stem cell transplantation (HSCT) is a curative therapy for a
range of hematologic malignancies and bone marrow failure syndromes. Following
transplantation, monitoring the proportion of donor-derived hematopoietic cells (termed
chimerism) is critical for clinical decision-making. Declining donor chimerism may
indicate graft rejection or impending disease relapse, prompting therapeutic
interventions such as donor lymphocyte infusion or immunosuppression
modification.[@Clark2015str]

The current standard for chimerism monitoring relies on short tandem repeat (STR)
analysis by capillary electrophoresis.[@Clark2015str] STR-based methods have
well-characterized limitations: sensitivity of 1--5% for minor component detection,
limited quantitative precision from the narrow dynamic range of fragment analysis, and
stutter artifacts that complicate interpretation. STR chimerism also requires a
dedicated laboratory workflow (separate sample handling, PCR, and capillary
electrophoresis) that runs independently of the laboratory's NGS-based testing and so
cannot share its economies of scale. Chimerism is
intrinsically a property of a cell population rather than of whole blood per se, and
clinical practice often incorporates lineage-specific testing of sorted cell subsets
when whole-blood signal is less informative or ambiguous.[@KharfanDabaja2021astct;
@Clark2025bjh; @Kakodkar2023alloseq]

Next-generation sequencing (NGS) approaches to chimerism monitoring offer the potential
to address these limitations through improved quantitative precision via digital allele
counting and lower limits of detection through deep sequencing. Several groups have
demonstrated the feasibility of SNP-based NGS chimerism using targeted panels ranging
from 24 to 202 markers.[@Aloisio2016amplicon; @Lee2019snp; @Vynck2021devyser;
@Kakodkar2023alloseq; @ZhangA2024comparison; @Blouin2024comparison; @Qama2026devyser]
Recent clinical validations have confirmed analytical sensitivity of 0.06--0.2% for
single-donor detection, with NGS detecting residual host DNA in the majority of samples
classified as full donor chimerism by STR.[@Qama2026devyser; @Blouin2024comparison]
Head-to-head comparisons have shown that NGS yields substantially more informative
markers than STR or quantitative PCR (qPCR) approaches, with corresponding improvements
in sensitivity and dynamic range.[@ZhangA2024comparison] The clinical value of this
improved sensitivity is supported by evidence that increasing mixed chimerism detected
by sensitive methods is predictive of leukemia relapse.[@ZhangR2024relapse] However, all
currently available NGS chimerism tools, including AlloSeq HCT (CareDx,
{{ tool_landscape.alloseq_n_markers }} SNPs), Devyser Chimerism (Thermo Fisher,
{{ tool_landscape.devyser_n_markers }} indels), NGStrack (GenDx), and ScisGo Chimerism
MD (Scisco Genetics), are proprietary commercial products, each tightly coupled to its
respective assay kit.[@Vynck2021devyser; @Blouin2024comparison; @Kakodkar2023alloseq;
@Qama2026devyser]

This creates a missed opportunity. Many clinical laboratories already run targeted NGS
panels for haematological malignancy characterisation, pharmacogenomics, or other
diagnostic purposes that include biallelic polymorphic markers, often for sample
identification or quality control. The markers most useful for chimerism are common ones
(high minor-allele frequency, so host and donor are more often genotypically different)
that sit in low linkage disequilibrium with each other (so each adds independent
information); panels designed for sample identification tend to have both properties
already. For example, the IDT rhAmpSeq Sample ID panel
includes {{ panel_specs.n_markers_panel }} SNPs designed for sample tracking, and
similar marker sets are incorporated into numerous clinical capture and amplicon panels.
These markers are sequenced at high depth as part of routine clinical workflows, and the
resulting data is already available. Lee et al. demonstrated this principle by using 121
SNPs embedded in a myeloid neoplasm panel for chimerism monitoring, but the analysis
required custom scripting with no reusable tool.[@Lee2019snp]

The statistical machinery for quantifying DNA mixtures from allelic read depths is well
established in the forensic genetics literature, including maximum likelihood estimation
for biallelic SNP mixtures[@CrysupWoerner2022; @Woerner2024demixtify] and
characterisation of panel informativity and per-marker amplification
bias.[@Vynck2022markers; @Vynck2023bias; @Vynck2025fabcase] The clinical chimerism
setting offers a key simplification: donor and host genotypes are known from
pre-transplant samples, so the tool does not have to reconstruct unknown contributor
profiles and can instead concentrate on the low-fraction signal that matters clinically.

Here we present allomix, an open-source tool that performs donor chimerism monitoring
from any set of biallelic markers already present in a laboratory's clinical NGS
workflow. Rather than requiring a dedicated chimerism assay, allomix repurposes the
polymorphic markers that laboratories are already sequencing. It answers two distinct
clinical questions with two complementary tests that read different markers: a magnitude
estimate (how much donor is present, with a confidence interval) and a residual-host
presence test (whether any host remains, even below the level the magnitude estimate can
quantify). Because the tool is meant to run inside routine laboratory operations, it
also carries a built-in set of quality and sample-integrity checks that flag sample
swaps, contamination from a genome other than the host or donor (for example another
patient co-pooled on the same flowcell), and unexpected donor-host relatedness from the
same marker data. We describe the approach in plain terms, present in silico validation
together with a demonstration on a public dataset of real titrated DNA mixtures, and
discuss the practical considerations for integrating chimerism monitoring into existing
NGS pipelines. The scope of this paper is the tool, its estimation method, and an
analytical characterization of how it behaves under controlled noise, with the
real-mixture analysis as an independent check on reads we did not simulate. Clinical
validation against STR chimerism on patient cohorts is a separate study and is
deliberately out of scope here. Because allomix is panel-agnostic, there is no single
assay to validate once: clinical performance depends on the specific marker set,
sequencing chemistry, and specimen type a laboratory uses, so validation is necessarily
done per laboratory and per panel, in the same way any laboratory-developed test is
validated locally before use. A single validation on one panel would not transfer to the
others the tool is meant to serve.
