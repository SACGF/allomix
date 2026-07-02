## Discussion

allomix lets a laboratory add post-HSCT chimerism monitoring to an NGS panel it already
runs, by repurposing polymorphic markers that are sequenced anyway, without a dedicated
assay or proprietary software. It does this with two complementary readouts (how much
donor is present, and whether any host remains) and a set of sample-integrity checks
built from the same marker data.

### Accuracy and sensitivity

The <1% mean absolute error across simulated depths (Supplementary Table S4), together
with the recovery of known fractions on the real SRP434573 titrated mixtures (Results),
is competitive with published performance for commercial tools: Kakodkar et al. reported
an analytical LoD of 0.3% for AlloSeq
HCT,[@Kakodkar2023alloseq] Pedini et al. and Vynck et al. showed comparable precision
for the Devyser system,[@Pedini2021devyser; @Vynck2021devyser] and Blouin et al.
reported R^2 = 0.9987 across the full chimerism range for ScisGo.[@Blouin2024comparison]
Qama et al. validated the Devyser assay at a LoD of 0.06% with high STR
concordance (R^2 = 0.998), and found that NGS detected residual host DNA (>0.1%) in 85%
of samples with >95% donor chimerism, against 5% by
STR.[@Qama2026devyser] That gap suggests STR-based definitions of full donor chimerism
can miss residual host haematopoiesis, and it is the clinical reason a sensitive method,
and a dedicated residual-host test, are worth having. Per-marker bias correction in
allomix is a precision refinement rather than an accuracy fix: it sharpens the
confidence intervals but barely moves the point estimate (Methods, Results), consistent
with bias having its largest relative effect at the extreme expected frequencies of
boundary fractions.[@Vynck2023bias]

### Overdispersion is the real limiter at clinical depth

The likelihood of Crysup and Woerner was evaluated at read depths of 2--100, where
sampling noise dominates. Clinical targeted panels run at 500--2,000x or higher, where
per-marker systematic effects (amplification efficiency, GC content, capture-probe
affinity) become the dominant source of variance, which is why allomix adds per-marker
bias correction and overdispersion modelling that the original derivation does not.
Overdispersion is also the dominant control on the LoD at high depth: as
depth grows, the per-marker variance stops falling and the effective depth is capped, so
the LoD saturates at a floor rather than continuing to improve as 1/√n (Supplementary
Figures S7, S8). At the overdispersion fitted from our real mixtures (heterozygous-class
$\rho \approx$ {{ overdispersion_lod_headline.real_rho_het_median }}) the in silico LoD is
roughly {{ overdispersion_lod_headline.fold_real_rho_vs_binomial }}-fold higher than the
pure-binomial value (rising to
{{ overdispersion_lod_headline.lod_at_real_rho_pct }}%, near the ~1% measured directly on
the same mixtures), so beyond that regime adding depth buys little and the panel's
overdispersion, not its coverage, sets the achievable LoD. Marker count is a separate
hard limit: depth and overdispersion act on the markers a panel already has, and no
amount of depth compensates for too few informative markers, so panel size floors the
achievable LoD independently of coverage (Results, relatedness). The noise-ablation
analysis makes the same point on point-estimate accuracy: adding overdispersion to the
full noise model raises RMSE more than amplification bias, depth non-uniformity, or
sequencing error individually ({{ supp_synthetic.ablation_rmse_overdispersion_pct }}%
versus {{ supp_synthetic.ablation_rmse_full_pct }}% for the otherwise-identical binomial
model; Supplementary Figure S4). This is what explains the gap between an
idealised in silico LoD and real-panel performance, and it is the same overdispersion
that raises the goodness-of-fit statistic on the real SRP434573 mixtures, though the QC
layer promotes a fit to REVIEW only when that misfit is large rather than merely
significant at depth (Supplementary Methods S10). For the same reason allomix does not
weight marker contributions by per-base or per-genotype quality (Supplementary Methods
S11).

Across the simulated depths the 95% intervals cover close to nominal
({{ depth_1000.ci_coverage_pct }}--{{ depth_200.ci_coverage_pct }}%) with overdispersion
fit separately for the donor-heterozygous and donor-homozygous marker classes (the default
per-marker-type model). A residual calibration gap remains on the real
titrated ladder, where the intervals mildly undercover the nominal fraction. The driver
is the co-pooled donor-homozygous contamination background, a dose-dependent excess of
host-allele reads that the per-marker correction reduces but does not fully remove,
together with real heterozygous-marker scatter that is not an exact beta-binomial,
rather than a shortfall of the model in controlled data. The realised mixing itself is
faithful: the donor-homozygous markers alone return the nominal fraction. Empirical
recalibration from training data is one route to closing that gap.

### Outlier-resistant estimation at the limit of detection

The one-sided outlier-resistant refit (Methods) matters most at the limit of detection.
When the true host fraction is around 1%, the handful of markers actually carrying that
signal sit off a donor-dominated fit and read as outliers, so a symmetric outlier rule
would discard exactly the markers of interest and collapse the estimate to zero. Removing
only the markers whose deviation points away from host presence, while protecting those
pointing toward it, is what recovers the donor LoD against an aberration-bearing recipient
background while leaving aberration-free samples unchanged (Results, copy-number stress
test).

### Marker independence and linkage disequilibrium

The likelihood is a product over markers, which assumes the markers are statistically
independent. Linkage disequilibrium (LD) between nearby markers would break that
assumption: correlated markers carry less information than the same number of
independent ones, so the point estimate stays unbiased but the profile-likelihood
intervals are too narrow. In practice this is governed by panel design.
Sample-identification and ancestry SNP sets, the panels allomix is meant to repurpose,
are chosen to sit far apart across the genome in approximate linkage equilibrium
precisely so each marker is an independent identity bit, and both panels used here (the
76-SNP rhAmpSeq sample-ID set and the {{ srp434573.panel_n_snps | commas }}-SNP MIP
panel) are of that design. LD only becomes a concern on panels with deliberately
clustered markers (several SNPs per amplicon, or dense tiling), which is a panel choice
a laboratory controls. A laboratory whose panel does carry clustered markers can restrict
allomix to an LD-pruned subset with a simple BED file or marker list, keeping one marker
per linked cluster, so the independence assumption holds without changing the panel.
Residual correlation does not fail silently: it inflates
per-marker variance in the same way overdispersion does, so it is partly absorbed by the
per-class overdispersion model and surfaces in the goodness-of-fit statistic. A marker
set dense enough for LD to matter would raise that statistic far enough to cross the
effect-size threshold and be caught as an overdispersed fit (Supplementary Methods S10)
rather than reported as a confidently wrong interval. LD can also be
turned to advantage: Kim et al. used phased SNPs in tight LD on the same read as a
concordance filter, requiring the linked alleles to agree before a read counts, which
suppresses sequencing error and pushes the detectable host fraction lower.[@Kim2021ld]
That haplotype-aware extension is a sensitivity upgrade for the residual-host presence
test rather than a correctness fix, and it depends on multiple phased markers per
amplicon, so we leave it to future work.

### Two complementary tests and the clinical relapse question

The clinical motivation for sensitivity below the STR threshold is early detection of
relapse, where residual or returning host DNA appears at fractions a magnitude estimate
cannot yet quantify. allomix addresses this with a separate residual-host presence test
that reads donor-homozygous markers where the host carries the donor-absent allele,
asking only whether that host-only allele exceeds the background-artifact floor. We
validated it as a capability (in silico, and on SRP434573 down to 1% host) but we do not
present it as a demonstrated relapse-detection result: the clinical thresholds that
would turn a positive presence call into an action, and its operating characteristics on
real patient samples, still have to be established. The design intent is that the two
tests are read together, the magnitude estimate for trajectory and the presence test for
the low-fraction tail, with the tool warning when the presence test detects host below
what the magnitude estimate resolves. On real data the test's sensitivity is set by how
well the per-marker error background is known: allomix can take a per-site,
per-direction empirical error table (Methods) so each marker's background is its own
measured REF-to-ALT or ALT-to-REF rate rather than a global average, which is the route
to pushing the detectable host fraction lower on a given panel. We used the global rate
here because the simulator and the SRP434573 reads do not exercise a calibrated per-site
table, so its benefit is a real-cohort question still to be measured.

### A safety suite for routine deployment

For a clinical laboratory, catching sample mix-ups matters, because a quantitative
chimerism number is misleading if it is computed on the wrong sample. allomix turns the
same marker data into a set of integrity checks: a relatedness check that flags sample
swaps from unexpected kinship, a consensus-site consistency test that catches a
wrong-patient VCF the fraction estimate cannot see, and an in-data contamination estimate
that separates a real low-fraction signal from co-pooled or foreign DNA. These are kept
distinct by reading marker sets defined by genotype geometry rather than by
re-thresholding one statistic (Methods), so they act as independent checks rather than
competing interpretations of the same number, all from data the assay already produces.

### Repurposing existing panels

The central advantage of allomix is that it works with markers laboratories already
sequence. Sample-identification SNPs, pharmacogenomic markers, and other polymorphic
loci included for quality control or diagnosis can serve double duty for chimerism,
eliminating a separate dedicated assay. Lee et al. demonstrated the principle with 121
SNPs in a myeloid panel but released no reusable tool;[@Lee2019snp] allomix generalises
this into a deployable tool. Vynck et al. showed
that three informative markers are enough to make a quantification identifiable, with
accuracy improving as markers are added, and that panels of about 20 markers with MAFs
near 0.5 give a >95% chance of at least three informative markers even for sibling pairs;
their FABCASE tool can assess panel sufficiency prospectively for a specific donor-host
pair.[@Vynck2022markers; @Vynck2025fabcase] Three markers is the floor for producing an
estimate at all, not a target for sensitivity: the limit of detection keeps falling as
informative markers are added (Figure 4), so reaching the low fractions that matter
clinically needs a panel on the order of 100 markers, within the roughly 24 to 202
markers commercial chimerism assays carry (Table 1), which yields the tens of informative
markers the low fractions require rather than a handful.

### Multi-donor chimerism

Multi-donor transplants (cord blood, sequential transplants) are increasingly common and
need simultaneous quantification of multiple donor fractions. Blouin et al. validated
multi-donor chimerism with ScisGo on clinical samples, reporting 0.5% sensitivity for
double-donor detection and successful triple-donor quantification, though with reduced
informative-marker counts.[@Blouin2024comparison] Our in silico sibling multi-donor
validation showed <2% per-donor MAE with correct ranking of asymmetric fractions;
sibling donors are the hardest case, and unrelated multi-donor settings would yield more
informative markers and better per-donor precision.

### Clinical workflow and cellular composition

By accepting standard VCFs, allomix is not tied to a single upstream alignment and
variant-calling pipeline, though it does impose specific requirements on how those VCFs
are produced. The requirement is the two-phase calling in Methods: genotypes from joint
germline calling (any joint germline caller, for example GATK), and admixture allele
depths from a forced bcftools pileup at the panel sites rather than from that germline
caller, so the minority alternative reads carrying the low-fraction signal are retained
rather than stripped. The germline-calling phase can use whichever joint caller a
laboratory already runs, but the admixture pileup must be configured as a separate step
rather than reusing one joint-calling pass for every sample. We have not validated
alternative callers against this workflow, so a laboratory substituting its own would
confirm the low-fraction reads survive, as documented for the pileup step in Methods.

Chimerism is a property of a cell population, not of whole blood per se. A donor fraction
measured in unfractionated blood is a lineage-abundance-weighted average, so a meaningful
change confined to a small lineage may show only a fractional whole-blood change, which
is why lineage-specific testing on sorted subsets is recommended in some clinical
contexts.[@Clark2025bjh; @KharfanDabaja2021astct; @Kakodkar2023alloseq] allomix is
lineage-agnostic by design: it estimates the donor fraction in whatever DNA the input VCF
represents, so the same validation applies to cell-sorted inputs. Where whole-blood
signal is adequate, sensitivity better than STR can reduce the need for sorting; where a
change is confined to a small lineage, sorting may still be required. Specimen choice is a
clinical decision outside the tool, and the tool's analytical sensitivity is necessary
but not sufficient for clinical sensitivity in any given scenario.

### Limitations and future directions

The validation here is primarily in silico, by the scope of this paper. The contribution
is the tool and its analytical characterization, with one demonstration on real reads
(the SRP434573 titrated mixtures), which recovered known fractions from 10% down to 1%
and resolved a three-person mixture. That is a co-pooled research panel with its own
contamination floor, an independent check on reads we did not simulate rather than a
substitute for clinical validation.

That SRP434573 floor is consistent with the source study's own analysis: working from
the UMI-tagged reads, Chu identified index misassignment as a low-level cross-sample
contamination in this dataset and controlled it with a per-sample reads-per-UMI
threshold.[@Chu2024mipsnp] The public FASTQs carry no UMI bases, so allomix runs on raw
depth and cannot apply that filter, which makes this dataset close to a worst case for a
pooled panel: a patterned-flowcell instrument where index hopping is most pronounced,
seven co-pooled individuals, and no UMI correction available. allomix instead separates
the floor from true signal in the data, by the carrier-dose response, and reaches a
residual floor near 0.5% host, comparable to what the UMI-based pipeline achieved for
minor-contributor identification. This is why allomix does not require UMIs: where they
are present they remain a useful orthogonal way to push the floor lower, but the
index-hopping contribution can also be held down at the indexing-design level (unique
dual indexes) or, as here, subtracted at each marker by the carrier-dose correction
above, which estimates the per-carrier contamination rate from the run's own
consensus-homozygous markers and removes it before the fit.

Clinical validation against STR chimerism, with
controlled cell-line dilution series, is a separate study. We are planning such a
validation at SA Pathology, running allomix in tandem with the accredited STR chimerism
assay across a patient cohort, to be reported separately. Because allomix is
panel-agnostic, that validation does not transfer between panels: a laboratory adopting
the tool validates it on its own panel and specimen types, as for any
laboratory-developed test, rather than relying on a single validation done once for the
tool. Blouin et al. describe a practical framework for such validation,
including run-level metrics and sample-level acceptance criteria, that is a useful model
for future allomix studies.[@Blouin2024comparison]

The residual-host presence test is a per-sample test. Under serial monitoring, applying
its per-timepoint p < 0.05 threshold repeatedly across a patient's timeline does not
control a family-wise or trend-level error rate, so the per-sample specificity does not
carry over unchanged to a multi-timepoint decision. A monitoring-context decision rule,
such as confirmation on a repeat draw or a trend test across timepoints, is a
clinical-validation-stage choice we leave to the adopting laboratory.

The LoD reported here is a best-case analytical figure, not a validated
assay limit (a real assay's LoD can only be higher). It is the CLSI EP17-A2
95%-detection criterion[@CLSIEP17A2] applied to simulated data under a noise model
calibrated from empirical panel data, with reads drawn from a binomial. Table 1 lays out
the full evidence ladder, from the ~1% real-data LoD on the SRP434573 mixtures to the
sub-0.2% analytical figures in near-binomial simulation, and places each against the
commercial assays. The gap between the analytical figure and real-panel performance is
expected to be driven largely by overdispersion (above), and wetlab validation will set
the floor allomix actually delivers in routine use. As a partial
check against real reads, we
sub-sampled the high-depth SRP434573 mixtures across depth and panel size and measured
the LoD on those reads directly (Figure 3): it falls with depth and panel size as the
simulated curves do and tracks them within the dilution grid, a confirmation of the
simulated trend rather than a replacement for the wetlab limit. Finally, the in silico work
characterises the estimator on bulk DNA mixtures and does not address the cellular
composition of the source specimen; clinical validation should include matched
whole-blood and lineage-sorted samples.

{# TODO: Add results from clinical validation once available #}
{# Planned validation studies: #}
{# - Concordance with STR chimerism on retrospective patient samples #}
{# - Controlled dilution series using cell lines or DNA mixtures #}
{# - Multi-donor detection with cord blood transplant samples #}
{# - Comparison of LOD with commercial tools on matched samples #}

Future development priorities include formal analytical validation following AMP
guidelines and longitudinal monitoring features including trend analysis and alerting
for clinically significant chimerism changes.

## Conclusions

allomix lets laboratories repurpose polymorphic markers already present in their
clinical NGS panels for donor chimerism monitoring after HSCT, reporting both how much
donor is present and whether any host remains, with built-in sample-integrity checks and
<1% mean absolute error in silico, without a dedicated assay, additional reagents, or
proprietary software. Two lines of clinical validation follow. We are running our own
validation study at SA Pathology, in tandem with the accredited STR chimerism assay, to
be reported separately; and because allomix is panel-agnostic, each adopting laboratory
validates it on its own panel and specimen types as a laboratory-developed test before
clinical use. We encourage both, and will support laboratories running their own
validation; where a laboratory can release titrated mixtures we will run allomix on them
and return the results, as we have done here with the public SRP434573 data.
