## Discussion

allomix lets a laboratory add post-HSCT chimerism monitoring to an NGS panel it already
runs, by repurposing polymorphic markers that are sequenced anyway, without a dedicated
assay or proprietary software. It does this with two complementary readouts (how much
donor is present, and whether any host remains) and a set of sample-integrity checks
built from the same marker data.

### Accuracy and sensitivity

The <1% mean absolute error across all depths is competitive with published performance
for commercial tools: Kakodkar et al. reported 0.3--1.5% MAE for AlloSeq
HCT,[@Kakodkar2023alloseq] Pedini et al. and Vynck et al. showed comparable precision
for the Devyser system,[@Pedini2021devyser; @Vynck2021devyser] and Blouin et al.
reported R^2 = 0.9987 across the full chimerism range for ScisGo.[@Blouin2024comparison]
Qama et al. validated the Devyser assay at a limit of detection of 0.06% with high STR
concordance (R^2 = 0.998), and found that NGS detected residual host DNA (>0.1%) in 85%
of samples called full donor chimerism (>95%) by STR, against 5% by
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
Overdispersion is also the dominant control on the limit of detection at high depth: as
depth grows, the per-marker variance stops falling and the effective depth is capped, so
the LoD saturates at a floor rather than continuing to improve as 1/√n (Supplementary
Figures S7, S8). At a fitted overdispersion consistent with our panel the in silico LoD
is roughly {{ overdispersion_lod_headline.fold_rho100_vs_binomial }}-fold higher than
the pure-binomial value, so beyond that regime adding depth buys little and the panel's
overdispersion, not its coverage, sets the achievable LoD. Marker count is a separate
hard limit: depth and overdispersion act on the markers a panel already has, and no
amount of depth compensates for too few informative markers, so panel size floors the
achievable LoD independently of coverage (Results, relatedness). The noise-ablation
analysis makes the same point on point-estimate accuracy: adding overdispersion to the
full noise model raises RMSE more than amplification bias, depth non-uniformity, or
sequencing error individually ({{ supp_synthetic.ablation_rmse_overdispersion_pct }}%
versus {{ supp_synthetic.ablation_rmse_full_pct }}% for the otherwise-identical binomial
model; Supplementary Figure S4). This is the honest explanation for the gap between an
idealised in silico LoD and real-panel performance, and it is the same overdispersion
that drives the goodness-of-fit REVIEW flags on the real SRP434573 mixtures.

Across the simulated depths the 95% intervals cover close to nominal
({{ depth_1000.ci_coverage_pct }}--{{ depth_200.ci_coverage_pct }}%) once overdispersion
is fit separately for the donor-heterozygous and donor-homozygous marker classes rather
than as a single shared parameter. A residual calibration gap remains on the real
titrated ladder, where the intervals mildly undercover the nominal fraction. The driver
is the co-pooled donor-homozygous contamination background, a dose-dependent excess of
host-allele reads that the per-marker correction reduces but does not fully remove,
together with real heterozygous-marker scatter that is not an exact beta-binomial,
rather than a shortfall of the model in controlled data. The realised mixing itself is
faithful: the donor-homozygous markers alone return the nominal fraction. Empirical
recalibration from training data is one route to closing that gap.

### Outlier-resistant estimation at the limit of detection

The relapsing recipient clone often carries copy-number changes, so a few markers sit
far off the diploid fit and an outlier-resistant refit is needed to recover the estimate
(Results, copy-number stress test). The trim is deliberately one-sided: it protects
markers whose deviation points toward host presence and removes only those pointing
away. The reason is specific to the low-fraction regime. When the true host fraction is
around 1%, the handful of markers actually carrying that signal sit off a
donor-dominated fit and read as outliers, so a symmetric outlier rule would discard
exactly the markers of interest and collapse the estimate to zero. At the limit of
detection, keeping a few artifacts is the safer trade than throwing away a real
low-fraction host signal, and this asymmetry is what recovers the donor LoD against an
aberration-bearing recipient background while leaving aberration-free samples unchanged.

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
a laboratory controls. Residual correlation does not fail silently: it inflates
per-marker variance in the same way overdispersion does, so it is partly absorbed by the
per-class overdispersion model and surfaces in the goodness-of-fit check that drives the
REVIEW flag, the same mechanism that flagged the overdispersed SRP434573 mixtures. A
marker set dense enough for LD to matter would therefore tend to be caught as an
overdispersed fit rather than reported as a confidently wrong interval. LD can also be
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
asking only whether that host-only allele exceeds the sequencing-error background. We
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

A quantitative chimerism number is dangerous if it is computed on the wrong sample, so
allomix turns the same marker data into a set of integrity checks: a relatedness check
that flags sample swaps from unexpected kinship, a consensus-site consistency test that
catches a wrong-patient VCF the fraction estimate cannot see, and an in-data
contamination estimate that separates a real low-fraction signal from co-pooled or
foreign DNA. These three low-fraction signals (residual host, contamination, gross swap)
are kept distinct by reading disjoint marker sets defined by genotype geometry rather
than by re-thresholding one statistic, which is what makes them independent rather than
competing interpretations of the same number. For a clinical laboratory, aggressive
detection of sample mix-ups is a desirable property, and building these checks from data
already in hand adds no extra cost at the bench.

### Repurposing existing panels

The central advantage of allomix is that it works with markers laboratories already
sequence. Sample-identification SNPs, pharmacogenomic markers, and other polymorphic
loci included for quality control or diagnosis can serve double duty for chimerism,
eliminating a separate dedicated assay. Lee et al. demonstrated the principle with 121
SNPs in a myeloid panel but required custom scripting with no reusable
tool;[@Lee2019snp] allomix generalises this into a deployable tool. Vynck et al. showed
that as few as three informative markers suffice for quantification, with accuracy
improving as markers are added, and that panels of about 20 markers with MAFs near 0.5
give a >95% chance of at least three informative markers even for sibling pairs; their
FABCASE tool can assess panel sufficiency prospectively for a specific donor-host
pair.[@Vynck2022markers; @Vynck2025fabcase] Sample-ID marker sets with tens of
polymorphic markers are therefore expected to be adequate for most clinical scenarios.

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

By accepting standard VCFs, allomix decouples chimerism analysis from upstream alignment
and variant calling. The one requirement is the two-phase calling in Methods: genotypes
from joint germline calling (for example with GATK), and admixture allele depths from a
forced bcftools pileup at the panel sites rather than from that germline caller, so the
minority alternative reads carrying the low-fraction signal are retained rather than
stripped. Both phases use standard tools, but the admixture pileup must be configured as
a separate step rather than reusing one joint-calling pass for every sample.

Chimerism is a property of a cell population, not of whole blood per se. A donor
fraction measured in unfractionated blood is a lineage-abundance-weighted average, so a
meaningful change confined to a small lineage may show only a fractional whole-blood
change. Lineage-specific testing on sorted subsets is recommended in some clinical
contexts, though no quantitative intervention threshold has been standardised for
sub-STR detection.[@Clark2025bjh; @KharfanDabaja2021astct; @Kakodkar2023alloseq] allomix
is lineage-agnostic by design: it estimates the donor fraction in whatever DNA the input
VCF represents, and the same validation applies to sorted inputs. Where whole-blood
signal is adequate, sensitivity better than STR can reduce the need for cell sorting and
its associated cost and turnaround; where a clinically important change is confined to a
small lineage, sorting is still required, and the same analysis applies to the sorted
input. Specimen choice is a clinical and workflow decision outside the tool; the tool's
analytical sensitivity is necessary but not sufficient for adequate clinical sensitivity
in any given monitoring scenario.

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
controlled cell-line dilution series, is a separate study. Because allomix is
panel-agnostic, that validation does not transfer between panels: a laboratory adopting
the tool validates it on its own panel and specimen types, as for any
laboratory-developed test, rather than relying on a single validation done once for the
tool. A validation on one panel would not establish performance on the others the tool
is meant to serve. Blouin et al. describe a practical framework for such validation,
including run-level metrics and sample-level acceptance criteria, that is a useful model
for future allomix studies.[@Blouin2024comparison]

The limit of detection reported here is a best-case analytical figure, not a validated
assay limit (a real assay's LoD can only be higher). It is the CLSI EP17-A2
95%-detection criterion[@CLSIEP17A2] applied to simulated data under a noise model
calibrated from empirical panel data, with reads drawn from a binomial. The in silico
LoD at 1,000x and 100 markers (magnitude estimate
{{ lod_headline.unrelated_lod_1000x_100markers_pct }}% unrelated,
{{ lod_headline.sibling_lod_1000x_100markers_pct }}% sibling; residual-host detection
{{ presence_lod_curve_headline.presence_unrelated_lod_1000x_100markers_pct }}%
unrelated, {{ presence_lod_curve_headline.presence_sibling_lod_1000x_100markers_pct }}%
sibling) sits within the range reported for commercial NGS chimerism assays (0.06% for
Devyser,[@Qama2026devyser] 0.3% for AlloSeq HCT,[@Kakodkar2023alloseq] 0.2--0.5% for
ScisGo[@Blouin2024comparison]), but those vendor figures come from dilution series on
real DNA, so the comparison is between a best-case analytical figure and reported wetlab
performance, not a head-to-head benchmark on matched samples. The gap is expected to be
driven largely by overdispersion (above), and wetlab validation will set the floor
allomix actually delivers in routine use. As a partial check against real reads, we
sub-sampled the high-depth SRP434573 mixtures across depth and panel size and measured
the LoD on those reads directly (Figure 5): it falls with depth and panel size as the
simulated curves do and tracks them within the dilution grid, but because these are
pseudo-replicates from one library rather than independent dilutions it confirms the
simulated trend rather than replacing the wetlab limit. Finally, the in silico work
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
proprietary software. Clinical validation on patient cohorts is the planned next study,
carried out per panel by each adopting laboratory before clinical use.
