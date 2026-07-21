## Materials and Methods

The statistical formulae (the likelihood, the error model, the confidence-interval and
presence-test calculations, and the multi-donor search) are given in full in the
Supplementary Methods; this section describes each mechanism in operational terms.

### Overview

allomix is implemented in Python (version 3.10 or later) and operates on standard
Variant Call Format (VCF) files. The workflow has four stages: parse the VCFs and decide
which markers are informative, estimate the donor fraction, optionally correct
per-marker amplification bias, and run a panel of quality and sample-integrity checks.
The tool works with any set of biallelic markers (SNPs or indels) and makes no
assumption about the specific panel.

The magnitude estimate rests on the mixture model. At each marker, recipient and donor
genotypes are already known, so the allele reads expected in a mixed sample are a
weighted blend of the two known genotypes. Near full donor chimerism, the usual clinical
situation, the small residual recipient shifts the read counts a little away from the donor's
expected pattern, and the less recipient that remains the smaller that shift.
allomix inverts that forward model: it finds the single donor fraction whose weighted
blend best matches the observed allele counts across all informative markers at once. Box
1 walks through this by counting reads at two markers (Marker 1 and Marker 2).

### Input requirements and two-phase upstream calling

allomix requires three sets of VCF files: recipient pre-transplant genotyping,
donor pre-transplant genotyping, and post-transplant admixture sample(s). VCF parsing
uses the cyvcf2 library.[@PedersenQuinlan2017cyvcf2] For each biallelic site the tool
reads the genotype (GT), allele depth (AD), total depth (DP), and genotype quality (GQ)
fields.

The two things allomix needs from upstream calling are best produced by different tools,
and keeping them separate matters for sensitivity. Recipient and donor genotypes (GT) come
from joint calling of the reference samples with a germline caller such as
GATK,[@McKenna2010gatk] where high-confidence filtering is an asset. Admixture allele
depths (AD) come from a separate forced pileup: bcftools mpileup[@Danecek2021bcftools] at
the panel sites followed by bcftools call constrained to the recipient/donor alleles, so every
panel site reports a reference and alternative count whether or not the alternative was
seen. The admixture is deliberately not taken from GATK: HaplotypeCaller in GVCF mode is
a local-reassembly caller that at homozygous-reference blocks keeps only reads supporting
the called allele, discarding the minority alternative reads that carry the low-fraction
donor signal before they reach the AD field. This is not specific to GATK. Any caller
that applies low-level artifact filtering, including somatic callers tuned to suppress
low-frequency noise, removes exactly the minority reads we are trying to measure, whereas
a forced pileup keeps them. The full rationale, including an empirical check on the
rhAmpSeq panel where zero alternative reads survived in joint-called admixture AD, is
documented on the allomix GitHub page.

### Which markers are informative

A marker is informative only when recipient and donor have different genotypes: those markers
carry the mixing-fraction signal and drive the donor-fraction estimate. Markers where
recipient and donor match carry no fraction information and are set aside, but they are not
wasted. The consensus-homozygous matches (where recipient and donor are homozygous for the
same allele) feed the contamination and sample-swap checks instead, which are independent
of the fraction estimate and are described under Quality control and sample-integrity
checks.

Among informative markers, the size of the signal depends on the genotype contrast, and
allomix sorts every informative marker into one of six genotype-contrast types, tracked
per donor, following Vynck et al.[@Vynck2023bias] The table below lists the six types,
together with the two non-informative consensus classes, and which check (fraction
estimate, residual-recipient presence, contamination correction, or consensus QC) uses each.

| Type | Recipient | Donor | Contrast | Fraction estimate | Residual-recipient presence (S7) | Contamination correction (S8) | Consensus QC (S8) |
|:--|:-----|:--|:--|:--|:--|:--|:--|
| 0  | 0/0 | 1/1 | full | yes | yes | yes | no |
| 1  | 1/1 | 0/0 | full | yes | yes | yes | no |
| 10 | 0/1 | 0/0 | half | yes | yes | no  | no |
| 11 | 0/1 | 1/1 | half | yes | yes | no  | no |
| 20 | 0/0 | 0/1 | half | yes | no  | no  | no |
| 21 | 1/1 | 0/1 | half | yes | no  | no  | no |
| cons-hom | hom       | hom (same) | none | no | no | no | contamination + swap |
| cons-het | 0/1       | 0/1 | none | no | no | no | shared-het balance |

**Marker classes and their uses.** Among the informative types, 0 and 1 give the maximum
allelic contrast (the minority allele has a single possible source); the heterozygous
types give half the contrast, because the heterozygous party already contributes the
allele on its own and the other party only nudges the count (Box 1, Marker 2). Every
informative type feeds the donor-fraction estimate. The residual-recipient presence test needs
a clean donor-absent allele to count, so it uses only the donor-homozygous types (0, 1,
10, 11); types 20 and 21 leave no donor-absent allele. The optional co-pool contamination
correction acts only on the fully homozygous-contrast types (0, 1), where the recipient
(donor-absent) allele reads are otherwise a clean background. The last two rows are the
non-informative consensus classes: at consensus-homozygous markers the minority allele
can only be a background artifact or foreign DNA, which drives contamination and
sample-swap detection, and at consensus-heterozygous markers the alternative-allele
fraction should sit near 0.5 whatever the mixing fraction, so a systematic skew flags
contamination, allelic imbalance, or a sample mix-up.

By default a site is used only if recipient and donor genotype quality is at least 20 and the
admixture depth is at least 100, and at least three informative markers are required to
report an estimate. Sex and mitochondrial contigs (X, Y, M) are excluded by default,
because in a sex-mismatched donor and recipient pair the expected recipient and donor allele
dosage on the sex chromosomes departs from the autosomal diploid model the estimator
assumes. They can be re-enabled per run once recipient and donor sex are known to match, and
the informative sex-chromosome markers that were dropped are reported.

#### Box 1. Reading the residual recipient off the counts

Take two markers in one post-transplant sample from a patient near full donor chimerism,
the usual clinical situation, where the sample is mostly donor and the question is how
much recipient remains. At each marker, recipient and donor genotypes were fixed before transplant,
so we know which person contributes which allele.

**Marker 1 is fully informative: recipient is A/A, donor is G/G.** The donor contributes only
G and the recipient never contributes G, so every A read can only have come from the recipient. In a
near-full-donor sample the minority of A reads is the residual recipient fraction, and the
majority of G reads is the donor.

| Counts at Marker 1 | A | G | Reading |
|:---|--:|--:|:---|
| Residual-recipient reading (minority allele) | 30 | 970 | 30 A reads can only be recipient = **3% residual recipient** |
| Magnitude reading (majority allele) | 30 | 970 | 970/1000 G reads = 97% donor, i.e. **3% recipient remaining** |

The same marker answers two questions. The minority allele, the one only the recipient could
have produced, tells you *whether any recipient remains*, which is the clinically urgent
question near full donor chimerism. The majority shift gives the magnitude, *how much
donor* is present and so how much recipient. These are the two complementary tests allomix
runs.

**Marker 2 is partially informative: recipient is A/G, donor is G/G.** The recipient carries the
donor-absent allele A on only one of its two alleles, so a residual recipient puts A reads into
the sample at half its true fraction. The recipient signal is real but half the size of a
fully informative marker.

| Counts at Marker 2 | A | G | Reading |
|:---|--:|--:|:---|
| Residual-recipient reading (minority allele) | 15 | 985 | 15 A reads = 1.5% observed, doubled for the half-signal = **3% residual recipient** |

This has three consequences. A pure-donor sample still shows a low background of reads
carrying the recipient-only allele, produced by background artifacts (miscalled bases and, on
co-pooled runs, index hopping and low-level cross-sample contamination), which is why
there is a measured floor rather than a clean zero. One marker on its own is noisy, which is
why allomix pools dozens of markers into a single estimate with a confidence interval.
Markers sometimes disagree by more than simple read-counting would predict, which is
why the interval is widened to reflect that extra scatter.

### How the donor fraction is estimated

allomix builds on the biallelic mixture genotype likelihood of Crysup and
Woerner,[@CrysupWoerner2022] the framework that also underlies the forensic mixture tool
Demixtify.[@Woerner2024demixtify] It estimates the donor fraction by maximum likelihood:
it scans candidate fractions and keeps the one under which the observed allele counts
are most probable, given the known recipient and donor genotypes. Three mechanisms make that
estimate realistic for clinical data; the full formulae are in the Supplementary Methods.

First, a background-artifact term. Even a marker where one allele should be absent shows a
low background of that allele, chiefly from miscalled bases. allomix sets the expected
level of these stray reads so it does not read pure background as real signal. The default
is a fixed per-base rate (1%) spread across the three possible miscalls, giving a
per-direction floor near 0.33%. This default is a conservative fallback, not a number fit
from the run: where the data support it, each marker can instead use its own measured
per-site, per-direction rate, estimated from homozygous calls in a training cohort
(Per-site error calibration, below). The other co-pooled artifacts named above, index
hopping and low-level cross-sample contamination, are not folded into this floor. They are
handled separately by the contamination estimate and its optional per-marker correction
(Sample-integrity checks, below), because they scale with a co-pooled genome's allele dose
rather than sitting at a flat per-base level.

Second, the allele counts are modeled so that marker-to-marker scatter can exceed simple
read-sampling noise. Capture and amplification do not treat every marker identically, so
at the high depths of clinical panels the spread across markers is wider than plain
binomial counting predicts. allomix captures this extra scatter with an overdispersion
parameter fit from the data; on a noisy panel it widens the confidence intervals, and on
a clean panel it has almost no effect. This term, not read depth, is the dominant limit
on sensitivity at clinical coverage (Discussion). The scatter is not uniform across
markers: donor-heterozygous markers carry a background allele fraction near one half,
where capture and amplification scatter is largest, whereas donor-homozygous markers sit
near zero or one, where a low-fraction recipient or donor signal appears. allomix therefore
fits a separate overdispersion parameter for each class rather than one shared value,
which stops the symmetric heterozygous scatter from rectifying into a spurious
low-fraction signal (demonstrated on real data in Results). When a class has too few
informative markers to identify its own parameter (fewer than 30 here), allomix falls
back to a single shared parameter for that sample.

Third, the fraction and the overdispersion are fitted together. allomix tries many
candidate donor fractions spread across the full 0% to 100% range, finds the best
overdispersion for each marker class at each candidate, and keeps the fraction whose
combined fit is best, then fine-tunes around that best-fitting fraction. The reported
confidence interval is a profile-likelihood interval: the range of fractions that the
data cannot confidently rule out. The donor fraction is a proportion bounded at 0% and
100%, and the clinically important samples sit right against the upper bound (near full
donor, with a small residual recipient), so the interval is constructed to respect those hard
boundaries rather than running past them, which a simpler symmetric interval would do
(Supplementary Methods).[@Wilks1938]

allomix runs this likelihood in the reverse direction from its original use: Crysup and
Woerner estimate unknown contributor genotypes at a known mixture fraction, whereas
allomix estimates the mixture fraction from contributor genotypes already known from the
pre-transplant samples, the simplification the clinical setting allows.

### The residual-recipient presence test

The magnitude estimate above answers how much recipient remains (equivalently, how much
donor). Near full donor chimerism, where the recipient fraction is near zero, its confidence
interval widens and it cannot cleanly separate a small residual recipient from zero, yet that
is exactly where the clinical concern (early relapse) sits. A separate presence
test answers the narrower question "is any recipient left?", and is built to stay sensitive at
and near zero, where the magnitude estimate loses precision. The two tests emphasise
different markers: the magnitude estimate uses all informative markers, while the presence
test reads only the subset that gives a one-sided recipient signal, so the presence-test
markers are a subset of the magnitude estimate's markers rather than a disjoint set.

The presence test looks only at markers where the donor is homozygous and the recipient
carries an allele the donor does not (Box 1, Marker 1, residual-recipient reading). At such a
marker, in a pure-donor sample, that recipient-only allele should appear at no more than the
background-artifact floor. Any consistent excess across these markers is residual recipient
DNA. allomix combines the recipient-only allele counts across all such markers into a
one-sided test against the per-marker error background and reports a p-value, a
recipient-fraction estimate, and a confidence interval. Because the honest answer when no
recipient is present sits at the zero boundary, the test uses the same boundary-aware
statistics as the confidence interval above (Supplementary Methods). The result is a
presence/absence call designed to detect residual recipient below the level the magnitude
estimate reliably quantifies.

### Per-marker bias correction

Capture and amplicon panels have small systematic per-marker biases that pull observed
allele frequencies away from their true values.[@Vynck2023bias] allomix can optionally
correct these using a per-marker bias table measured at heterozygous sites (true
frequency 0.5), applied on the log-odds scale so it does not overcorrect at the extreme
frequencies that dominate low-fraction samples (formula and table construction in
Supplementary Methods S6). As shown in Results, it narrows the confidence intervals but
barely moves the point estimate, by design.

### Per-site error calibration

The default error term is a single symmetric per-base rate, but the background that
matters for low-fraction detection is direction-specific: REF-to-ALT and ALT-to-REF
substitution rates differ (oxidation, strand bias, flanking context), and the dominant
direction sets the floor a residual-recipient signal must clear. allomix can therefore read an
optional per-site, per-direction empirical error table (the `estimate-errors` command,
which pools reads at homozygous calls: the minority-allele rate at homozygous-reference
calls estimates the REF-to-ALT rate and at homozygous-alternative calls the ALT-to-REF
rate), and each marker then uses its own measured rate in the matching direction in both
the magnitude likelihood and the presence-test background. The in silico results below
use the global rate, because the simulator draws errors from that same uniform model, so
a per-site table cannot change simulated data; it sharpens detection only on real reads,
where error rates genuinely vary by site and direction. Where a set of reference
individuals shares a panel, run and chemistry, the substitution background is shared and
one table can be pooled across them, which fills both substitution directions at more
sites than a per-patient table can. Supplementary Figure S15 quantifies what the choice
of error model buys on real reads, comparing the flat default, per-patient tables and a
pooled table. Pooling applies to the substitution background only; run-specific index
hopping and cross-sample contamination are handled by the separate contamination
machinery described below.

### Multi-donor estimation

The single-donor model extends directly to two donors plus recipient: instead of one
fraction, allomix searches the pair of donor fractions (constrained to sum to no more
than one, with the recipient as the remainder) over a triangular grid, then refines locally
and reports a profile-likelihood interval for each donor (Supplementary Methods). A
marker is informative if the recipient differs from any donor, and informative counts are
tracked per donor, since related donors can leave few markers that separate one donor
from the other. The recipient-plus-two-donor scope is a practical default set by the common
clinical case and the growing cost of searching a higher-dimensional fraction simplex,
not a hard limit of the likelihood, which generalizes to more components; settings such
as sequential transplants or some cord-blood mixtures that exceed two donors are untested
here.

### Outlier-resistant refit for off-fit markers

A few markers can sit far off the fit, most often because the recipient's malignant clone
carries a copy-number variant (CNV) there, which is common in the hematological cancers
for which most transplants are performed (see Recipient copy-number variants below). allomix
optionally re-estimates the fraction after
down-weighting such markers, using a median-based outlier-resistant scale. The trim is
deliberately one-sided: it removes markers whose deviation points *away* from recipient
presence, and it protects markers whose deviation points *toward* recipient presence. This
matters at low recipient fraction, where the handful of markers carrying the real low-level
recipient signal sit off a donor-dominated fit and a symmetric outlier rule would discard
exactly those markers, collapsing a real low-fraction recipient signal to zero. The
trade is explicit: at the limit of detection (LoD) we would rather keep a few artifacts than
throw away a real low-fraction recipient signal. The refit engages only when the number of
off-fit markers exceeds chance expectation, so clean samples are unchanged, and it is
floored so it does not over-trim small panels. When too many markers are excluded, the
sample is flagged for review rather than reported as a confident estimate.

### Quality control and sample-integrity checks

allomix reports a three-level verdict per sample: PASS, REVIEW (an estimate is produced
but a reliability flag should be checked), or FAIL (the result is not usable). Basic
checks cover marker sufficiency, sequencing depth, coverage uniformity (whether enough
markers reach a set share of the sample's mean depth, catching an uneven capture that an
acceptable mean would otherwise hide), and confidence-interval width. A goodness-of-fit
check compares the marker-to-marker scatter against what the model expects, computed both
before and after the outlier-resistant refit so the refit cannot hide a bad fit by
discarding its own outliers. Because a chi-squared test is overpowered at panel depth (it
reaches significance for a misfit too small to change the result), it is promoted to
REVIEW on effect size, not bare significance (Supplementary Methods S10). The presence
test is cross-checked against the magnitude estimate, raising a soft warning when
residual recipient is detected at a level below what the magnitude estimate can quantify. The verdict and
the presence-test p-values are analytical outputs on the data quality and the fit; they
are not a validated clinical call, and turning them into clinical action requires
PASS/REVIEW/FAIL and presence-test thresholds each laboratory sets and validates on its
own panel.

The sample-integrity checks are separated by genotype geometry rather than by
re-thresholding one number, each reading a marker class the magnitude estimate never uses
(full tests and estimators in Supplementary Methods S8, S10):

- **Residual recipient** at donor-homozygous markers where the recipient carries the donor-absent
  allele (the presence test above).
- **Contamination by a non-recipient, non-donor genome** at consensus-homozygous markers,
  where recipient and every donor are homozygous for the same allele so the minor allele can
  only be a background artifact or foreign DNA (for example another patient co-pooled on
  the same flowcell via index hopping). allomix estimates it as the excess minor-allele
  signal over a data-internal error floor, and separates real foreign DNA from a uniform
  error elevation by its dose-response (foreign reads scale with co-pooled allele dose).
  When that dose-response is significant, an optional per-marker correction (off by
  default) subtracts the dose-predicted contamination before the magnitude fit; because
  it only ever removes recipient-allele reads it is an asymmetric downward adjustment, so it is
  gated rather than applied unconditionally. An optional VCF header field recording the
  flowcell and lane adds a pure-metadata index-hopping flag, kept separate from the
  in-data estimate.[@Costello2018indexswap]
- **A sample swap or contaminating genome** at those same consensus-homozygous markers
  when the minority allele is individually significant rather than a faint background,
  catching a wrong-patient VCF the magnitude estimate cannot see. As with goodness-of-fit,
  promotion to REVIEW is gated on effect size, since a genuine swap mismatches roughly
  half the consensus sites.
- **Allele imbalance at shared-heterozygous markers**, where recipient and every donor are
  heterozygous so the alternative-allele fraction should sit near 0.5 whatever the mixing
  fraction; a systematic skew flags contamination, copy-number or allelic imbalance, or a
  sample mix-up.

Two further checks guard sample identity and marker quality. A relatedness check
estimates the kinship between each pair of input samples directly from the genotypes,
using an allele-frequency-free coefficient in the style of
somalier,[@Pedersen2020somalier] which suits a panel-agnostic tool with unknown panel
allele frequencies. Its logic is asymmetric: losing an expected close relationship is
treated as a swap signature and fails, a declared identical pair always fails (a
syngeneic donor is unmeasurable), and a milder relatedness shift only triggers review.
Separately, a read-level artifact filter screens donor-homozygous markers for soft-clip,
read-position, and strand-bias artifacts using bcftools mpileup annotations, judging
strand bias by effect size rather than significance (at high depth a real allele's mild
strand skew is highly significant but harmless, whereas an artifact is extreme regardless
of depth). It runs as a screen on the pileup output rather than delegated to a caller's
built-in filtering, for the same reason the admixture is piled up rather than called: a
caller strips the low-fraction signal along with the artifacts before we ever see the
reads. The filter auto-disables on single-strand amplicon panels. Finally, a
genotype/allele-depth consistency guard drops reference-sample markers whose called
genotype contradicts their own allele depths, so a marginally rescued heterozygous call
does not feed systematic bias into the estimate; it is applied only to the reference
samples, never to the admixture.

### Simulation framework

For validation, allomix includes a simulator that blends two genotype VCFs at a specified
donor fraction and layers four calibrated sources of measurement noise: per-marker
amplification bias, non-uniform depth, sequencing error, and locus dropout (provenance
and values in Supplementary Table S1, distributions against empirical data in Figures
S1-S3). One caveat applies to every in silico result: the simulator draws reads from the
same model the estimator fits, so accuracy and LoD are best-case figures showing the
estimator recovers truth when the data match its model, not that the model matches real
sequencing. The main runs draw from a binomial, so overdispersion (the dominant real-data
noise at clinical depth; Discussion) and co-pooled contamination are absent by
construction; the real-read SRP434573 analysis (Results) is the independent check. An
overdispersed beta-binomial draw is available for the overdispersion characterisation
(Figures S7, S8) only. For longitudinal monitoring we simulated a six-timepoint post-HSCT
trajectory (day +14 to +365, donor fraction 15% to 97% with a 3-percentage-point dip at
day +180) at 500x with five replicates. Full generative details, including the
copy-number and sibling-segregation extensions, are in Supplementary Methods.

### Limit of detection

Limit of detection and limit of blank (LoB) follow the CLSI guideline
EP17-A2,[@CLSIEP17A2] as in published evaluations of comparable NGS chimerism
assays.[@Vynck2021devyser] LoB is the 95th percentile of the estimated donor fraction on
a pure-recipient blank, and LoD is the lowest true fraction at which at least 95% of
replicates exceed LoB, read from a logistic fit to the dilution series. We computed LoB
and LoD separately for each donor/recipient pair and report the median across pairs with a
10th-to-90th-percentile band. The sweep covered two relatedness levels (unrelated, full
sibling), five depths (100x to 2,000x), six panel sizes (25 to 400 markers, nested so a
smaller panel is a strict prefix of a larger one), and seven true donor fractions (0 to
5%). A nested design (multiple donor/recipient pairs, each with multiple sequencing
replicates that vary only read-sampling noise) keeps identity-by-descent variation
between pairs in the reported band rather than leaking it into the central estimate;
details and the full grid are in the Supplementary Methods. Because the simulator draws
reads from a binomial, this LoD is a near-binomial analytical best case rather than a
validated assay limit; a real assay's LoD can only be higher (Results, Discussion).

### Recipient copy-number variants

HSCT recipients are usually hematological malignancy patients whose residual or
relapsing neoplastic clone often carries CNVs, which break the assumption that every
locus is diploid. To test sensitivity to this, the simulator applies a per-marker
recipient CNV before read sampling, modelling the recipient as a mixture of
normal diploid cells and a CNV-bearing clone, with the expected allele fraction a
copy-number-weighted average rather than the diploid mean. Three CNV types are
modelled (copy-neutral loss of heterozygosity (cnLoH), single-copy deletion, and single-copy
gain) at a controllable burden, applied only to the admixture sample to match the
two-phase workflow in which recipient genotypes come from a clean reference. We computed
the EP17-A2 LoD in two directions (recipient clone as the minor component, the relapse
direction; and donor as the minor component against a CNV-bearing recipient
background) at {{ cnv_loh_headline.ref_markers }} markers and
{{ cnv_loh_headline.ref_depth }}x depth, with and without the one-sided
outlier-resistant refit, sweeping burden over 0, 10%, 25%, and 50%.

### Real-data validation dataset

To test allomix on real reads, we used a public dataset of titrated DNA mixtures (SRA
study SRP434573, BioProject PRJNA960854), companion data to a doctoral thesis (Chu,
Huazhong University of Science and Technology, 2024) with no associated journal article.
These are artificial mixtures of DNA from unrelated individuals, not post-transplant
patient samples, so they exercise the estimator on real reads at known mixing fractions
but do not carry the biology of clinical chimerism. The dataset is
{{ srp434573.n_individuals }} unrelated individuals captured on a
{{ srp434573.panel_n_snps | commas }}-SNP molecular inversion probe (MIP)
sample-identification panel on an {{ srp434573.platform }}, with reads merged to
single-end inserts (each marker read from a single strand). At a median on-target depth
of {{ srp434573.depth_median | commas }}x this is the high-depth biallelic
sample-identification SNP marker class allomix is built for, so it exercises the estimator
on its intended input despite not being transplant material; the reconstructed panel
intervals and inferred recipient/donor labels (below) are the caveats on it, not the marker
chemistry. Pairwise two-person mixtures and one three-person mixture were prepared at
known major:minor ratios, giving minor-contributor fractions of
{{ srp434573.dilution_ladder }}.

We assigned the minor (titrated) contributor to the recipient role and the major
contributor to the donor. This mirrors the usual clinical situation, where a patient is
near full donor chimerism and the recipient is the small residual fraction the assay is
trying to measure, so the monitored quantity is the recipient fraction and each dilution
series reads as a declining-chimerism (relapse) trajectory. The assignment is a labelling
convention based on mixing fraction, not a biological recipient/donor identity the thesis
provides, and it affects only the recipient/donor labelling, not the genotyping or the
recovered fractions.
Reads were aligned to hg38 with bwa-mem;[@LiDurbin2009bwa; @Li2013bwamem] because every
MIP amplicon shares start and end coordinates, duplicate marking was not applied (it
would discard almost all coverage). Because the thesis publishes no panel coordinates,
panel intervals were reconstructed from the aligned reads as high-depth amplicon
clusters ({{ srp434573.n_intervals | commas }} intervals:
{{ srp434573.n_intervals_autosomal | commas }} autosomal and
{{ srp434573.n_intervals_chrx }} on chrX), and genotypes and admixture allele depths
were generated with the two-phase pipeline above. allomix was run with default
parameters; on this single-strand panel the strand-bias
artifact filter auto-disables. The dataset covers the unrelated-donor case only (no
related or sibling donors) and carries a co-pooled contamination floor characterised in
Results.

### Real-data limit of detection by subsampling

To obtain a LoD on real reads rather than simulation, we reduced the depth and panel size
of the high-depth SRP434573 mixtures until the LoD moved into the dilution series'
measurable range. From each two-person mixture we sub-sampled reads to a target mean depth (100x
to 2,000x) by a single per-sample binomial keep rate applied uniformly across markers,
which is the statistical analogue of FASTQ read subsampling and preserves the real
locus-to-locus depth variation, then re-applied the depth filter so low-depth loci drop
out. Panel size was set by nested random subsets of the informative markers (prefixes of
one permutation per mixture and replicate, so each curve is monotonic in panel size), and
multiple sub-sampling replicates were drawn per depth-and-panel cell. For each cell we
read the LoD off the same logistic fit as the simulated sweep, separately for the
magnitude estimate (a sample is detected when its 95% confidence interval for the recipient
fraction excludes zero) and the residual-recipient presence test (detected when the
presence-test p-value is below 0.05). Both use blank-free per-sample detection rules, so
neither needs an EP17 blank sample, which this dataset does not provide as used here.
Because these
are pseudo-replicates sub-sampled from one library rather than independent low-depth
libraries, the result shows whether the real-data LoD tracks the simulated curves within
the dilution grid rather than serving as an independent wet-lab limit.

### Software availability

allomix is implemented in Python with dependencies on cyvcf2, NumPy,[@Harris2020numpy]
SciPy,[@Virtanen2020scipy] and Jinja2. It is available under the MIT license at
https://github.com/SACGF/allomix, installable via pip (`pip install allomix`). The
command-line interface provides `detect` for single-sample or multi-timepoint analysis,
`timeline` for consolidated multi-timepoint reporting, `estimate-bias` for panel bias
calibration, and `estimate-errors` for per-site error calibration. The repository also
includes a reference Snakemake[@Molder2021snakemake] workflow that runs the two-phase
upstream calling (GATK joint genotyping of recipient and donor, then forced bcftools pileup
of the admixture samples at the panel sites) and the optional bwa-mem alignment step
ahead of it. Each analysis can emit a tab-separated summary, a structured JSON record,
and a single-file HTML or PDF report that presents the estimate, its confidence interval,
the quality-control verdict, and the method explanations for clinical review. The
HTML and PDF reports are rendered from templates a laboratory can override to match its
own reporting style (logo, layout, and wording); the PDF is for attaching to a record or
laboratory information system. The repository documentation includes a
guide for qualifying an existing clinical panel for chimerism use, covering
informative-marker sufficiency, per-marker characterization, inclusion thresholds, and
building the bias and error correction tables.
