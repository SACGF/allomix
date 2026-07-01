## Materials and Methods

The statistical formulas (the likelihood, the error model, the confidence-interval and
presence-test calculations, and the multi-donor search) are given in full in the
Supplementary Methods; this section describes each mechanism in operational terms.

### Overview

allomix is implemented in Python (version 3.10 or later) and operates on standard
Variant Call Format (VCF) files. The workflow has four stages: parse the VCFs and decide
which markers are informative, estimate the donor fraction, optionally correct
per-marker amplification bias, and run a panel of quality and sample-integrity checks.
The tool works with any set of biallelic markers (SNPs or indels) and makes no
assumption about the specific panel.

The method rests on a single idea, the mixture model. At each marker, host and donor
genotypes are already known, so the allele reads expected in a mixed sample are a
weighted blend of the two known genotypes: the more donor in the sample, the more the
read counts shift from the host's expected pattern toward the donor's. allomix inverts
that forward model: it finds the single donor fraction whose weighted blend best matches
the observed allele counts across all informative markers at once. Box 1 walks through
this by counting reads at two markers (Marker 1 and Marker 2).

### Input requirements and two-phase upstream calling

allomix requires three sets of VCF files: host (recipient) pre-transplant genotyping,
donor pre-transplant genotyping, and post-transplant admixture sample(s). VCF parsing
uses the cyvcf2 library.[@PedersenQuinlan2017cyvcf2] For each biallelic site the tool
reads the genotype (GT), allele depth (AD), total depth (DP), and genotype quality (GQ)
fields.

The two things allomix needs from upstream calling are best produced by different tools,
and keeping them separate matters for sensitivity. Host and donor genotypes (GT) come
from joint calling of the reference samples with a germline caller such as
GATK,[@McKenna2010gatk] which gives high-confidence germline calls. Admixture allele
depths (AD) come from a separate forced pileup: bcftools mpileup[@Danecek2021bcftools]
at the panel sites followed by bcftools call constrained to the host/donor reference and
alternative alleles, so every panel site reports a reference and alternative count
whether or not the alternative was seen. Admixture samples are deliberately not taken
from GATK. HaplotypeCaller in GVCF mode is a local-reassembly caller that, at
homozygous-reference blocks, keeps only reads supporting the called allele, so the
minority alternative reads that carry the low-fraction donor signal are discarded before
joint calling and never reach the AD field. Forcing a pileup at the known panel sites
keeps them. This is not specific to GATK: any caller that applies low-level artifact
filtering, including somatic callers tuned to suppress low-frequency noise, is built to
remove exactly the minority reads that carry the low-fraction signal we are trying to
measure. A forced pileup applies no such filtering. The full rationale, including an
empirical check on the rhAmpSeq panel where zero alternative reads survived in
joint-called admixture AD, is documented on the allomix GitHub page.

### Which markers are informative

A marker is informative only when host and donor have different genotypes: those markers
carry the mixing-fraction signal and drive the donor-fraction estimate. Markers where
host and donor match carry no fraction information and are set aside, but they are not
wasted. The consensus-homozygous matches (where host and donor are homozygous for the
same allele) feed the contamination and sample-swap checks instead, which are independent
of the fraction estimate and are described under Quality control and sample-integrity
checks. Among informative markers, the size of the signal
depends on the genotype contrast. When one party is homozygous for an allele the other
never carries (for example host A/A, donor G/G), the marker is fully informative: the
minority allele can come from only one of the two people. When one party is
heterozygous, the marker is partially informative: that party already contributes the
allele on its own, so the other party only nudges the count and the per-marker signal is
roughly half the size (Box 1, Marker 2). allomix classifies every informative marker
into one of six genotype-contrast types and tracks them per donor; the full taxonomy,
following Vynck et al.,[@Vynck2023bias] is in the Supplementary Methods. By default a
site is used only if host and donor genotype quality is at least 20 and the admixture
depth is at least 100, and at least three informative markers are required to report an
estimate.

#### Box 1. Reading the residual host off the counts

Take two markers in one post-transplant sample from a patient near full donor chimerism,
the usual clinical situation, where the sample is mostly donor and the question is how
much host remains. At each marker, host and donor genotypes were fixed before transplant,
so we know which person can make which allele.

**Marker 1 is fully informative: host is A/A, donor is G/G.** The donor makes nothing but
G and the host never makes G, so every A read can only have come from the host. In a
near-full-donor sample the minority of A reads is the residual host fraction, and the
majority of G reads is the donor.

| Counts at Marker 1 | A | G | Reading |
|:---|--:|--:|:---|
| Residual-host reading (minority allele) | 30 | 970 | 30 A reads can only be host = **3% residual host** |
| Magnitude reading (majority allele) | 30 | 970 | 970/1000 G reads = **97% donor** |

The same marker answers two questions. The minority allele, the one only the host could
have produced, tells you *whether any host remains*, which is the clinically urgent
question near full donor chimerism. The majority shift tells you *how much donor* is
present. These are the two complementary tests allomix runs.

**Marker 2 is partially informative: host is A/G, donor is G/G.** The host carries the
donor-absent allele A on only one of its two alleles, so a residual host puts A reads into
the sample at half its true fraction. The host signal is real but half the size of a
fully informative marker.

| Counts at Marker 2 | A | G | Reading |
|:---|--:|--:|:---|
| Residual-host reading (minority allele) | 15 | 985 | 15 A reads = 1.5% observed, doubled for the half-signal = **3% residual host** |

This has three consequences. A pure-donor sample still shows a low background of reads
carrying the host-only allele, produced by background artifacts (miscalled bases and, on
co-pooled runs, index hopping and low-level cross-sample contamination), which is why
there is a measured floor rather than a clean zero. One marker on its own is noisy, which is
why allomix pools dozens of markers into a single estimate with a confidence interval.
And markers sometimes disagree by more than simple read-counting would predict, which is
why the interval is widened to reflect that extra scatter.

### How the donor fraction is estimated

allomix builds on the biallelic mixture genotype likelihood of Crysup and
Woerner,[@CrysupWoerner2022] the framework that also underlies the forensic mixture tool
Demixtify.[@Woerner2024demixtify] It estimates the donor fraction by maximum likelihood:
it scans candidate fractions and keeps the one under which the observed allele counts
are most probable, given the known host and donor genotypes. Three mechanisms make that
estimate realistic for clinical data; the full formulas are in the Supplementary Methods.

First, a background-artifact term. Even a marker where one allele should be absent shows a
low background of that allele, from miscalled bases and, on co-pooled runs, index hopping
and low-level cross-sample contamination. allomix models this floor with a fixed per-base
rate (default 1%) spread across the possible miscalls, which sets the expected level of
stray reads and keeps the estimate from reading pure background as real signal.

Second, the allele counts are modeled so that marker-to-marker scatter can exceed simple
read-sampling noise. Capture and amplification do not treat every marker identically, so
at the high depths of clinical panels the spread across markers is wider than plain
binomial counting predicts. allomix captures this extra scatter with an overdispersion
parameter fit from the data; on a noisy panel it widens the confidence intervals, and on
a clean panel it has almost no effect. This term, not read depth, is the dominant limit
on sensitivity at clinical coverage (Results, Discussion). The scatter is not uniform
across markers: at markers where the donor is heterozygous the background allele
fraction sits near one half, where capture and amplification scatter is largest and
symmetric, whereas at markers where the donor is homozygous the background sits near
zero or one, which is where a low-fraction host or donor signal appears. Because the two
classes carry different amounts of scatter, allomix fits a separate overdispersion
parameter for each (donor-heterozygous and donor-homozygous) rather than one shared
value. At a low mixture the near-0.5 heterozygous class is close to pure noise, so its
fitted overdispersion is small and its effective weight drops toward zero, which keeps its
symmetric scatter from rectifying into a spurious low-fraction signal; the
donor-homozygous class, where the real low-fraction signal sits, keeps its weight. When a
class has too few informative markers to identify its own parameter (fewer than 30 here),
allomix falls back to a single shared parameter for that sample.

Third, the fraction and the overdispersion are fit together. allomix tries many
candidate donor fractions spread across the full 0% to 100% range, finds the best
overdispersion for each marker class at each candidate, and keeps the fraction whose
combined fit is best, then fine-tunes around that best-fitting fraction. The reported
confidence interval is a profile-likelihood interval: the range of fractions that the
data cannot confidently rule out. The donor fraction is a proportion bounded at 0% and
100%, and the clinically important samples sit right against the upper bound (near full
donor, with a small residual host), so the interval is constructed to respect those hard
boundaries rather than running past them, which a simpler symmetric interval would do
(Supplementary Methods).[@Wilks1938]

allomix uses the same likelihood for the reverse problem: Crysup and Woerner estimate
unknown contributor genotypes at a known mixture fraction, whereas allomix estimates the
mixture fraction from contributor genotypes already known from the pre-transplant samples,
the simplification the clinical setting allows.

### The residual-host presence test

The magnitude estimate above answers "how much donor?" Near full donor chimerism its
confidence interval widens and it cannot cleanly separate a small residual host from zero,
yet that is exactly where the clinical concern (early relapse) sits. A separate presence
test answers the narrower question "is any host left?", and is built to stay sensitive at
and near zero, where the magnitude estimate loses precision. The two tests emphasise
different markers: the magnitude estimate uses all informative markers, while the presence
test reads only the subset that gives a one-sided host signal, so the presence-test
markers are a subset of the magnitude estimate's rather than a disjoint set.

The presence test looks only at markers where the donor is homozygous and the host
carries an allele the donor does not (Box 1, Marker 1, residual-host reading). At such a
marker, in a pure-donor sample, that host-only allele should appear at no more than the
background-artifact floor. Any consistent excess across these markers is residual host
DNA. allomix combines the host-only allele counts across all such markers into a
one-sided test against the per-marker error background and reports a p-value, a
host-fraction estimate, and a confidence interval. Because the honest answer when no
host is present sits at the zero boundary, the test uses the same boundary-aware
statistics as the confidence interval above (Supplementary Methods). The result is a
presence/absence call designed to detect residual host below the level the magnitude
estimate reliably quantifies.

### Per-marker bias correction

Capture and amplicon panels have small systematic per-marker biases that pull observed
allele frequencies away from their true values.[@Vynck2023bias] allomix can correct
these using a per-marker bias table measured at heterozygous sites, where the true
frequency is known to be 0.5. The correction is applied proportionally (on the log-odds
scale) rather than as a flat additive shift, because an additive shift overcorrects at the
extreme expected frequencies that dominate low-fraction samples; the formula is in the
Supplementary Methods. The table can be built from reference samples called the same way
as the admixture, or from a patient cohort at markers where host and every donor are
heterozygous. Correction is optional and is skipped when no table is supplied. As shown
in Results, it sharpens the confidence intervals but barely moves the point estimate, by
design.

### Per-site error calibration

The default error term is a single symmetric per-base rate, but the background that
matters for low-fraction detection is direction-specific: REF-to-ALT and ALT-to-REF
substitution rates differ (oxidation, strand bias, flanking context), and the dominant
direction sets the floor a residual-host signal must clear. allomix can therefore read a
per-site, per-direction empirical error table estimated from a cohort of reference
samples called the same way as the admixture (the `estimate-errors` command pools reads
at homozygous calls: the minority-allele rate at homozygous-reference calls estimates
the REF-to-ALT rate, and at homozygous-alternative calls the ALT-to-REF rate, with a
floor so a zero observed rate cannot make a single stray read produce an infinite
penalty). When supplied, each marker uses its own measured rate in the matching
direction in both the magnitude likelihood and the residual-host presence-test
background, falling back to the global symmetric rate where a direction was not
measured. The table is optional and skipped when none is supplied. The in silico results
below use the global rate, because the simulator draws errors from that same uniform
model, so a per-site table cannot change simulated data; it sharpens detection only on
real reads, where error rates genuinely vary by site and direction.

### Multi-donor estimation

The single-donor model extends directly to two donors plus host: instead of one
fraction, allomix searches the pair of donor fractions (constrained to sum to no more
than one, with the host as the remainder) over a triangular grid, then refines locally
and reports a profile-likelihood interval for each donor (Supplementary Methods). A
marker is informative if the host differs from any donor, and informative counts are
tracked per donor, since related donors can leave few markers that separate one donor
from the other.

### Outlier-resistant refit for aberrant markers

A few markers can sit far off the fit, most often because the recipient's malignant clone
carries a copy-number change there, which is common in the haematological cancers that
bring patients to transplant (see Recipient copy-number aberrations below). allomix
optionally re-estimates the fraction after
down-weighting such markers, using a median-based outlier-resistant scale. The trim is
deliberately one-sided: it removes markers whose deviation points *away* from host
presence, and it protects markers whose deviation points *toward* host presence. This
matters at low host fraction, where the handful of markers carrying the real low-level
host signal sit off a donor-dominated fit and a symmetric outlier rule would discard
exactly those markers, collapsing a real low-fraction host signal to zero. The
trade is explicit: at the LoD we would rather keep a few artifacts than
throw away a real low-fraction host signal. The refit engages only when the number of
off-fit markers exceeds chance expectation, so clean samples are unchanged, and it is
floored so it does not over-trim small panels. When too many markers are excluded, the
sample is flagged for review rather than reported as a confident estimate.

### Quality control and sample-integrity checks

allomix is meant to run inside routine laboratory operations, so it reports a
three-level verdict per sample: PASS, REVIEW (an estimate is produced but a reliability
flag should be checked), or FAIL (the result is not usable). Basic checks cover
marker sufficiency, sequencing depth, and confidence-interval width. A goodness-of-fit
check compares the marker-to-marker scatter against what the model expects and is
computed both before and after the outlier-resistant refit, gated on the worse of the
two, so the refit cannot hide a genuinely bad fit by discarding its own outliers. The
presence test is cross-checked against the magnitude estimate, and a warning is raised
when residual host is detected below the level the magnitude estimate resolves; this
stays a soft warning because its behaviour on real samples is still being characterised.

Three of these checks look for low-fraction signals that the magnitude estimate alone
cannot see, and they are separated by genotype geometry rather than by re-thresholding
one number. The residual-host markers (donor-homozygous, host carrying the donor-absent
allele) are distinct from the consensus-homozygous markers that the contamination and
swap checks both read, and the consensus-homozygous markers are never used by the
magnitude estimate at all:

- **Residual host** shows up at donor-homozygous markers where the host carries the
  donor-absent allele (the presence test above).
- **Contamination by a non-host, non-donor genome** shows up at markers where host and
  every donor are homozygous for the same allele, so the minor allele cannot come from any
  expected contributor and can only be a background artifact or foreign DNA (for example
  another patient co-pooled on the same flowcell via index hopping). allomix estimates
  contamination as the excess minor-allele signal over a data-internal error floor (the
  low percentile of per-site minor fractions), using the median across sites so that a
  few gross miscall sites do not dominate, and capping clearly miscalled sites; the full
  estimator is in the Supplementary Methods. It also separates the magnitude of
  contamination from its mechanism by distinguishing a dose-response signature (foreign
  reads scaling with co-pooled allele dose) from a uniform error elevation. When that
  dose-response is statistically significant, an optional per-marker correction (off by
  default) subtracts the dose-predicted contamination from each donor-homozygous
  host-allele count before the magnitude fit: the per-carrier rate is calibrated per
  run, gated on the consensus-homozygous dose-response significance, and applied at the
  level of the informative donor-homozygous markers themselves, with the flat error
  floor left to the per-site error model. Because the correction only ever removes
  host-allele reads, it is an asymmetric downward adjustment and is therefore gated
  rather than applied unconditionally (full estimator in the Supplementary Methods). An
  optional VCF header field recording the flowcell and lane lets allomix additionally
  flag, from metadata alone, when an admixture sample shares a sequencing run with the
  host (an index-hopping risk), kept separate from the in-data contamination
  estimate.[@Costello2018indexswap]
- **A sample swap or contaminating genome** shows up when those same
  consensus-homozygous markers carry a minority allele that is individually significant
  rather than a faint background. allomix runs a consistency test across these sites
  that catches a wrong-patient VCF the magnitude estimate cannot see, because the
  magnitude estimate only ever looks at informative markers and never at the consensus
  sites.

Two further checks guard sample identity and marker quality. A relatedness check
estimates the kinship between each pair of input samples directly from the genotypes,
using an allele-frequency-free coefficient in the style of
somalier,[@Pedersen2020somalier] which suits a panel-agnostic tool with unknown panel
allele frequencies. Its logic is asymmetric to mirror real failure modes: losing an
expected close relationship is treated as a swap signature and fails, a declared
identical pair always fails (a syngeneic donor is unmeasurable), and a milder
relatedness shift only triggers review. Separately, a read-level artifact filter screens
donor-homozygous markers for soft-clip, read-position, and strand-bias artifacts using
bcftools mpileup annotations, judging strand bias by effect size rather than statistical
significance (at high depth a real allele's mild strand skew is highly significant but
harmless, whereas an artifact is extreme regardless of depth). These checks are run as a
screen on the pileup output rather than delegated to a variant caller's built-in read
filtering, for the same reason the admixture is piled up rather than called: a caller
applies its filters before we ever see the reads and would strip the low-fraction signal
along with the artifacts, whereas screening here flags a suspect marker without discarding
minority reads elsewhere. The filter auto-disables on single-strand amplicon panels, where
almost every marker is read from one strand.
Finally, a genotype/allele-depth consistency guard drops reference-sample markers whose
called genotype contradicts their own allele depths, so a marginally rescued
heterozygous call in a small joint call does not feed systematic bias into the estimate;
this guard is applied only to the reference samples, never to the admixture, whose
allele balance is not expected at 0, 0.5, or 1.

### Simulation framework

For validation, allomix includes a simulator that generates synthetic chimeric VCFs by
blending two genotype VCFs at a specified donor fraction. The simulated fraction is the
donor proportion in the analysed DNA. One caveat applies to every in silico result that
follows. The simulator draws reads from the same mixture and error model the estimator
fits, so the estimator is being tested partly against its own assumptions. In silico
accuracy and LoD are therefore best-case figures that show the estimator
recovers the truth when the data match its model, not evidence that the model matches
real sequencing. The independent check is the real-read SRP434573 analysis (Results),
produced on a panel, platform, and noise process the simulator did not generate. In
particular, overdispersion is the part the simulator cannot claim to have solved: the
main runs draw reads from a binomial, so any noise the shared model omits, of which
overdispersion is the dominant one at clinical depth (Discussion), is absent from the
simulated data by construction rather than handled by the estimator. The simulator
likewise does not generate index hopping or co-pooled cross-sample contamination, so the
background artifacts those produce on real data (Results) are not present in the simulated
results.

On top of the shared mixture and error model, the simulator layers four sources of
measurement noise calibrated from empirical panel data:

1. **Per-marker amplification bias**, drawn from a heavy-tailed distribution calibrated
   from {{ panel_empirical.n_het_total | commas }} heterozygous observations across
   {{ panel_empirical.n_bias_markers | fmt('g') }} markers in
   {{ panel_empirical.n_vcfs | fmt('g') }} joint-called VCFs from a
   {{ panel_specs.n_markers_panel }}-SNP rhAmpSeq sample-identification panel
   (per-marker bias SD {{ panel_empirical.sd_bias }}; Supplementary Table S1).
2. **Non-uniform depth**, drawn from a log-normal distribution matching the empirical
   depth coefficient of variation of {{ panel_empirical.mean_sample_depth_cv }}
   (Supplementary Table S1).
3. **Sequencing errors**, using the same per-base error model as the likelihood, so the
   simulator and estimator share a consistent generative model.
4. **Locus dropout**, with each marker producing zero reads at the empirical no-call
   rate of {{ panel_empirical.mean_nocall_pct }}%.

Reads are drawn from a binomial by default; the simulator can instead draw from an
overdispersed (beta-binomial) distribution at a chosen concentration, which is used to
characterise how the LoD depends on overdispersion (Supplementary Figures
S7, S8) but not in the main validation. To evaluate longitudinal monitoring we simulated
a six-timepoint post-HSCT trajectory (day +14 to day +365) with true donor fractions
from 15% (early engraftment) to 97% (full donor), including a clinically relevant
3-percentage-point dip at day +180, at 500x depth with five independent replicates. Full
generative details, including the copy-number and sibling-segregation extensions below,
are in the Supplementary Methods.

### Limit of detection

Limit of detection (LoD) and limit of blank (LoB) follow the CLSI guideline
EP17-A2,[@CLSIEP17A2] as in published evaluations of comparable NGS chimerism
assays.[@Vynck2021devyser] LoB is the 95th percentile of the estimated donor fraction on
a pure-host blank, and LoD is the lowest true fraction at which at least 95% of
replicates exceed LoB, read from a logistic fit to the dilution series. We computed LoB
and LoD separately for each donor/host pair and report the median across pairs with a
10th-to-90th-percentile band. The sweep covered two relatedness levels (unrelated, full
sibling), five depths (100x to 2,000x), six panel sizes (25 to 400 markers, nested so a
smaller panel is a strict prefix of a larger one), and seven true donor fractions (0 to
5%). A nested design (multiple donor/host pairs, each with multiple sequencing
replicates that vary only read-sampling noise) keeps identity-by-descent variation
between pairs in the reported band rather than leaking it into the central estimate;
details and the full grid are in the Supplementary Methods. Because the simulator draws
reads from a binomial, this LoD is a near-binomial analytical best case rather than a
validated assay limit; a real assay's LoD can only be higher (Results, Discussion).

### Recipient copy-number aberrations

HSCT recipients are usually haematological-malignancy patients whose residual or
relapsing clone often carries copy-number changes, which break the assumption that every
locus is diploid. To test sensitivity to this, the simulator applies a per-marker
recipient aberration before read sampling, modelling the recipient as a mixture of
normal diploid cells and an aberrant clone, with the expected allele fraction a
copy-number-weighted average rather than the diploid mean. Three aberration types are
modelled (copy-neutral loss of heterozygosity, single-copy deletion, and single-copy
gain) at a controllable burden, applied only to the admixture sample to match the
two-phase workflow in which recipient genotypes come from a clean reference. We computed
the EP17-A2 LoD in two directions (recipient clone as the minor component, the relapse
direction; and donor as the minor component against an aberration-bearing recipient
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
single-end inserts (each marker read from a single strand). Pairwise two-person mixtures
and one three-person mixture were prepared at known major:minor ratios, giving
minor-contributor fractions of {{ srp434573.dilution_ladder }}.

We assigned the minor (titrated) contributor to the host (recipient) role and the major
contributor to the donor. This mirrors the usual clinical situation, where a patient is
near full donor chimerism and the recipient is the small residual fraction the assay is
trying to measure, so the monitored quantity is the host fraction and each dilution
series reads as a declining-chimerism (relapse) trajectory. These samples are not
biological host and donor, and the thesis does not label which contributor is which; the
assignment is purely a labelling convention based on mixing fraction and affects only the
host/donor labelling, not the genotyping or the recovered fractions.
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

To obtain a LoD on real reads rather than simulation, we subsampled the
high-depth SRP434573 mixtures until their LoD rose into the dilution series' measurable
window. From each two-person mixture we sub-sampled reads to a target mean depth (100x
to 2,000x) by a single per-sample binomial keep rate applied uniformly across markers,
which is the statistical analogue of FASTQ read subsampling and preserves the real
locus-to-locus depth variation, then re-applied the depth filter so low-depth loci drop
out. Panel size was set by nested random subsets of the informative markers (prefixes of
one permutation per mixture and replicate, so each curve is monotonic in panel size), and
multiple sub-sampling replicates were drawn per depth-and-panel cell. For each cell we
read the LoD off the same logistic fit as the simulated sweep, separately for the
magnitude estimate (a sample is detected when its 95% confidence interval for the host
fraction excludes zero) and the residual-host presence test (detected when the
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
upstream calling (GATK joint genotyping of host and donor, then forced bcftools pileup
of the admixture samples at the panel sites) and the optional bwa-mem alignment step
ahead of it. Each analysis can emit a tab-separated summary, a structured JSON record,
and a single-file HTML report that presents the estimate, its confidence interval, the
quality-control verdict, and collapsible method explanations for clinical review. The
HTML report is rendered from templates a laboratory can override to match its own
reporting style (logo, layout, and wording). The repository documentation includes a
guide for qualifying an existing clinical panel for chimerism use, covering
informative-marker sufficiency, per-marker characterization, inclusion thresholds, and
building the bias and error correction tables.
