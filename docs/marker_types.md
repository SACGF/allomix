# Marker types and what they are used for

allomix classifies every marker by comparing the host and donor genotypes, then
uses each class for a different job: estimating the donor fraction, testing
whether any host remains, correcting for co-pool contamination, or checking
sample integrity. This page explains that classification and which check uses
each marker class. For where the classification lives in the code, see
[`genotype.classify_markers`](architecture.md); for choosing a panel with enough
informative markers, see the [Panel Guide](panel_guide.md).

## Informative vs non-informative

A marker is **informative** only when the host and donor have different
genotypes. Those markers carry the mixing-fraction signal and drive the
donor-fraction estimate. Markers where host and donor match carry no fraction
information and are set aside, but they are not wasted: the consensus matches
feed the contamination and sample-swap checks instead, which are independent of
the fraction estimate.

Informative-marker counts are **pair-specific**. The same panel gives different
counts for different host/donor pairs, and unrelated pairs share more informative
markers than siblings. Every `allomix detect` run reports how many input markers
were informative and how many were used in the fit.

## The eight marker classes

Among informative markers, the size of the signal depends on the genotype
contrast. allomix sorts every informative marker into one of six genotype-contrast
types (tracked per donor, following Vynck et al.), plus two non-informative
consensus classes. The table lists all eight and which check uses each.

| Type | Host | Donor | Contrast | Fraction estimate | Residual-host presence | Contamination correction | Consensus QC |
|:--|:--|:--|:--|:--|:--|:--|:--|
| 0  | 0/0 | 1/1 | full | yes | yes | yes | no |
| 1  | 1/1 | 0/0 | full | yes | yes | yes | no |
| 10 | 0/1 | 0/0 | half | yes | yes | no  | no |
| 11 | 0/1 | 1/1 | half | yes | yes | no  | no |
| 20 | 0/0 | 0/1 | half | yes | no  | no  | no |
| 21 | 1/1 | 0/1 | half | yes | no  | no  | no |
| cons-hom | hom | hom (same) | none | no | no | no | contamination + swap |
| cons-het | 0/1 | 0/1 | none | no | no | no | shared-het balance |

**Full-contrast types (0, 1).** Host and donor are homozygous for opposite
alleles, so the minority allele has a single possible source. These give the
maximum allelic contrast.

**Half-contrast types (10, 11, 20, 21).** One party is heterozygous. That party
already contributes the allele on its own, so the other party only nudges the
count, and the signal is half the size of a full-contrast marker.

**Which check uses which type.**

- **Donor-fraction estimate:** every informative type (0, 1, 10, 11, 20, 21).
- **Residual-host presence test:** only the donor-homozygous types (0, 1, 10,
  11), because the test needs a clean donor-absent allele to count. Types 20 and
  21 leave no donor-absent allele.
- **Co-pool contamination correction:** only the fully homozygous-contrast types
  (0, 1), where the host (donor-absent) allele reads are otherwise a clean
  background.
- **Consensus QC:** the two non-informative classes. At consensus-homozygous
  markers the minority allele can only be a background artifact or foreign DNA,
  which drives contamination and sample-swap detection. At consensus-heterozygous
  markers the alternative-allele fraction should sit near 0.5 whatever the mixing
  fraction, so a systematic skew flags contamination, allelic imbalance, or a
  sample mix-up.

By default a site is used only if host and donor genotype quality is at least 20
and the admixture depth is at least 100, and at least three informative markers
are required to report an estimate.

## Worked example: reading the residual host off the counts

Take two markers in one post-transplant sample from a patient near full donor
chimerism, the usual clinical situation, where the sample is mostly donor and the
question is how much host remains. At each marker, host and donor genotypes were
fixed before transplant, so we know which person can make which allele.

**Marker 1 is fully informative (type 0): host is A/A, donor is G/G.** The donor
makes nothing but G and the host never makes G, so every A read can only have come
from the host. In a near-full-donor sample the minority of A reads is the residual
host fraction, and the majority of G reads is the donor.

| Counts at Marker 1 | A | G | Reading |
|:---|--:|--:|:---|
| Residual-host reading (minority allele) | 30 | 970 | 30 A reads can only be host = **3% residual host** |
| Magnitude reading (majority allele) | 30 | 970 | 970/1000 G reads = 97% donor, i.e. **3% host remaining** |

The same marker answers two questions. The minority allele, the one only the host
could have produced, tells you *whether any host remains*, which is the clinically
urgent question near full donor chimerism. The majority shift gives the magnitude,
*how much donor* is present and so how much host. These are the two complementary
tests allomix runs.

**Marker 2 is partially informative (type 11): host is A/G, donor is G/G.** The
host carries the donor-absent allele A on only one of its two alleles, so a
residual host puts A reads into the sample at half its true fraction. The host
signal is real but half the size of a fully informative marker.

| Counts at Marker 2 | A | G | Reading |
|:---|--:|--:|:---|
| Residual-host reading (minority allele) | 15 | 985 | 15 A reads = 1.5% observed, doubled for the half-signal = **3% residual host** |

This has three consequences. A pure-donor sample still shows a low background of
reads carrying the host-only allele, produced by background artifacts (miscalled
bases and, on co-pooled runs, index hopping and low-level cross-sample
contamination), which is why there is a measured floor rather than a clean zero.
One marker on its own is noisy, which is why allomix pools dozens of markers into
a single estimate with a confidence interval. And markers sometimes disagree by
more than simple read-counting would predict, which is why the interval is widened
to reflect that extra scatter.
