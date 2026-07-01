---
title: "allomix: donor chimerism monitoring from existing clinical NGS panels"
author:
  - David M Lawrence
  - Xiaowen Wu
  - Hamish S Scott
  - Yasir Kusay
  - Chung Hoow Kok
  - David M Ross
  - Anna L Brown
  - Wendy T Parker
---

## Abstract

Monitoring the fraction of donor-derived cells (chimerism) after allogeneic
hematopoietic stem cell transplantation (HSCT) guides decisions about relapse and graft
rejection, but it is usually run as a separate dedicated assay. Many clinical
next-generation sequencing (NGS) panels already include biallelic polymorphic markers
for sample identification or quality control, sequenced at high depth as part of routine
testing. We present allomix, an open-source tool that turns those markers into a
quantitative chimerism readout, so a laboratory can add monitoring to a panel it already
runs without new reagents or a proprietary kit. allomix takes standard VCF files and
reports the donor fraction, with single-donor or multi-donor (host plus up to two
donors) support, profile-likelihood confidence intervals, and a built-in suite of
quality and sample-integrity checks. It runs two complementary tests: a
maximum-likelihood estimate of how much donor is present, and a separate residual-host
presence test that asks whether any host signal remains below the quantification limit.
Each gives its own sensitivity figure. On a public dataset of real titrated DNA
mixtures, allomix recovered known host fractions accurately from 10% down to 1% host,
tracked a dilution ladder read as a declining-chimerism monitoring trajectory, and
resolved a three-person mixture, and the residual-host presence test returned a positive
call at all {{ srp434573.presence_n_total | dp(0) }} mixtures with no false-positive call
at the {{ srp434573.zero_host_n | dp(0) }} pure-donor (true-0%-host) controls the same
dataset provides. Measured by subsampling those mixtures, the real-data limit of
detection was about {{ subsample_lod_headline.mle_lod_1000x_100markers_pct | dp(0) }}%
host at the reference operating point (100 markers, 1,000x depth); this is the headline
sensitivity figure. The method floor with the dataset's confounds removed sits well below
that, from two independent estimates that agree: a semi-synthetic remix of the same real
reads (co-pooled contamination confound removed) tracked host fractions down to about
0.1%, and, as a best-case analytical ceiling from near-binomial simulated data, the in
silico limit of detection (the lowest fraction recovered in at least 95% of replicates,
CLSI EP17-A2; unrelated donor, same operating point) was
{{ lod_headline.unrelated_lod_1000x_100markers_pct }}% for the magnitude estimate and
{{ presence_lod_curve_headline.presence_unrelated_lod_1000x_100markers_pct }}% residual
host for the presence test. That the real-reads and pure-simulation floors agree
indicates the ~1% real-data ceiling reflects this dataset's contamination and dilution
limits rather than the estimator; these in silico figures are analytical bounds, not
wetlab assay limits, so a real assay's limit can only be higher. They sit within the
range reported for commercial NGS chimerism kits, but those vendor numbers are
clinically-validated dilution-series limits, not a head-to-head comparison. Across
simulated panels at depths from 50x to 1,000x, mean absolute error stayed below 1% (from
{{ depth_1000.mean_abs_error_pct | dp(2) }}% at 1,000x to
{{ depth_50.mean_abs_error_pct | dp(2) }}% at 50x). At the lowest 0.5% titration the
point estimates scatter against a co-pooled contamination floor in that dataset, which the
presence call there cannot be separated from. This paper presents the tool with its analytical
characterization and a real-data demonstration; clinical validation against STR
chimerism remains to be done and, because allomix is panel-agnostic, is necessarily
carried out by each laboratory on its own panel rather than once for the tool.
