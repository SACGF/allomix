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
Each gives its own sensitivity figure. At 100 markers and 1,000x depth the in silico
limit of detection (the lowest fraction recovered in at least 95% of replicates, CLSI
EP17-A2; unrelated donor) was {{ lod_headline.unrelated_lod_1000x_100markers_pct }}% for
the magnitude estimate and
{{ presence_lod_curve_headline.presence_unrelated_lod_1000x_100markers_pct }}% residual
host for the presence test, within the range reported for commercial NGS chimerism kits,
with the caveat that these are best-case figures from near-binomial simulated data
rather than wetlab assay limits. Across simulated panels at depths from 50x to 1,000x,
mean absolute error stayed below 1% (from {{ depth_1000.mean_abs_error_pct | dp(2) }}%
at 1,000x to {{ depth_50.mean_abs_error_pct | dp(2) }}% at 50x), and a simulated
six-timepoint engraftment trajectory was tracked accurately including a small dip. On a
public dataset of real titrated DNA mixtures, allomix recovered known host fractions
accurately from 10% down to 1% and resolved a three-person mixture, and the residual-host
presence test returned a positive call at all {{ srp434573.presence_n_total | dp(0) }}
mixtures with no false-positive call at the {{ srp434573.zero_host_n | dp(0) }} pure-donor
(true-0%-host) controls the same dataset provides. At the lowest 0.5% titration the
point estimates scatter against a co-pooled contamination floor in that dataset, which the
presence call there cannot be separated from. This paper presents the tool with its analytical
characterization and a real-data demonstration; clinical validation against STR
chimerism remains to be done and, because allomix is panel-agnostic, is necessarily
carried out by each laboratory on its own panel rather than once for the tool.
