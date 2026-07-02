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
runs without new reagents or a proprietary kit.

allomix takes standard VCF files and reports the donor fraction, with single-donor or
multi-donor (host plus up to two donors) support, profile-likelihood confidence
intervals, and a built-in suite of quality and sample-integrity checks. It runs two
complementary tests: a maximum-likelihood estimate of how much donor is present, and a
separate residual-host presence test that asks whether any host remains below the
quantification limit.

On a public dataset of real titrated DNA mixtures, allomix recovered known host fractions
from 10% down to 1% host, resolved a three-person mixture, and its residual-host presence
test returned a positive call at all {{ srp434573.presence_n_total | dp(0) }} mixtures
with no false positive at the {{ srp434573.zero_host_n | dp(0) }} pure-donor
(true-0%-host) controls. Measured by subsampling those mixtures, the real-data limit of
detection was about {{ subsample_lod_headline.mle_lod_1000x_100markers_pct | dp(0) }}%
host at a reference operating point of 100 markers and 1,000x depth. Once the dataset's
co-pooled contamination confound is removed, two independent estimates put the method's
own floor near 0.1% to 0.2% host: a semi-synthetic remix of the same real reads and a
near-binomial simulation agree, which places the ~1% ceiling on this dataset's
contamination and dilution limits rather than the estimator. These lower figures are
analytical bounds, not wetlab assay limits, so a real assay's limit can only be higher.
Across simulated panels at depths from 50x to 1,000x, mean absolute error stayed below 1%
(from {{ depth_1000.mean_abs_error_pct | dp(2) }}% at 1,000x to
{{ depth_50.mean_abs_error_pct | dp(2) }}% at 50x). This paper presents the tool with its
analytical characterization and a real-data demonstration; clinical validation against
STR chimerism remains to be done and, because allomix is panel-agnostic, is carried out
by each laboratory on its own panel.
