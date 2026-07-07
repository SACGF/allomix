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

Monitoring the fraction of donor-derived cells in the recipient (chimerism) after allogeneic
hematopoietic stem cell transplantation guides decisions about relapse and graft
rejection. Standard chimerism testing is usually run as a separate assay. Many clinical sequencing panels
already include biallelic polymorphic markers for sample identification, sequenced at
high depth during routine testing. We present allomix, an open-source tool that turns
those markers into a quantitative chimerism readout, letting a laboratory add monitoring
to a panel it already runs without new sequencing workflows.

allomix takes standard VCF files and reports the donor fraction for single-donor or
multi-donor (recipient plus up to two donors) mixtures, with profile-likelihood confidence
intervals and quality checks. It runs a maximum-likelihood magnitude estimate alongside a
separate test for residual recipient below the quantification limit.

On public titrated DNA mixtures, allomix recovered known recipient fractions from 10% to 1%,
resolved a three-person mixture, and called residual recipient at all
{{ srp434573.presence_n_total | dp(0) }} mixtures with no false positive in the
{{ srp434573.zero_host_n | dp(0) }} pure-donor controls. Subsampling gave a real-data
limit of detection near {{ subsample_lod_headline.mle_lod_1000x_100markers_pct | dp(0) }}%
recipient at 100 markers and 1,000x depth; semi-synthetic and near-binomial analyses place the
estimator's floor near 0.1% to 0.2%. Across simulated panels from 50x to 1,000x depth,
mean absolute error stayed below 1%. These are analytical bounds; clinical validation against
standard chimerism methods is left to each adopting laboratory.
