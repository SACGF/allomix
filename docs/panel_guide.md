# Qualifying Your Own Panel for allomix

allomix is panel-agnostic: it works on whatever bi-allelic markers (SNPs or
indels) are present in your input VCFs. The common case is not designing a panel
from scratch, it is repurposing a panel your lab already runs (a sample-ID,
forensic, ancestry, or pharmacogenomics SNP panel) for post-transplant chimerism
monitoring. This guide covers how to check that such a panel is fit for the job
and how to build the correction tables allomix uses.

The tooling for this already exists but is spread across `scripts/`, CLI
subcommands, and two other docs. This guide stitches those together in the order
you would actually use them. It links into the
[Joint Calling Guide](joint_calling.md) and [Bias Estimation Guide](estimate_bias.md)
rather than repeating them.

Steps, in order:

1. [Do I have enough informative markers?](#1-do-i-have-enough-informative-markers)
2. [Define the panel sites](#2-define-the-panel-sites)
3. [Characterize the panel on a reference cohort](#3-characterize-the-panel-on-a-reference-cohort)
4. [Set marker inclusion thresholds](#4-set-marker-inclusion-thresholds)
5. [Build the correction tables](#5-build-the-correction-tables)
6. [Appendix: selecting a panel from scratch](#appendix-selecting-a-panel-from-scratch)
7. [Known limitations](#known-limitations)

## 1. Do I have enough informative markers?

A marker is informative for a given transplant only if it distinguishes host from
donor, so informative-marker counts are **pair-specific**: the same panel gives
different counts for different host/donor pairs, and unrelated pairs share more
informative markers than siblings.

The literature gives a usable rule of thumb. Vynck et al. (2022, Clin Chim Acta
532:123-129) showed that as few as three informative markers can support
quantification, with accuracy improving as markers are added, and that a panel of
about 20 markers with minor allele frequencies (MAFs) near 0.5 gives more than a
95% chance of at least three informative markers even for sibling pairs. Sample-ID
marker sets with tens of markers therefore have comfortable headroom for unrelated
donors, and the margin narrows (but usually holds) for siblings. If you need to
know whether a specific panel will be sufficient for a specific donor-host pair
before sequencing, the FABCASE tool assesses that prospectively (Vynck 2025, Int J
Lab Hematol 47:690-697).

Where allomix reports the count: every `allomix detect` run reports how many
input markers were informative and how many were used in the fit (the HTML report
states it as "N of M input markers were informative; K used in the fit", and the
same `n_informative` / per-donor counts are in the structured JSON). For a
multi-donor run the per-donor informative counts are reported separately. Use
these to confirm each real case clears the rule-of-thumb floor, and watch for
pairs that sit near three informative markers, where precision will be weakest.

Handling of markers in approximate linkage disequilibrium is discussed in
[issue #9](https://github.com/SACGF/allomix/issues/9); for panel design the
relevant point is spacing markers so they are close to independent (see the
[appendix](#appendix-selecting-a-panel-from-scratch)).

## 2. Define the panel sites

allomix needs the panel loci as genomic positions so the upstream pipeline can
genotype them. Which tool you use depends on what you already have.

**You have a BED of panel positions.** Turn it into a sites-only VCF with
`scripts/build_panel_vcf.py`, which reads REF/ALT from a source VCF and writes a
VCF suitable for forced genotyping:

```bash
python scripts/build_panel_vcf.py panel.bed source.vcf.gz panel.sites.vcf.gz --bgzip bgzip --tabix tabix
```

**You want to filter a population VCF to your capture region.** For the GATK
force-output path, `scripts/build_force_output_panel.sh` filters a population VCF
(for example gnomAD v4.1 sites) to positions overlapping your capture BED that are
PASS, biallelic SNPs, above an allele-frequency threshold:

```bash
scripts/build_force_output_panel.sh capture.bed gnomad.sites.vcf.gz 0.05 panel_alleles.vcf.gz
```

**You have BAMs but no usable BED.** Amplicon and MIP panels sometimes publish
only a kit name or an awkward probe list, or nothing. Because these assays
concentrate reads into tight high-depth footprints, the panel can be recovered
from coverage alone with `scripts/recover_panel_bed.py`. It keeps positions
covered at `--min-depth` in at least `--min-samples` BAMs and merges adjacent kept
positions into one interval per amplicon:

```bash
python scripts/recover_panel_bed.py --bam-glob 'output/bam/*.bam' --out panel.bed --min-samples 50
```

Set `--min-samples` relative to your cohort size (a clear majority of samples is a
reasonable starting point). Lower it to tolerate more per-sample dropout at the
cost of admitting more off-target positions. The SRP434573 public dataset is a
worked example: its MIP panel ships no coordinates, and
`paper/scripts/build_srp434573_panel_bed.sh` (a thin wrapper around
`recover_panel_bed.py`) recovers 1052 hg38 intervals from 64 BAMs. See
`paper/public_data/SRP434573/README.md` for that write-up.

Once you have the sites, the [Joint Calling Guide](joint_calling.md) covers what
to do with them: GATK joint calling of the host/donor reference samples for
genotypes, and forced `bcftools mpileup` at the same sites for the admixture
samples.

## 3. Characterize the panel on a reference cohort

Run a set of reference samples (any cohort genotyped through your normal pipeline;
the samples used for bias training work well) and measure how each marker behaves.
`scripts/measure_panel_bias.py` reports per-marker characteristics from
joint-called genotyping VCFs:

```bash
python scripts/measure_panel_bias.py vcf_list.txt --output output/panel_stats
```

It writes two files:

- `output/panel_stats_per_marker.tsv` — one row per marker, with `call_rate`,
  `n_nocall`, `mean_depth`, `depth_cv`, `median_bias` / `mean_bias` (het VAF
  deviation from 0.5), `het_ratio_vs_hwe` (an allele-dropout signal), and the het
  count available for bias estimation.
- `output/panel_stats_facts.csv` — one row of panel-wide summaries
  (`mean_nocall_rate`, `markers_gt5pct_nocall`, `mean_depth`, `p95_abs_bias`,
  `max_abs_bias`, `markers_low_het`, an `ado_estimate`, and so on).

A worked output is committed in `paper/empirical_results/` (the 76-marker
rhAmpSeq sample-ID panel run through this script), which is a useful reference for
what healthy per-marker numbers look like.

Before using samples for bias training, screen them with
`scripts/qc_bias_samples.py`, which flags samples with too many no-calls, an
implausible heterozygosity rate, or a skewed VAF balance:

```bash
python scripts/qc_bias_samples.py joint.vcf.gz --output-samples pass_samples.txt --output-metrics qc_metrics.tsv
```

This is the same screen described in
[the Bias Estimation Guide, step 3](estimate_bias.md#3-sample-level-qc); the pass
list feeds the bias estimation in step 5 below.

## 4. Set marker inclusion thresholds

The characterization in step 3 produces statistics, not decisions. This section
suggests how to turn them into a marker-inclusion list. **These are starting
points, not validated cutoffs.** Tune them against your own cohort, and expect a
sample-ID panel already validated for identity work to need little pruning.

allomix has no built-in per-marker inclusion command yet (see
[Known limitations](#known-limitations)); for now this is a manual pass over
`panel_stats_per_marker.tsv`. Candidate criteria for dropping a marker:

- **Low call rate / high dropout.** Markers that frequently fail to genotype add
  no signal and inflate noise. A reasonable starting cutoff is dropping markers
  with `call_rate` below about 0.90 (equivalently no-call rate above ~0.10, the
  default `qc_bias_samples.py` uses at the sample level). The
  `markers_gt5pct_nocall` panel-level count tells you how many markers are even
  near this line.
- **Extreme amplification bias.** A large, consistent het-VAF offset from 0.5
  (`mean_bias` / `median_bias`) means the marker systematically over- or
  under-calls one allele. Bias correction (step 5) handles moderate bias, so you
  do not need to drop markers just because they are biased. Reserve exclusion for
  the extreme tail (for example `|mean_bias|` well beyond the panel `p95_abs_bias`),
  where the correction is least reliable.
- **Allele dropout.** A `het_ratio_vs_hwe` well below 1 (a het deficit relative to
  Hardy-Weinberg expectation) suggests one allele drops out at that marker, which
  biases fractions in a way bias correction does not fully fix. Treat a strong het
  deficit as a reason to drop.
- **Unstable depth.** A very high `depth_cv` at a marker means its depth swings
  widely across samples, which weakens the per-marker weighting. This is usually a
  secondary criterion; combine it with the above rather than dropping on depth
  variability alone.

Keep the excluded-marker list under version control with a one-line reason per
marker, so the panel definition is auditable. Re-running step 3 after pruning
confirms the panel-level summaries (mean call rate, bias percentiles) improved as
expected.

## 5. Build the correction tables

With a qualified marker set, build the tables allomix uses to correct systematic
effects. These are calibration steps you run once per panel (and refresh if the
assay changes).

- **Per-marker bias** (`allomix estimate-bias`): the multiplicative het-VAF
  correction characterized in step 3. Full workflow, input modes, and output
  format are in the [Bias Estimation Guide](estimate_bias.md). Use the
  QC-passed sample list from step 3 as the training set.
- **Per-site error rates** (`allomix estimate-errors`): empirical per-direction
  error rates from force-called hom-ref background positions, used by the
  host-presence detection test. Run it on the same reference cohort.
- **Co-pooled contamination** (`allomix build-contamination-table`): a per-marker
  table of dose-predicted co-pooled contamination, for labs that pool samples and
  see low-level cross-contamination (see
  [issue #12](https://github.com/SACGF/allomix/issues/12)). Build it from a joint
  VCF containing host, donor(s), and the co-pooled samples.

The [Joint Calling Guide](joint_calling.md) describes how the upstream Snakemake
pipeline automates the genotyping these tables depend on. The `estimate-bias`
subcommand is documented in [CLI usage](cli.md#estimate-bias); `estimate-errors`
and `build-contamination-table` are run the same way (`allomix <command> --help`
lists their options).

## Appendix: selecting a panel from scratch

If you are choosing markers rather than qualifying an existing panel, the design
goals are well established and the same across sample-ID panels:

- **MAF near 0.5.** Markers with allele frequencies near 0.5 maximize the chance
  that a random host/donor pair differs at the marker, which is what makes it
  informative. This is why the ~20-marker rule of thumb in step 1 assumes MAFs
  near 0.5.
- **Approximate linkage equilibrium.** Space markers far apart across the genome
  so their genotypes are close to independent; clustered markers double-count the
  same signal. See [issue #9](https://github.com/SACGF/allomix/issues/9) on
  handling residual linkage disequilibrium.
- **Genome distribution.** Spreading markers across chromosomes reduces the chance
  that a single copy-number or loss-of-heterozygosity event in the recipient
  distorts many markers at once.

A practical shortcut is to adopt the sites list distributed with
[Somalier](https://github.com/brentp/somalier), a widely used relatedness and
sample-swap QC tool. Its released sites VCFs (a curated set of about 17,700 common
SNPs, provided for hg19, GRCh37, hg38, and a CHM13/T2T build) were selected for
allele frequencies near 0.5 and spread across the genome, so they already match
the design goals above. Using the same sites has a second payoff: you can run
Somalier itself on the sequenced host, donor, and admixture samples for
relatedness and sample-swap checks, a useful independent cross-check on
donor/host identity before chimerism analysis. Pull the sites file for your build
from the Somalier [releases page](https://github.com/brentp/somalier/releases) and
feed it into the site-definition step above.

For prospective per-pair sufficiency (will this panel resolve this specific
donor-host pair), use FABCASE (Vynck 2025) rather than re-deriving the statistics.
Designing a new panel is a larger exercise than this guide covers; the references
above are the right starting points.

## Known limitations

Two things a guide can point at but the tool does not yet automate:

- **No built-in per-marker inclusion command.** Step 4 is currently a manual pass
  over the per-marker statistics. A convenience command (`allomix panel-qc`) that
  applies recommended cutoffs and emits a marker-inclusion report is a candidate
  follow-up, not yet implemented.
- **No integrated panel-sufficiency check.** allomix reports informative-marker
  counts per run but does not assess prospective per-pair sufficiency. Use FABCASE
  and the ~20-marker rule of thumb for that.
