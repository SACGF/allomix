# Bias Estimation with allomix estimate-bias

Capture and amplification panels introduce systematic per-marker shifts in observed VAF. A marker that should read 0.5 in a heterozygous sample might consistently read 0.45 or 0.55 due to differential amplification efficiency. These biases are consistent across samples sequenced with the same panel, so they can be estimated once from a training set and then corrected during chimerism estimation.

Bias correction is optional but recommended, particularly at low donor fractions (<5%) where small systematic errors in expected allele frequencies have the greatest impact on accuracy.

## When to run

Run `estimate-bias` once per panel, using a set of genotyping VCFs from samples with known heterozygous calls across the marker set. Re-run if you change the capture panel, library preparation protocol, or sequencer platform.

## What VCFs to use

The estimation relies on markers called heterozygous (`GT 0/1`) and assumes true heterozygotes sit at VAF = 0.5, so any deviation is attributed to panel bias. This means the training samples must have clean germline heterozygotes: samples where LOH, copy number alterations, or clonal dominance could shift het VAF away from 0.5 for biological reasons will corrupt the estimate.

**Use donor VCFs.** Donors are healthy individuals whose blood gives reliable germline heterozygotes at the marker loci. You are already genotyping them as part of the transplant workup, so donor VCFs accumulate naturally over time without extra effort.

Host VCFs from haematology patients are generally not suitable as training data. Pre-transplant blood from haem patients can have LOH, chromosomal amplifications, or clonal hematopoiesis that shifts het VAF away from 0.5 at individual markers. If host samples come from a germline source (buccal swab, skin, nail clipping), they can be included.

Do not use admixture timepoint VCFs: post-HSCT blood is a mixture of donor and host, so markers are not truly heterozygous.

More samples give more reliable estimates. With fewer than ~10 donors, estimates at low-frequency markers will be noisy. Aim for 20+ samples if available. The median estimator is robust to sporadic outliers (a donor with an incidental CNV at one marker will not dominate), but systematic bias from a bad sample type will not average out.

The VCFs do not need to be joint-called for this step. Independent per-sample VCFs are fine.

## Command

Two input modes are supported. Use whichever matches how your donor VCFs are organised.

### Multiple per-sample VCFs

```bash
allomix estimate-bias \
    --vcfs sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz \
    --output bias_table.tsv
```

Shell globbing works:

```bash
allomix estimate-bias \
    --vcfs /path/to/genotyping/*.vcf.gz \
    --output bias_table.tsv
```

The first sample in each VCF is used. This mode suits workflows where each genotyping run produces a separate per-sample VCF.

### Named samples from a joint-called VCF

```bash
allomix estimate-bias \
    --vcf joint_called.vcf.gz \
    --samples DONOR_001 DONOR_002 DONOR_003 \
    --output bias_table.tsv
```

Use this when your donor genotyping VCFs were joint-called alongside other samples and all donors are in a single multi-sample VCF. The sample names must match those in the VCF header.

### Options

| Option | Default | Description |
|---|---|---|
| `--vcfs` | | Per-sample VCFs, one per file (mutually exclusive with `--vcf`) |
| `--vcf` | | Joint-called multi-sample VCF (use with `--samples`) |
| `--samples` | | Sample names to extract from `--vcf` |
| `--output` / `-o` | `bias_table.tsv` | Output path for the bias table TSV |
| `--min-het` | 1 | Minimum heterozygous observations required to include a marker |

**`--min-het`**: Markers with fewer than this many het observations are excluded from the table. With a small training set (e.g. 5 donors), the default of 1 means a bias estimate can be based on a single observation, which is unreliable. With larger training sets, raising this to 5 or 10 filters out markers with sparse data.

## Output format

The bias table is a TSV with one row per marker:

```
chrom   pos     ref  alt  bias       n_het
chr1    925952  G    A    -0.012341  23
chr1    931279  A    G     0.008712  19
...
```

- **bias**: median(observed\_het\_VAF - 0.5) across training samples. Positive means the ALT allele is preferentially amplified; negative means the REF allele is.
- **n\_het**: number of heterozygous observations used to estimate bias at this marker.

## Using the bias table

Pass the bias table to `monitor` or `timeline` with `--bias-table`:

```bash
allomix monitor \
    --vcf patient001_joint.vcf.gz \
    --host-sample HOST_001 \
    --donor-sample DONOR_001 \
    --sample TP1_20240101 \
    --bias-table bias_table.tsv \
    --output results.tsv
```

Markers in the admixture VCF that are not present in the bias table are analysed without correction. Markers in the bias table that are not present in the VCF are silently ignored.

To run without correction even when a bias table is specified, use `--no-bias-correction`.
