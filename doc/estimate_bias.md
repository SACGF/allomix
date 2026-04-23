# Bias Estimation with allomix estimate-bias

Capture and amplification panels introduce systematic per-marker shifts in observed VAF. A marker that should read 0.5 in a heterozygous sample might consistently read 0.45 or 0.55 due to differential amplification efficiency. These biases are consistent across samples sequenced with the same panel, so they can be estimated once from a training set and then corrected during chimerism estimation.

Bias correction is optional but recommended, particularly at low donor fractions (<5%) where small systematic errors in expected allele frequencies have the greatest impact on accuracy.

## When to run

Run `estimate-bias` once per panel, using a set of genotyping VCFs from samples with known heterozygous calls across the marker set. Typical sources:

- A batch of host or donor genotyping VCFs from patients processed with the same panel
- Any VCFs produced by the same sequencing protocol where true heterozygotes can be identified from the genotype calls

Re-run if you change the capture panel, library preparation protocol, or sequencer platform.

## What VCFs to use

Use **genotyping VCFs** (host or donor samples), not admixture timepoint VCFs. The estimation relies on markers called heterozygous (`GT 0/1`). Admixture samples are rarely truly heterozygous at informative loci.

More samples give more reliable bias estimates. With fewer than ~10 samples, estimates at low-frequency markers will be noisy. Aim for 20+ samples if available.

The VCFs do not need to be joint-called for this step. Independent per-sample VCFs are fine.

## Command

```bash
allomix estimate-bias \
    --vcfs sample1.vcf.gz sample2.vcf.gz sample3.vcf.gz \
    --output bias_table.tsv
```

Multiple VCFs can be listed directly or using shell globbing:

```bash
allomix estimate-bias \
    --vcfs /path/to/genotyping/*.vcf.gz \
    --output bias_table.tsv
```

### Options

| Option | Default | Description |
|---|---|---|
| `--vcfs` | required | One or more VCF files to use as training data |
| `--output` / `-o` | `bias_table.tsv` | Output path for the bias table TSV |
| `--min-het` | 1 | Minimum heterozygous observations required to include a marker |

**`--min-het`**: Markers with fewer than this many het observations are excluded from the table. With a small training set (e.g. 5 VCFs), the default of 1 means a bias estimate can be based on a single observation, which is unreliable. With larger training sets, raising this to 5 or 10 filters out markers with sparse data.

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
