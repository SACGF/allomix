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

If you do not have enough donor VCFs yet, an alternative path is to joint-call a cohort from archived BAMs on the same panel and filter for clean samples. See [Building a training cohort from BAMs](#building-a-training-cohort-from-bams) below.

## Building a training cohort from BAMs

If you do not have enough clean donor VCFs but do have archived BAMs from the same panel, you can build a training set end-to-end. This suits labs where donor VCFs accumulate slowly, or where historical data was only kept as BAMs. Because these BAMs are typically from diseased-blood samples, they carry the CNV and clonal hematopoiesis risks described above, so the QC stage is not optional.

### 1. Assemble a sample CSV

Select BAMs from a panel version that captures your marker set. Exclude older runs that predate coverage at those sites, and aim for 100-200 samples so each marker accumulates enough hets.

```bash
# Example: idt_haem BAMs from 2024 onward
ls /tau/data/clinical_hg38/idt_haem/Haem_2[456]_*/1_BAM/*.hg38.bam \
  | shuf -n 200 > sample_bams.txt

# Convert to the samples.csv format the pipeline expects
{
  echo "sample_id,bam_filename"
  awk '{
    n = split($0, parts, "/")
    name = parts[n]
    sub(/\.hg38\.bam$/, "", name)
    print name "," $0
  }' sample_bams.txt
} > samples_bias_training.csv

# Check sample_ids are unique (empty output = all unique)
cut -d, -f1 samples_bias_training.csv | tail -n +2 | sort | uniq -d
```

If duplicates appear, prepend the run folder to the sample_id to make them unique.

### 2. Joint-call at marker sites

Run `pipeline/Snakefile` with `intervals` set to the marker bed. Restricting to the bed keeps the job cheap (a handful of sites vs the whole panel) and produces hom-ref calls with AD at every site for every sample, which is what bias estimation needs.

```bash
snakemake -s pipeline/Snakefile \
  --config ref=/path/to/hg38.fa \
           samples_csv=samples_bias_training.csv \
           intervals=/tau/ngs_pipelines/hg38_reference_files/capture_kits/idt_rhampseq_sid/v1/idt_rhampseq_sid_SNPsQC.bed \
           output_dir=output/bias_training \
  --cores $(nproc)
```

Output: `output/bias_training/samples_bias_training.idt_rhampseq_sid_SNPsQC.vcf.gz`. The VCF filename prefix defaults to the samples CSV basename; override with `--config output_prefix=foo` if you prefer a different name.

The sample column names in the final VCF come from each BAM's `@RG SM:` tag, not the `sample_id` in the CSV. Extract the real names once the job completes:

```bash
bcftools query -l output/bias_training/samples_bias_training.idt_rhampseq_sid_SNPsQC.vcf.gz
```

### 3. Sample-level QC

`scripts/qc_bias_samples.py` computes per-sample metrics and flags samples unsuitable for bias training:

- **No-call rate** across the marker sites. Catches BAMs from a panel version that does not cover the markers (everything will be `./.`).
- **Het rate**. Unrelated individuals at common SNPs should land around 0.3-0.5. Far outside that range suggests a problem.
- **Mean `|VAF - 0.5|` at het calls**. Catches samples with heavy CNV or LOH that skew allele balance across many sites.

```bash
python scripts/qc_bias_samples.py \
  output/bias_training/samples_bias_training.idt_rhampseq_sid_SNPsQC.vcf.gz \
  --output-samples output/bias_training/pass_samples.txt \
  --output-metrics output/bias_training/qc_metrics.tsv
```

Review `qc_metrics.tsv` before trusting the pass list. Defaults are lenient starting points:

| Option | Default | Description |
|---|---|---|
| `--min-dp` | 100 | Minimum AD-sum depth for a call to count |
| `--max-nocall-rate` | 0.10 | Exclude samples with no-call rate above this |
| `--min-het-rate` | 0.15 | Exclude samples with het rate below this |
| `--max-het-rate` | 0.60 | Exclude samples with het rate above this |
| `--max-mean-vaf-dev` | 0.15 | Exclude samples where mean \|VAF-0.5\| at hets exceeds this |
| `--min-hets-for-vaf-dev` | 10 | Minimum hets required to apply the VAF-deviation check |

The `--max-mean-vaf-dev` check includes panel bias itself in the signal (every sample is pushed the same way at biased sites), so if real panel bias is large a blanket threshold may over-reject. If that happens, an iterative approach works: run `estimate-bias` on all samples first to get rough per-site biases, subtract them, then flag samples with residual skew.

### 4. Run estimate-bias

Feed the pass list into the joint-VCF mode of `estimate-bias`:

```bash
allomix estimate-bias \
  --vcf output/bias_training/samples_bias_training.idt_rhampseq_sid_SNPsQC.vcf.gz \
  --samples $(cat output/bias_training/pass_samples.txt) \
  --output output/bias_training/bias_table.tsv \
  --min-het 30
```

With 100+ samples passing QC, raise `--min-het` to around 30: every site should accumulate plenty of hets, and a stricter threshold drops any locus with unexpectedly thin coverage.

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
