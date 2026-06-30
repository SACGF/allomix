# Bias Estimation with allomix estimate-bias

Capture and amplification panels introduce systematic per-marker shifts in observed VAF. A marker that should read 0.5 in a heterozygous sample might consistently read 0.45 or 0.55 due to differential amplification efficiency. These biases are consistent across samples sequenced with the same panel, so they can be estimated once from a training set and then corrected during chimerism estimation.

Bias correction is optional. It targets per-marker amplification bias on panels where that bias is large; it does not address overdispersion, locus dropout, or cross-contamination, and on a well-behaved panel the uncorrected estimate is often already accurate. Validate that a table helps on your own data before relying on it.

## How the correction is applied

Each marker's bias is the median deviation of observed heterozygous VAF from 0.5, so it is measured where the expected reference weight is 0.5. allomix applies it multiplicatively, in logit space, rather than as a flat additive shift:

```
w_corrected = expit(logit(w) - logit(0.5 + bias))
```

At a heterozygous expected weight (0.5) this reproduces the measured deviation; at an informative marker whose expected VAF is near 0 or 1 (the common case at low chimerism) it is a small proportional nudge instead of a large additive jump. An earlier additive form (`w - bias`) overcorrected at those extreme VAFs and degraded the fit (issue #20).

## Two things to get right: caller consistency and marker coverage

**Caller consistency.** Per-marker bias is caller-specific: realignment, BAQ, and indel handling differ between callers, so a bias measured under one caller does not transfer to data called another way. The recommended two-phase pipeline genotypes host and donor with GATK joint calling but produces the admixture AD with `bcftools mpileup`. A bias table estimated from the GATK panel VCFs applied to mpileup admixture data makes results worse, not better. Estimate bias from data called the **same way as the admixture AD** (mpileup at the panel sites), or use the both-het mode below, which reads the admixture VCFs directly.

**Marker coverage.** A marker is only correctable where it was measured as a heterozygote. The markers that are informative for a given host/donor pair are, by construction, homozygous in both contributors (e.g. host 0/0, donor 1/1), so their bias cannot be measured from that pair at all. The bias for those markers has to come from other samples in which the marker is heterozygous. This is why bias correction needs a table built across a cohort, and cannot be estimated inline from a single host/donor pair.

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

If duplicates appear, either prepend the run folder to the sample_id to make them unique, or drop the duplicates and top up to your target count by sampling a few more BAMs from `all_recent_bams.txt`.

### 2. Joint-call at marker sites

Run `pipeline/Snakefile` with `intervals` set to the marker bed. Restricting to the bed keeps the job cheap (a handful of sites vs the whole panel) and produces hom-ref calls with AD at every site for every sample, which is what bias estimation needs.

For a one-off run, pass config on the command line:

```bash
snakemake -s pipeline/Snakefile \
  --config ref=/path/to/hg38.fa \
           samples_csv=samples_bias_training.csv \
           intervals=/tau/ngs_pipelines/hg38_reference_files/capture_kits/idt_rhampseq_sid/v1/idt_rhampseq_sid_SNPsQC.bed \
           output_dir=output/bias_training \
  --cores $(nproc)
```

For a repeatable setup (recommended when you expect to re-run on new cohorts), put the config in a small yaml file and reference it with `--configfile`. For example `pipeline/bias_training.yaml`:

```yaml
ref: "/path/to/hg38.fa"
samples_csv: "output/samples_bias_training.csv"
intervals: "/tau/ngs_pipelines/hg38_reference_files/capture_kits/idt_rhampseq_sid/v1/idt_rhampseq_sid_SNPsQC.bed"
output_dir: "output/bias_training"
```

Tool paths (for example a non-PATH `gatk: "/tau/tools/gatk-4.1.3.0/gatk"`) go in `pipeline/tools.yaml`, which the Snakefile loads automatically, not in this per-run file (issue #30).

Then run:

```bash
snakemake -s pipeline/Snakefile --configfile pipeline/bias_training.yaml --cores $(nproc)
```

Output: `output/bias_training/samples_bias_training.idt_rhampseq_sid_SNPsQC.vcf.gz`. The VCF filename prefix defaults to the samples CSV basename; override with `output_prefix: foo` in the yaml (or `--config output_prefix=foo` inline) if you prefer a different name.

The sample column names in the final VCF come from each BAM's `@RG SM:` tag, not the `sample_id` in the CSV. Extract the real names once the job completes:

```bash
bcftools query -l output/bias_training/samples_bias_training.idt_rhampseq_sid_SNPsQC.vcf.gz
```

Expect the biallelic-site count reported by the QC step to be equal to or slightly less than your panel size (e.g. 71/76). Sites where no sample in the cohort carries an ALT allele will not appear as biallelic in the joint VCF. That is expected and not a sign of a problem.

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

### Both-het mode (caller-consistent, from admixture VCFs)

When you do not have reference samples called the same way as the admixture AD (the two-phase pipeline does not produce mpileup'd pure references), estimate bias from the admixture VCFs themselves, at markers where the host and every donor are heterozygous. There the true admixture VAF is 0.5 regardless of the mixing fraction, so the admixture AD gives the per-marker bias directly, from the same caller being analysed.

```bash
allomix estimate-bias --both-het \
    --vcf patient_genotypes.vcf.gz \
    --host-sample HOST --donor-sample DONOR \
    --admix-vcfs patient.admix.vcf.gz \
    --output bias_table.tsv
```

A pair's both-het markers are non-informative for that **same** pair, so a single patient's table only helps **other** patients. Build the table across a cohort (run per patient and accumulate, or pass several `--admix-vcfs`), then apply the pooled table to patients whose informative markers it covers. Coverage grows with cohort size, since a marker that is both-het in one pair is often informative in another. The standard het-site mode (above), run on a reference cohort processed the same way as the admixture, gives broader coverage per sample and is preferred when such references exist.

### Options

| Option | Default | Description |
|---|---|---|
| `--vcfs` | | Per-sample VCFs, one per file (het-site mode; mutually exclusive with `--vcf`) |
| `--vcf` | | Joint-called multi-sample VCF (het-site mode with `--samples`, or genotype source for `--both-het`) |
| `--samples` | | Sample names to extract from `--vcf` (het-site mode) |
| `--both-het` | off | Estimate from admixture samples at host+donor both-het markers (caller-consistent cohort table) |
| `--host-sample` | | Host sample name in `--vcf` (`--both-het` mode) |
| `--donor-sample` | | Donor sample name in `--vcf` (`--both-het` mode; repeat for multi-donor) |
| `--admix-vcfs` | | Admixture VCFs supplying the both-het observations (`--both-het` mode) |
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
