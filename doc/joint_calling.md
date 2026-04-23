# Joint Calling for allomix

allomix requires a single multi-sample VCF produced by GATK joint calling. This document explains why, and provides a ready-to-use Snakemake pipeline to produce one from BAM files.

## Why joint calling?

allomix detects low-fraction donor chimerism by examining allele depths (AD fields) at markers where host and donor have different genotypes. When a patient is 98% host / 2% donor, most markers will be called homozygous-reference in the admixture sample. The question is whether the AD field at those sites carries a two-element value like `498,2` (ref, alt) or a single-element value like `500`.

If each sample is called independently, the variant caller has no reason to report an ALT allele at a site that looks homozygous-reference. The AD field will be single-element or absent, and allomix loses that marker entirely.

Joint calling solves this. When all samples (host, donor, admixture timepoints) are called together, ALT alleles discovered in the donor are propagated to every sample in the callset. The admixture sample gets a proper two-element AD even at sites called hom-ref, preserving the small ALT read counts that carry the chimerism signal.

This matters most below ~5% donor fraction. Above that, enough markers show heterozygous calls in the admixture sample that independent calling can still work, but joint calling gives better results at all fractions.

## Which variant caller?

Use GATK HaplotypeCaller in `-ERC GVCF` mode. The joint-calling workflow requires GVCFs (genomic VCFs that include reference blocks), and HaplotypeCaller is the standard tool for producing these.

## Pipeline overview

The Snakemake pipeline in `pipeline/` automates three GATK steps:

```
BAMs (per sample)
  |
  v
HaplotypeCaller -ERC GVCF     (parallel, one job per sample)
  |
  v
CombineGVCFs                   (merge all GVCFs)
  |
  v
GenotypeGVCFs                  (joint genotype -> final VCF)
  |
  v
joint_called.vcf.gz            (input to allomix)
```

The pipeline uses `CombineGVCFs` rather than `GenomicsDBImport` to merge GVCFs. `CombineGVCFs` does not require an intervals file, which keeps the pipeline simple for small targeted panels. For cohorts larger than ~100 samples, consider switching to `GenomicsDBImport` (which requires an intervals BED but scales better).

## Prerequisites

- **GATK 4.x** on `$PATH` (tested with 4.4+)
- **Snakemake** (`pip install snakemake`)
- Indexed reference genome (`.fa` + `.fai` + `.dict`)
- Indexed BAM files (`.bam` + `.bai`)

## Input format

A CSV file with two columns:

```csv
sample_id,bam_filename
HOST_001,/path/to/host.hg38.bam
DONOR_001,/path/to/donor.hg38.bam
TP1_20240101,/path/to/timepoint1.hg38.bam
TP2_20240201,/path/to/timepoint2.hg38.bam
```

All samples for a single patient (host, donor(s), and all admixture timepoints) should be listed together. The `sample_id` values will become the sample names in the output VCF.

## Configuration

Edit `pipeline/config.yaml`:

```yaml
ref: "/path/to/hg38.fa"
samples_csv: "output/test_samples.csv"
output_dir: "output/joint_call"
```

Or pass config values on the command line:

```bash
snakemake -s pipeline/Snakefile \
    --config ref=/path/to/hg38.fa samples_csv=my_samples.csv \
    --cores 8
```

## Running the pipeline

```bash
# Full run (HaplotypeCaller parallelised across available cores)
snakemake -s pipeline/Snakefile --configfile pipeline/config.yaml --cores 8

# Dry run (show what would execute without running anything)
snakemake -s pipeline/Snakefile --configfile pipeline/config.yaml --cores 8 -n

# Clean all output
snakemake -s pipeline/Snakefile --configfile pipeline/config.yaml clean
```

Per-sample GVCF calling is the slow step but parallelises well. For a small targeted panel (e.g. 76 SNPs within a larger capture panel), each HaplotypeCaller job typically finishes in a few minutes.

## Output

| File | Description |
|---|---|
| `output/joint_call/joint_called.vcf.gz` | Final multi-sample joint-called VCF, ready for `allomix monitor` or `allomix timeline` |
| `output/joint_call/gvcfs/*.g.vcf.gz` | Per-sample intermediate GVCFs |
| `output/joint_call/combined.g.vcf.gz` | Merged GVCF (intermediate) |
| `output/joint_call/logs/` | GATK log files for each step |

## Adding new timepoints

When a new post-HSCT sample arrives:

1. Add its BAM to the samples CSV
2. Re-run the pipeline. Snakemake will only call HaplotypeCaller on the new sample (existing GVCFs are cached), then re-run `CombineGVCFs` and `GenotypeGVCFs`
3. Re-run allomix on the updated `joint_called.vcf.gz`

## Restricting to specific loci

By default the pipeline calls variants across all regions covered by the BAM files. For targeted panels this is fine since the BAMs only contain reads at panel loci. If you want to restrict calling to specific positions (e.g. only the 76 IDT Sample ID SNPs within a larger panel), add an intervals BED file to the HaplotypeCaller rule:

```bash
snakemake -s pipeline/Snakefile \
    --config ref=/path/to/hg38.fa samples_csv=my_samples.csv intervals=/path/to/targets.bed \
    --cores 8
```

This is not currently implemented in the pipeline but is a straightforward addition to the HaplotypeCaller shell command (`-L {config[intervals]}`).
