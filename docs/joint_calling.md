# Genotyping Pipeline for allomix

allomix needs two distinct things from upstream:

1. Reliable genotypes (`GT`) for the **host and donor** reference samples, so we can classify each panel marker as informative or not.
2. Raw per-allele read counts (`AD`) for the **admixture samples**, so the chimerism MLE and the low-fraction host-presence detector can count REF/ALT reads at panel sites.

These two needs have different best-in-class tools. The pipeline in `pipeline/Snakefile` runs them in two phases against per-patient CSVs (all CSVs in `pipeline/sample_csvs/` are processed in one DAG by default).

## Why not GATK joint calling for everything?

The original version of this pipeline used GATK joint calling for all samples (host, donor, every admixture timepoint) on the theory that joint calling propagates donor-discovered ALT alleles into the admixture sample's `AD`, preserving the rare ALT reads that carry the low-fraction signal.

That theory is wrong. `HaplotypeCaller -ERC GVCF` is a local-reassembly caller: at hom-ref blocks it only tracks reads supporting the called (reference) allele. Minority ALT reads that fall inside a hom-ref block are not recorded in the GVCF, so `CombineGVCFs` and `GenotypeGVCFs` never see them. Joint calling propagates the **site** (so the admixture sample gets a row at every panel position) but does not propagate the **reads** (so the `AD` at that row is ref-only).

We verified this empirically on the rhAmpSeq SID panel: across ~9 million reads at admixture-sample hom-ref calls in joint-called VCFs, zero ALT reads were retained in `FORMAT/AD`. The low-fraction signal the joint-calling step was supposed to preserve had been stripped before the joint call ever happened.

No combination of GATK flags fixes this. `-ERC BP_RESOLUTION`, `GenotypeGVCFs --include-non-variant-sites`, `-A DepthPerAlleleBySample`, and `GenotypeGVCFs --force-output-intervals panel.vcf` all affect site emission or annotation but none of them recover minority ALT reads from a hom-ref block, because those reads were never written to the GVCF. Conpair, somalier, and demixtify all extract `AD` from raw pileups for exactly this reason.

## Two-phase architecture

```
per-patient CSV (sample_id, bam_filename, sample_type)
  |
  +-- HOST + DONOR rows ------------+
  |                                 v
  |          GATK HaplotypeCaller -ERC GVCF (per sample, parallel)
  |          GATK CombineGVCFs
  |          GATK GenotypeGVCFs
  |                                 |
  |                                 v
  |                      <patient>.vcf.gz   (panel sites + host/donor GTs)
  |                                 |
  |   (phase 1b, optional)          |
  |   bcftools mpileup -a AD,DP on HOST+DONOR BAMs
  |     at panel sites + amplicon midpoints
  |     |                           |
  |     v                           |
  |   allomix estimate-errors --genotype-vcf panel --homref-vcf midpoints
  |     |                           |
  |     v                           |
  |   <patient>.error_table.tsv   (per-site, both-direction background)
  |                                 |
  +-- ADMIX rows                    |
                |                   |
                v                   v
            bcftools mpileup -a FORMAT/AD,FORMAT/DP -R panel
            bcftools call  -m -C alleles -T panel.targets
                |
                v
            bcftools merge
                |
                v
            <patient>.admix.vcf.gz   (raw AD at every panel site)
```

Phase 1 (GATK) is used only for what it is good at: producing high-confidence germline genotypes. `AD` from phase 1 is never read for chimerism work, so the AD-stripping behaviour is irrelevant there.

Phase 2 (bcftools mpileup) handles every admixture timepoint. `bcftools mpileup -a FORMAT/AD,FORMAT/DP` writes raw REF/ALT base counts directly from the pileup, with no local reassembly to filter minority reads. `bcftools call -m -C alleles -T panel.targets.tsv.gz` then constrains genotyping to the phase-1 panel's REF/ALT pair at every panel position, so the output VCF has a `GT` + `AD` row for every panel site in every admix sample regardless of whether the ALT was observed.

Phase 1b (optional, `build_error_table: true`, on by default) builds a per-patient sequencing-error background for the host-presence detector. The detector asks whether the donor-absent allele appears above the per-site error rate, so at low host fraction the answer is dominated by how clean that background is: the flat `--error-rate` default (0.01 -> 0.33% per direction) sits above the real per-site rate (~0.05%) and over-attributes genuine host reads to error, underestimating the fraction at the bottom of a dilution series. Phase 1b measures the background directly.

It cannot use GATK for the same reason phase 2 does not: the GVCF reassembly strips minority ALT reads at hom-ref blocks, which is the exact ref->alt signal an error rate needs. So phase 1b runs `bcftools mpileup -a FORMAT/AD,FORMAT/DP` on the HOST and DONOR BAMs (the pure, uncontaminated references) at two sets of positions:

- the panel sites (same constrained `call -m -A -C alleles -T panel.targets` as phase 2): a reference sample that is hom-ref there gives the ref->alt rate, one that is hom-alt gives alt->ref.
- one force-called background position per amplicon, the BED-interval midpoint (`call -m -A`, no panel allele): almost always hom-ref in every individual, so the few stray ALT reads there are clean ref->alt error. The variant-only joint call emits no all-hom-ref sites, so without these the ref->alt direction would be undersampled.

`allomix estimate-errors --genotype-vcf <panel> --homref-vcf <midpoints> --samples HOST DONOR...` pools both into one per-site, both-direction table (`<patient>.error_table.tsv`). A midpoint that is actually polymorphic in some individual is dropped by the estimator's `--max-vaf-homref` guard, so no external allele-frequency resource is needed. Pass the table to `allomix monitor --error-table`. Phase 1b needs the HOST/DONOR BAMs only, so it runs for every patient regardless of whether admix timepoints exist yet, and it needs `allomix` on `PATH` (set `allomix:` in config to point at it).

The two phases live in one Snakefile and share one DAG. Snakemake skips phase-1 work that already exists when only new admix timepoints are added.

## Why force a genotype at every panel site

When a `panel_alleles_vcf` is configured, phase 1 runs `GenotypeGVCFs --force-output-intervals panel.vcf` so the output VCF has a host/donor genotype at every panel position, not just the positions GATK called as variant.

Note this is `--force-output-intervals`, not `--alleles`. GenotypeGVCFs has no `--alleles` option (that argument lives on HaplotypeCaller, for the old GENOTYPE_GIVEN_ALLELES behaviour). `--force-output-intervals` is the GenotypeGVCFs-native mechanism for "emit a genotype at these given sites even if non-variant in the samples". It is mutually exclusive with `--include-non-variant-sites` (the broader all-sites variant); since force-output already genotypes the panel sites whether or not they are variant, we use it on its own.

The reason for forcing is not the hom-ref/hom-ref sites it adds. Those sites are uninformative for chimerism (no allele distinguishes host from donor) and the estimator masks them anyway. GATK joint calling already emits every informative site without forcing: a site where host and donor differ has at least one non-ref allele, so it is variant in the joint call and both samples get genotyped there regardless.

The payoff is the marginal informative marker that does not clear the calling QUAL threshold (`stand-call-conf`, default 30) in a small two-sample joint call. With only host + donor in the call, a real het with modest evidence can fall below QUAL 30 and drop out of the VCF entirely. Forced output pins its genotype back. A large pooled joint call would clear the threshold on the strength of many samples; forcing recovers the same markers without pooling unrelated patients (see "Why one CSV per patient?").

At the >1000x depth of our rhAmpSeq SID deployment this rarely bites, since a true het clears QUAL 30 easily, so for us the force-output is mostly an organizational convenience: a constant set of panel rows per patient (with the uninformative hom-ref ones as harmless filler) sitting alongside whatever else GATK discovers in the wider `intervals` BED. It matters more for lower-depth panels. Either way it is the cheap fix (it only changes the `genotype_gvcfs` step), so we keep it on.

The panel VCF and the `intervals` BED are additive: HaplotypeCaller and GenotypeGVCFs use `-L <intervals.bed>` for discovery scope, and `--force-output-intervals <panel.vcf.gz>` pins the panel sites in addition. The total per-patient site count therefore varies (wide-BED discoveries change between patients) while the panel sites are constant.

To confirm every panel site came through, intersect the patient VCF against the panel and count:

```bash
bcftools isec -n=2 -w1 output/genotypes/<patient>.<bed_tag>.vcf.gz panel.vcf.gz | grep -vc '^#'
```

A short count points to a genuine coverage gap at a panel site (no reference block in the combined GVCF for `--force-output-intervals` to act on), which is a QC finding worth surfacing rather than hiding behind a forced `DP=0` row.

## Why not a somatic variant caller?

Somatic callers (Mutect2 and similar) are designed to *decide* whether a low-fraction event is a real mutation versus an artifact, and are tuned to reject low-fraction events as noise. The sub-1% donor reads we want to measure are exactly what those filters discard.

allomix is not discovering variants. It is quantifying a mixture of two known germline genotypes at a fixed panel. The caller's only job is to emit honest two-element `AD` values at known sites. That is what raw pileup gives us; that is not what a somatic caller is built for.

The sensitivity below 1% comes from the statistical model aggregating evidence across all markers (see `paper/methods.md`), not from a per-site detection threshold.

## Why one CSV per patient?

Joint calling across patients adds samples to the VCF with no informative markers for each other (no shared ancestry between unrelated patient pairs), inflating run time and output size for no benefit. One CSV per patient also keeps the host/donor pairing explicit in the filename, lets each patient run independently, and keeps the phase-1 panel small (just that patient's two reference samples).

## Input format

A CSV per patient with three columns:

```csv
sample_id,bam_filename,sample_type
HOST_001,/path/to/host.hg38.bam,HOST
DONOR_001,/path/to/donor.hg38.bam,DONOR
TP1_20240101,/path/to/timepoint1.hg38.bam,ADMIX
TP2_20240201,/path/to/timepoint2.hg38.bam,ADMIX
```

`sample_type` must be `HOST`, `DONOR`, or `ADMIX` (case-insensitive). The pipeline errors clearly if the column is missing. Multiple `DONOR` rows are allowed (two-donor transplant) and are joint-called together in phase 1.

If a patient has only pre-transplant baseline samples (no ADMIX rows yet), the pipeline runs phase 1 only and produces just the panel VCF. Adding ADMIX rows later triggers only phase-2 work.

## Prerequisites

- **GATK 4.x** on `$PATH` (tested with 4.4+)
- **bcftools 1.x** on `$PATH` (tested with 1.18+; needs `mpileup` and `call -C alleles`)
- **Snakemake** (`pip install snakemake`)
- Indexed reference genome (`.fa` + `.fai` + `.dict`)
- Indexed BAM files (`.bam` + `.bai`)

## Configuration

Config is split in two (issue #30):

- `pipeline/tools.yaml` — per-machine tool paths (`gatk`, `bcftools`, `samtools`, `bwa`, `tabix`, `bgzip`, `allomix`) and resources (`java_mem`, `align_threads`). Set up once per machine. The Snakefile loads it automatically, so you rarely touch it after the first run. Verify a new machine with `snakemake -s pipeline/Snakefile validate_tools`.
- `pipeline/config.yaml` — per-run settings. Edit this (or pass values on the command line) for each cohort. Key options:

```yaml
ref: "/path/to/hg38.fa"
samples_csv: "patient_4MO.csv"
output_dir: "output/genotypes"
intervals: "/path/to/targets.bed"   # REQUIRED capture-panel BED

# Phase 2 filters (defaults shown)
max_depth: 100000      # set well above panel coverage so >1000x is not capped
min_mapq: 20
min_baseq: 20
read_filters: "UNMAP,SECONDARY,QCFAIL,DUP"
```

The phase-2 filters approximate sensible GATK-equivalent defaults. They are tunable per patient if you need stricter or looser thresholds. A `--configfile` you pass (and any `--config key=val`) overrides both default files.

## Running

```bash
# All patients in pipeline/sample_csvs/ in one DAG
snakemake -s pipeline/Snakefile --configfile pipeline/config.yaml --cores 16

# Single patient (override the directory with one CSV path)
snakemake -s pipeline/Snakefile \
    --config ref=/path/to/hg38.fa samples_csv=patient_4MO.csv \
    --cores 8

# Point at a different directory of patient CSVs
snakemake -s pipeline/Snakefile \
    --config samples_csv_dir=/path/to/csvs \
    --cores 16

# Dry run (show DAG without executing)
snakemake -s pipeline/Snakefile --configfile pipeline/config.yaml -n

# Clean all output
snakemake -s pipeline/Snakefile --configfile pipeline/config.yaml clean
```

Phase-1 HaplotypeCaller and phase-2 pileup are both per-sample and parallelise well across patients in a single Snakemake invocation.

## Output

| File | Description |
|---|---|
| `output/genotypes/<patient>.vcf.gz` | Phase 1: GATK joint-called VCF for HOST + DONOR. Source of host/donor `GT`. |
| `output/genotypes/<patient>.admix.vcf.gz` | Phase 2: multi-sample admix VCF with raw pileup `AD` at every panel site. Source of admix `AD`. |
| `output/genotypes/gvcfs/*.g.vcf.gz` | Phase 1 per-sample GVCFs (intermediate, shared across patients by sample ID) |
| `output/genotypes/<patient>/combined.g.vcf.gz` | Phase 1 per-patient combined GVCF (intermediate) |
| `output/genotypes/<patient>/admix/per_sample/*.vcf.gz` | Phase 2 per-admix-sample VCFs (intermediate) |
| `output/genotypes/<patient>/admix/targets.tsv.gz` | Phase 1 panel sites in `bcftools call -C alleles` format |
| `output/genotypes/logs/` | Per-rule log files |

allomix then reads `<patient>.vcf.gz` for the host/donor genotypes (`--genotype-vcf`) and `<patient>.admix.vcf.gz` for the admix allele depths (`--admix-vcf`).

## Adding new timepoints

1. Add the new BAM as an `ADMIX` row in the patient CSV
2. Re-run the pipeline. Snakemake skips phase 1 (already done) and only runs `pileup_admix` for the new sample plus `merge_admix`
3. Re-run allomix against the updated admix VCF

## Restricting to specific loci

For targeted panels the BAMs typically only contain reads at panel loci, so no intervals file is needed. To restrict phase 1 (GATK) explicitly:

```yaml
intervals: "/path/to/targets.bed"
```

or pass `--config intervals=/path/to/targets.bed` on the command line. Phase 2 is automatically restricted to phase-1 sites by the constrained-allele call step, so no separate phase-2 intervals are required.
