# SRP434573 — public mixture dataset for allomix validation

Real targeted-sequencing data of known DNA mixtures, used to validate allomix
against ground-truth fractions on real reads (not just in-silico). Sourced and
vetted in [SACGF/allomix issue #16](https://github.com/SACGF/allomix/issues/16).

## What this is

- **SRA study `SRP434573`** / BioProject `PRJNA960854`.
- Companion work: **Chu Xufeng (褚旭峰), PhD thesis, Huazhong University of
  Science & Technology (HUST), Tongji Medical College, 2024**, supervisor Xiong
  Bo. DOI `10.27157/d.cnki.ghzku.2024.000769`. There is **no journal paper** —
  cite the thesis / accession directly. Full provenance is in issue #16.
- **1062-autosomal-SNP MIP capture panel + UMI**, Illumina HiSeq 3000, reads
  pear-merged to single-end ~96 bp inserts. Raw on-panel depth ~1000–1900x.
- **7 unrelated individuals** (4 male M1–M4, 3 female F1–F3). Pairwise
  two-person mixtures plus one three-person mixture, at known ratios.

This covers the **unrelated-donor** case. HSCT donors are often HLA-matched
siblings; there is no related-donor series here (issue #16 notes CEPH-1463 as a
possible related-donor analog for later).

## Ground-truth decoding (from the thesis)

Mixture sample aliases are `1_<N>_<X>-<Y>` = a major:minor ratio of **1:N**, so
the **minor contributor fraction = 1/(1+N)**:

| ratio | 1:9 | 1:19 | 1:39 | 1:79 | 1:99 | 1:199 |
|---|---|---|---|---|---|---|
| minor % | 10 | 5 | 2.5 | 1.25 | 1 | 0.5 |

### Role mapping (minor = HOST)

The **minor (titrated) contributor is labelled HOST** and the **major
(background) contributor DONOR**. This mirrors our common clinical case: the
residual / recurring patient (host) is the small fraction detected against a
donor-dominated graft, so the titration series (minor from 10% down to 0.5%)
exercises allomix exactly as a relapse / declining-chimerism series would. The
monitored fraction (the `known_minor_pct` column of `manifest.tsv`) is therefore
the **host** fraction.

File names are stable composition labels `mix_<minor>_into_<major>` (spike-in
into background), independent of the role mapping, so they do not change if the
mapping is flipped.

- `single1..45` and bare `F1 F2 F3 M1 M2 M3 M4` are pure single-source runs.
  Only the seven bare individuals appear in mixtures, so only those are used as
  host/donor genotype references. The 45 `single*` runs are **not used** here.
- `-degraded` variants are degraded-DNA versions of the `M1-M2` mixtures.
- `1_3_5_F2-M1-M2` is the one **three-person** mixture (1:3:5 of F2:M1:M2), i.e.
  the "host + up to 2 donors = 3 genomes" case. Following the minor = HOST rule,
  the monitored minority is **host = F2 (1/9 ≈ 11%)** and the two donors are
  **M1 (3/9 ≈ 33%)** and **M2 (5/9 ≈ 56%)**.

### One open assumption (does not affect joint calling)

The thesis confirms the **1:N** convention but does **not** state which of the
two names in `X-Y` is the minor. We take the **first** as the minor (= HOST),
because the data is structured as one contributor titrated across ratios and
backgrounds (e.g. `M3` is first against F1/F2/F3/M4 — one spike-in titrated into
four backgrounds). This is an inference, flagged as issue #16 question 1.

It does **not** change the joint-calling genotypes: HOST and DONOR are genotyped
identically in phase 1. It only sets which contributor allomix later treats as
host vs donor. If the authors confirm the opposite ordering, flip
`MINOR_IS_FIRST` in `make_sample_csvs.py` and regenerate.

## Files here

| File | What |
|---|---|
| `fetch_ena_metadata.sh` | Pulls `ena_runs.tsv` + `download_fastqs.sh` from the ENA API. |
| `ena_runs.tsv` | ENA run → sample alias table (source of the mapping). |
| `download_fastqs.sh` | ENA "download all" wget script (FASTQs). |
| `make_sample_csvs.py` | Regenerates the CSVs + manifest from `ena_runs.tsv`. |
| `sample_csvs/*.csv` | One per-patient CSV per mixture (the pipeline input). |
| `manifest.tsv` | Ground truth: patient, sample, run, role, known minor %. |
| `config.yaml` | allomix joint-calling config for this dataset. |
| `genotypes/*.vcf.gz` | Committed snapshot of the joint-called genotype + admix VCFs (see below). |
| `genotypes/*.error_table.tsv` | Optional per-patient host-presence error tables (issue #23), present once built TAU-side (see "Error tables" below). |
| `SRP434573.midpoints.bed` | One background position per amplicon (BED midpoints), the force-called hom-ref sites the error tables learn the ref->alt rate from. |

11 patient CSVs: 10 two-person mixture pairs + 1 three-person mixture. Each maps
a contributor pair to one "patient": the two pure individuals as HOST/DONOR
references, and that pair's dilution series as the ADMIX timepoints.

Regenerate (e.g. to change the BAM directory):

```bash
python make_sample_csvs.py --bam-dir /tau/data/chimerism/SRP434573/bam
```

## Committed genotype snapshot (default for the paper build)

`genotypes/` holds a 2.6 MB snapshot of the joint-called outputs: per-mixture
`<mix>.SRP434573.vcf.gz` (host/donor genotypes) and `<mix>.admix.vcf.gz` (raw AD
at panel sites), with `.tbi` indexes. These are the only files the paper's
real-data section consumes, so the build reproduces from a fresh checkout with
just `allomix` and `cyvcf2` (issue #21): no FASTQ download, no alignment, no
joint calling.

The paper Snakefile reads this directory directly. `paper/scripts/run_srp434573_allomix.py`
and `paper/scripts/probe_contam_median_srp434573.py` prefer a freshly joint-called
`output/genotypes/SRP434573` if one is present (the full from-scratch run below),
and otherwise fall back to this committed snapshot. The bulky joint-calling
intermediates (gVCFs, per-mixture work dirs) are not committed; they are not
needed downstream.

The full pipeline below is optional and only needed to regenerate the snapshot
from raw reads.

## Semi-synthetic sub-0.5% mixtures (issue #5)

The real titration bottoms out at a 0.5% minor (host) fraction. To see allomix
behaviour below that, each two-person pair's two pure reference BAMs are blended
with `samtools view --subsample` (via `scripts/mix_bams.sh`) at host fractions
0.1-0.5% (5 independent subsample seeds each), then joint-called the same way.
These points are **semi-synthetic**: real reads, real panel noise, real
GATK/bcftools path, but an artificial mixing ratio. They are always labelled as
such in the figures so a synthetic fraction is never presented as a measured one.
The 0.5% synthetic point doubles as a cross-check against the real 0.5%.

The committed snapshot lives in `genotypes_synthetic/` (per-pair
`<pair>.synthetic.SRP434573.vcf.gz` genotypes + `<pair>.synthetic.admix.vcf.gz`
raw AD). The paper build consumes it through
`paper/scripts/run_srp434573_allomix.py` (which reuses each pair's real
`genotypes/<pair>.error_table.tsv`, since the host/donor individuals are
unchanged) and `paper/scripts/generate_srp434573_synthetic_facts.py`. When the
snapshot is absent (fresh checkout before generation), the synthetic run is
skipped and the facts/figure degrade to an `n_points=0` stub, so the build stays
green.

### Regenerating (TAU-side, where the BAMs are)

`samtools view --subsample` needs the aligned BAMs, which live on `/tau`. Run the
driver, then the normal pipeline over the synthetic CSVs, then copy the VCFs into
the committed snapshot:

```bash
# 1. Subsample+merge the pure reference BAMs and write the synthetic sample CSVs.
python scripts/make_semisynthetic_srp434573.py \
    --bam-dir /tau/data/chimerism/SRP434573/synthetic_bam
# (run with --dry-run first to inspect the mix_bams commands without touching BAMs)

# 2. Joint-call + pileup the synthetic mixtures with the normal pipeline.
snakemake -s pipeline/Snakefile \
    --configfile paper/public_data/SRP434573/config.yaml \
    --config samples_csv_dir=output/semisynthetic_csv \
             output_dir=output/genotypes/SRP434573_synthetic \
    --cores 16

# 3. Copy the per-pair genotype + admix VCFs into the committed snapshot.
mkdir -p paper/public_data/SRP434573/genotypes_synthetic
cp output/genotypes/SRP434573_synthetic/*.synthetic.*.vcf.gz* \
   paper/public_data/SRP434573/genotypes_synthetic/
```

The driver prints these exact follow-on commands. Convention note: `mix_bams.sh`
treats its `DONOR_BAM` argument as the minor (titrated) contributor, while
SRP434573/allomix treats the **host** as the minor monitored fraction, so the
driver passes each pair's allomix-host BAM as `mix_bams`' DONOR argument. Getting
that backwards would invert every fraction.

## Prerequisites before running (full from-scratch pipeline, optional)

### 0. (Re)generate the ENA manifest + download list

`ena_runs.tsv` and `download_fastqs.sh` are committed, but both come straight
from the ENA API and can be regenerated for reproducibility:

```bash
./fetch_ena_metadata.sh           # SRP434573 by default
python make_sample_csvs.py        # rebuild CSVs + manifest from ena_runs.tsv
```

### 1. Download

```bash
cd /tau/data/chimerism/SRP434573    # or wherever you keep it
bash /path/to/download_fastqs.sh    # ~693 MB of FASTQs
```

### 2. Alignment (FASTQ → BAM)

The pipeline aligns the FASTQs for you: `config.yaml` sets `fastq_dir` /
`bam_dir`, so Phase 0 runs `bwa mem` on each `<fastq_dir>/<run>.fastq.gz` to
`<bam_dir>/<run>.hg38.bam` (single-end, read-group SM = sample_id, sorted and
indexed). The only manual step is **bwa-indexing the reference once**:

```bash
bwa index /path/to/hg38.fa
```

The aligner is configurable like the other tools (`bwa:` / `samtools:` paths,
`align_threads:`). BAMs are **not** duplicate-marked on purpose: this is an
amplicon/MIP panel where every amplicon read shares start/end coordinates, so
dup-marking would discard almost all coverage. The `bam_dir` / `bam_suffix` in
`config.yaml` must match the CSV `bam_filename` paths (the
`make_sample_csvs.py` defaults); regenerate the CSVs with `--bam-dir` if you
move the BAMs.

UMI collapsing is **not possible** from these FASTQs — the deposited reads carry
no UMI bases (issue #16 verified this empirically), so this is raw-depth only.

### 3. Panel BED (`intervals:`)

GATK needs a capture-panel BED. The thesis publishes no coordinates, so the
panel is rebuilt straight from the aligned BAMs: this is a MIP/amplicon assay,
so every captured locus piles thousands of reads into one tight ~95 bp footprint
and the panel self-recovers as high-depth clusters (off-target background sits
far below). `SRP434573.bed` (committed here, hg38) was built with:

```bash
python paper/scripts/build_srp434573_panel_bed.py \
    --bam-glob 'output/bam/*.bam' \
    --out paper/public_data/SRP434573/SRP434573.bed
```

It recovers 1052 intervals (1025 autosomal + 27 on chrX), in agreement with the
issue #16 laptop probe (~1053) and the stated ~1062 SNPs. A position is kept
when it is covered at >=100x (MAPQ/BASEQ >=20) in at least 50 of the 64 runs,
and adjacent kept positions are merged into one interval per amplicon; every
recovered cluster is amplicon-shaped (80-100 bp). Lower `--min-samples` to
tolerate more per-sample dropout (e.g. 40 -> 1092 intervals). The 27 chrX
clusters are genuine ~95 bp amplicons captured across the runs even though the
thesis describes the panel as autosomal; pass `--bam-glob` plus a grep/filter,
or drop them downstream, if a strictly autosomal BED is wanted. The `ref:` build
in `config.yaml` must match this BED (hg38).

## Running

```bash
# from repo root
snakemake -s pipeline/Snakefile \
    --configfile paper/public_data/SRP434573/config.yaml --cores 16

# dry run (inspect the DAG without executing)
snakemake -s pipeline/Snakefile \
    --configfile paper/public_data/SRP434573/config.yaml -n
```

Outputs land in `output/genotypes/SRP434573/`: per-patient `*.vcf.gz`
(host/donor genotypes) and `*.admix.vcf.gz` (raw AD at panel sites). Feed those
to allomix and compare its estimate of the monitored minority (the **host**
fraction, `known_minor_pct` in `manifest.tsv`) against the known values.

To refresh the committed snapshot from a fresh run, copy the two VCFs (plus
`.tbi`) for each mixture into `genotypes/`:

```bash
cp output/genotypes/SRP434573/*.SRP434573.vcf.gz* \
   output/genotypes/SRP434573/*.admix.vcf.gz* \
   paper/public_data/SRP434573/genotypes/
```

### Error tables (issue #23)

The presence-test fraction reads too low at the lowest dilutions (0.5%, 1%)
because, with no per-site error table, the detector falls back to the flat
`--error-rate` default, which is larger than the real signal at those fractions
(root cause in issue #22). The committed genotype VCFs cannot supply a clean
ref->alt background: GATK joint calling emitted variant sites only, so there are
no all-hom-ref sites, and GATK strips the minority ALT reads at hom-ref blocks
anyway (`doc/joint_calling.md`).

The fix is the pipeline's phase 1b: a raw `bcftools mpileup` of the HOST/DONOR
BAMs at the panel sites plus the amplicon midpoints (`SRP434573.midpoints.bed`),
fed to `allomix estimate-errors`. This runs TAU-side (the BAMs are on `/tau`).
`build_error_table: true` is the default, so the standard pipeline run above
already produces `output/genotypes/SRP434573/<mix>.error_table.tsv`. Copy them
into the snapshot:

```bash
cp output/genotypes/SRP434573/*.error_table.tsv \
   paper/public_data/SRP434573/genotypes/
```

Once committed, `run_srp434573_allomix.py` passes each table to `allomix monitor
--error-table` automatically (it falls back to the default model when absent), so
the paper figures pick up the data-derived background on the next build. The
tables hold only aggregated per-site rates and read counts (no patient
identifiers or genotypes), consistent with the data-access rule.

The admix VCF also carries the index-hopping metadata (issue #12) in its header:
one `##allomixRunUnit` line per admix sample with a recoverable sequencing run
unit (flowcell+lane), flagged when it shares one with the host. allomix reads
these back and surfaces `run_unit` / `index_hop_risk` columns. For this dataset
no lines are written: SRA renamed the reads `<accession>.<n>` and stripped the
`@RG PU` tag, so the flowcell is unrecoverable; the (optional) metadata is simply
absent and the flag degrades to "cannot determine" rather than a false negative.
On real lab BAMs (with PU tags or Illumina read names) the flowcell resolves and
the flag is meaningful. Inspect it with `bcftools view -h <patient>.admix.vcf.gz
| grep allomixRunUnit`.

See `doc/joint_calling.md` for the two-phase rationale.
