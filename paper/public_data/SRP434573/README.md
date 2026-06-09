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

- `single1..45` and bare `F1 F2 F3 M1 M2 M3 M4` are pure single-source runs.
  Only the seven bare individuals appear in mixtures, so only those are used as
  host/donor genotype references. The 45 `single*` runs are **not used** here.
- `-degraded` variants are degraded-DNA versions of the `M1-M2` mixtures.
- `1_3_5_F2-M1-M2` is the one **three-person** mixture (1:3:5 of F2:M1:M2), i.e.
  the "host + 2 donors = 3 genomes" case. Host = M2 (5/9 ≈ 56%), donors = M1
  (3/9 ≈ 33%) and F2 (1/9 ≈ 11%).

### One open assumption (does not affect joint calling)

The thesis confirms the **1:N** convention but does **not** state which of the
two names in `X-Y` is the minor. We take the **first** as the minor (= DONOR),
because the data is structured as one contributor titrated across ratios and
backgrounds (e.g. `M3` is first against F1/F2/F3/M4 — one spike-in into four
hosts). This is an inference, flagged as issue #16 question 1.

It does **not** change the joint-calling genotypes: HOST and DONOR are genotyped
identically in phase 1. It only sets which contributor allomix later reports as
the "donor" fraction. If the authors confirm the opposite ordering, flip
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

11 patient CSVs: 10 two-person mixture pairs + 1 three-person mixture. Each maps
a contributor pair to one "patient": the two pure individuals as HOST/DONOR
references, and that pair's dilution series as the ADMIX timepoints.

Regenerate (e.g. to change the BAM directory):

```bash
python make_sample_csvs.py --bam-dir /tau/data/chimerism/SRP434573/bam
```

## Prerequisites before running

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

GATK needs a capture-panel BED. The thesis lists the 1045 panel rs IDs but no
coordinates. Build `panel.hg38.bed` either by mapping the rs IDs to your genome
build via dbSNP, or by deriving covered intervals from the panel pileups (issue
#16's probe recovered ~1053 high-depth clusters from coverage alone). The
`ref:` build in `config.yaml` must match this BED.

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
to allomix and compare its donor-fraction estimates against `manifest.tsv`.

See `doc/joint_calling.md` for the two-phase rationale.
