# Building the paper

The paper build is orchestrated by Snakemake (`paper/Snakefile`). All seven
validation and figure scripts in `paper/scripts/` are independent and run in
parallel, then `vibepaper build` renders the final document from the facts CSVs
they produce in `output/facts/`.

> A full build runs the heavy validation simulations (LoD sweeps, presence
> sweep, subsample LoD, real-data runs, calibration batches) and takes several
> hours. When you only need to validate a change, run the single rule or script
> that produces the figure you touched.

## Environment

The paper dependencies (snakemake >=8) require Python 3.11+, even though the
core tool runs on 3.10+. Pin the venv to Python 3.13. Do not use
`--python '>=3.11'`: it resolves to the newest available interpreter (e.g.
3.14), and snakemake's transitive dependency `immutables==0.21` ships no wheel
for 3.14, so it falls back to compiling a C extension from source and fails
without a compiler installed.

```bash
uv venv --python 3.13 && source .venv/bin/activate
uv pip install -e ".[dev,scripts,paper]"
```

## Running the build

```bash
snakemake -s paper/Snakefile --cores $(nproc)              # full parallel build
snakemake -s paper/Snakefile --cores $(nproc) --forceall   # force rerun everything
snakemake -s paper/Snakefile --cores 1 paper               # render paper only (facts must exist)
snakemake -s paper/Snakefile clean                         # remove generated output
```

Snakemake tracks file timestamps, so editing a script or its input data reruns
only the affected rule and the downstream paper build.

## Quick build (for previewing)

```bash
snakemake -s paper/Snakefile --cores $(nproc) --config quick=1
```

Quick mode shrinks the loops in the heavy validation rules (far fewer pairs,
replicates, depths, and panel sizes). The estimates are then low-iteration and
noisy, so **every figure is stamped with a "QUICK BUILD" watermark** (see
`paper/scripts/paper_quick.py`); do not use quick-build figures for publication.
Setting the environment variable `ALLOMIX_PAPER_QUICK=1` has the same effect.

## Faster LoD sweep

The LoD sweep (`run_lod_validation.py`) is roughly 99% of build time. Pass
`--config fast_grid=1` to route it through the opt-in vectorized grid estimator
(about 6.5x faster, max `lod_pct` deviation 0.0011 pp vs the exact estimator).
The exact estimator is the default, so **omit `fast_grid` for the final
publication build**. It composes with `quick=1` (quick shrinks the grid;
fast_grid swaps the estimator):

```bash
snakemake -s paper/Snakefile --cores $(nproc) --config fast_grid=1            # full grid, fast estimator
snakemake -s paper/Snakefile --cores $(nproc) --config quick=1 fast_grid=1    # quick + fast (fastest, not for publication)
```

## Output formats and system dependencies

The build produces a Word document (`output/allomix_paper_<date>.docx`) and
rendered Markdown (`.md`). pandoc is bundled by vibepaper (via
`pypandoc-binary`), so no system pandoc is needed.

It also produces a **PDF when weasyprint's system libraries (pango, cairo,
gdk-pixbuf) are present**, and skips the PDF without failing the build when they
are not (the rule uses `vibepaper build --md --pdf-if-available`). To get the
PDF, install those libraries:

```bash
# Debian/Ubuntu (other distros: install the equivalent pango/cairo/gdk-pixbuf packages)
sudo apt-get install -y libpango-1.0-0 libpangocairo-1.0-0 libgdk-pixbuf-2.0-0 libcairo2
```

`vibepaper build --is-pdf-available` reports whether the PDF toolchain is
usable.

## Appendix: regenerating the real-data inputs (SRP434573)

Not part of a normal build. Skip this unless you are changing how the real-data
inputs themselves are produced.

The real-data figures come from the public SRP434573 titrated-mixture dataset. A
snapshot of the joint-called genotype and admix VCFs is committed under
`paper/public_data/SRP434573/genotypes`, and the Snakefile reads it directly, so
a fresh checkout builds those figures with no FASTQ download, no alignment, no
joint calling, and no access to the internal BAMs.

### Reuse vs rebuild

One rule governs all of it: **a from-scratch run always wins, a committed
snapshot is the fallback, and a missing snapshot degrades rather than fails.**
Nothing is ever silently stale, and the build stays green on a bare checkout.

Four things resolve this way. The first two choose input data, the third chooses
a model, the fourth chooses which stored results a figure summarises.

| What resolves | Preferred (from scratch) | Fallback (committed) | If neither |
|---|---|---|---|
| Genotype + admix VCFs | `output/genotypes/SRP434573` | `paper/public_data/SRP434573/genotypes` | build fails, the snapshot is required |
| Semi-synthetic sub-0.5% VCFs | `output/genotypes/SRP434573_synthetic` | `.../genotypes_synthetic` | synthetic run skipped, facts degrade to an `n_points=0` stub |
| Error table used by a run | `pooled.error_table.tsv` | that mixture's `<mix>.error_table.tsv` | flat `--error-rate` default |
| Error-table arm results (Figure S15) | `output/error_table_arms/<arm>` | `pooled` from the ordinary build output, `flat` and `per_mixture` from `.../error_table_arms` | that arm is dropped from the figure |

The resolvers are `resolve_srp434573_genotypes_dir` and
`resolve_srp434573_synthetic_dir` in `paper/scripts/srp434573_common.py`,
`resolve_error_table` in `paper/scripts/run_srp434573_allomix.py`, and `arm_dir`
in `paper/scripts/generate_error_table_arms_facts.py`.

Two consequences worth knowing:

- **Snapshots are committed only for what the build cannot otherwise produce.**
  The pooled arm is not committed, because the pipeline defaults to the pooled
  table and so the ordinary build already writes exactly that arm; committing a
  copy would duplicate an artifact and let the two drift. Only the flat and
  per-mixture arms are stored.
- **Presence of a fresh directory is the whole trigger.** A stale
  `output/genotypes/SRP434573` from an old run silently outranks a newer
  committed snapshot. Delete it if you want the snapshot used.

To confirm which source a build actually used, check the paths in the Snakemake
job's `input:` line, or the log line each script prints naming the directory it
read.

Regenerating either snapshot needs the aligned BAMs, which are not public. The
full procedure (download, alignment, panel BED recovery, joint calling, error
tables, and the semi-synthetic blending) is documented in
[`paper/public_data/SRP434573/README.md`](../paper/public_data/SRP434573/README.md).
The run itself is:

```bash
# from the repo root
snakemake -s pipeline/Snakefile \
    --configfile paper/public_data/SRP434573/config.yaml \
                 <your_data_paths.yaml> \
    --cores 16
```

Three separate kinds of config are in play, and mixing them up is the usual
first failure:

| What | Holds | Where |
|---|---|---|
| Run config | `intervals`, `samples_csv_dir`, `output_dir`, mpileup settings | `paper/public_data/SRP434573/config.yaml` |
| Data paths | `ref`, `fastq_dir`, `bam_dir` | `<your_data_paths.yaml>`, yours to write |
| Tool paths | `gatk`, `bcftools`, `samtools`, `bwa`, `tabix`, `bgzip`, resource limits | `pipeline/tools.yaml`, loaded automatically |

`paper/public_data/SRP434573/config.yaml` ships placeholders for the three data
paths, so it cannot run unedited. Supply real values either by editing that file
or, better, by keeping them in a separate `<your_data_paths.yaml>` layered on top.
A tools file will **not** do this job: it holds executables only, so passing one
here leaves `ref` at its placeholder and the run fails with a
`MissingInputException` naming `/path/to/hg38.fa`.

Note that both paths go on a **single** `--configfile`. The option takes a list,
so a repeated `--configfile` flag does not layer: the last one silently replaces
all earlier ones. Within the one flag, later values override earlier ones, so the
data-paths config goes second to override the placeholders. Check what actually
got loaded on the `Config file(s):` line of the snakemake startup banner.

Outputs land in `output/genotypes/SRP434573/`. Nothing downstream reads them
there by default, so refresh the committed snapshot to pick them up:

```bash
cp output/genotypes/SRP434573/*.SRP434573.vcf.gz* \
   output/genotypes/SRP434573/*.admix.vcf.gz* \
   output/genotypes/SRP434573/*.error_table.tsv \
   paper/public_data/SRP434573/genotypes/
```

That is 56 files: for each of the 11 mixtures a genotype VCF, an admix VCF, their
two `.tbi` indexes and an error table, plus one `pooled.error_table.tsv` spanning
all seven reference individuals. If the pooled table is missing, the run was
built without `build_pooled_error_table` and the presence test falls back to the
thinner per-mixture tables.

If you keep tool paths somewhere other than `pipeline/tools.yaml`, add that file
to the same `--configfile` list. Check the tools on a machine with:

```bash
snakemake -s pipeline/Snakefile validate_tools \
    --configfile paper/public_data/SRP434573/config.yaml
```

Two things about that invocation are easy to get wrong. The target must come
**before** `--configfile`: that option takes a list, so a target written after it
is swallowed as another config path and the run dies with a `FileNotFoundError`
naming the target. And `validate_tools` still needs a run config even though it
only checks executables, because the Snakefile discovers the sample CSVs at parse
time and aborts before any rule runs if `samples_csv_dir` is unset.
