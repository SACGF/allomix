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
