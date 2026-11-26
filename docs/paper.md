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

The paper dependencies (snakemake >=8) require Python 3.11+, the same floor
as the core tool. Pin the venv to Python 3.13. Do not use
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

## Wording that depends on facts

The paper needs vibepaper with claims sidecars (pinned to a GitHub commit in
`pyproject.toml` until a PyPI release has them). Every number in the paper is a
`{{ namespace.field }}` fact, and so is every word that makes a claim about one.
Direction, magnitude, approximation and a/an come from vibepaper word filters
(`approx`, `number_word`, `article`, `fold`, ...).

The prose holds no checks: it is plain Markdown plus `{{ facts }}`, so co-authors
can edit it freely. Each section's checks live beside it in a claims sidecar
(`paper/results.md` → `paper/results.claims.toml`):

```toml
[[claim]]                      # an assert tied to the wording it checks
quote   = "All ... of ... two-person mixtures gate in"
assert  = "srp434573.correction_n_gated == srp434573.correction_n_mixtures"
message = "'All ... of ... two-person mixtures gate in'"

[[ok]]                         # a hedge, definition or cited figure
quote  = "most plausibly index hopping"
reason = "hedge: the mechanism is attributed by the source study, not measured here"

[[check]]                      # an invariant with no wording
assert  = "qc_gating.n_review_clinical == srp434573.n_review"
message = "the gating comparison's REVIEW count is the one quoted earlier"
```

- **A `quote` matches one sentence** of the section (or one table row), ignoring
  case, `*` and line wrapping; `…` or `...` stands for a `{{ fact }}` or any text.
  Quote the wording the check is about, at least three words. A `[[claim]]` covers
  the claim words inside its quote; `quote` can be a list when a claim rests on
  several sentences. The build fails if an assert is false and prints the fact
  values with the sidecar line and the matched sentence.
- **When you edit a sentence a quote matches**, `vibepaper lint` reports the quote
  as no longer matching. Re-check the claim against the facts, then update the quote.
  `vibepaper claims show` prints the prose with each section's checks under the
  sentences they match, for review.
- **Ordering and totality** ("all", "none", "rises", "lower") need a `[[claim]]`
  whose assert is the condition the wording depends on.
- **Judgement words** ("close to", "near", "comparable to", "essentially unchanged",
  "unchanged") have one tolerance each, defined under `[claims]` in `paper.toml`
  with the reason for it, and are used as tests: `assert = "a is near(b)"`.
  Do not write a new hand-made tolerance into an assert; if a word needs a different
  tolerance, it is a different word.
- **Hedges and non-claim uses** of a claim word ("most plausibly", "at least one",
  "estimated versus known" as an axis label) get an `[[ok]]` quoting the phrase, with
  a reason. lint reports an `[[ok]]` that no longer accepts anything.
- **Hand-typed numbers.** A result belongs in a fact. Design values (titration
  levels, depths, the nominal 95%, p < 0.05) are listed under `[lint]
  allowed_numbers`; a one-off constant (the heterozygous allele fraction 0.5, the
  chi-squared critical value) or a figure quoted from a cited study gets an `[[ok]]`
  quoting it.
- **Citations after a full stop** (`text.[@key] Next`) stop lint splitting the
  sentences there, so a quote may need to be longer to match only one. Prefer one
  short `[[ok]]` per cited figure over one long quote.
- A derived value a filter or assert needs (a ratio, a count of exceptions, a
  min/max across conditions) is emitted by the analysis script as a fact. It is
  never computed in the template for display.

### Typed facts and provenance

Producers write facts through `paper/scripts/facts_io.py` (`write_fact_csv`), which
writes a vertical `field,value,type` CSV with vibepaper's `write_facts` and a
`<name>.provenance.json` sidecar with the allomix version, date and command. Types
come from the naming convention: `*_pct` is a `percent` (already multiplied by
100), `n_*` / `*_n` a `count`, `*_rate` / `detect_rate_*` / `*_detect_0p<N>` a
`fraction`, `*_ratio` / `*_fold` / `fold_*` a `ratio`; everything else is untyped.
Follow the convention when you add a fact, and format a typed percent with
`{{ x | pct(N) }}` rather than `{{ x | dp(N) }}%`. Scripts that read a facts CSV use
`read_fact_csv`. Multi-row data tables keep their own `csv` writers, with
`lineterminator="\n"`.

The Methods cite the allomix version from `srp434573._prov`, and `vibepaper lint`
flags a paragraph that mixes facts produced by different allomix versions.

To retype a facts directory produced before this convention without rerunning the
slow simulations, run `python paper/scripts/convert_facts_to_typed.py` from the
repository root (or pass another facts directory). It converts one-row CSVs in
place, rewrites CRLF data tables as LF, and writes no provenance sidecars (the
cheap summary rules, which list `facts_io.py` as an input, rerun and write theirs).

### Frozen facts

`paper/frozen_facts.csv` holds values that are deliberately not measured by this
pipeline, each with its source and type: the commercial-kit LoDs (`percent`) and
marker counts (`count`) quoted in Table 1 and the Discussion. vibepaper applies it on
top of the facts CSVs. Add a cited constant with
`vibepaper freeze ns.field --value V --type percent --reason "..."`, format it like a
measured fact (`{{ tool_landscape.devyser_lod | pct(2) }}`), and commit the file.

### Checks

```bash
vibepaper lint --numbers --check   # claim words, hand-typed numbers, fact types (CI runs this)
vibepaper lint --unused            # facts the paper never quotes; skip with [lint] unused_ok
vibepaper check --check            # every sidecar assert against the facts that exist (CI runs this)
vibepaper claims show              # the prose with each sentence's checks under it
snakemake -s paper/Snakefile --cores 1 paper_check   # lint + check, as one rule
```

`vibepaper check` reports an assert as MISSING when its facts are not produced yet,
so it works on a partial facts directory. A `quick=1` build cannot render the paper,
because the low-iteration LoD grid fails the asserts that pin the analytical floor;
use `paper_check` there to see which.

After re-running any analysis, run `vibepaper diff` before re-reading the paper. It
names the facts that changed and the unchecked sentences that use them, with a word
diff of each changed paragraph. Facts are not committed, so to compare two runs keep
a copy of the old `output/facts` and use `vibepaper diff --against-dir <old facts>`;
paste `--format md` output into a re-run PR.

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
