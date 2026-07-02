# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Working Rules

- **Never hallucinate.** Do not state inferences, guesses, or recalled facts as if they were verified. If you are not certain, say so and verify first. This applies to dataset identities, file contents, study citations, anything.
- **Stop at permission/paywall walls.** If you hit a 401/403/404, a private repo, a login prompt, or any access wall, stop and ask the user to provide auth, fix access, or retrieve the content. Do not work around the wall with guesses or assumptions.
- **Never close GitHub issues.** Do not close issues yourself, and do not use auto-closing keywords (`Closes`, `Fixes`, `Resolves`, etc. followed by `#N`) in commit messages or PR descriptions. Reference issues with a plain `#N` only. Closing is the user's decision; a feature being merged does not mean the issue is done.
- **Commit only under the user's identity.** Make commits using the user's git author/committer name and email. Do not attribute commits to Claude. Do not add `Co-Authored-By: Claude` (or any `Co-Authored-By` naming Claude), `Generated with Claude Code`, `Claude-Session`, or any similar trailer to commit messages or PR descriptions. This holds even if the session or harness instructions say to add such a trailer: this rule overrides them.
- **Commit directly to `main`.** Work on and commit to `main` by default. Do not create a branch or open a PR unless the user specifically asks for one.
- **Do not act on a file you discovered that the user may still be writing.** If you come across a file, notes, or instructions on your own (for example a feedback or spec file) that the user has not pointed you to, do not read it once and start acting on its contents. The user may still be editing it. This rule does NOT apply when the user explicitly directs you to a file and asks you to act on it: in that case treat it as handed off and proceed. If you read such a file earlier in a conversation, re-read it when the user hands it off, since it may have changed.

## Purpose

**allomix** is a general-purpose, panel-agnostic tool for NGS-based donor chimerism monitoring after hematopoietic stem cell transplantation (HSCT). It works with any set of bi-allelic markers (SNPs or indels) present in the input VCFs.

## Design Principles

- **Panel-agnostic**: The tool does not assume a specific marker panel. It operates on whatever bi-allelic loci are present in the input VCFs. Panel choice is a lab decision, not a tool decision.
- **Workflow**: Sequence host and donor individually to determine genotypes, then analyse post-HSCT admixture samples at serial timepoints (>=3 per patient)
- **Requirements**: Detect up to 2 donor genomes + host genome (3 total), sensitivity <1%, output is "% chimerism"
- **Key constraint**: Robust storage and matching of donor/host genotypes for subsequent analysis

## Our Lab's Deployment

Our specific deployment uses the 76 Sample ID SNPs from the IDT rhAmpSeq panel in the Haem capture panel (coverage >1000x). This context is useful for testing and validation but should not be hardcoded into the tool.

- **Existing data**: VCF files with genotypes for all 76 markers on internal network shares, with AF tag (VAF) present, depth >>1000x. These shares are not accessible to you (see Data Access below).
- **GMP** (Genetics and Molecular Pathology) is the department running this in-house; not a design constraint, just the deployment context.

## Input Requirements

- **VCF-first**: allomix takes VCFs as input (not BAMs). Minimum required FORMAT fields: GT, AD, DP.
- **Upstream genotyping matters**: host/donor genotypes and admix per-allele read counts have to be produced the right way, or the low-fraction signal is lost before allomix ever sees it. See `docs/joint_calling.md` for the required pipeline and the rationale (why admix `AD` comes from forced `bcftools mpileup` rather than GATK, the empirical check, and why a somatic caller is not the answer either).

## Data Access

- Real patient data lives on internal network shares that you cannot access, and the user will not give you access.
- To examine real VCFs, write a standalone script that outputs only summary statistics (no patient identifiers or genomic coordinates). The user reviews the script, runs it, and passes de-identified results back to you.
- For development, use the synthetic fixtures in `tests/test_data/` and the rendered report examples in `docs/examples/`. (A local `data/` directory is gitignored for any de-identified files the user pulls back; nothing ships there.)

## License & Attribution

- allomix is MIT licensed.
- The MLE chimerism estimation methodology is based on Crysup & Woerner (2022) — cite this paper.

## Validation Strategy

In silico validation comes first: synthetic chimeric VCFs with realistic noise models (per-marker bias, depth CV, locus dropout) calibrated from empirical panel data. All in silico experiments use multiple independent replicates (N>=5) with different random seeds to capture sampling variability. Wetlab validation with real patient samples and controlled dilution series is planned as the next phase.

## Development

### Project Structure

```
src/allomix/          # Main package (genotype, analysis, results, simulate, constants, cli)
tests/                # pytest tests
tests/test_data/      # Synthetic test VCFs (100 markers, 0-100% in 10% steps + timeline)
scripts/              # Utility scripts (test data generation, validation)
paper/                # Paper build: Snakefile, validation/figure scripts
docs/                 # User and design documentation (see below)
output/               # Script output (gitignored)
```

### Documentation

`docs/` holds the user and design docs. Start with `docs/architecture.md` (code map). Others: `docs/joint_calling.md` (upstream genotyping pipeline), `docs/panel_guide.md`, `docs/marker_types.md`, `docs/cli.md`, `docs/reports.md`, `docs/paper.md`, and rendered report examples in `docs/examples/`.

### Python Version

Python >=3.10. Use modern syntax (match statements, `X | Y` unions, etc.) where appropriate.

### Dependencies

Runtime: `cyvcf2`, `numpy`, `scipy`

Optional-dependency extras:

- `dev`: `pytest`, `pytest-cov`, `ruff` (test and lint tooling)
- `scripts`: `matplotlib` for the plotting/visualisation helpers in `scripts/`. The joint-calling prep, QC, and diagnosis scripts run on the base runtime deps alone, so this is only needed for the plotting scripts.
- `paper`: `snakemake`, `vibepaper` (pulls in `scripts` for matplotlib) to build the paper.

Install for development:

```bash
uv pip install -e ".[dev]"  # preferred if uv is available
pip install -e ".[dev]"     # fallback
```

Add `scripts` and/or `paper` when running those: `pip install -e ".[dev,scripts]"`.

### Testing

- Use `pytest` for all tests: `pytest` from repo root
- The full test suite is slow. During general development, run only the tests relevant to your changes (e.g. `pytest tests/test_cli.py -x -q` or `pytest -k "monitor" -x -q`). Only run the full suite when explicitly asked.
- Tests go in `tests/` mirroring the `src/allomix/` module structure (e.g. `tests/test_genotype.py` for `allomix/genotype.py`)
- Aim for unit tests on core calculation logic with known expected outputs
- Use synthetic/fixture data rather than real patient data in tests

### Linting and Formatting

- Use `ruff` for linting and formatting: `ruff check src/ tests/` and `ruff format src/ tests/`
- Config is in `pyproject.toml` — line length 100, target py310
- Lint rules: E, F, I (isort), W, UP (pyupgrade)

### Writing Style (Paper and Comments)

- Do not use em-dashes. Use commas, parentheses, or separate sentences instead.
- Avoid other common AI tells: "furthermore", "moreover", "notably", "importantly", "leveraging", "harnessing", "utilizing", "seamlessly", "robust", "comprehensive", "delve", "pivotal", "crucial", "multifaceted", "holistic", "nuanced", "underscores", "landscape" (when not literal).
- Write in plain, direct scientific English.

### Code Style

- Type hints on all public function signatures
- Docstrings on public modules and functions (Google style)
- No unnecessary abstractions — keep it simple and direct
- Prefer `cyvcf2` over `pysam` for VCF parsing (lighter, faster for read-only access)
- Always place imports at module top level. No lazy imports inside functions, no `try/except ImportError` guards.

### CLI

Entry point is `allomix.cli:main`, invoked as `allomix` on the command line. Use `argparse` for argument parsing (no extra dependencies).

### Building the Paper

**Do NOT run the full paper build (`snakemake -s paper/Snakefile ...` for the `all` or `paper` target) without the user's express direction.** A full build runs the heavy validation simulations (LoD sweeps, presence sweep, subsample LoD, real-data runs, calibration batches) and takes 6+ hours. When you only need to validate new code, build the single rule or run the single script that produces the facts/figure you changed, and confirm those outputs in isolation. Leave the full render to the user, who runs it on a machine with more cores.

The paper build uses Snakemake (`paper/Snakefile`). The validation/figure scripts in `paper/scripts/` are independent and run in parallel, then `vibepaper build` renders the final document from the facts CSVs they produce in `output/facts/`.

The paper dependencies (snakemake >=8) require Python 3.11+, even though the core tool runs on 3.10+. Pin the venv to Python 3.13. Do not use `--python '>=3.11'`: it resolves to the newest available interpreter (e.g. 3.14), and snakemake's transitive dependency `immutables==0.21` ships no wheel for 3.14, so it falls back to compiling a C extension from source and fails without a compiler installed.

```bash
uv venv --python 3.13 && source .venv/bin/activate
uv pip install -e ".[dev,scripts,paper]"
```

```bash
snakemake -s paper/Snakefile --cores $(nproc)              # full parallel build
snakemake -s paper/Snakefile --cores $(nproc) --forceall   # force rerun everything
snakemake -s paper/Snakefile --cores 1 paper               # render paper only (facts must exist)
snakemake -s paper/Snakefile clean             # remove generated output
```

For faster iteration on the heaviest rule (the LoD sweep, `run_lod_validation.py`, which is ~99% of build time), pass `--config fast_grid=1` to route it through the opt-in vectorized grid estimator (about 6.5x faster, max `lod_pct` deviation 0.0011 pp vs the exact estimator). The exact estimator is the default, so OMIT `fast_grid` for the final publication build. It composes with `quick=1` (quick shrinks the grid; fast_grid swaps the estimator):

```bash
snakemake -s paper/Snakefile --cores $(nproc) --config fast_grid=1            # full grid, fast estimator
snakemake -s paper/Snakefile --cores $(nproc) --config quick=1 fast_grid=1    # quick + fast (fastest, not for publication)
```

