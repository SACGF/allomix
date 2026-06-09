# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Working Rules

- **Never hallucinate.** Do not state inferences, guesses, or recalled facts as if they were verified. If you are not certain, say so and verify first. This applies to dataset identities, file contents, study citations, anything.
- **Stop at permission/paywall walls.** If you hit a 401/403/404, a private repo, a login prompt, or any access wall, stop and ask the user to provide auth, fix access, or retrieve the content. Do not work around the wall with guesses or assumptions.
- **Never close GitHub issues.** Do not close issues yourself, and do not use auto-closing keywords (`Closes`, `Fixes`, `Resolves`, etc. followed by `#N`) in commit messages or PR descriptions. Reference issues with a plain `#N` only. Closing is the user's decision; a feature being merged does not mean the issue is done.

## Purpose

**allomix** is a general-purpose, panel-agnostic tool for NGS-based donor chimerism monitoring after hematopoietic stem cell transplantation (HSCT). It works with any set of bi-allelic markers (SNPs or indels) present in the input VCFs.

## Design Principles

- **Panel-agnostic**: The tool does not assume a specific marker panel. It operates on whatever bi-allelic loci are present in the input VCFs. Panel choice is a lab decision, not a tool decision.
- **Workflow**: Sequence host and donor individually to determine genotypes, then analyse post-HSCT admixture samples at serial timepoints (>=3 per patient)
- **Requirements**: Detect up to 2 donor genomes + host genome (3 total), sensitivity <1%, output is "% chimerism"
- **Key constraint**: Robust storage and matching of donor/host genotypes for subsequent analysis

## Our Lab's Deployment

Our specific deployment uses the 76 Sample ID SNPs from the IDT rhAmpSeq panel in the Haem capture panel (coverage >1000x). This context is useful for testing and validation but should not be hardcoded into the tool.

- **Existing data**: VCF files with genotypes for all 76 markers on TAU server at `/tau/data/clinical_hg38/idt_rhampseq_sid/` with AF tag (VAF) present, depth >>1000x
- **GMP** (Genetics and Molecular Pathology) is the department running this in-house; not a design constraint, just the deployment context.

## Input Requirements

- **VCF-first**: allomix takes VCFs as input (not BAMs). Minimum required FORMAT fields: GT, AD, DP.
- **Two-phase upstream**: Host/donor `GT` should come from GATK joint calling of those reference samples. Admix sample `AD` should come from forced `bcftools mpileup` at the panel sites, not from GATK — `HaplotypeCaller -ERC GVCF` strips minority ALT reads at hom-ref blocks, which is exactly the low-fraction signal we need. See `doc/joint_calling.md` for the full rationale (including the empirical check and why a somatic caller is also not the right answer).

## Data Access

- Patient data on `/tau` is NOT directly accessible. To examine real VCFs, write a standalone script that outputs only summary statistics (no patient identifiers or genomic coordinates) and ask the user to run it.
- De-identified example VCFs are in `data/` for development use.

## License & Attribution

- allomix is MIT licensed.
- The MLE chimerism estimation methodology is based on Crysup & Woerner (2022) — cite this paper.
- Code is independently reimplemented. Do NOT copy code from AGPL (Demixtify), non-commercial (Conpair, chimerism_smmip), or unlicensed repos. The math is published science and freely reimplementable.
- MIT-licensed repos (All-FIT, FABCASE, somalier) can be referenced/adapted.

## Validation Strategy

In silico validation comes first: synthetic chimeric VCFs with realistic noise models (per-marker bias, depth CV, locus dropout) calibrated from empirical panel data. All in silico experiments use multiple independent replicates (N>=5) with different random seeds to capture sampling variability. Wetlab validation with real patient samples and controlled dilution series is planned as the next phase.

## Background Materials

- `claude/` — Planning documents, decision records, and reference tool analysis
- Commercial product evaluations (CareDx AlloSeq HCT, Thermo Fisher Devyser) were reviewed during planning; specs are captured in `claude/step4_reference_tool_analysis.md` and the README comparison table

## Development

### Project Structure

```
src/allomix/          # Main package (genotype, chimerism, qc, report, simulate, cli)
tests/                # pytest tests
tests/test_data/      # Synthetic test VCFs (100 markers, 0-100% in 10% steps + timeline)
scripts/              # Utility scripts (test data generation, validation)
data/                 # De-identified example VCFs from real pipeline
claude/               # Planning documents, decision records, research notes
output/               # Script output (gitignored)
```

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

The paper build uses Snakemake (`Snakefile` in repo root). All 7 validation/figure scripts in `paper/scripts/` are independent and run in parallel, then `vibepaper build` renders the final document from the facts CSVs they produce in `output/facts/`.

```bash
snakemake -s paper/Snakefile --cores $(nproc)              # full parallel build
snakemake -s paper/Snakefile --cores $(nproc) --forceall   # force rerun everything
snakemake -s paper/Snakefile --cores 1 paper               # render paper only (facts must exist)
snakemake -s paper/Snakefile clean             # remove generated output
```

### Versioning

Version is defined in both `pyproject.toml` and `src/allomix/__init__.py`. Keep them in sync.
