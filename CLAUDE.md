# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

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
- **GMP constraint**: New workflow must be conducted entirely within GMP to avoid donor/host genome mix-up during data analysis (identified as a problem in external review of current workflow)

## Files

- `email.txt` — Internal email thread establishing project requirements, data locations, and questions
- `Brochure -One Lambda Devyser Chimerism.pdf` / `One_Lambda_Devyser_Chimerism_brochure.txt` — Thermo Fisher/Devyser product using 24 indel markers, LOD 0.05%, Advyser software
- `MAR0094_Rev5_AlloSeq_HCT_sales_brochure-DIGITAL.pdf` / `AlloSeq_HCT_brochure.txt` — CareDx product using 202 bi-allelic SNPs, LOD 0.22%, HCT software
- `CareDx AlloSeq HCT Demo 18.3.26 WP notes.docx` — Notes from CareDx commercial software demo

## Commercial Products Evaluated

| Product | Vendor | Markers | Type | LOD | Software |
|---|---|---|---|---|---|
| AlloSeq HCT | CareDx | 202 SNPs (bi-allelic, 22 autosomes) | SNP | 0.22% | HCT Software (web-based) |
| Devyser Chimerism for NGS | One Lambda / Thermo Fisher | 24 indel markers (17 chromosomes) | Indel | 0.05% | Advyser (desktop) |

Both support up to 3 genomes (host + 2 donors) and provide timeline visualization of chimerism over time.

## In-House Approach (rhAmpSeq Sample ID)

- Uses IDT xGen Human ID Hybridization Capture Panel: 76 SNPs with 229 probes, discrimination power >1 in 5 million
- Panel info: https://sg.idtdna.com/pages/products/next-generation-sequencing/workflow/xgen-ngs-hybridization-capture/pre-designed-hyb-cap-panels/human-id-hyb-panel
- VCF data already generated via GATK per-sample calling from rhAmpSeq amplicon data
- Suggestion to use VariantGrid (VG) as database for donor/host genotypes matched to subsequent post-HSCT samples

## Development

### Project Structure

```
src/allomix/          # Main package
tests/                # pytest tests
docs/                 # Documentation
claude/               # Planning documents and research notes
```

### Python Version

Python >=3.10. Use modern syntax (match statements, `X | Y` unions, etc.) where appropriate.

### Dependencies

Runtime: `cyvcf2`, `numpy`, `scipy`

Dev: `pytest`, `pytest-cov`, `ruff`

Install for development:

```bash
pip install -e ".[dev]"
```

### Testing

- Use `pytest` for all tests: `pytest` from repo root
- Tests go in `tests/` mirroring the `src/allomix/` module structure (e.g. `tests/test_genotype.py` for `allomix/genotype.py`)
- Aim for unit tests on core calculation logic with known expected outputs
- Use synthetic/fixture data rather than real patient data in tests

### Linting and Formatting

- Use `ruff` for linting and formatting: `ruff check src/ tests/` and `ruff format src/ tests/`
- Config is in `pyproject.toml` — line length 100, target py310
- Lint rules: E, F, I (isort), W, UP (pyupgrade)

### Code Style

- Type hints on all public function signatures
- Docstrings on public modules and functions (Google style)
- No unnecessary abstractions — keep it simple and direct
- Prefer `cyvcf2` over `pysam` for VCF parsing (lighter, faster for read-only access)

### CLI

Entry point is `allomix.cli:main`, invoked as `allomix` on the command line. Use `argparse` for argument parsing (no extra dependencies).

### Versioning

Version is defined in both `pyproject.toml` and `src/allomix/__init__.py`. Keep them in sync.
