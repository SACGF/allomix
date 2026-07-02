# Changelog

## [Unreleased]

## [0.4.0] - 2026-07-02

### Added

- **PDF report output** #35
- **`panel-qc` subcommand** #37
- **Variant-caller mismatch warnings** #42
- **Per-sample uniformity and shared het balance QC** #38
- **Run command stored in reports**

### Changed

- **Renamed `monitor` subcommand to `detect`**
- **`--admix-vcf` now repeatable** across VCFs
- **Renamed CLI args** (`--panel-vcf` -> `--genotype-vcf`, `--samples` -> `--sample`)
- **QC review flag made less sensitive** #40

## [0.3.0] - 2026-06-29

### Added

- **Standalone HTML report** #27
- **Split HET/HOM marker overdispersion** #33

### Changed

- **Snakemake config split** into separate tool configuration (where GATK, samtools,
  etc. are installed, set once per machine) and per-run configuration (#30). Tool
  paths no longer need to be re-specified for every run.

[Unreleased]: https://github.com/SACGF/allomix/compare/v0.4.0...HEAD
[0.4.0]: https://github.com/SACGF/allomix/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/SACGF/allomix/compare/v0.2.0...v0.3.0
