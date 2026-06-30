# Changelog

## [Unreleased]

## [0.3.0] - 2026-06-29

### Added

- **Standalone HTML report** #27
- **Split HET/HOM marker overdispersion** #33

### Changed

- **Snakemake config split** into separate tool configuration (where GATK, samtools,
  etc. are installed, set once per machine) and per-run configuration (#30). Tool
  paths no longer need to be re-specified for every run.

[Unreleased]: https://github.com/SACGF/allomix/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/SACGF/allomix/compare/v0.2.0...v0.3.0
