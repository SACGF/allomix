# Changelog

## [Unreleased]

### Added

- **Contig classification module** `allomix.contigs` (#50). `classify_contig(chrom, pos)`
  labels a marker as autosome, chrX/chrY PAR or non-PAR, MT, or other (alt, decoy,
  unplaced, random, HLA). The pseudoautosomal mask is the exact union of the GRCh37
  and GRCh38 NCBI intervals, so no genome build is needed; a test cross-checks the
  constants against `bioutils.par` when the fork is installed.
- `MarkerGenotypes.n_par_excluded` and `MarkerGenotypes.n_other_contig_excluded`
  counters, and matching stderr notes from `detect` when either is non-zero.

### Changed

- **Pseudoautosomal (PAR) markers are now always excluded** from the informative
  set, whether or not `--use-sex-chroms` is given, and counted separately from the
  non-PAR sex-chromosome exclusion (#50).
- **Markers on non-primary contigs** (alt, decoy, unplaced, random) are now always
  excluded and counted (#50). Existing fixtures and public data contain none, so
  results on them are unchanged.
- `is_sex_chrom` moved from `allomix.genotype` to `allomix.contigs` (still importable
  from `genotype`). It now also returns True for non-primary contigs, consistent
  with its "not an autosome" meaning.

## [0.4.2] - 2026-07-02

### Changed

- **README:** removed the commercial-product comparison table (the sourced version
  lives in the paper).

## [0.4.1] - 2026-07-02

### Fixed

- **Relative documentation links in README** now use absolute GitHub URLs so they
  resolve on the PyPI project page.
- **Simulator reproducibility across Python versions.** Binomial allele-count draws
  now always use numpy's generator instead of `random.Random.binomialvariate`
  (which only exists on Python 3.12+), so simulated data is identical across
  interpreter versions.

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

[Unreleased]: https://github.com/SACGF/allomix/compare/v0.4.2...HEAD
[0.4.2]: https://github.com/SACGF/allomix/compare/v0.4.1...v0.4.2
[0.4.1]: https://github.com/SACGF/allomix/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/SACGF/allomix/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/SACGF/allomix/compare/v0.2.0...v0.3.0
