# Changelog

## [Unreleased]

### Added

- **Sex-matched chrX routing** (#46). Non-PAR chrX markers now enter the
  estimate when the recipient and donor are the same sex: through the normal
  diploid path for a female pair, and at homozygous/homozygous sites only for
  a male pair, where a het call in either reference sample (a genotyping error
  on a hemizygous chromosome) drops the marker. For a sex-mismatched or
  unresolved pair chrX stays excluded, and the stderr note now says why.
  `assess_sex` runs before `classify_markers` so the pair status can drive the
  routing; `classify_markers` takes `pair_status`, and `genotype.admit_contig`
  is the pure routing table.
- `--contig-policy {sex_aware,autosomes_only}` on `detect` and `timeline`
  (default `sex_aware`); the value is recorded in the JSON `params`.
- `MarkerGenotypes.n_chrx_male_het_dropped` / `QCReport.n_chrx_male_het_dropped`,
  in the TSV (new last column) and JSON, in the HTML footer line when non-zero,
  and as a stderr note.
- **Bias-fitting sex guard.** `estimate_biases` takes `sample_sexes` and skips
  het observations at non-PAR chrX/chrY for male samples;
  `estimate_biases_both_het` takes `host_sex` / `donor_sexes` and uses a non-PAR
  chrX both-het site only when every party is female. `allomix estimate-bias`
  (and inline `--estimate-bias`) infers the sexes with the same inference
  `detect` uses and prints a one-line summary. `estimate-errors` uses only
  homozygous sites and needs no guard.
- `allomix.sex_types`, a leaf module holding the `Sex` and `PairStatus` enums
  (re-exported by `allomix.qc.sex`).
- Validation: `scripts/validate_sexchrom_routing.py` (in-silico, N seeds) and
  the `srp434573_chrx` paper rule (`paper/scripts/run_srp434573_chrx.py`,
  facts `output/facts/srp434573_chrx.csv`: the five same-sex public titrations
  with chrX on vs off). Autosomal results on the existing fixtures and public
  mixtures are unchanged.
- **Sex inference for reference samples** (`allomix.qc.sex`, #50). Recipient and
  donor sex is inferred from non-PAR chrX heterozygosity (binomial likelihood
  ratio of a female model, het rate from the sample's own autosomes, against a
  male spurious-het model), with honest `ambiguous` and `unavailable` states.
  Attached to the result as `result.sex` with a host/donor pair status
  (`matched_female` / `matched_male` / `mismatched` / `unknown`). Thresholds are
  starting values from the public SRP434573 individuals, to be finalised from
  `scripts/sex_calibration_summary.py` on the internal cohort.
- `--donor-sex` (one per `--donor-sample`, like `--donor-relationship`), and
  `--recipient-sex` is now parsed (`F`/`female`, `M`/`male`, any case) as well
  as displayed. A declared sex resolves an ambiguous or unavailable inference;
  a confident inference that contradicts the declaration is a QC **FAIL**
  (sample mix-up).
- TSV/JSON columns `host_sex`, `donor_sex`, `sex_pair`, `sex_source`,
  `donor_sex_source`, `n_chrx_used`, `n_par_excluded`, `n_other_contig_excluded`,
  `n_informative_sex_chrom_excluded`; a `sex` object in the JSON analysis
  payload with the full per-sample inference. The HTML header shows recipient
  and donor sex as "declared / inferred" and the footer shows the pair status
  and chrX marker count in place of the old included/excluded line.
- `genotype.ContigPolicy` (`SEX_AWARE`, the default; `AUTOSOMES_ONLY`; and the
  diagnostic `ALL_PRIMARY`).
- **Contig classification module** `allomix.contigs` (#50). `classify_contig(chrom, pos)`
  labels a marker as autosome, chrX/chrY PAR or non-PAR, MT, or other (alt, decoy,
  unplaced, random, HLA). The pseudoautosomal mask is the exact union of the GRCh37
  and GRCh38 NCBI intervals, so no genome build is needed; a test cross-checks the
  constants against `bioutils.par` when the fork is installed.
- `MarkerGenotypes.n_par_excluded` and `MarkerGenotypes.n_other_contig_excluded`
  counters, and matching stderr notes from `detect` when either is non-zero.

### Changed

- **`--use-sex-chroms` is retired** (#50). It is now a hidden option that exits
  with an error explaining that sex chromosomes are handled automatically from
  the inferred and declared sex (`--recipient-sex` / `--donor-sex`). It will be
  removed entirely in a later release. Its replacement is the sex-aware chrX
  routing above (#46); `--contig-policy autosomes_only` keeps chrX out.
- `classify_markers` and `analyse_sample` take `contig_policy: ContigPolicy`
  instead of `use_sex_chroms: bool`; `analyse_sample` also takes
  `declared_host_sex` and `declared_donor_sexes`.
- Version bumped to 0.5.0.
- **Pseudoautosomal (PAR) markers are now always excluded** from the informative
  set under every contig policy, and counted separately from the non-PAR
  sex-chromosome exclusion (#50).
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
