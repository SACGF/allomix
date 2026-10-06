# allomix architecture and code map

A reading guide for the `src/allomix/` package. The goal is to let a reviewer
walk the code in dependency order and know what each module owns before opening
it. For the upstream pipeline rationale, see `docs/joint_calling.md`.

## What the tool does

allomix estimates donor chimerism (the donor fraction of a DNA mixture) after
stem-cell transplant, from per-sample VCFs at a set of bi-allelic markers. Host
and donor genotypes come from GATK joint calling; the admixture sample's allele
depths come from a forced `bcftools mpileup` at the panel sites. allomix
classifies each marker, fits the donor fraction by maximum likelihood, and also
runs a separate detection test for whether the host is present at all.

## Data flow

```
  host VCF      donor VCF(s)     admix VCF
     |              |               |
     +----- genotype.parse_vcf -----+        cyvcf2 -> MarkerData
                    |
            genotype.classify_markers        host vs donor -> InformativeMarker
                    |                         (Vynck marker types; PAR and non-primary
                    |                          contigs dropped; non-PAR chrX routed on the
                    |                          qc.sex pair status, chrY/MT dropped)
                    v
            analysis.analyse_sample  ------------------------------+
                    |                                              |
        +-----------+-----------+                                  |
        |                       |                                  |
  chimerism.estimate_*    host_presence.host_presence_test          host_presence.donor_hom_markers
  (donor fraction MLE)    (is host present at all?)          (per-marker detail,
        |                       |                             artifact-flagged)
        +-----------+-----------+                                  |
                    |                                              |
       qc.sex_mismatch.sex_mismatch_check                          |
       (sex-mismatched pair only: independent                      |
        chrX copy-number / chrY depth estimate,                    |
        compared with the MLE, never blended)                      |
                    |                                              |
              qc.assess_quality                          scripts/ diagnostic plots
                    |
            report.to_tsv / to_json / timeline_json
```

`analysis.analyse_sample` is the single per-sample entry point. Both the CLI
(`cli._run_single_sample`) and the `scripts/` diagnostics call it, so the
classify -> estimate -> presence -> select path is defined in exactly one place.

## Package layout

The package groups modules by pipeline stage. Top-level modules are the
cross-cutting pieces (constants, marker types, the result dataclasses, the
orchestrator, the CLI, the simulator); the four subpackages each own one stage:

- `calibration/` -- per-marker tables estimated from cohort data: `bias`,
  `error_rates`, `contamination_table`.
- `estimate/` -- the MLE model: `likelihood` and the `chimerism` estimators.
- `qc/` -- quality / identity checks, detection, and metadata: `qc`,
  `host_presence`, `sample_contamination`, `relatedness`, `sex`, `sex_mismatch`,
  `runmeta`. Its
  `__init__` is empty on purpose (`qc.qc` <-> `results` would otherwise form a
  partial-initialisation cycle).
- `report/` -- output formatting (`report`) and the `html/` rendering subpackage
  (HTML reports, and the PDF renderer in `html/pdf.py`, which shares the same
  templates via WeasyPrint).

## Modules (dependency order)

| Module | Owns | Key public surface |
| --- | --- | --- |
| `contigs.py` | Contig classification as a pure function of `(chrom, pos)`: autosome, chrX/chrY split into PAR and non-PAR, MT, or `OTHER` (alt, decoy, unplaced, random, HLA). Owns the PAR mask (exact union of GRCh37 and GRCh38, so no genome build is needed) and the name-only `is_sex_chrom` predicate. No allomix imports. | `ContigClass`, `classify_contig`, `is_sex_chrom`, `PAR_REGIONS` |
| `genotype.py` | VCF parsing (cyvcf2) and marker classification (the eight marker classes are described in [Marker types](marker_types.md)). The canonical home of `MarkerKey`/`marker_key`. `ContigPolicy` decides which contig classes may enter the informative set (`SEX_AWARE` default, routing non-PAR chrX on the host/donor sex pair; `AUTOSOMES_ONLY`; `ALL_PRIMARY` for diagnostic views); `admit_contig` is the pure routing table. | `parse_vcf`, `classify_markers`, `admit_contig`, `ContigPolicy`, `MarkerData`, `InformativeMarker`, `MarkerGenotypes`, `marker_type`, `MarkerKey` |
| `sex_types.py` | The `Sex` and `PairStatus` enums, in a leaf module (no allomix imports) so `genotype` can route on the pair status while `qc.sex` infers it from `genotype.MarkerData`. Re-exported by `qc.sex`. | `Sex`, `PairStatus` |
| `estimate/chimerism.py` | The donor-fraction MLE: beta-binomial likelihood, grid + Brent (single donor) / Nelder-Mead (multi), profile-likelihood CIs. | `estimate_single_donor_bb`, `estimate_multi_donor`, `ChimerismResult`, `MultiDonorResult`, `detection_limit` |
| `calibration/bias.py` | Per-marker amplification-bias table (median het-VAF deviation), used to shift the expected REF weight in the MLE. | `estimate_biases`, `save_bias_table`, `load_bias_table` |
| `calibration/error_rates.py` | Per-site, per-direction empirical error table (panel of normals). Same key shape as `bias`. | `estimate_error_rates`, `save_error_table`, `load_error_table` |
| `qc/host_presence.py` | Host-presence detection at donor-homozygous markers, plus the read-level artifact filter. Independent of the fraction MLE. | `host_presence_test`, `donor_hom_markers`, `DonorHomMarker`, `HostPresenceResult`, `ArtifactThresholds` |
| `qc/sex.py` | Sex inference per reference sample from non-PAR chrX heterozygosity (binomial LR, female model from the sample's autosomal het rate vs a male spurious-het model), reconciliation with a declared sex, and the host/donor pair status. Run by `analyse_sample` before `classify_markers` (the pair status drives the chrX routing) and attached to the result before QC, like relatedness; a declared-vs-inferred conflict is a QC FAIL. With forced-pileup depths, chrY relative depth graded with the lab thresholds resolves an ambiguous or unavailable chrX call to male (never the reverse, never overruling a confident call). Also holds the depth-ratio helpers `sex_mismatch` uses. | `infer_sex`, `assess_sex`, `pair_status`, `parse_declared_sex`, `grade_chry_depth`, `chry_depth_ratio`, `Sex`, `PairStatus`, `SexInference`, `SexResult` |
| `qc/sex_mismatch.py` | Sex-chromosome cross-check of the donor fraction for a sex-mismatched single-donor pair (#48). Two readouts with their own n and CI: the chrX copy-number-weighted allele-fraction fit (non-PAR chrX markers where host and donor differ, male party homozygous; beta-binomial with profiled rho, grid + Brent, profile CI; the same bias and error handling as the main estimator) and, when forced-pileup depths are supplied, the experimental chrY depth ratio normalised by the male reference sample. Reports `basis` (`chrY-depth` / `chrX-cn` / none) and concordance with the MLE; never blended into `donor_pct`. Attached by `analyse_sample` as `result.sex_mismatch`; a discordant pair is a soft warning in `qc.assess_quality` (status unchanged). | `sex_mismatch_check`, `fit_chrx_fraction`, `expected_ref_fraction`, `chry_male_fraction`, `region_depths_from_vcf`, `is_concordant`, `SexMismatchResult` |
| `qc/qc.py` | Quality verdict: marker counts, beta-binomial goodness-of-fit, identity checks (relatedness, sex, swap, sex-chromosome cross-check), PASS/REVIEW/FAIL with reasons. | `assess_quality`, `QCReport` |
| `analysis.py` | The shared single-sample pipeline that ties classify -> estimate -> presence -> QC together. | `analyse_sample`, `AdmixtureSampleAnalysis` |
| `report/report.py` | Output formatting (TSV, JSON, timeline JSON) for single- and multi-donor results. | `to_tsv`, `to_json`, `timeline_json` |
| `cli.py` | Argument parsing and the `detect` / `timeline` / `estimate-bias` / `estimate-errors` commands. Thin: parses input, calls `analyse_sample`, formats output. | `main` |
| `simulate.py` | Standalone synthetic-VCF generator for in-silico validation. Dependency-light, plain-text VCF I/O, so its parser is `parse_text_vcf` (not `genotype.parse_vcf`) and it keeps its own `alt_dose`. | `blend_vcfs`, `build_joint_vcf`, `parse_text_vcf` |

## Two analysis paths

allomix answers two different questions, kept deliberately separate:

1. **How much donor?** `chimerism.estimate_single_donor_bb` /
   `estimate_multi_donor` fit the donor fraction by maximum likelihood over all
   informative markers (beta-binomial, with optional bias and per-site error
   tables). This is the headline `donor_pct`.
2. **Is the host present at all?** `host_presence.host_presence_test` is a one-sided
   detection test at the markers where every donor is homozygous and the host
   carries the donor-absent allele. That allele sits at the sequencing-error
   background in a pure-donor sample, so its pooled read counts give a p-value
   and a separate low-level host-fraction estimate, more sensitive than the MLE
   CI near full donor.

## Marker keys and tables

Every marker is keyed by `(chrom, pos, ref, alt)`. This shape is defined once as
`genotype.MarkerKey` (built by `genotype.marker_key`) and imported by `bias`,
`error_rates`, and `host_presence`, so the bias table, the error table, and the
detector all join on the same key.

## Contig classes and the PAR mask

`contigs.classify_contig(chrom, pos)` is the only place allomix reads a
coordinate against reference annotation. `classify_markers` drops markers on
`OTHER` contigs (mapping hazards) and in the pseudoautosomal regions (ambiguous
copy number) before anything else, under every `ContigPolicy`, counting them
into `MarkerGenotypes.n_other_contig_excluded` and `n_par_excluded`; the CLI
reports both on stderr when non-zero. The remaining non-PAR chrX, chrY and MT
markers are routed by `genotype.admit_contig` from the `ContigPolicy` and the
host/donor sex pair status that `qc.sex.assess_sex` computed first (it reads
the same contig classification to pick the non-PAR chrX sites it infers sex
from). Under the default `SEX_AWARE` policy:

| class | pair | result |
|---|---|---|
| non-PAR chrX | `matched_female` | used through the normal diploid path (`n_chrx_used`) |
| non-PAR chrX | `matched_male` | used where host and every donor are homozygous (`n_chrx_used`); a het call in any party drops the marker (`n_chrx_male_het_dropped`), since a het on a hemizygous chromosome is a genotyping error |
| non-PAR chrX | `mismatched`, `unknown`, or no pair status | excluded; informative ones counted in `n_informative_sex_chrom_excluded` |
| non-PAR chrY, MT | any | excluded from the estimate; informative ones counted in `n_informative_sex_chrom_excluded` |

For a sex-matched pair the parties share the chrX copy number, so the diploid
dosage arithmetic is exact (for hom/hom contrasts the ploidy cancels), and
downstream code (robust refit, per-type overdispersion, host-presence test)
sees chrX markers as ordinary markers. For a mismatched pair the parties
differ in copy number and the diploid model does not hold, which is why the
gate exists. Those excluded chrX markers are not wasted: for a mismatched
single-donor pair `qc.sex_mismatch` fits them on their own under the
copy-number-weighted expectation (and reads chrY depth when a midpoint pileup
is supplied) as an independent cross-check of the autosomal estimate.
`AUTOSOMES_ONLY` (`--contig-policy autosomes_only`) excludes every
non-autosome regardless of pair status; the diagnostic `ALL_PRIMARY` admits
everything primary. Identity QC (`relatedness`, `sample_contamination`,
`shared_het_balance`) stays autosomal under every policy. The calibration side
has a matching guard: `calibration.bias.estimate_biases` skips non-PAR chrX /
chrY het observations from male samples, and `estimate_biases_both_het`
trusts a non-PAR chrX both-het site only when every party is female.

The PAR intervals are vendored constants: the exact union of the GRCh37 and
GRCh38 NCBI values (overlapping intervals merged, so PAR1 is one interval and
each PAR2 is two), 1-based inclusive, with the per-build source values in the
module comment. Because PAR markers are always excluded, over-excluding the
positions that are PAR in only one build (the PAR1 edges, and the GRCh37 chrX
PAR2 interval, about 330 kb of non-PAR Xq28 in GRCh38 coordinates) is
acceptable and removes any need to know the genome build. The transcription is
guarded by
`tests/test_contigs.py::test_constants_match_bioutils_union`, which compares the
constants with `bioutils.par` (converting its interbase, 0-based half-open
intervals to 1-based inclusive) and is skipped unless that module is
importable. It lives on a fork branch (proposed upstream as biocommons/bioutils
PR #88), not in any released bioutils; to run the check:

```bash
pip install "git+https://github.com/davmlaw/bioutils@add-par-regions"
```

## Diagnostics in `scripts/`

The `scripts/` directory holds standalone, regenerable diagnostic and
data-generation tools (not part of the installed package). The host-presence
plots (`plot_host_presence_per_marker.py`, `host_presence_manhattan.py`) consume
`analysis.analyse_sample` and the public `host_presence.donor_hom_markers`, so the
markers and pooled lines they draw match the `detect` batch exactly, including
the sex-chromosome and artifact-filter handling. See `docs/scripts.md`.
