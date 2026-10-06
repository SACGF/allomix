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
                    |                          contigs dropped; non-PAR X/Y/MT dropped
                    |                          under ContigPolicy.AUTOSOMES)
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
  `host_presence`, `sample_contamination`, `relatedness`, `sex`, `runmeta`. Its
  `__init__` is empty on purpose (`qc.qc` <-> `results` would otherwise form a
  partial-initialisation cycle).
- `report/` -- output formatting (`report`) and the `html/` rendering subpackage
  (HTML reports, and the PDF renderer in `html/pdf.py`, which shares the same
  templates via WeasyPrint).

## Modules (dependency order)

| Module | Owns | Key public surface |
| --- | --- | --- |
| `contigs.py` | Contig classification as a pure function of `(chrom, pos)`: autosome, chrX/chrY split into PAR and non-PAR, MT, or `OTHER` (alt, decoy, unplaced, random, HLA). Owns the PAR mask (exact union of GRCh37 and GRCh38, so no genome build is needed) and the name-only `is_sex_chrom` predicate. No allomix imports. | `ContigClass`, `classify_contig`, `is_sex_chrom`, `PAR_REGIONS` |
| `genotype.py` | VCF parsing (cyvcf2) and marker classification (the eight marker classes are described in [Marker types](marker_types.md)). The canonical home of `MarkerKey`/`marker_key`. `ContigPolicy` decides which contig classes may enter the informative set (`AUTOSOMES` default; `ALL_PRIMARY` for diagnostic views). | `parse_vcf`, `classify_markers`, `ContigPolicy`, `MarkerData`, `InformativeMarker`, `MarkerGenotypes`, `marker_type`, `MarkerKey` |
| `estimate/chimerism.py` | The donor-fraction MLE: beta-binomial likelihood, grid + Brent (single donor) / Nelder-Mead (multi), profile-likelihood CIs. | `estimate_single_donor_bb`, `estimate_multi_donor`, `ChimerismResult`, `MultiDonorResult`, `detection_limit` |
| `calibration/bias.py` | Per-marker amplification-bias table (median het-VAF deviation), used to shift the expected REF weight in the MLE. | `estimate_biases`, `save_bias_table`, `load_bias_table` |
| `calibration/error_rates.py` | Per-site, per-direction empirical error table (panel of normals). Same key shape as `bias`. | `estimate_error_rates`, `save_error_table`, `load_error_table` |
| `qc/host_presence.py` | Host-presence detection at donor-homozygous markers, plus the read-level artifact filter. Independent of the fraction MLE. | `host_presence_test`, `donor_hom_markers`, `DonorHomMarker`, `HostPresenceResult`, `ArtifactThresholds` |
| `qc/sex.py` | Sex inference per reference sample from non-PAR chrX heterozygosity (binomial LR, female model from the sample's autosomal het rate vs a male spurious-het model), reconciliation with a declared sex, and the host/donor pair status. Attached to the result by `analyse_sample` before QC, like relatedness; a declared-vs-inferred conflict is a QC FAIL. chrY depth is reserved for a later phase. | `infer_sex`, `assess_sex`, `pair_status`, `parse_declared_sex`, `Sex`, `PairStatus`, `SexInference`, `SexResult` |
| `qc/qc.py` | Quality verdict: marker counts, beta-binomial goodness-of-fit, identity checks (relatedness, sex, swap), PASS/REVIEW/FAIL with reasons. | `assess_quality`, `QCReport` |
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
reports both on stderr when non-zero. Non-PAR chrX, chrY and MT markers are
then dropped under the default `AUTOSOMES` policy (informative ones counted in
`n_informative_sex_chrom_excluded`) and kept under the diagnostic
`ALL_PRIMARY` policy (`n_chrx_used`). `qc.sex` reads the same classification
to pick the non-PAR chrX sites it infers sex from; the sex-aware chrX routing
of #46 will sit between these two steps.

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
importable. It lives on a fork branch, not in any released bioutils; to run the
check:

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
