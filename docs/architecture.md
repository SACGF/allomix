# allomix architecture and code map

A reading guide for the `src/allomix/` package. The goal is to let a reviewer
walk the code in dependency order and know what each module owns before opening
it. For the project history and the rationale behind each design decision, see
`claude/allomix_overall_plan.md`; for the upstream pipeline rationale, see
`docs/joint_calling.md`.

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
                    |                         (Vynck marker types; sex chroms dropped)
                    v
            analysis.analyse_sample  ------------------------------+
                    |                                              |
        +-----------+-----------+                                  |
        |                       |                                  |
  chimerism.estimate_*    detect.host_presence_test          detect.donor_hom_markers
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

## Modules (dependency order)

| Module | Owns | Key public surface |
| --- | --- | --- |
| `genotype.py` | VCF parsing (cyvcf2) and marker classification. The canonical home of `MarkerKey`/`marker_key`. | `parse_vcf`, `classify_markers`, `MarkerData`, `InformativeMarker`, `MarkerGenotypes`, `marker_type`, `MarkerKey` |
| `chimerism.py` | The donor-fraction MLE: beta-binomial likelihood, grid + Brent (single donor) / Nelder-Mead (multi), profile-likelihood CIs. | `estimate_single_donor_bb`, `estimate_multi_donor`, `ChimerismResult`, `MultiDonorResult`, `detection_limit` |
| `bias.py` | Per-marker amplification-bias table (median het-VAF deviation), used to shift the expected REF weight in the MLE. | `estimate_biases`, `save_bias_table`, `load_bias_table` |
| `error_rates.py` | Per-site, per-direction empirical error table (panel of normals). Same key shape as `bias`. | `estimate_error_rates`, `save_error_table`, `load_error_table` |
| `detect.py` | Host-presence detection at donor-homozygous markers, plus the read-level artifact filter. Independent of the fraction MLE. | `host_presence_test`, `donor_hom_markers`, `DonorHomMarker`, `HostPresenceResult`, `ArtifactThresholds` |
| `qc.py` | Quality verdict: marker counts, beta-binomial goodness-of-fit, PASS/REVIEW/FAIL with reasons. | `assess_quality`, `QCReport` |
| `analysis.py` | The shared single-sample pipeline that ties classify -> estimate -> presence -> QC together. | `analyse_sample`, `SampleAnalysis` |
| `report.py` | Output formatting (TSV, JSON, timeline JSON) for single- and multi-donor results. | `to_tsv`, `to_json`, `timeline_json` |
| `cli.py` | Argument parsing and the `detect` / `timeline` / `estimate-bias` / `estimate-errors` commands. Thin: parses input, calls `analyse_sample`, formats output. | `main` |
| `simulate.py` | Standalone synthetic-VCF generator for in-silico validation. Dependency-light, plain-text VCF I/O, so its parser is `parse_text_vcf` (not `genotype.parse_vcf`) and it keeps its own `alt_dose`. | `blend_vcfs`, `build_joint_vcf`, `parse_text_vcf` |

## Two analysis paths

allomix answers two different questions, kept deliberately separate:

1. **How much donor?** `chimerism.estimate_single_donor_bb` /
   `estimate_multi_donor` fit the donor fraction by maximum likelihood over all
   informative markers (beta-binomial, with optional bias and per-site error
   tables). This is the headline `donor_pct`.
2. **Is the host present at all?** `detect.host_presence_test` is a one-sided
   detection test at the markers where every donor is homozygous and the host
   carries the donor-absent allele. That allele sits at the sequencing-error
   background in a pure-donor sample, so its pooled read counts give a p-value
   and a separate low-level host-fraction estimate, more sensitive than the MLE
   CI near full donor. See `claude/20_host_presence_detection_plan.md`.

## Marker keys and tables

Every marker is keyed by `(chrom, pos, ref, alt)`. This shape is defined once as
`genotype.MarkerKey` (built by `genotype.marker_key`) and imported by `bias`,
`error_rates`, and `detect`, so the bias table, the error table, and the
detector all join on the same key.

## Diagnostics in `scripts/`

The `scripts/` directory holds standalone, regenerable diagnostic and
data-generation tools (not part of the installed package). The host-presence
plots (`plot_host_presence_per_marker.py`, `host_presence_manhattan.py`) consume
`analysis.analyse_sample` and the public `detect.donor_hom_markers`, so the
markers and pooled lines they draw match the `detect` batch exactly, including
the sex-chromosome and artifact-filter handling. See `docs/scripts.md`.
