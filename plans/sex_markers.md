# Sex-chromosome markers: #50, #46, #48 as one piece of work

Plan for the three open sex-marker issues, done together. Decisions from the
2026-10-06 review are folded in; the remaining open questions are at the
bottom.

Evidence key, as used in the issues: **[data]** = observed in the repository,
the committed public data, or the lab's pipeline code during planning,
**[derived]** = closed-form from the model, **[speculative]** = expected payoff,
not yet measured.

---

## Decisions taken (2026-10-06 review)

| topic | decision |
|---|---|
| PAR data | Vendor the PAR boundaries into allomix with provenance; cross-check test against `bioutils.par` from the fork when installed. The module is not in any released bioutils (0.6.1, Nov 2024); it lives on the fork branch `davmlaw/bioutils@add-par-regions` (interbase coordinates, flat JSON like `_data/assemblies`) and is proposed upstream as biocommons/bioutils PR #88 (opened 2026-10-06). Switch to a `bioutils>=` pin once a release carries it. |
| Genome build | Not needed. PAR mask is the union of the GRCh37 and GRCh38 intervals (see "Why no genome build" below). Build inference and `--genome-build` are dropped from the plan. |
| chrY content | The rhAmpSeq SID BED carries three chrY sex-typing amplicons (SRY, ZFY, AMELY intron 3), all non-PAR. The idt_haem kit has none. |
| Male-male chrX | Included in Phase 2 (drop het calls, route hom/hom). |
| PAR markers | Excluded everywhere, counted and reported as skipped. |
| Declared sex | Add `--donor-sex`; declared-vs-inferred conflict is a QC **FAIL** (sample mix-up). |
| `--use-sex-chroms` | Keep as a hidden option that exits with an explanatory error; note in CHANGELOG that it is removed in a later release. |
| Calibration | `scripts/sex_calibration_summary.py` (written) runs on the internal data and prints de-identified summaries. |
| #48 scope | Single donor first. |
| Alt / decoy / unplaced contigs | Classified `OTHER` and excluded from everything (mapping hazards). |
| Paper | Add the chrX-on and sex-mismatch arms to the paper pipeline; decide what to publish after seeing the output. |

## Why do them together

The three issues share one interface and one set of fixtures:

- All three need a marker's contig classified as autosome / chrX non-PAR /
  chrX PAR / chrY non-PAR / chrY PAR / chrMT / other (#50 item 1).
- All three need a per-sample sex call with an honest "ambiguous" state (#50
  item 3). #46 gates chrX inclusion on it. #48 gates its whole module on it.
  The metadata-disagreement check (#50 item 5) is the same class of signal #48
  produces.
- #46 and #48 are two consumers of the same copy-number arithmetic: the
  mixture's expected allele fraction at a locus is
  `sum_i f_i * cn_i * allele_frac_i / sum_i f_i * cn_i`. When every contributor
  has the same copy number it collapses to the autosomal `/2` model, which is
  exactly why sex-matched chrX (#46) can use the existing estimator unchanged
  and sex-mismatched chrX (#48) cannot. The simulator already implements this
  formula as `simulate.cn_weighted_vaf` for host CNV clones. [data]
- Retiring `--use-sex-chroms` (#50 item 6) only makes sense once the
  replacement behaviour (#46) exists.
- Test fixtures, simulator support, reporting columns, HTML header fields and
  the paper methods paragraph are each touched once instead of three times.

The work is still sequenced into phases that land independently, each with its
own tests, so the estimator-affecting part (#46) can be revalidated on its own
and #48 can trail without blocking anything.

## What the data says (planning findings)

### The public SRP434573 cohort separates the sexes cleanly on chrX het rate [data]

The committed genotypes in `paper/public_data/SRP434573/genotypes/` cover 7
unrelated individuals of known sex (README: M1-M4 male, F1-F3 female), 27 chrX
amplicons, hg38, GATK diploid calls everywhere (no haploid GTs on chrX).
Called non-PAR chrX sites per reference sample, DP >= 20, no GQ filter:

| sample | sex | chrX called | chrX het | het rate | autosomal het rate |
|---|---|---:|---:|---:|---:|
| F1 | F | 24 | 17 | 0.71 | 0.50 |
| F2 | F | 27 | 15 | 0.56 | 0.51 |
| F3 | F | 25 | 11 | 0.44 | 0.52 |
| M1 | M | 26 | 1 | 0.04 | 0.52 |
| M2 | M | 26 | 1 | 0.04 | 0.53 |
| M3 | M | 26 | 1 | 0.04 | 0.51 |
| M4 | M | 24 | 2 | 0.08 | 0.50 |

With a GQ >= 20 filter added (the calibration script's default), the spurious
male hets drop to zero in M1-M3 and one in M4, and the female rates are
unchanged. [data] So the GQ filter already does most of the work the
bias-fitting guard in Phase 2 exists for, and the guard is belt and braces.

### All 27 public chrX amplicons are non-PAR [data]

Positions run from 8.5 Mb to 131.3 Mb on hg38; PAR1 ends at 2,781,479 and PAR2
starts at 155,701,383. On the one panel we can validate against the PAR mask
is a no-op, so PAR handling has to be covered by synthetic fixtures.

### chrX capture efficiency is panel-specific, so depth is not a clean sex signal [data]

Median chrX depth over median autosomal depth in the admix samples:

| mixture | sexes | panel version | chrX / autosome depth |
|---|---|---|---:|
| F1 into F3 | F, F | v2 | 1.45 to 1.67 |
| M3 into F1 | M minor, F major | v2 | 1.45 to 1.64 |
| M1 into M2 | M, M | v1 | 0.52 to 0.55 |

Female diploid chrX comes out at 1.5x autosomal depth on this panel, not 1.0x,
and the two panel versions differ in chrX amplicon count. A chrX depth ratio
therefore only means something relative to a same-panel, same-sex reference.
This puts sex inference on heterozygosity (primary) with chrY depth as the
secondary signal, and it means #48's chrX depth readout is weaker than its
allele-fraction readout.

### The deployment SID BED has chrY sex-typing amplicons, and the lab already imputes sex from them [data]

`idt_rhampseq_sid_SNPsQC.bed` carries SRY (chrY:2787394), ZFY
(chrY:2979984) and AMELY intron 3 (chrY:6869956), each a 2 bp "amplicon"
interval rather than a SNP. All three are outside PAR1 on GRCh38 (SRY is 5.9 kb
past its end). The lab's amplicon pipeline (NGS-pipelines,
`combine_target_region_mean_depth_for_amplicon_and_impute_sex.py`) imputes sex
for rhAmpSeq SID runs from mean depth at five sex-typing regions (chrX_AMELX,
chrX_ZFX, SRY, ZFY, AMELY):

```
rel = max(SRY, ZFY, AMELY) / median(all non-sex-region depths)
Undetermined if dropout rate > 0.1
M  if rel > 0.5;  M* if rel > 0.15;  F* if rel > 0.05;  else F
```

allomix's chrY-depth secondary should start from the same statistic and
thresholds, so the two tools agree on the same sample. Two consequences:

- GATK emits no VCF record at an invariant chrY site, so chrY depth cannot come
  from the genotype VCF. In the allomix pipeline it is available for HOST/DONOR
  from the forced midpoint pileup (`refs/midpoints.vcf.gz`, derived from the
  `intervals` BED), but the ADMIX pileup is at panel SNP sites only. Phase 3
  needs a pipeline change to pile up the admix BAMs at the midpoints too.
- The GATK `DP` in the committed public genotype VCFs is capped at 50 (GATK's
  downsampling), so depth for calibration must come from BAMs, which is why
  the calibration script has a `samtools depth` path. [data]

### Which panel the clinical data is on [data]

Two data sources appear in the repository, and they are different assays:

- `/tau/data/clinical_hg38/idt_rhampseq_sid/`: the rhAmpSeq Sample ID
  **amplicon** assay (NGS-pipelines scheme `idt_rhampseq_sid_*`, general scheme
  `GERMLINE_AMPLICON`). Source of the 210 joint-called VCFs behind
  `paper/empirical_results/` (depth CV 0.43, bias SD 0.018) and of the sex
  imputation above.
- `/tau/data/clinical_hg38/idt_haem/`: the **hybrid-capture** haem panel with
  the SID SNP sites added as probes. Source of the bias-training cohort
  (`docs/estimate_bias.md`), the clinical validation batches (pipeline run with
  `intervals = union_sid_haem_vendor_probes.bed`, hence the
  `.union_sid_haem_vendor_probes.vcf.gz` suffix in the diagnostic scripts), and
  the "Internal QC plot for haem run" work.

`paper/supplementary.md` (Table S1 text) says the panel characterisation VCFs
came from the SID SNPs "incorporated into a custom hematology capture panel",
while `paper/empirical_results/README.md` says they came from the
`idt_rhampseq_sid` directory, i.e. the amplicon assay. One of those is wrong;
see open question 3. For this plan it matters because the sex-typing amplicons
are SID amplicon targets, and whether the haem capture gives any depth at those
chrY positions is unknown (open question 2).

### Current code shape [data]

- `genotype.is_sex_chrom` is a flat predicate over `{X, Y, M, MT}`
  (`src/allomix/genotype.py:187`). Used by `classify_markers` (drop unless
  `use_sex_chroms`), `qc/relatedness.py`, `qc/sample_contamination.py`, and
  `scripts/plot_informative_karyogram.py`.
- `--use-sex-chroms` is plumbed through `cli._run_single_sample`,
  `analysis.analyse_sample`, `classify_markers`, the HTML parameter footer, and
  three diagnostic scripts that pass `use_sex_chroms=True` to show chrX/Y in
  genomic views.
- `--recipient-sex` is free text rendered in the HTML header. No `--donor-sex`.
- Test fixtures in `tests/test_data/` are autosome-only.
- allomix 0.4.2 is on PyPI. PyPI rejects direct-URL requirements, so a
  git-pinned fork cannot appear in `dependencies` or in any extra; hence
  vendoring.

### Why no genome build (answer to review question 11)

allomix has never needed a genome build because nothing in it interprets a
coordinate: markers are joined across samples by `(chrom, pos, ref, alt)` as an
opaque key, and every statistic is per marker. The sex work is the first place
a position has to be read against reference annotation, and only for one
purpose: deciding whether a chrX or chrY site is pseudoautosomal. PAR
boundaries differ between GRCh37 and GRCh38 (X PAR1 ends at 2,699,520 vs
2,781,479; X PAR2 starts at 154,931,044 vs 155,701,383; the Y boundaries differ
more), so a build-naive mask would be wrong in the band between the two
builds' boundaries.

Since the decision is to exclude PAR markers anyway, the fix is to over-exclude:
mask the **exact union** of both builds' PAR intervals (the two PAR1 intervals
overlap and merge; the two PAR2 intervals do not overlap on either chromosome,
so they stay as two intervals each, and the gap between them is non-PAR in both
builds). A site is then excluded if it is PAR in at least one build. The cost
is the positions that are PAR in one build only: on X, 10001-60000 and
2699521-2781479 at the PAR1 boundary, and the GRCh37 PAR2 interval
154931044-155260560, which is about 330 kb of non-PAR Xq28 in GRCh38
coordinates (the GRCh38 PAR2 interval lies beyond GRCh37's chrX length, so it
costs nothing there). Neither panel we know has a marker within 5 Mb of any of
these. Outside the union, PAR status is identical in both builds, so no build
is needed. The classification is a pure function of `(chrom, pos)`, like
everything else in allomix, and the build inference, the `--genome-build` flag,
and the fixture-header problem all go away. (Phase 0 implementation note: a
min-to-max span instead of the exact union would also have dropped the 770 kb
gap on X and 2.1 Mb on Y that are non-PAR in both builds, so the exact union is
what is implemented.)

---

## Phase 0: contig classification and PAR (#50 items 1, 2, 4)

New module `src/allomix/contigs.py`.

- `ContigClass` enum: `AUTOSOME`, `X_NONPAR`, `X_PAR`, `Y_NONPAR`, `Y_PAR`,
  `MT`, `OTHER`. `OTHER` is any contig that is not chr1-22/X/Y/M (alt, decoy,
  unplaced, random): excluded from every analysis, counted and reported.
- `classify_contig(chrom: str, pos: int) -> ContigClass`. Handles the `chr`
  prefix as `is_sex_chrom` does today.
- PAR intervals as module constants: the union of GRCh37 and GRCh38, 1-based
  inclusive, with the per-build source values and NCBI provenance in the
  comment. A test imports `bioutils.par` when installed
  (`pytest.importorskip`) and asserts the union of its two builds equals the
  constants, so the transcription hazard is covered without a runtime
  dependency on the fork. `docs/` records how to install the fork for that
  test.
- `is_sex_chrom(chrom)` stays, documented as "not an autosome", for the
  identity-QC callers that only have a contig name and only need that. It gains
  `OTHER` contigs (they are not autosomes either).
- `classify_markers` drops `OTHER` contigs unconditionally with a new counter
  (`n_other_contig_excluded`). Any marker on an alt or decoy contig in today's
  fixtures or public data would change results; a check during implementation
  confirms there are none (expected: none).

Tests: `tests/test_contigs.py` with boundary cases on both builds' values (the
GRCh37 X/Y PAR1 asymmetry, PAR2 ends, union behaviour in the band between the
builds, `chr`-prefix handling, alt/decoy/random names).

## Phase 1: sex inference, metadata, reporting, flag retirement (#50 items 3, 5, 6)

New module `src/allomix/qc/sex.py`.

- `Sex` enum: `FEMALE`, `MALE`, `AMBIGUOUS`, `UNAVAILABLE`. `UNAVAILABLE`
  means the panel has no usable non-PAR chrX content; `AMBIGUOUS` means there
  was content and it did not resolve.
- `SexInference` dataclass: `sex`, `n_x_sites`, `n_x_het`, `x_het_rate`,
  `log10_lr` (female vs male), `chry_rel_depth | None`, `n_y_sites`, plus
  `declared: Sex | None` and `source` (`inferred`, `declared`,
  `inferred+declared`, or `conflict`).
- `infer_sex(markers: list[MarkerData], min_dp, min_gq) -> SexInference`.
  Primary signal: het count among called non-PAR chrX sites passing the depth
  and GQ filters, compared by likelihood ratio between a female model (het
  probability from the sample's own autosomal het rate, the panel's design het
  rate, scaled by a constant to be calibrated) and a male model (het
  probability = a residual spurious-het rate). Call when `|log10 LR|` exceeds
  a threshold; otherwise `AMBIGUOUS`. Minimum site count below which the result
  is `UNAVAILABLE`. Secondary: when non-PAR chrY depth is available (see Phase 3
  on where it comes from), the lab statistic `max(chrY) / median(non-sex)` with
  the lab thresholds (0.5 / 0.15 / 0.05). Used to break an ambiguous het-rate
  call, never to overrule a confident one.
  Thresholds start from the 7 public individuals and are finalised from the
  calibration script's output on the internal cohort.
- Sex is inferred for host and each donor from their reference genotypes. It is
  not inferred for the admixture sample, which is a mixture.
- `pair_status(host_sex, donor_sexes) -> MATCHED_FEMALE | MATCHED_MALE |
  MISMATCHED | UNKNOWN`. `UNKNOWN` whenever any party is `AMBIGUOUS` or
  `UNAVAILABLE`.
- Declared sex: `--recipient-sex` gains parsing (`F`, `female`, `M`, `male`,
  case-insensitive; anything else is kept as display text, ignored by logic,
  with a warning) and a new `--donor-sex` (append, one per `--donor-sample`,
  like `--donor-relationship`). Declared vs inferred disagreement on a
  confident inference is a QC **FAIL** with an explicit sample mix-up message,
  in the identity block of `qc.assess_quality`. Declared sex resolves an
  `AMBIGUOUS` inference for routing purposes (recorded as `source=declared`).
- Pipeline: `analysis.analyse_sample` computes `SexInference` per reference
  sample once and attaches `result.sex` (host + donors + pair status) before
  QC, mirroring how relatedness is attached.
- Reporting:
  - TSV/JSON: `host_sex`, `donor_sex` (per donor, joined like the relatedness
    columns), `sex_pair`, `sex_source`, `n_chrx_used`, `n_par_excluded`,
    `n_other_contig_excluded`, `n_informative_sex_chrom_excluded` (kept).
  - HTML: header shows recipient and donor sex as "declared / inferred";
    footer parameter line replaces "Sex chromosomes: included / excluded" with
    the pair status and chrX marker count; skipped PAR / other-contig counts
    shown when non-zero.
  - stderr message in `cli._run_single_sample` reworded; it no longer suggests
    a flag. PAR and other-contig exclusions are reported when non-zero.
- Retire `--use-sex-chroms`: hidden argparse option that exits with an error
  explaining the automatic behaviour and pointing at `--recipient-sex` /
  `--donor-sex`. CHANGELOG notes the option is removed entirely in a later
  release. `analyse_sample` and `classify_markers` lose the boolean and gain a
  `contig_policy` argument (default: the Phase 2 rules). The three diagnostic
  scripts that pass `use_sex_chroms=True` for genomic views get an explicit
  "include everything, diagnostic only" policy value. Version bump to 0.5.0.
- Docs: `docs/cli.md`, `docs/reports.md`, `docs/architecture.md` (new modules
  in the table), `docs/panel_guide.md` (new section on sex-chromosome content,
  including the sex-typing amplicons), `CHANGELOG.md`, and
  `paper/methods.md:103-107`, which currently describes the flag.

Tests: `tests/test_sex.py` on synthetic `MarkerData` (clean female, clean
male, male with spurious hets, too few sites, no chrX), declared-vs-inferred
conflict producing FAIL in `tests/test_qc.py`, CLI flag parsing and the
retired-flag error in `tests/test_cli.py`, report columns in
`tests/test_report.py` and `tests/test_html_report.py`.

## Phase 2: chrX for sex-matched pairs (#46)

Routing rule inside `classify_markers`, per shared marker, by `ContigClass`:

| class | rule |
|---|---|
| `AUTOSOME` | unchanged |
| `X_PAR`, `Y_PAR` | excluded, counted (`n_par_excluded`), reported when non-zero |
| `X_NONPAR`, pair `MATCHED_FEMALE` | included through the normal diploid path. Dosage math is exact. |
| `X_NONPAR`, pair `MATCHED_MALE` | included only where host and every donor are homozygous. A het call on a hemizygous chromosome is a genotyping error by definition, so het sites are dropped and counted (`n_chrx_male_het_dropped`). For hom/hom contrasts the ploidy cancels and the diploid math is exact. [derived] |
| `X_NONPAR`, pair `MISMATCHED` or `UNKNOWN` | excluded, counted in `n_informative_sex_chrom_excluded` as today |
| `Y_NONPAR`, `MT` | always excluded from the MLE |
| `OTHER` | always excluded, counted (`n_other_contig_excluded`) |

Everything else in the informative-marker path (robust refit, per-type
overdispersion, host-presence test) sees chrX markers as ordinary markers,
which is the point. Identity QC (`relatedness`, `sample_contamination`,
`shared_het_balance`) stays autosomal; a third-party contaminant of unknown
sex has unknown chrX copy number, so those checks are cleaner without it.

Bias-fitting guard:

- `calibration.bias.estimate_biases` takes per-sample sex (inferred inside
  `cmd_estimate_bias` from the same parsed markers) and skips non-PAR chrX/chrY
  het observations from male samples.
- `estimate_biases_both_het` requires every party to be female at non-PAR chrX
  sites.
- Check `calibration.error_rates` and `panel-qc` for the same exposure during
  implementation; hom sites on male chrX are legitimately pure, so error-rate
  estimation is probably fine, but verify rather than assume.

Simulator support (needed for fixtures and validation):

- `simulate` gains chrX-aware genotype generation: a `host_sex`/`donor_sex`
  pair, a count of non-PAR chrX markers, a spurious-het rate for male chrX,
  and GATK-style diploid encoding of hemizygous calls. A few PAR markers at
  real PAR coordinates.
- `blend_vcfs` uses `cn_weighted_vaf` with per-contributor chrX copy number
  (1 or 2) so sex-mismatched chrX mixes simulate correctly. Autosomal output is
  byte-identical to today (regression test).

Fixtures: `tests/test_data/sexchrom/` with small FF, MM and MF host/donor/admix
sets at a couple of fractions, plus a PAR marker and an alt-contig marker,
generated by a new `scripts/` generator (or an extension of
`generate_test_data.py`). Existing fixtures are untouched so existing expected
values hold.

Validation (single Snakemake rule or standalone script, not the full build):

1. In silico, N >= 5 seeds: FF pair with chrX on vs off at matched autosomal
   marker count, measuring LoD and bias. Expectation: improvement
   proportional to the added informative markers, no bias shift. [speculative]
2. In silico: MF pair with chrX forced in, to document the bias the gate
   prevents.
3. In silico: MM pair with spurious hets, with and without the het drop.
4. Real data: SRP434573 same-sex titrations (F1 into F3, F3 into F2, F2 into
   F1; M1 into M2, M3 into M4) against truth, chrX on vs off, as an extra arm
   of `paper/scripts/run_srp434573_allomix.py` with its own facts file and
   figure. What goes into the paper text is decided after seeing it.

Acceptance: autosomal-only results are byte-identical to current output on the
existing fixtures and on the SRP434573 runs (the gate must not perturb anything
when no chrX is routed in). The full test suite and the paper's quick build
pass before the behaviour change is considered landed.

## Phase 3: sex-mismatch cross-check (#48)

New module `src/allomix/qc/sex_mismatch.py`. Active only when `pair_status`
is `MISMATCHED` with confident sexes on both sides. Single donor; multi-donor
pairs skip the module with a note.

Readouts, in priority order, each with its own n and CI:

1. **chrY depth ratio** (needs non-PAR chrY depth in admix and in the pair's
   male reference sample). `R = max(chrY region depth) / median(non-sex region
   depth)`, the lab statistic, in the admix sample, normalised by the same
   ratio in the male reference sample of the pair, which is sequenced on the
   same panel. Estimated male fraction `= R_admix / R_male_ref`, CI by
   bootstrap over chrY regions (three on the SID panel, so the CI will be wide
   and honest). Self-normalising per panel, which the depth table above shows
   is necessary. [speculative]
   Data path: the allomix pipeline already piles up HOST/DONOR BAMs at the
   `intervals` midpoints (`refs/midpoints.vcf.gz`); the admix BAMs need the
   same midpoint pileup added (new rule, output alongside the admix panel VCF),
   and allomix `detect` gains an optional `--admix-depth-vcf` / per-sample
   depth input to read it. Without that input the readout is `NA`.
2. **chrX allele-fraction estimate** at non-PAR chrX markers informative for
   the pair, using the copy-number-weighted expectation with host and donor
   chrX copy number 1 or 2. Fitted by a one-dimensional beta-binomial MLE
   (grid + Brent, profile CI) in this module, reusing the per-marker
   likelihood pieces from `estimate.likelihood` where they take a precomputed
   expected fraction, and otherwise standing alone. It does not touch
   `PLOIDY`, `ref_dose`, or `estimate_single_donor_bb`. The signal is
   asymmetric: a male minor contributor at a female-homozygous site gives
   `f/(2-f)` (half the autosomal signal); a female hom-alt minor at a male
   hom-ref site gives `2f/(1+f)` (double). [derived]
3. chrX depth ratio is not pursued: a 10% male contribution moves female chrX
   depth by 5%, well inside the depth CV seen in the public data.

Concordance: compare the sex-chromosome estimate with the autosomal MLE
(CI overlap, or absolute difference within a tolerance in percentage points
chosen during validation). Discordance is a soft QC warning in the identity
block. It was planned to promote to REVIEW, but on the public mismatched
titrations the chrX readout's CI missed the truth in 1 of 27 (5 to 8 pp wide
from 11 to 20 markers), so, as for the host-presence disagreement, it does not
change the status until its real-sample behaviour is mapped. The headline
`donor_pct` is never blended.

Reporting: TSV/JSON columns `sexchrom_frac`, `sexchrom_ci`, `sexchrom_basis`
(`chrY-depth`, `chrX-cn`, `NA`), `sexchrom_n`, `sexchrom_concordant`; HTML
panel next to the host-presence block.

Validation: SRP434573 has five sex-mismatched titrations with truth at 0.5 to
10% minor (M3 into F1/F2/F3, F2 into M1/M2), which exercise readout 2 in both
directions, added to the paper pipeline as a facts file and figure. Plus
in-silico MF mixes from the Phase 2 simulator. Readout 1 has no data until the
haem-run chrY depth question is settled (open question 2) or the GIAB trio arm
in `plans/wetlab_validation_design.md` (kid <- mum is sex-mismatched) is
sequenced.

## Sequencing and size

| phase | lands as | estimator affected | needs revalidation | status |
|---|---|---|---|---|
| 0 | contigs module + tests | no (unless an `OTHER` contig marker exists in current data, see Phase 0) | byte-identity check | done, `0981963` |
| 1 | sex module, flags, reporting, flag retirement, docs, 0.5.0 | no | byte-identity check on fixtures | done, `d6057f8` |
| 2 | routing, bias guard, simulator, fixtures, validation rule | yes, for matched pairs | in-silico sweeps + SRP434573 same-sex arms | done, `9467b12` + `e6c3907` |
| 3 | sex-mismatch module, concordance, reporting, pipeline midpoint pileup for admix | no (independent estimate) | SRP434573 mismatched arms | done, `04881c8` (chrY readout untested) |

Before Phase 1 ships: run `scripts/sex_calibration_summary.py` on the internal
cohort (the rhAmpSeq SID joint VCFs with the lab's imputed sex as declared sex,
and the haem-run VCFs and BAMs) and set the thresholds from its output.

---

## Status (end of 2026-10-06)

All four phases are on `main`, version 0.5.0. Nothing is blocked on code;
the remaining items need internal data or a decision.

### Landed

| commit | phase | content |
|---|---|---|
| `0981963` | 0 | `allomix.contigs`: contig classes, build-free PAR mask (exact union of GRCh37 and GRCh38), PAR and non-primary contigs always excluded and counted |
| `9467b12` | 2 (simulator) | sex-aware chrX genotypes, copy-number-weighted blending, `tests/test_data/sexchrom/` FF / MM / MF fixtures |
| `d6057f8` | 1 | `allomix.qc.sex`: chrX het-rate inference, `--donor-sex`, parsed `--recipient-sex`, FAIL on declared-vs-inferred conflict, `--use-sex-chroms` retired (hidden, errors), report columns and HTML header, 0.5.0 |
| `46ae855`, `fa9e963` | 0 | cross-check test against `bioutils.par` in its interbase convention; upstream PR noted |
| `e6c3907` | 2 | `ContigPolicy.SEX_AWARE` routing (FF through the diploid path, MM hom/hom only with het drop), bias-fitting guard, `--contig-policy`, `paper/scripts/run_srp434573_chrx.py` + rule `srp434573_chrx`, `scripts/validate_sexchrom_routing.py` |
| `04881c8` | 3 | `allomix.qc.sex_mismatch`: chrX copy-number readout and (experimental, untested) chrY depth readout, `--admix-depth-vcf` / `--ref-depth-vcf`, pipeline rule `pileup_admix_bg`, HTML panel, `paper/scripts/run_srp434573_sexmismatch.py` + rule, `scripts/validate_sex_mismatch.py` |

Also written: `scripts/sex_calibration_summary.py` (de-identified calibration
summary for the internal cohort) and `scripts/build_union_bed.py` (rebuilds
the union intervals BED). The PAR data source is proposed upstream as
biocommons/bioutils PR #88 (fork branch `davmlaw/bioutils@add-par-regions`,
interbase coordinates, flat JSON like `_data/assemblies`); allomix keeps the
vendored constants until a release carries the module.

### Validation results

Phase 2, same-sex pairs. On the five same-sex SRP434573 titrations chrX adds 9
to 17 informative markers (1.5 to 3%) and moves estimates by at most 0.1 pp,
CI width ratio 0.985, QC status unchanged. In silico (5 seeds, 60 autosomal +
20 chrX): the FF gate adds no bias and tightens the CI in proportion to the
marker gain; an MF pair with chrX forced in is biased by +1.6 pp at 10% host;
the MM het drop removes a -0.9 to -2.3 pp bias. Autosomal output is
byte-identical to before.

Phase 3, mismatched pairs. On the five mismatched SRP434573 titrations the
chrX readout is concordant with the MLE in 26 of 27; chrX bias +0.3 to +1.3 pp
host with CIs 5 to 8 pp wide (MLE 0.7 to 0.8 pp) from 11 to 20 markers. The
one discordant case (F2 into M2 at 5% host) has a chrX CI that misses the
truth. In silico the readout's bias is within 0.25 pp at 1 to 10% host in both
directions.

### Departures from the plan as first written

- **No genome build.** The PAR mask is the exact union of both builds'
  intervals, so classification is a pure function of `(chrom, pos)`.
- **Discordance is a soft warning, not REVIEW.** Because of the one real
  titration whose chrX CI misses the truth, the cross-check does not change QC
  status until its real-sample dispersion is understood, as for the
  host-presence disagreement warning. Revisit with `SEXCHROM_CONCORDANCE_PP`.
- **`Sex` and `PairStatus` live in `allomix.sex_types`** (a leaf module) to
  avoid a `genotype` / `qc.sex` import cycle; `qc.sex` re-exports them.
- **The bias guard is opt-in for library callers** (`sample_sexes=None` keeps
  old behaviour); the CLI always passes inferred sexes.

### Remaining

Needs internal data (user, other machine):

1. Run `scripts/sex_calibration_summary.py` on the rhAmpSeq SID joint VCFs
   (with the lab's imputed sex as declared sex) and on haem-run VCFs and BAMs.
   Set `MALE_X_SPURIOUS_HET`, `SEX_LR_THRESHOLD` and `MIN_X_SITES` from the
   output. Known sensitivity: a male with 2 spurious hets in 20 chrX sites and
   an autosomal het rate of 0.38 is `ambiguous` under the current constants
   (log10 LR -1.1 vs threshold 2.0); on the public panel the GQ filter leaves
   0 to 1 hets per male, and declared sex resolves ambiguity by design.
2. The same run's `--bam-list --sex-bed` path answers open question 2 (chrY
   depth on haem runs).

Code, once the above is known:

3. Wire chrY relative depth into sex inference as the secondary signal
   (`SexInference.chry_rel_depth` is reserved and always None today).
4. Exercise the chrY depth readout of Phase 3 on real data, or on the GIAB
   trio arm of `plans/wetlab_validation_design.md`.
5. Run the paper rules `srp434573_chrx` and `srp434573_sexmismatch` in a full
   build and decide what goes in the text (decision 10).
6. Full test suite run before the next release (only phase subsets were run).
7. Remove the hidden `--use-sex-chroms` option in a later release.

## Open questions

1. **Which SID BED went into the union BED.** `scripts/build_union_bed.py` is
   written; the two inputs are believed to be
   `capture_kits/idt_rhampseq_sid/v1/idt_rhampseq_sid_SNPsQC.bed` and the
   idt_haem kit's vendor probe BED. Confirm the paths, and whether the SID side
   was the SNPsQC BED or the vendor target BED
   (`vendor_files/rhampseq-sample-id-panel-target-hg38.bed`).
2. **Is there depth at the SID sex-typing amplicons on haem-run BAMs?** The
   SID SNPs were added to the haem capture as probes; whether SRY / ZFY / AMELY
   / chrX_AMELX / chrX_ZFX were too is unknown. The calibration script's
   `--bam-list --sex-bed` path answers it directly. If there is no depth,
   the chrY readout of Phase 3 only works on rhAmpSeq SID runs.
3. **Paper Table S1 provenance.** `paper/supplementary.md` says the
   characterisation VCFs are from the SID SNPs inside the haem capture;
   `paper/empirical_results/README.md` says they are from the
   `idt_rhampseq_sid` amplicon directory. Which is right? If the amplicon
   assay, the paper text needs correcting, and the simulation defaults (depth
   CV 0.43, bias SD 0.018) describe the amplicon assay rather than the haem
   capture the clinical runs use.
4. **Are chrX_AMELX and chrX_ZFX in the SID BED too?** The pasted chrY rows
   suggest the X counterparts exist (the lab script uses them). If so they are
   non-SNP chrX regions that must be excluded from het-rate inference and from
   the MLE, which the 2 bp interval width and absence of a GATK record will
   usually do on their own, but the calibration output will show if one leaks
   through as a "site".
5. **Concordance tolerance.** `SEXCHROM_CONCORDANCE_PP` (2 pp) is provisional
   and mostly moot on real data, where the chrX CIs are 5 to 8 pp wide and
   concordance is driven by CI overlap. Decide after more mismatched pairs
   whether the chrX fit needs a dispersion floor and whether discordance
   should promote to REVIEW.
