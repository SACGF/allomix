# Sex-chromosome markers: #50, #46, #48 as one piece of work

Plan for the three open sex-marker issues, done together. Open questions for
review are at the bottom.

Evidence key, as used in the issues: **[data]** = observed in the repository or
the committed public data during planning, **[derived]** = closed-form from the
model, **[speculative]** = expected payoff, not yet measured.

---

## Why do them together

The three issues share one interface and one set of fixtures:

- All three need a marker's contig classified as autosome / chrX non-PAR /
  chrX PAR / chrY non-PAR / chrY PAR / chrMT, which needs PAR intervals, which
  need a genome build (#50 items 1 and 2).
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
  replacement behaviour (#46) exists. Done separately, #50 would ship a removal
  with nothing to point to, or a narrowed flag that #46 then removes again.
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
Counting called non-PAR chrX sites per reference sample:

| sample | sex | chrX called | chrX het | het rate | autosomal het rate |
|---|---|---:|---:|---:|---:|
| F1 | F | 24 | 17 | 0.71 | 0.50 |
| F2 | F | 27 | 15 | 0.56 | 0.51 |
| F3 | F | 25 | 11 | 0.44 | 0.52 |
| M1 | M | 26 | 1 | 0.04 | 0.52 |
| M2 | M | 26 | 1 | 0.04 | 0.53 |
| M3 | M | 26 | 1 | 0.04 | 0.51 |
| M4 | M | 24 | 2 | 0.08 | 0.50 |

Two things follow. The het-rate signal is strong even with ~25 sites. And every
male carries one or two spurious chrX het calls, which is the bias-fitting and
dosage hazard the issues predicted, now observed.

### All 27 public chrX amplicons are non-PAR [data]

Positions run from 8.5 Mb to 131.3 Mb on hg38; PAR1 ends at 2,781,479 and PAR2
starts at 155,701,383. So on the one panel we can validate against, the PAR
mask is a no-op, and PAR handling has to be covered by synthetic fixtures.

### chrX capture efficiency is panel-specific, so depth is not a clean sex signal [data]

Median chrX depth over median autosomal depth in the admix samples:

| mixture | sexes | panel version | chrX / autosome depth |
|---|---|---|---:|
| F1 into F3 | F, F | v2 | 1.45 to 1.67 |
| M3 into F1 | M minor, F major | v2 | 1.45 to 1.64 |
| M1 into M2 | M, M | v1 | 0.52 to 0.55 |

Female diploid chrX comes out at 1.5x autosomal depth on this panel, not 1.0x,
and the two panel versions differ in chrX amplicon count (25 or 26 vs 19). A
chrX depth ratio therefore only means something relative to a same-panel,
same-sex reference. This pushes sex inference onto heterozygosity (primary) with
depth as a secondary signal only, and it means #48's chrX depth readout is
weaker than its allele-fraction readout.

### No panel in hand carries chrY [data]

SRP434573 has zero chrY amplicons. The deployment panel's chrY content is
unknown. The `scripts/build_panel_vcf.py` docstring says BED sites absent from
the joint VCF are "typically sex-chromosome or non-SNP markers that the joint
call never saw an alt at", which hints the deployment BED does carry
sex-chromosome targets without a cohort SNP. If any are chrY, their depth is
still available from the forced admix `mpileup` even with no polymorphism, but
only if those sites survive into the admix target list. See open question 2.

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
- Nothing reads `##contig` lengths today. `qc/runmeta.read_run_units` is the
  existing pattern for reading header lines via `VCF.raw_header`.
- Test fixtures in `tests/test_data/` are autosome-only and their `##contig`
  lines all carry chr1's length; `simulate.write_genotype_vcf` writes contigs
  with no length at all. Build inference from the header will need fixtures
  with real lengths, and a defined fallback when lengths are absent.
- The PAR data lives in a bioutils fork
  (`davmlaw/bioutils@add-par-regions`, module `bioutils.par`, 1-based
  inclusive, generated from NCBI `par_align.gff`, GRCh37 and GRCh38). Not yet
  proposed upstream; no upstream PR or issue exists. bioutils is already an
  optional `scripts` dependency but is not installed in the current venv.
- allomix 0.4.2 is on PyPI. PyPI rejects any direct-URL requirement, so a
  git-pinned fork cannot appear in `dependencies` or in any extra.

---

## Phase 0: contig classification, genome build, PAR (#50 items 1, 2, 4)

New module `src/allomix/contigs.py`.

- `GenomeBuild` enum: `GRCH37`, `GRCH38`, `UNKNOWN`.
- `ContigClass` enum: `AUTOSOME`, `X_NONPAR`, `X_PAR`, `Y_NONPAR`, `Y_PAR`,
  `MT`, `X_UNKNOWN_PAR`, `Y_UNKNOWN_PAR`. The last two are what chrX/chrY
  positions become when the build is unknown, so "we cannot tell if this is
  PAR" is a first-class state rather than a guess.
- `classify_contig(chrom: str, pos: int, build: GenomeBuild) -> ContigClass`.
  Handles the `chr` prefix as `is_sex_chrom` does today.
- PAR intervals for both builds as module constants, 1-based inclusive, with
  the NCBI provenance in the comment. 16 numbers. A test imports
  `bioutils.par` when it is installed (`pytest.importorskip`) and asserts the
  constants equal `get_par_map` for both builds, so the transcription hazard
  the fork was built to avoid is covered without a runtime dependency on the
  fork. See open question 1.
- `infer_build(vcf_path) -> GenomeBuild` from `##contig` lengths: chr1 is
  248,956,422 in GRCh38 and 249,250,621 in GRCh37; chrX is 156,040,895 vs
  155,270,560. Use chr1 and chrX together and return `UNKNOWN` if they
  disagree or are absent. Verify `cyvcf2.VCF` exposes `seqnames`/`seqlens`
  before relying on them; fall back to parsing `raw_header` otherwise.
- CLI: `--genome-build {GRCh37,GRCh38}` override on `detect`, `timeline`,
  `estimate-bias`, `estimate-errors`, `panel-qc`. Inferred by default; a
  mismatch between the flag and the header is an error, not a warning.
- `is_sex_chrom(chrom)` stays, documented as "not an autosome", for the
  identity-QC callers that only have a contig name and only need that.
  Alt/decoy/unplaced contigs keep today's behaviour (treated as autosomal);
  flagged in open question 9.
- Graceful degradation: `UNKNOWN` build makes every chrX/chrY marker
  `*_UNKNOWN_PAR`, which every downstream consumer treats as "exclude", so a
  VCF with no contig lengths behaves exactly as today.

Tests: `tests/test_contigs.py` with boundary cases on both builds (the
GRCh37 X/Y PAR1 asymmetry, PAR2 ends, `chr`-prefix handling), build inference
from real-length headers, from the current fixtures' wrong lengths (should be
`UNKNOWN`, since chr1 and chrX disagree), and from a header with no lengths.

## Phase 1: sex inference, metadata, reporting, flag retirement (#50 items 3, 5, 6)

New module `src/allomix/qc/sex.py`.

- `Sex` enum: `FEMALE`, `MALE`, `AMBIGUOUS`, `UNAVAILABLE`. `UNAVAILABLE`
  means the panel has no usable non-PAR chrX content (or the build is
  unknown); `AMBIGUOUS` means there was content and it did not resolve.
- `SexInference` dataclass: `sex`, `n_x_sites`, `n_x_het`, `x_het_rate`,
  `log10_lr` (female vs male), `chry_depth_ratio | None`, `n_y_sites`,
  plus `declared: Sex | None` and `source` (`inferred`, `declared`,
  `inferred+declared`, or `conflict`).
- `infer_sex(markers: list[MarkerData], build, min_dp) -> SexInference`.
  Primary signal: het count among called non-PAR chrX sites with
  `DP >= min_dp`, compared by likelihood ratio between a female model (het
  probability taken from the sample's own autosomal het rate, which is the
  panel's design het rate, scaled by a constant to be calibrated) and a male
  model (het probability = a spurious-het rate, about 0.04 from the table
  above). Call when `|log10 LR|` exceeds a threshold; otherwise `AMBIGUOUS`.
  Minimum site count below which the result is `UNAVAILABLE`. Secondary: if
  non-PAR chrY sites exist, `median DP(chrY) / median DP(autosomes)`; a
  female should sit near zero on any panel, so an absolute threshold is
  defensible here where it is not for chrX. Used only to break an ambiguous
  het-rate call, never to overrule a confident one. [speculative] until seen
  on data with chrY.
  Thresholds are calibrated on the 7 public individuals, which is thin. A
  de-identified summary script for the internal cohort is proposed in open
  question 7.
- Sex is inferred for host and each donor from their reference genotypes.
  It is not inferred for the admixture sample, which is a mixture.
- `pair_status(host_sex, donor_sexes) -> MATCHED_FEMALE | MATCHED_MALE |
  MISMATCHED | UNKNOWN`. `UNKNOWN` whenever any party is `AMBIGUOUS` or
  `UNAVAILABLE`.
- Declared sex: `--recipient-sex` gains parsing (`F`, `female`, `M`, `male`,
  case-insensitive; anything else is kept as display text and ignored by
  logic, with a warning) and a new `--donor-sex` (append, one per
  `--donor-sample`, like `--donor-relationship`). Declared vs inferred
  disagreement on a confident inference is a QC warning that promotes to
  REVIEW, alongside the relatedness mismatch and swap checks in the identity
  block of `qc.assess_quality`. Declared sex may resolve an `AMBIGUOUS`
  inference for routing purposes (recorded as `source=declared`). See open
  question 5.
- Pipeline: `analysis.analyse_sample` computes `SexInference` per reference
  sample once and attaches `result.sex` (host + donors + pair status) before
  QC, mirroring how relatedness is attached.
- Reporting:
  - TSV/JSON: `host_sex`, `donor_sex` (per donor, joined like the relatedness
    columns), `sex_pair`, `sex_source`, `n_chrx_used`,
    `n_informative_sex_chrom_excluded` (kept for compatibility).
  - HTML: header shows recipient and donor sex as "declared / inferred";
    footer parameter line replaces "Sex chromosomes: included / excluded" with
    the pair status and chrX marker count.
  - stderr message in `cli._run_single_sample` reworded; it no longer
    suggests a flag.
- Retire `--use-sex-chroms`. argparse keeps the name as a hidden option that
  exits with an explanatory error pointing at the automatic behaviour, so
  existing pipelines fail loudly rather than silently changing meaning.
  `analyse_sample` and `classify_markers` lose the boolean and gain a
  `contig_policy` argument (default: the rules in Phase 2). The three
  diagnostic scripts that pass `use_sex_chroms=True` for genomic views get an
  explicit "include everything, diagnostic only" policy value. Version bump
  to 0.5.0 for the breaking CLI change. See open question 6.
- Docs: `docs/cli.md`, `docs/reports.md`, `docs/architecture.md` (two new
  modules in the table), `docs/panel_guide.md` (new section on sex-chromosome
  content), `CHANGELOG.md`, and `paper/methods.md:103-107`, which currently
  describes the flag.

Tests: `tests/test_sex.py` on synthetic `MarkerData` (clean female, clean
male, male with two spurious hets, too few sites, unknown build), declared-vs-
inferred conflict producing REVIEW in `tests/test_qc.py`, CLI flag parsing and
the retired-flag error in `tests/test_cli.py`, report columns in
`tests/test_report.py` and `tests/test_html_report.py`.

## Phase 2: chrX for sex-matched pairs (#46)

Routing rule inside `classify_markers`, per shared marker, by `ContigClass`:

| class | rule |
|---|---|
| `AUTOSOME` | unchanged |
| `X_PAR`, `Y_PAR` | excluded from the MLE, counted separately (`n_par_excluded`). Diploid in both sexes, but on GRCh37 reads split between the X and Y PAR copies and on GRCh38 the Y PAR is usually hard-masked, so genotype quality there depends on the reference preparation. Safer out. See open question 4. |
| `X_NONPAR`, pair `MATCHED_FEMALE` | included through the normal diploid path. Dosage math is exact. |
| `X_NONPAR`, pair `MATCHED_MALE` | included only where host and every donor are homozygous. A het call on a hemizygous chromosome is a genotyping error by definition, so het sites are dropped and counted (`n_chrx_male_het_dropped`). For hom/hom contrasts the ploidy cancels and the diploid math is exact. [derived] See open question 3. |
| `X_NONPAR`, pair `MISMATCHED` or `UNKNOWN` | excluded, counted in `n_informative_sex_chrom_excluded` as today |
| `Y_NONPAR`, `MT`, `*_UNKNOWN_PAR` | always excluded from the MLE |

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
  sites (a both-het call in a male is an error, not a VAF-0.5 anchor).
- Check `calibration.error_rates` and `panel-qc` for the same exposure during
  implementation; hom sites on male chrX are legitimately pure, so error-rate
  estimation is probably fine, but verify rather than assume.

Simulator support (needed for fixtures and validation):

- `simulate` gains chrX-aware genotype generation: a `host_sex`/`donor_sex`
  pair, a count of non-PAR chrX markers, a spurious-het rate for male chrX
  (default from the table above), and GATK-style diploid encoding of
  hemizygous calls. A few PAR markers at real PAR coordinates for both builds.
- `blend_vcfs` uses `cn_weighted_vaf` with per-contributor chrX copy number
  (1 or 2) so sex-mismatched chrX mixes simulate correctly. Autosomal output is
  byte-identical to today (regression test).
- `write_genotype_vcf` and `write_joint_vcf` write real GRCh38 contig lengths
  so build inference works on simulated data.

Fixtures: `tests/test_data/sexchrom/` with small FF, MM and MF host/donor/admix
sets at a couple of fractions, generated by a new `scripts/` generator (or an
extension of `generate_test_data.py`). Existing fixtures are untouched so
existing expected values hold.

Validation (single Snakemake rule or standalone script, not the full build):

1. In silico, N >= 5 seeds: FF pair with chrX on vs off at matched autosomal
   marker count, measuring LoD and bias. Expectation: improvement
   proportional to the added informative markers, no bias shift. [speculative]
2. In silico: MF pair with chrX forced in, to document the bias the gate
   prevents. This is the regression that justifies the gate.
3. In silico: MM pair with spurious hets, with and without the het drop.
4. Real data: SRP434573 same-sex titrations (F1 into F3, F3 into F2, F2 into
   F1; M1 into M2, M3 into M4) against truth, chrX on vs off. Fits
   `paper/scripts/run_srp434573_allomix.py` as an extra arm. Whether it
   becomes a paper fact or stays a validation output is open question 10.

Acceptance: autosomal-only results are byte-identical to current output on the
existing fixtures and on the SRP434573 runs (the gate must not perturb anything
when no chrX is routed in). The full test suite and the paper's quick build
pass before the behaviour change is considered landed.

## Phase 3: sex-mismatch cross-check (#48)

New module `src/allomix/qc/sex_mismatch.py`. Active only when `pair_status`
is `MISMATCHED` with confident sexes on both sides. Single donor first;
multi-donor is open question 8.

Readouts, in priority order, each with its own n and CI:

1. **chrY depth ratio** (needs non-PAR chrY sites). `R = median DP(chrY) /
   median DP(autosomes)` in the admix sample, normalised by the same ratio in
   the male reference sample of the pair (host or donor), which is sequenced on
   the same panel. Estimated male fraction `= R_admix / R_male_ref`, CI by
   bootstrap over chrY sites. Self-normalising per panel, which the depth table
   above shows is necessary. Untestable on any data in hand; ships gated and
   labelled experimental. [speculative]
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
chosen during validation). Discordance is a QC warning in the identity block
and promotes to REVIEW. The headline `donor_pct` is never blended.

Reporting: TSV/JSON columns `sexchrom_frac`, `sexchrom_ci`, `sexchrom_basis`
(`chrY-depth`, `chrX-cn`, `NA`), `sexchrom_n`, `sexchrom_concordant`; HTML
panel next to the host-presence block.

Validation: SRP434573 has five sex-mismatched titrations with truth at 0.5 to
10% minor (M3 into F1/F2/F3, F2 into M1/M2), which exercise readout 2 in both
directions. Plus in-silico MF mixes from the Phase 2 simulator. Readout 1 has
no data until the deployment panel's chrY content is known (open question 2)
or the GIAB trio arm in `plans/wetlab_validation_design.md` (kid <- mum is
sex-mismatched) is sequenced.

## Sequencing and size

| phase | lands as | estimator affected | needs revalidation |
|---|---|---|---|
| 0 | contigs module + tests | no | no |
| 1 | sex module, flags, reporting, flag retirement, docs, 0.5.0 | no (autosomal output unchanged) | byte-identity check on fixtures |
| 2 | routing, bias guard, simulator, fixtures, validation rule | yes, for matched pairs | yes: in-silico sweeps + SRP434573 same-sex arms |
| 3 | sex-mismatch module, concordance, reporting | no (independent estimate) | SRP434573 mismatched arms |

Phases 0 and 1 are small and mostly plumbing. Phase 2 is the largest and the
only one that changes what the validated estimator sees. Phase 3 can trail.
Each phase is committed to `main` with its tests; nothing in a later phase is
needed for an earlier one to be correct.

---

## Open questions

1. **PAR data dependency.** PyPI rejects direct-URL requirements, so the fork
   cannot be a runtime or extras dependency while allomix stays on PyPI. I
   recommend vendoring the 16 PAR boundary numbers into `allomix/contigs.py`
   with provenance, and a test that cross-checks them against `bioutils.par`
   whenever the fork is installed. If bioutils releases the module upstream,
   switch to `bioutils>=<that version>` as a runtime dependency. Alternative:
   drop the PyPI constraint and pin the fork by git URL. Which?
2. **Deployment panel chrY content.** This decides whether #48's primary
   readout is testable before wetlab, and whether the sex-inference chrY
   secondary ever runs. One command on the internal BED answers it:
   `cut -f1 <panel.bed> | sort | uniq -c`. Related: does the forced admix
   `mpileup` target list include BED sites that have no record in the panel
   VCF (the "typically sex-chromosome" sites `build_panel_vcf.py` skips)? If
   not, chrY depth never reaches allomix even if the BED has chrY.
3. **Male-male chrX in Phase 2 or later?** The issue suggested female-female
   first. The male-male case adds one filter (drop het calls) and the
   hom/hom math is exact, so I recommend including it. The cost is one more
   fixture and validation arm (SRP434573 has M1 into M2 and M3 into M4).
4. **PAR markers: exclude everywhere (recommended) or treat as autosomal?**
   Diploid in both sexes, but genotype quality depends on how the reference
   handled the Y PAR. Excluding costs nothing on the public panel (zero PAR
   amplicons) and is the conservative default. Counted and reported either
   way.
5. **Declared sex semantics.** (a) Add `--donor-sex`? (b) May a declared sex
   resolve an ambiguous inference for chrX routing (recommended: yes,
   recorded as `source=declared`)? (c) Declared-vs-inferred conflict:
   REVIEW (recommended, consistent with the relatedness mismatch) or FAIL?
6. **`--use-sex-chroms` retirement style.** Hidden option that errors with an
   explanation (recommended), or silent removal (argparse "unrecognized
   arguments")? Either way a 0.5.0 bump.
7. **Threshold calibration on internal data.** The het-rate thresholds come
   from 7 public individuals. I can write a standalone script that, for each
   reference sample on the internal shares, prints only: declared sex (if
   known), number of called non-PAR chrX sites, number het, median chrX and
   autosomal depth, and chrY site count and depth if any. No identifiers or
   coordinates. Worth running before Phase 1 ships, or acceptable to ship on
   the public calibration and tune after?
8. **#48 multi-donor scope.** Single donor only in the first cut
   (recommended); with two donors the sex-chromosome estimate is a sum over
   whichever contributors differ from the host and is harder to read.
9. **Alt, decoy and unplaced contigs** are treated as autosomal today and would
   stay that way. Markers on them would be odd on any sample-ID panel. Leave
   as is, or classify them `OTHER` and exclude?
10. **Paper impact.** Should the SRP434573 chrX-on arm (Phase 2) and the
    sex-mismatch readout (Phase 3) become paper facts and figures, or stay as
    validation outputs under `output/` until the wetlab phase? Adding them
    means new rules in `paper/Snakefile` and text in Methods and Results.
11. **Build inference fallback when the header has no contig lengths.** The
    plan treats that as `UNKNOWN` and excludes chrX, which preserves today's
    behaviour for such files. The alternative is to require `--genome-build`
    explicitly whenever the panel has chrX content and the header is silent.
    Recommended: `UNKNOWN` plus a stderr note saying chrX was excluded for
    want of a build, so nothing breaks but the cost is visible.
