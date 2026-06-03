# Donor Chimerism Tool — Overall Build Plan

Master plan for allomix, a general-purpose, panel-agnostic NGS chimerism monitoring tool.

Completed steps are summarised to the decision + rationale; the full implementation detail lives in the code, tests, `doc/`, and the per-step detail docs (`claude/*_plan.md`). Open work keeps enough context to be actioned directly.

---

## Step 1: Project Setup ✅
allomix on PyPI (v0.0.1); package under `src/allomix/`; `monitor` + `timeline` CLI. Context in `CLAUDE.md` / `README.md`.

## Step 2: Input Format — VCF, not BAM ✅
**Decision: VCF-only, no BAM in v1.** VCF AD (ref/alt/depth) provides everything the estimator needs, and the upstream pipeline can be shaped to emit what we want. Min FORMAT fields: GT, AD, DP. Detail: `claude/step2_bam_vs_vcf_decision.md`.

## Step 3: Synthetic Test Data ✅
`src/allomix/simulate.py` + `scripts/generate_*.py`. VCF blending with binomial-sampled allele counts; per-marker capture bias (`marker_bias_sd`); optional overdispersion `rho` (Step 21). Plain-text VCF I/O so the simulator is dependency-light.

## Step 4: Reference Tools ✅
**Why our approach:** MLE likelihood from Crysup & Woerner 2022 (Demixtify Formula 5, simplified by known genotypes); per-marker amplification-bias correction (Vynck); grid-search + Brent refinement (All-FIT / Conpair). **License: reimplement the published math independently, cite the paper; do not copy AGPL (Demixtify) / non-commercial (Conpair) code.** Detail: `claude/historical/step4_reference_tool_analysis.md`.

## Step 5: Implementation Plan ✅
6 modules (genotype, chimerism, bias, qc, report, cli); MLE algorithm + test plan. Detail: `claude/step5_implementation_plan.md`.

## Step 6: Core Algorithm (single-donor) ✅
MLE estimation end-to-end: VCF parse → Vynck marker classification → grid search + Brent + profile-likelihood CI → QC → TSV/JSON. CLI commands in `README.md`.

## Step 7: Multi-Donor Support ✅
host + 2 donors. Sibling-trio simulation (Mendelian segregation), markers informative for ANY donor, triangular grid → Nelder-Mead → per-donor profile CIs. CLI auto-detects donor count. Validated on sibling donors. Detail: `claude/multi_donor_plan.md`.

## Step 8: Bias Correction ✅
`bias.py`: per-marker bias = median(VAF_het − 0.5) from a training cohort, applied to the expected REF weight in the MLE. ~15% MAE / ~25% max-error reduction at 2000x. CLI: `estimate-bias` + `monitor --bias-table`.

## Step 9: In-Silico Validation ✅
Simulator models four noise sources, all calibrated from empirical panel characterisation (210 VCFs / 18,047 samples, 76-SNP rhAmpSeq): per-marker amplification bias (SD ~0.018), non-uniform depth (CV 0.43), sequencing error (ε=0.01), locus dropout (1.6%). Sub-2% MAE across depths (50–1000x) and relatedness (unrelated→sibling). **CI coverage was low (25–58%) because the plain binomial does not model systematic bias / non-uniform depth — this motivated the beta-binomial work (Step 13) and the per-site error/dropout steps (14, 15).**

## Step 10: VariantGrid Integration 🔲
VG stores donor/host genotypes, exports VCF, ingests allomix JSON per patient, renders the timeline chart. JSON schema + API integration TBD.

## Step 11: Real Sample Validation ✅
Joint-called VCFs for the idt_rhampseq_sid panel produced on /tau. Batch runners drive allomix across the patient list. **Key interpretive note: allomix reports bulk-DNA chimerism (a cell-type-weighted average), so it is not apples-to-apples with sorted-cell flow (CD45/CD3/CD13) and tracks CD13 myeloid more closely than CD45 in lineage-disparate samples (e.g. RCAR).** The current two-phase-pipeline results are run9 → run10 (Steps 23, 27); these supersede the old all-GATK `validation_run_new_bias2`. Next phase: controlled dilution series for quantitative accuracy.

## Step 12: Per-Marker Likelihood Context Refactor ✅
Shipped as `PanelCalibration` (`chimerism.py`): one dataclass bundling per-marker biases + per-direction error rates, built once per estimator call and threaded through the likelihood closures, replacing the separate `marker_biases`/`marker_errors` kwargs. No behaviour change.

## Step 13: Beta-Binomial Goodness-of-Fit ✅
The MLE already fit overdispersion `ρ`, but `qc.py` standardised residuals by binomial variance so `gof_pval` was ~0 even on good fits. Fixed: GoF now uses beta-binomial variance with df corrected for fitted params, plus an error-adjusted `p_alt` at saturated (f≈0/1) markers so a ~1% error residual no longer blows up the chi-sq. Detail: `claude/beta_binomial_plan.md`.

## Step 14: Empirical Per-Site Error Rates ✅ (mechanism) — operational table-build is Step 28
Replaced the global `--error-rate 0.01` with empirically measured per-site, per-direction rates (ALT-read rate at hom-ref sites, REF-read rate at hom-alt). Shipped: `error_rates.py` (`estimate_error_rates`, save/load), the asymmetric REF/ALT error model in the beta-binomial likelihood, the `estimate-errors` CLI subcommand, `--error-table` on `monitor`/`timeline`, and tests. This is the per-site background the host-presence detector wants (Step 20) and the principled fix for the TP53-style artifact (Step 27/28). **Building and adopting a real panel-of-normals table on the /tau data is the operational follow-up tracked in Step 28** (run9/run10 still ran `global-fallback`).

## Step 15: Per-Site Dropout Rate 🔲
Integrate per-site no-call rate (already available from the bias cohort) into the likelihood so flaky sites are downweighted. Detail: `claude/per_site_dropout_plan.md`.

## Step 16: GQ-Weighted Marker Contributions 🔲
Replace the hard `--min-gq 20` cutoff with a per-marker likelihood weight so borderline genotypes stay informative but downweighted. Small expected gain; cheap once the per-marker likelihood is being touched. Detail: `claude/16_gq_weighted_markers_plan.md`.

## Step 17: Per-Base Quality-Aware Likelihood 🔲 (skeptical, may not ship)
Weight each read's likelihood contribution by its base quality (Conpair-style). **Decision gate: revisit only if Steps 14–16 land and real-data CIs are still wider than we want.** It is the most invasive remaining algorithm step (upstream `bcftools mpileup -a FORMAT/QS`, new module, simulator + CLI changes) for a likely small gain on a Q30+ panel, and Step 14 attacks the same slack far more cheaply. Detail: `claude/17_bq_aware_plan.md`.

## Step 18: Publication 🟡 IN PROGRESS
Method paper (vibepaper), target JMD Technical Advance. In-silico validation, multi-donor, and simulation-calibration supplementary figures done. Cite Crysup & Woerner 2022 + Vynck.

Remaining:
- [ ] Add bias-stability figure (`fig_bias_stability.png`) to `results.md` (validates the fixed-bias-per-marker assumption).
- [ ] Decide whether the ablation study (Fig S4) adds a no-overdispersion (binomial) baseline.
- [ ] **Rewrite the joint-calling references for the two-phase pipeline** (`methods.md:11`, `discussion.md:29` still carry the wrong "joint calling preserves admix AD" claim; `supplementary.md` is fine). Gate on the real-data results being final.
- [ ] Real-sample validation (Step 11 dilution series) still needed before submission.

## Step 19: Intronic Shoulder Marker Evaluation 🔲
Our capture panel's depth extends past exon boundaries into flanking introns, a source of extra high-heterozygosity (near-0.5 MAF) markers → tighter estimate, lower LOD. **Open risk to check: read-end mapping bias.** These SNPs sit near read ends, so ALT reads carry a terminal mismatch and are preferentially soft-clipped / MAPQ-penalised, dropping ALT from AD and skewing VAF toward REF. This is allele-asymmetric and NOT caught by the (allele-blind) depth filter. Analysis (summary-stats script against /tau, no coordinates/IDs out): per-marker median het VAF, depth, dropout vs intron offset — does het VAF leave 0.5 only after depth has already fallen below QC, or while depth is still healthy? Latter → need an explicit allele-balance filter. From a 2026-05-27 design discussion. (Note: this is the same read-end/soft-clip failure mode the Step 27 artifact filter keys on.)

## Step 20: Host-Presence Detection at Donor-Homozygous Markers ✅ (route A) / 🔲 (route B)
A detection test for "is the host present at all?", separate from the fraction MLE, for low-level host re-occurrence (relapse). Uses markers where the donor is homozygous and the host carries the donor-absent allele: that allele sits at the sequencing-error background in pure donor, so its pooled counts give a one-sided test + an LRT host-fraction estimate against that background — freed from the MLE's single shared `ρ` and global error rate, which blunt it at low fraction.

Route A (reported alongside the MLE) is **done**: `src/allomix/detect.py` (`host_presence_test`), on by default in `monitor`, results in `batch.tsv` (`host_present_p`, `host_f_est`, CI, `host_err_source`, `host_artifact_filtered`). The achievable limit is error-floor-bound until the per-site background (Step 14 / Step 28) is supplied (`host_err_source=global-fallback` today). Route B (unified two-component likelihood) is still TODO. Detail: `claude/20_host_presence_detection_plan.md`.

## Step 21: Calibrate Simulator Overdispersion for Realistic LoD 🟡 IN PROGRESS
**Key finding (2026-05-28):** the simulator drew from a pure binomial, so its in-silico LoD (~0.13–0.32%) was optimistic by ~3–5x vs real run3 LoDs (~0.5–1%). Beta-binomial variance approaches `p(1-p)/(ρ+1)`, so effective depth caps near `ρ+1` reads — **overdispersion, not depth, is the dominant LoD control at clinical coverage.** Done: `rho` arg in the simulator; `plot_lod_saturation.py` + `run_overdispersion_lod.py` (Supp S7/S8); `scripts/diagnose_sample.py` prints per-sample fitted `rho`.

TODO:
- [ ] Calibrate `rho` from real per-sample fits and re-run the headline `lod_validation` so it reflects real overdispersion (expensive job; warn first).
- [ ] The simulator applies one global `rho` to every marker/allele, including the near-zero donor-absent allele where overdispersion is not physical. A marker-type/allele-aware model is needed before `rho` can validate host-presence detection (Step 20).
- [ ] Decide whether `discussion.md`'s headline LoD switches to the overdispersion-calibrated number.

## Step 22: Pileup / Two-VCF Model ✅
**Decision (2026-05-29): committed to pileup-only.** Why: GATK `HaplotypeCaller -ERC GVCF` strips minority ALT reads at hom-ref blocks (verified: 0 ALT across ~9M reads), destroying the low-fraction signal. The two-phase pipeline uses GATK only for HOST/DONOR `GT` and forced `bcftools mpileup` for ADMIX `AD`. Migration landed: `--vcf` removed from `monitor`/`timeline` (kept on `estimate-bias`/`estimate-errors`), fixtures rebuilt as panel/admix pairs. Full rationale in `doc/joint_calling.md` (including why a somatic caller is also wrong).

## Step 23: Widen Force-Output Panel ✅
Recovered marginal markers (2026-05-29). Full write-up + before/after numbers: `claude/2026-05-29_wider_panel_validation_notes.md`. What was needed beyond the wider panel:
1. gnomAD v4.1-derived panel build (`scripts/build_force_output_panel.sh`); recommended `output/union_sid_haem_gnomad_af05.vcf.gz` (258 sites).
2. `bcftools call -A` so the admix VCF keeps the panel ALT at hom-ref sites (~48 informative SNPs/patient otherwise lost in the join).
3. `-e 'ALT="."'` in the `panel_tsv` rule (force-output REF-only rows produce malformed PL under `-A`).
4. Skip indels in `parse_vcf` (pileup can't count indel reads like GATK reassembly).
5. GT/AD consistency check on host/donor (drops GATK miscalls where a called het has AD VAF <0.35 or >0.65).

Result (run9): n_informative up across the board, donor% matches flow on every sample.
- **Open (user owns):** share run9 gains + "any detection" results with the post-doc; LNAN host-presence p=0.16 is borderline — discuss.

## Step 24: Overdispersion / REVIEW Samples (NDAD, BHOA, PCAH) 🔲
These come up `gof_pval = 0.0000` (QC=REVIEW): the chimerism fraction matches flow, but residual per-marker variance exceeds beta-binomial expectation. The TP53 artifact (Step 27) was one contributor and is now filtered. **Likely resolved by the per-site error null (Step 28)**, which down-weights each noisy locus by its own measured background rather than inflating a global `rho`. If a gap remains after Step 28, refit `rho` via `scripts/diagnose_sample.py` or move to per-marker-type overdispersion (Step 21).

## Step 25: Host-Presence Visualisation and Per-Marker Diagnostics ✅ (2026-06-01)
Standalone host-presence diagnostic plots (`scripts/`, not the package), all artifact-filter-aware and mirroring the `monitor` run (sex-chroms off by default; Step 27 filter convention):
- `plot_presence_lod_curve.py` — detection probability vs spiked level (binomial vs beta-binomial).
- `plot_host_presence_per_marker.py` — dose-normalised implied host fraction per marker vs rank, with the pooled MLE line.
- `host_presence_manhattan.py` — genomic view, chromosome-banded, nearest-gene labels on upregulated markers (uses `output/refseq109_genes.bed`).
- `host_presence_markers_vcf.py` — per-sample VCF of donor-hom markers with per-marker INFO for VEP / driver-panel intersection.
- `plot_chimerism_comparison.py` draws the host-presence estimate + CI as a green/grey diamond beside each point (the `*_presence.png` plots); the old standalone `plot_host_presence.py` is retired.

**Convention settled:** donor % wherever a chimerism value is plotted; host fraction only in the per-marker diagnostics. **Refactor done (2026-06-02):** `allomix.analysis.analyse_sample` (returns `SampleAnalysis`) + public `detect.donor_hom_markers` give one shared "run a sample" path that `cli` and the plots consume (no private-internal imports; sex-chrom + artifact handling defined once). Also consolidated `marker_key` into `genotype.py` and renamed `simulate.parse_vcf`→`parse_text_vcf`. Reading guide: `doc/architecture.md`.

- [ ] **LoD-overlay full re-run (pending):** the presence overlay on `fig5_lod_curves.png` is a proof-of-concept. For the real run, pick one realistic per-site error rate (ideally Step 14's empirical value) and run `run_lod_validation.py` + `run_presence_lod.py` at it (presence LoD collapses to ~14% at 1% error, so the overlay needs matched error). Expensive; warn first.
- [ ] **Remaining refactor:** `host_presence_markers_vcf.py` (uses public `select_donor_hom_markers`) and `diagnose_sample.py` (imports `qc._error_adjusted_p_alt`) should migrate onto `analyse_sample`/public surfaces when next touched.

## Step 26: Sex-Chromosome Handling ✅ (2026-06-01)
X/Y/M allele dosage is wrong in sex-mismatched transplants. `classify_markers(use_sex_chroms=False)` (default) excludes them and reports `n_sex_chrom_excluded`; CLI `--use-sex-chroms` re-enables. Reportable runs exclude them; the Manhattan / markers-VCF diagnostics keep them visible for investigation, but `plot_host_presence_per_marker.py` now defaults to excluding them (with `--use-sex-chroms` to opt back in) so it mirrors the reported run. Cost in run9: 5/6/7 chrX markers for NDAD/BHOA/PCAH. Re-enable per run once host+donor sex are confirmed matched (sex being added to the project xls; user owns).

## Step 27: TP53/17p "Clonal LOH" Signal — REFUTED as an Alignment Artifact ❌ (2026-06-01)

The host-presence diagnostics flagged chr17:7676483 (TP53 intron 3) as a several-fold host-allele spike in NDAD/BHOA/PCAH, which looked like it might be a clonal-LOH relapse signal. It is **not** — it is a single-base alignment artifact. How that was settled:

- **Segmental check (decisive).** It is a lone single-base spike; immediate neighbours are flat at high depth (chr17:7676301 at 182 bp, rs1042522 at 7676154 at 329 bp, both ~0%). No CN-LOH/UPD/amplification has a footprint under 182 bp, so both the relapse and somatic-rescue (UPD) readings are excluded.
- **Read geometry.** The host-allele reads are strand-skewed (e.g. 2:34, 4:76), soft-clipped (SCBZ −5 to −11), read-position-biased (RPBZ 4–11) vs ~0 at clean neighbours — misalignment. Structural cause: a 16-bp deletion at chr17:7676325 (low-complexity PIN3 region) makes spanning reads soft-clip and dump spurious bases at 7676483.
- **Not host-specific.** The artifact VAF is ~2–3% essentially constant across the whole host-fraction range, *including pure-donor samples with zero host* (PNOL 3.04%, GBRI 1.93%) and RCAR at 58% host (2.37%). A host-derived cause would scale with host fraction; it doesn't. Control: LNAN (host+donor both genuinely het) shows the real allele at 56% with balanced strands.
- **Why the donor-pileup "proof" was wrong.** A clean true-donor pileup only rules out a *universal reference* artifact, not a per-library alignment artifact in the admix pileup. The admix reads must be judged on their own, and they fail every read-quality test.

**Handling (done).** Read-level artifact filter in `allomix.detect` (`ArtifactThresholds` + `_is_artifact_marker`), on by default in `host_presence_test`, togglable via `--no-artifact-filter`. Drops donor-hom markers whose donor-absent allele shows extreme strand skew (minor strand <10%, by effect size not p-value — a p-value over-drops at high depth), soft-clip bias (|SCBZ|>3), or read-position bias (|RPBZ|>6). `parse_vcf` captures DP4/RPBZ/SCBZ/BQBZ; `host_artifact_filtered` is in the TSV. Auto-drops 7676483 + its intron-3 neighbours, no hardcoded blacklist, generalises to new panels. **This is the cheap, control-free stopgap; the principled fix is the per-site empirical null (Step 28), which is fraction-preserving.**

**Effect (run9 → run10, `output/validation_run10/`).** Two host-presence detections were artifact-driven and correctly flip to not-detected: BHOA (p 1.3e-4 → 0.16) and PCAH-TP2 (`14_MO_IDH_APM5`, p 1.9e-13 → 1.0; 2 of its 3 markers were intron-3 artifacts, the only real one is in IDH1). Real detections preserved (QUDO-TP2, RCAR, BCOL, PCAH-TP1). Chimerism MLE unaffected (filter is host-presence-only).

**Caveat for any future write-up:** a deep-research literature pass found that LOH/CNV distorting chimerism is textbook (mostly STR), the feature version exists only narrowly (STR case series + NGS HLA-loss/6p assays), and per-marker LOH at non-HLA driver loci in an SNP chimerism panel *appears* under-explored — but that is moot here since the signal was an artifact, and the search was web-based ("appears novel" ≠ proven). **If revisited, do NOT cite the two claims that were refuted in verification: the "12% vs 2%" HLA-loss discordance figure, and "HLA-CLN is a general-LOH precedent."** The full original hypothesis/mechanism/literature write-up is in git history (pre-2026-06-01 cleanup) if needed.

## Step 28: Per-Site Empirical Error Null (Panel of Normals) 🔲 (the principled background fix)
The Step 27 bias filter is a stopgap (drops whole loci, heuristic thresholds). The principled fix is the per-site empirical error table (`error_rates.estimate_error_rates`), already consumable via `monitor --error-table` (`detect._resolve_e_per_marker`). run9/run10 ran without one (`host_err_source=global-fallback`), so a 2–3% background site like 7676483 trips against the global floor. The table is calibration not filtering: fraction-preserving, no per-read judgment, down-weights (not drops) noisy loci, and regenerates from controls when the panel changes.

- [ ] **Build the table** with `estimate-errors` on host-free samples piled through the SAME admix path (forced `bcftools mpileup`), NOT GATK (GATK reassembly hides the artifact). Host-free cohort = fully-reconstituted pure-donor timepoints PNOL (`6_MO`), GBRI (`30_MO`), QUDO-TP1 (`5_MO`); these emit the 7676483 background at full strength with zero host, so the table will learn it.
- [ ] **Wire it in** (`--error-table` on `monitor` / `run_csv_batch`); confirm `host_err_source` → `per-site`/`mixed` and the 7676483 calls collapse even with the bias filter off.
- [ ] **Re-check NDAD/BHOA GoF (Step 24)** — expect REVIEW to clear once artifact markers are down-weighted by their own background.
- [ ] Then the bias filter + any 7676483 blacklist are safety belts (no-controls fallback for new panels), not the mechanism.

## Step 29: Host CNV / CN-LoH in the Simulator (issue #13) ✅ (2026-06-03)
The recipient clone often carries somatic copy-number changes the diploid VAF model ignores. Added a copy-number-aware mixture to `simulate.py` (`HostAberration`, `cn_weighted_vaf`, `assign_cnv_aberrations` for cnloh/deletion/gain) and a both-directions minor-component LoD sweep (`paper/scripts/run_cnv_loh_validation.py` + `plot_cnv_loh.py`, Fig 8; Snakefile `cnv_loh_validation`/`cnv_loh_plot`). **Key finding:** relapse detection (recipient = minor component) is ~unaffected because the host-hom detector markers are CN-LoH-immune, but donor detection in mixed chimerism degrades badly with aberration burden (baseline ~0.5–1.1% → >20% undetectable by 25% burden; deletion worst, gain mildest). Fix shipped: median/MAD robust refit (`_robust_refit`; `estimate_*` `robust={off,auto,force}`, CLI `--robust`/`--robust-k` on monitor/timeline, default `auto`, threaded via `analysis.analyse_sample`). Gated so clean samples are byte-identical, floored at `ROBUST_MIN_MARKERS`, flagged REVIEW above `ROBUST_REVIEW_FRACTION`=0.15. Recovers gain (4.4%→0.7%) and low-burden deletion/CN-LoH LoD; cannot rescue ≥25% burden (flagged instead). Still TODO: per-marker allele-balance QC flag (overlaps Step 14), host-depth CN profile, robust-loss MLE, corrupting the host reference GT.

---

## Notes / gotchas (2026-06-01)

- **run10 is the current-code canonical validation batch** (`output/validation_run10/batch.tsv`), from `scripts/run_csv_batch.py` (CSV-driven, not Snakemake), filter on by default. run9 predates current `main`, so its `donor_pct`/`n_informative` differ slightly — code evolution, NOT the filter (which touches only host-presence).
- **Presence/comparison plots need `Donor` + `Chimerism result TP2` flow columns that `run_csv_batch` does not emit** (run9 had them joined from `Chimerism project patient list_run2.xlsx`). For run10 they were merged ad-hoc into `output/validation_run10/batch_flow.tsv`. GAP: no scripted join — either add one or have `run_csv_batch` optionally merge the xlsx columns.
- **Presence-plot regen (run10):** `plot_chimerism_comparison.py output/validation_run10/batch_flow.tsv --compare-tsv .../run2 .../run3 --labels run2 run3 run10 --flow-column "Chimerism result TP2" --label-code --output output/run2_run3_run10_presence.png`. `--label-code` is required or the x-axis shows full sample IDs instead of patient codes.
- **Only `output/run1_vs_run2.png` was ever sent to colleagues.** Every other comparison/presence/manhattan PNG is internal and regenerable.
- **Ad-hoc-script gotcha:** panel VCF sample column order is NOT consistently (host, donor) across patients — select by name from the CSV `sample_type` (as `run_csv_batch` does). A column-index script swaps some patients (e.g. BHOA is donor, host) and reports a ~99% "host fraction".
