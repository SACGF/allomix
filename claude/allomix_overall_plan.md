# Donor Chimerism Tool — Overall Plan

Master plan for allomix, a general-purpose, panel-agnostic NGS chimerism monitoring tool.

Completed steps are summarised to the decision + rationale; the full implementation detail lives in the code, tests, `doc/`, and the per-step detail docs (`claude/*_plan.md`). Open work keeps enough context to be actioned directly.

**How this file is ordered (2026-06-10 reprioritisation).** Tasks are listed in priority order, grouped into tiers. The `Step N` labels are stable identifiers (referenced from code comments, commits, and issues), NOT a running order, so they are deliberately out of numeric sequence. A `=== do version 2 of the paper here ===` marker shows where the line for the next colleague-facing paper draft falls. A `=== not worth doing ===` section at the bottom holds tasks judged not worth the effort, with enough detail kept to revisit if that judgement changes.

**The governing finding for this reprioritisation:** on real data the limiting factors at low fraction are overdispersion (Step 21) and, for co-pooled panels, the contamination floor and bad sites, NOT the per-marker refinements the older detailed plans chased. Per-marker bias correction does not move the estimate (Step 31), and the same logic applies to per-marker GQ and per-base-quality weighting (now in "not worth doing"). The principled lever that does attack the dominant background noise is the per-site empirical error null (Step 28).

---

## Remaining paper-edit work (from full re-read of the comments)
- ✅ **Contamination-floor figure for SRP434573 (#19). DONE.** Co-pooled floor as a dose-response
  boxplot, **Supplementary Figure S13**. `probe_contam_median_srp434573.py` now also writes a
  long-format `output/facts/srp_contam_persite.csv` (per pooled dp>=500 site: `n_carriers`,
  `n_alleles` het=1/hom=2 dose, `minor_frac`); the six `srp_contam.csv` medians are unchanged.
  New `plot_srp_contam.py` renders the carrier-count x-axis (0..5) boxplot to
  `figS13_srp_contam.png` (log y, per-bin median line, n_sites per box). Dave picked carrier-count
  over allele-dose after reviewing both; the allele-dose view (0..8+) stays behind a default-off
  `--allele-dose` flag, not in the build. Wired as Snakefile rules `srp_contam_facts` (now emits
  the per-site CSV too) + new `srp_contam_plot`, added to `ALL_FIGS`, captioned as S13 in
  `supplementary.md` (after S12), and cross-referenced from the Results contamination paragraph in
  `results.md`. Lint + `snakemake -n` dry run pass; full `vibepaper build` not yet run, so the S13
  caption Jinja is rendered-by-inspection only (all `srp_contam.*` keys already used in Results).
- ✅ **Alternative CI view for Figure 4 (#21). DONE.** `plot_srp434573.py:plot_logy()` (already
  existed) now writes `output/facts/figS12_srp434573_logy.png` (MLE + presence each with 95% CI vs
  known, grouped by mixture). Wired into the Snakefile as rule `srp434573_logy_plot`, added to
  `ALL_FIGS`, embedded as **Figure S12** in `supplementary.md`, and cross-referenced from the
  Figure 4 caption in `results.md`.
- ✅ **Quantify "markers used" on SRP434573 (#14). DONE.** The per-sample `n_used` count was already
  written to `srp434573_two_person.tsv`, so no rerun was needed. Added `markers_used_median`/`_min`/
  `_max` facts in `generate_srp434573_facts.py` (median 534, range 493-595) and cited them in the
  real-mixtures paragraph of `results.md`.
- **Real-data prominence in Methods (#13, partial/optional).** Results now leads with LOD and puts
  the real-mixture section third; the Methods "Real-data validation dataset" subsection is still
  last (after Simulation / LoD / CNV). Kept there as a logical build-up. If Dave wants real data
  more prominent in Methods too, move that subsection ahead of the CNV stress-test methods.


# Tier A: Paper version 2 (do now)

Version 2 is a colleague-facing draft, not an external submission, so it should be built almost entirely from results we already have. The work here is writing, figures, and one correctness fix, not new algorithm development. We can always do a version 3 once the Tier B accuracy work lands.

## Step 18: Paper v2 🟡 IN PROGRESS
Method paper (vibepaper), eventual target JMD Technical Advance. In-silico validation, multi-donor, and simulation-calibration supplementary figures are done. Cite Crysup & Woerner 2022 + Vynck.

The original v2 priority items are all done: joint-calling claim fixed (two-phase design in `methods.md`/`discussion.md`), real-data SRP434573 section added to `results.md`, bias-stability figure embedded (Fig S9), and the binomial baseline added to the ablation (Fig S4). Dave's 2026-06-10 review pass is also applied; the only remaining paper edits are the figures/facts in "Remaining paper-edit work" at the top of this file.

Note on LoD honesty: the expensive overdispersion-recalibrated LoD re-run (Step 21) is NOT required for v2. `discussion.md` already states the in-silico LoD reflects near-binomial sampling and explains the gap to real overdispersion (Supp S7/S8), which is an honest enough caveat for a colleague draft.

=== do version 2 of the paper here ===

Everything below waits for version 3 or later. The Tier B items change real-data numbers, so the natural v3 is "rerun the real-data section after Step 28 lands and the REVIEW samples clear."

---

# Tier B: Next, highest-value accuracy work (toward v3)

## Step 28: Per-Site Empirical Error Null (Panel of Normals) 🔲 (the principled background fix)
The single highest-value algorithm change. The Step 27 bias filter is a stopgap (drops whole loci, heuristic thresholds). The principled fix is the per-site empirical error table (`error_rates.estimate_error_rates`), already consumable via `monitor --error-table` (`detect._resolve_e_per_marker`). run9/run10 ran without one (`host_err_source=global-fallback`), so a 2–3% background site like 7676483 trips against the global floor. The table is calibration not filtering: fraction-preserving, no per-read judgment, down-weights (not drops) noisy loci, and regenerates from controls when the panel changes.

- [ ] **Build the table** with `estimate-errors` on host-free samples piled through the SAME admix path (forced `bcftools mpileup`), NOT GATK (GATK reassembly hides the artifact). Host-free cohort = fully-reconstituted pure-donor timepoints PNOL (`6_MO`), GBRI (`30_MO`), QUDO-TP1 (`5_MO`); these emit the 7676483 background at full strength with zero host, so the table will learn it.
- [ ] **Wire it in** (`--error-table` on `monitor` / `run_csv_batch`); confirm `host_err_source` → `per-site`/`mixed` and the 7676483 calls collapse even with the bias filter off.
- [ ] **Re-check NDAD/BHOA GoF (Step 24)** — expect REVIEW to clear once artifact markers are down-weighted by their own background.
- [ ] Then the bias filter + any 7676483 blacklist are safety belts (no-controls fallback for new panels), not the mechanism.

## Step 24: Overdispersion / REVIEW Samples (NDAD, BHOA, PCAH) 🔲
These come up `gof_pval = 0.0000` (QC=REVIEW): the chimerism fraction matches flow, but residual per-marker variance exceeds beta-binomial expectation. The TP53 artifact (Step 27) was one contributor and is now filtered. **Likely resolved by the per-site error null (Step 28)**, which down-weights each noisy locus by its own measured background rather than inflating a global `rho`. Mostly a verification check that follows directly from Step 28. If a gap remains, refit `rho` via `scripts/diagnose_sample.py` or move to per-marker-type overdispersion (Step 21).

## Step 31: Caller-mismatch warning between panel and admix VCFs 🔲 (cheap, self-contained)
A bias table estimated from GATK-called panel het sites and applied to `bcftools mpileup` admix data makes results worse, because per-marker bias is caller-specific (issue #11). Cheap guard: detect the likely source of each VCF (`DP4`/`I16` present = mpileup, absent with GATK annotations = GATK) and warn on a mismatch. Small, independent of the cohort work, prevents a real foot-gun. **[likely]**

Context: per-marker amplification bias on this panel is real and reproducible but is not the low-fraction limiting factor, and correcting it (now logit-space, issue #20) does not move the estimate. Global mean bias is ~0, the spread is already absorbed by the jointly-fit overdispersion `rho`, and the fully-informative low-fraction markers are homozygous in both contributors so their bias cannot be measured from a single pair. Per-marker bias correction is therefore deprioritized as a low-fraction accuracy lever (kept, harmless, small interior gain). **[data]**

## Step 33: Cohort / multi-sample phase 🔲
Several follow-ups block on a multi-sample entry point that `monitor` (single patient) and `estimate-bias` (single pair) do not provide. Consolidated in `cohort_phase_plan.md`. Payload 1 (bad-site blacklist) is the highest-value piece and is `[data]`-backed: bad loci directly contaminate the low-fraction signal.
- **Cohort-recurrence bad-site detection:** a panel-level blacklist for loci systematically inconsistent across the repeated cohort (Observation 4). **[likely]**
- **Batch / run-level contamination QC:** group admix samples by `##allomixRunUnit` and flag a whole flowcell lane when its samples share an elevated floor (Observation 3). **[likely]**
- **`estimate-bias --both-het` cohort entry point:** a pair's both-het markers only help *other* pairs, so the pooled table (the only way to reach the homozygous-in-both fully-informative markers) needs the multi-sample entry point (Observation 7). **[likely]**

## Step 19: Intronic Shoulder Marker Evaluation 🔲 (needs a /tau analysis first)
Our capture panel's depth extends past exon boundaries into flanking introns, a source of extra high-heterozygosity (near-0.5 MAF) markers → tighter estimate, lower LOD. **Open risk to check first: read-end mapping bias.** These SNPs sit near read ends, so ALT reads carry a terminal mismatch and are preferentially soft-clipped / MAPQ-penalised, dropping ALT from AD and skewing VAF toward REF. This is allele-asymmetric and NOT caught by the (allele-blind) depth filter. Analysis (summary-stats script against /tau, no coordinates/IDs out): per-marker median het VAF, depth, dropout vs intron offset — does het VAF leave 0.5 only after depth has already fallen below QC, or while depth is still healthy? Latter → need an explicit allele-balance filter. From a 2026-05-27 design discussion. (Same read-end/soft-clip failure mode the Step 27 artifact filter keys on.)

## Step 21: Calibrate Simulator Overdispersion for Realistic LoD 🟡 IN PROGRESS (paper polish, expensive)
**Key finding (2026-05-28):** the simulator drew from a pure binomial, so its in-silico LoD (~0.13–0.32%) was optimistic by ~3–5x vs real run3 LoDs (~0.5–1%). Beta-binomial variance approaches `p(1-p)/(ρ+1)`, so effective depth caps near `ρ+1` reads — **overdispersion, not depth, is the dominant LoD control at clinical coverage.** Done: `rho` arg in the simulator; `plot_lod_saturation.py` + `run_overdispersion_lod.py` (Supp S7/S8); `scripts/diagnose_sample.py` prints per-sample fitted `rho`. The honest in-text caveat is already in `discussion.md`, so this is paper polish, not a v2 blocker.

TODO:
- [ ] Calibrate `rho` from real per-sample fits and re-run the headline `lod_validation` so it reflects real overdispersion (expensive job; warn first).
- [ ] The simulator applies one global `rho` to every marker/allele, including the near-zero donor-absent allele where overdispersion is not physical. A marker-type/allele-aware model is needed before `rho` can validate host-presence detection (Step 20).
- [ ] Decide whether `discussion.md`'s headline LoD switches to the overdispersion-calibrated number.

## Step 30: Per-marker contamination term for the presence test 🔲 (strongest sub-1% lever, speculative)
The contamination floor currently feeds the presence background and reported LoD as a per-sample scalar (Observation 2, done). A donor-homozygous presence marker where a co-pooled genome carries the donor-absent allele is the one whose count is actually inflated; apportion the floor to those specific markers instead of taxing every donor-absent marker. With cohort genotypes this is identifiable; in deployment co-pooled genotypes are not visible, but the dose-response means the floor can still be apportioned. This is the strongest lever for sub-1% accuracy on co-pooled panels: on SRP434573 the floor (~0.06% on F2-M1, up to ~0.2% on M3-F3) competes directly with a 0.5-1% target, so even a correct per-sample floor bounds achievable low-end accuracy. **[speculative]**

Re-check before trusting the low end: at the lowest fractions the contamination subtraction can pull the presence estimate *below* the MLE host, inverting the presence-vs-MLE ordering. On `mix_F2_into_M1` 0.5% the presence estimate (0.051%) sits below the MLE host (0.19%). The `qc.py` host-presence-vs-MLE REVIEW warning fires on `mle_host < f_host_mle` (presence above MLE), so the floor shifts its operating point at <1%. Re-check that gate's behaviour before trusting it there. **[data]**

---

# Tier C: Lower priority / needs its own investigation

## Step 29: Recalibrate the robust-drop REVIEW for the low end 🔲
With the one-sided trim in place the drop fractions are much smaller, but a bare >15% drop count still does not separate "expected low-fraction trimming" from a real CNV/LoH/genotyping problem. Tie the threshold to the host fraction, or report the dominant *reason* for the drops (direction of trimmed residuals, clustering on a chromosome arm) instead of a count. Re-measure the post-fix drop distribution before changing the threshold. **[speculative]**

## Step 32: Validate contamination-flag thresholds on clean controls / other panels 🔲
`CONTAMINATION_WARN_FRACTION` (0.2%), `CONTAMINATION_REVIEW_FRACTION` (1%), and the empirical p10 floor were all set from SRP434573 alone. The p10-floor estimator assumes a useful fraction of no-carrier sites and a particular allele-frequency spectrum, and even on clean data it reports a non-zero floor (~0.15% at 1% error, 2000x) from binomial sampling spread (shrinks as ~1/sqrt(depth)). Measure the false-positive rate on genuinely uncontaminated high-depth samples and confirm the empirical floor behaves on panels with different marker allele frequencies. Validation work, not a code change, but it gates trusting any LoD floor built on the contamination estimate in deployment. **[speculative]**

## Step 15: Per-Site Dropout Rate 🔲 (cheap, marginal)
Down-weight chronically flaky sites in the likelihood: scale each marker's log-likelihood contribution by `w_s = 1 - d_s`, where `d_s` is the per-site no-call rate. The training data already exists (`scripts/measure_panel_bias.py` reports per-marker `nocall_rate` across the cohort; ~1.6% on the rhAmpSeq panel with a long tail). Implementation is a weighted pseudo-likelihood: a new `dropout.py` mirroring `bias.py` (estimate/save/load a TSV), a `marker_dropouts` kwarg threaded through `total_log_likelihood_bb` / `_multi_bb` and the estimators, an `estimate-dropout` subcommand, and a `--dropout-table` flag (missing sites default to `d_s = 0`, backwards-compatible). Composes multiplicatively with bias correction. This is a quality-of-fit (CI) improvement, not a point-estimate improvement, and dropout is only ~1.6% on the panel, so it is low value on its own. Worth doing only because it is cheap and the per-site downweighting machinery overlaps with Step 28's panel-of-normals scan (same host-free cohort, same "down-weight a locus by its own measured behaviour" idea); fold it in alongside Step 28 rather than as a standalone effort. (Detail plan deleted 2026-06-10; this paragraph carries the essence.)

## Step 34: Lower-priority / needs-its-own-investigation 🔲
- **Input-quality QC predicting poor fits.** The v1 library (F2-M1) fit worse than v2 (M3-F3) at matched fractions. A pre-estimate input-quality check (depth uniformity across the panel, per-site dropout) might predict which samples will fit poorly and warrant manual review before the estimate is attempted. **[speculative]**
- **chrX recovery for sex-matched pairs.** The panel carries ~27 chrX amplicons that we default-drop. For a host/donor pair of the same sex they are usable informative markers; inferring sex from the data and keeping chrX when it matches adds markers, which matters most at low fractions where every informative marker counts. **[speculative]**

---

## Notes / gotchas (2026-06-01)

- **run10 is the current-code canonical validation batch** (`output/validation_run10/batch.tsv`), from `scripts/run_csv_batch.py` (CSV-driven, not Snakemake), filter on by default. run9 predates current `main`, so its `donor_pct`/`n_informative` differ slightly — code evolution, NOT the filter (which touches only host-presence).
- **Presence/comparison plots need `Donor` + `Chimerism result TP2` flow columns that `run_csv_batch` does not emit** (run9 had them joined from `Chimerism project patient list_run2.xlsx`). For run10 they were merged ad-hoc into `output/validation_run10/batch_flow.tsv`. GAP: no scripted join — either add one or have `run_csv_batch` optionally merge the xlsx columns.
- **Presence-plot regen (run10):** `plot_chimerism_comparison.py output/validation_run10/batch_flow.tsv --compare-tsv .../run2 .../run3 --labels run2 run3 run10 --flow-column "Chimerism result TP2" --label-code --output output/run2_run3_run10_presence.png`. `--label-code` is required or the x-axis shows full sample IDs instead of patient codes.
- **Only `output/run1_vs_run2.png` was ever sent to colleagues.** Every other comparison/presence/manhattan PNG is internal and regenerable.
- **Ad-hoc-script gotcha:** panel VCF sample column order is NOT consistently (host, donor) across patients — select by name from the CSV `sample_type` (as `run_csv_batch` does). A column-index script swaps some patients (e.g. BHOA is donor, host) and reports a ~99% "host fraction".

---

# === not worth doing ===

These were judged not worth the effort during the 2026-06-10 reprioritisation. The detailed implementation plans were deleted; enough is kept here to rebuild them if the judgement changes. Both are per-marker weighting refinements with small expected gains, and the governing finding (overdispersion + contamination dominate the low-fraction limit, and per-marker bias correction does not move the estimate) applies to both.

## Step 16: GQ-Weighted Marker Contributions ❌ not worth doing
Idea: replace the hard `--min-gq 20` cutoff with a per-marker likelihood weight so borderline genotypes stay informative but downweighted, using the Phred posterior `P(GT correct) = 1 / (1 + 10^(-GQ/10))` over host + donor genotypes (`w_gq = prod_c P(GT_c correct)`), with a low hard floor (GQ 10). A new `gq_weight.py` helper, precomputed per-marker weights threaded through the aggregators/estimators, and `--gq-weighted` / `--gq-weight-scheme` / `--gq-floor` CLI flags. Composes multiplicatively with Step 15 dropout weighting.

Why not: the deployment panel runs >1000x and almost every host/donor call is GQ=99, so the weight is ~1 for nearly every marker and the method degenerates to the status quo. The only real lift is recovering a handful of GQ 10–20 markers on low-coverage host samples, a small effect that does not touch the overdispersion/contamination floor that actually limits low-fraction accuracy. Revisit only if real data shows we are routinely losing informative markers to the GQ 20 cutoff on weak host genotyping.

## Step 17: Per-Base Quality-Aware Likelihood ❌ not worth doing
Idea: weight each read's likelihood contribution by its base quality (Conpair-style), via a per-marker effective error rate `e = 10^(-meanPhred/10)` derived from `FORMAT/QS`. Tiers T1 (per-marker aggregated error) default, T2 (per-allele) behind a flag, T3 (per-read) out of scope.

Why not: it is the most invasive remaining algorithm step for the smallest expected gain. The production VCF carries no per-base quality, so it requires a new upstream pipeline stage (`bcftools mpileup -a FORMAT/QS` + `bcftools annotate` at the panel sites, a second BAM pass), a new `bq.py` module, simulator changes (O(depth) per marker), and CLI/report plumbing. On a Q30+ panel the per-marker effective error barely differs from the flat `--error-rate`, and the mean-phred-to-error conversion underestimates `e` under Jensen when BQs are heterogeneous. The dominant low-fraction noise is overdispersion and contamination, neither of which this addresses, and Step 28 attacks the background far more cheaply. The original decision gate ("revisit only if Steps 14–16 land and real-data CIs are still wider than we want") was never met, and Steps 15/16 are themselves now deprioritised. Revisit only if, after Step 28, real-data CIs remain wider than clinically acceptable and per-base quality is shown to be the residual driver.
