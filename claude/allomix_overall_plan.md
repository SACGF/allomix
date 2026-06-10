# Donor Chimerism Tool — Overall Plan

Master plan for allomix, a general-purpose, panel-agnostic NGS chimerism monitoring tool.

Completed steps are summarised to the decision + rationale; the full implementation detail lives in the code, tests, `doc/`, and the per-step detail docs (`claude/*_plan.md`). Open work keeps enough context to be actioned directly.

---

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

## Step 21: Calibrate Simulator Overdispersion for Realistic LoD 🟡 IN PROGRESS
**Key finding (2026-05-28):** the simulator drew from a pure binomial, so its in-silico LoD (~0.13–0.32%) was optimistic by ~3–5x vs real run3 LoDs (~0.5–1%). Beta-binomial variance approaches `p(1-p)/(ρ+1)`, so effective depth caps near `ρ+1` reads — **overdispersion, not depth, is the dominant LoD control at clinical coverage.** Done: `rho` arg in the simulator; `plot_lod_saturation.py` + `run_overdispersion_lod.py` (Supp S7/S8); `scripts/diagnose_sample.py` prints per-sample fitted `rho`.

TODO:
- [ ] Calibrate `rho` from real per-sample fits and re-run the headline `lod_validation` so it reflects real overdispersion (expensive job; warn first).
- [ ] The simulator applies one global `rho` to every marker/allele, including the near-zero donor-absent allele where overdispersion is not physical. A marker-type/allele-aware model is needed before `rho` can validate host-presence detection (Step 20).
- [ ] Decide whether `discussion.md`'s headline LoD switches to the overdispersion-calibrated number.

## Step 24: Overdispersion / REVIEW Samples (NDAD, BHOA, PCAH) 🔲
These come up `gof_pval = 0.0000` (QC=REVIEW): the chimerism fraction matches flow, but residual per-marker variance exceeds beta-binomial expectation. The TP53 artifact (Step 27) was one contributor and is now filtered. **Likely resolved by the per-site error null (Step 28)**, which down-weights each noisy locus by its own measured background rather than inflating a global `rho`. If a gap remains after Step 28, refit `rho` via `scripts/diagnose_sample.py` or move to per-marker-type overdispersion (Step 21).

## Step 28: Per-Site Empirical Error Null (Panel of Normals) 🔲 (the principled background fix)
The Step 27 bias filter is a stopgap (drops whole loci, heuristic thresholds). The principled fix is the per-site empirical error table (`error_rates.estimate_error_rates`), already consumable via `monitor --error-table` (`detect._resolve_e_per_marker`). run9/run10 ran without one (`host_err_source=global-fallback`), so a 2–3% background site like 7676483 trips against the global floor. The table is calibration not filtering: fraction-preserving, no per-read judgment, down-weights (not drops) noisy loci, and regenerates from controls when the panel changes.

- [ ] **Build the table** with `estimate-errors` on host-free samples piled through the SAME admix path (forced `bcftools mpileup`), NOT GATK (GATK reassembly hides the artifact). Host-free cohort = fully-reconstituted pure-donor timepoints PNOL (`6_MO`), GBRI (`30_MO`), QUDO-TP1 (`5_MO`); these emit the 7676483 background at full strength with zero host, so the table will learn it.
- [ ] **Wire it in** (`--error-table` on `monitor` / `run_csv_batch`); confirm `host_err_source` → `per-site`/`mixed` and the 7676483 calls collapse even with the bias filter off.
- [ ] **Re-check NDAD/BHOA GoF (Step 24)** — expect REVIEW to clear once artifact markers are down-weighted by their own background.
- [ ] Then the bias filter + any 7676483 blacklist are safety belts (no-controls fallback for new panels), not the mechanism.

## SRP434573 run follow-ups (issue #16)

Open items from the real-data run on SRP434573 (HiSeq 3000 MIP panel, ~1050 single-strand sites, known two/three-person mixtures titrated 10% -> 0.5%, minor contributor mapped to HOST so it reads as a declining-chimerism series). Confidence marked **[data]** (measured here), **[likely]** (well-supported inference), **[speculative]** (needs its own check). The wave-1 fixes from this run (one-sided robust trim, pre-trim GoF, contamination floor into LoD + presence background) are already implemented and validated in silico.

## Step 29: Recalibrate the robust-drop REVIEW for the low end 🔲
With the one-sided trim in place the drop fractions are much smaller, but a bare >15% drop count still does not separate "expected low-fraction trimming" from a real CNV/LoH/genotyping problem. Tie the threshold to the host fraction, or report the dominant *reason* for the drops (direction of trimmed residuals, clustering on a chromosome arm) instead of a count. Re-measure the post-fix drop distribution before changing the threshold. **[speculative]**

## Step 30: Per-marker contamination term for the presence test 🔲
The contamination floor currently feeds the presence background and reported LoD as a per-sample scalar (Observation 2, done). A donor-homozygous presence marker where a co-pooled genome carries the donor-absent allele is the one whose count is actually inflated; apportion the floor to those specific markers instead of taxing every donor-absent marker. With cohort genotypes this is identifiable; in deployment co-pooled genotypes are not visible, but the dose-response means the floor can still be apportioned. This is the strongest lever for sub-1% accuracy on co-pooled panels: on SRP434573 the floor (~0.06% on F2-M1, up to ~0.2% on M3-F3) competes directly with a 0.5-1% target, so even a correct per-sample floor bounds achievable low-end accuracy. **[speculative]**

Re-check before trusting the low end: at the lowest fractions the contamination subtraction can pull the presence estimate *below* the MLE host, inverting the presence-vs-MLE ordering. On `mix_F2_into_M1` 0.5% the presence estimate (0.051%) sits below the MLE host (0.19%). The `qc.py` host-presence-vs-MLE REVIEW warning fires on `mle_host < f_host_mle` (presence above MLE), so the floor shifts its operating point at <1%. Re-check that gate's behaviour before trusting it there. **[data]**

## Step 31: Caller-mismatch warning between panel and admix VCFs 🔲
A bias table estimated from GATK-called panel het sites and applied to `bcftools mpileup` admix data makes results worse, because per-marker bias is caller-specific (issue #11). Cheap guard: detect the likely source of each VCF (`DP4`/`I16` present = mpileup, absent with GATK annotations = GATK) and warn on a mismatch. Small, self-contained, independent of the cohort work. **[likely]**

Context: per-marker amplification bias on this panel is real and reproducible but is not the low-fraction limiting factor, and correcting it (now logit-space, issue #20) does not move the estimate. Global mean bias is ~0, the spread is already absorbed by the jointly-fit overdispersion `rho`, and the fully-informative low-fraction markers are homozygous in both contributors so their bias cannot be measured from a single pair. Per-marker bias correction is therefore deprioritized as a low-fraction accuracy lever (kept, harmless, small interior gain). **[data]**

## Step 32: Validate contamination-flag thresholds on clean controls / other panels 🔲
`CONTAMINATION_WARN_FRACTION` (0.2%), `CONTAMINATION_REVIEW_FRACTION` (1%), and the empirical p10 floor were all set from SRP434573 alone. The p10-floor estimator assumes a useful fraction of no-carrier sites and a particular allele-frequency spectrum, and even on clean data it reports a non-zero floor (~0.15% at 1% error, 2000x) from binomial sampling spread (shrinks as ~1/sqrt(depth)). Measure the false-positive rate on genuinely uncontaminated high-depth samples and confirm the empirical floor behaves on panels with different marker allele frequencies. Validation work, not a code change, but it gates trusting any LoD floor built on the contamination estimate in deployment. **[speculative]**

## Step 33: Cohort / multi-sample phase 🔲
Several follow-ups block on a multi-sample entry point that `monitor` (single patient) and `estimate-bias` (single pair) do not provide. Consolidated in `cohort_phase_plan.md`:
- **Batch / run-level contamination QC:** group admix samples by `##allomixRunUnit` and flag a whole flowcell lane when its samples share an elevated floor (Observation 3). **[likely]**
- **Cohort-recurrence bad-site detection:** a panel-level blacklist for loci systematically inconsistent across the repeated cohort (Observation 4). **[likely]**
- **`estimate-bias --both-het` cohort entry point:** a pair's both-het markers only help *other* pairs, so the pooled table (which is the only way to reach the homozygous-in-both fully-informative markers) needs the multi-sample entry point (Observation 7). **[likely]**

## Step 34: Lower-priority / needs-its-own-investigation 🔲
- **Input-quality QC predicting poor fits.** The v1 library (F2-M1) fit worse than v2 (M3-F3) at matched fractions. A pre-estimate input-quality check (depth uniformity across the panel, per-site dropout) might predict which samples will fit poorly and warrant manual review before the estimate is attempted. **[speculative]**
- **chrX recovery for sex-matched pairs.** The panel carries ~27 chrX amplicons that we default-drop. For a host/donor pair of the same sex they are usable informative markers; inferring sex from the data and keeping chrX when it matches adds markers, which matters most at low fractions where every informative marker counts. **[speculative]**

## Notes / gotchas (2026-06-01)

- **run10 is the current-code canonical validation batch** (`output/validation_run10/batch.tsv`), from `scripts/run_csv_batch.py` (CSV-driven, not Snakemake), filter on by default. run9 predates current `main`, so its `donor_pct`/`n_informative` differ slightly — code evolution, NOT the filter (which touches only host-presence).
- **Presence/comparison plots need `Donor` + `Chimerism result TP2` flow columns that `run_csv_batch` does not emit** (run9 had them joined from `Chimerism project patient list_run2.xlsx`). For run10 they were merged ad-hoc into `output/validation_run10/batch_flow.tsv`. GAP: no scripted join — either add one or have `run_csv_batch` optionally merge the xlsx columns.
- **Presence-plot regen (run10):** `plot_chimerism_comparison.py output/validation_run10/batch_flow.tsv --compare-tsv .../run2 .../run3 --labels run2 run3 run10 --flow-column "Chimerism result TP2" --label-code --output output/run2_run3_run10_presence.png`. `--label-code` is required or the x-axis shows full sample IDs instead of patient codes.
- **Only `output/run1_vs_run2.png` was ever sent to colleagues.** Every other comparison/presence/manhattan PNG is internal and regenerable.
- **Ad-hoc-script gotcha:** panel VCF sample column order is NOT consistently (host, donor) across patients — select by name from the CSV `sample_type` (as `run_csv_batch` does). A column-index script swaps some patients (e.g. BHOA is donor, host) and reports a ~99% "host fraction".
