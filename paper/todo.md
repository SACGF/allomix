# To do

Open work on the paper. Delete an item when it is done. How to work on the repo is in
`CLAUDE.md` and `docs/paper.md`.

## Needs a full pipeline run

- [ ] **Regenerate the full facts and provenance sidecars** (#56 item 1).
  Locally, lint notes the missing `lod_headline`, `presence_lod_curve_headline`,
  `subsample_lod_headline`, `supp_synthetic`, `sim_calibration`, `panel_specs`,
  `panel_empirical` and `bias_*` facts, and `vibepaper check` reports 15 asserts as
  MISSING until they exist.
- [ ] **Confirm the two `near` asserts** (Discussion and Figure S8). They compare the LoD
  at the fitted ρ with the measured subsample LoD
  (`subsample_lod_headline.mle_lod_1000x_100markers_pct`), replacing a hand-typed "~1%".
  The local overdispersion value is 0.82%; only a full build shows whether the pair is
  within a factor of two.

## Pipeline

- [ ] **`bias_with_bias` depends on rule order.** `generate_paper_facts.py` (rule
  `paper_facts`) loads `output/bias_comparison/vcfs/true_biases.tsv` if it exists, but
  that file is written by `compare_bias_correction.py` (rule `bias_comparison`) and is
  not an input of `paper_facts`. So `bias_with_bias.*` changes with which rule ran first
  (a quick build gave `ci_coverage_pct` 100.0 in parallel, 91.7 when rerun after
  `bias_comparison`). Declare the dependency, or have `paper_facts` write its own table.
- [ ] **`panel_empirical` reads its one-row CSV with `csv.DictReader`**, not
  `vibepaper.facts.read_facts`, because vibepaper 0.13's `read_facts` upcasts every
  integer in a horizontal CSV to float (`1732` → `1732.0`). Switch to `read_facts` once
  vibepaper fixes it (report upstream), or make
  `paper/empirical_results/panel_empirical.csv` vertical.

## Unchecked wording

- [ ] Results, Figure 2: "the estimate floors near the dataset's co-pooled contamination
  level" has no fact for the 0.5%-rung estimate, so nothing checks it (lint doesn't flag
  "near").

## For the author

- [ ] **Kakodkar 0.3%:** `discussion.md` says "Kakodkar et al. reported an analytical LoD
  of 0.3% for AlloSeq HCT", but the Table 1 footnote calls 0.3% the vendor-stated figure.
  Check against the source.
- [ ] **Unused facts:** decide whether to cite or list in `[lint] unused_ok`:
  `lod_saturation.*` (Figure S7 fit parameters) and
  `srp434573_synthetic.{real,syn}_05_mle_med_pct` / `syn_minus_real_05_mle_pct`.
  `val_host_*` / `validation_summary` are in `unused_ok` but nothing reads them: candidates
  for removal from `generate_paper_facts.py`.
- [ ] **Word limits:** add a `[limits]` table to `paper.toml` once the target journal is
  chosen (`vibepaper wordcount` prints the counts).
