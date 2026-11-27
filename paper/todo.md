# To do

Open work on the paper. Delete an item when it is done. How to work on the repo is in
`CLAUDE.md` and `docs/paper.md`.

## Needs a full pipeline run

- [ ] **One non-quick build for real PASS/FAIL on every assert** (#56 item 1). A clean
  `quick=1 fast_grid=1` build from a fresh clone (2026-11-27) produces every facts file
  with its provenance sidecar: `vibepaper check` gives 85 asserts, 72 pass, 13 fail,
  0 missing. The 13 fail on quick-mode numbers (coarse LoD grid, few replicates), so
  only a full build shows whether they pass.
- [ ] **Confirm the two `near` asserts** (Discussion and Figure S8). They compare the LoD
  at the fitted ρ with the measured subsample LoD
  (`subsample_lod_headline.mle_lod_1000x_100markers_pct`), replacing a hand-typed "~1%".
  The local overdispersion value is 0.82%; only a full build shows whether the pair is
  within a factor of two.

## Pipeline

- [ ] **`panel_empirical` reads its one-row CSV with `csv.DictReader`**, not
  `vibepaper.facts.read_facts`, because vibepaper 0.13's `read_facts` upcasts every
  integer in a horizontal CSV to float (`1732` → `1732.0`). Switch to `read_facts` once
  vibepaper fixes it (SACGF/vibepaper#56), or make
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
