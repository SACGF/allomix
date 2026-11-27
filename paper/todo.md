# To do

Open work on the paper. Delete an item when it is done. How to work on the repo is in
`CLAUDE.md` and `docs/paper.md`.

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
- [ ] Results, real-data LoD: "in the three titrated to 0.5% recipient, both the magnitude
  and presence-test LoD reached at or below that level" has no assert, and the magnitude
  half does not hold: the `to_0.5pct` MLE LoD bottoms out at 0.51% (2000x/400 markers,
  `subsample_lod_summary.csv`); only the presence test reaches 0.5% (censored). The
  headline facts read only `mixture_set == "all"`, so a fix needs `to_0.5pct` facts from
  `run_subsample_lod.py` (a `subsample_lod` rerun) or wording without the number.

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
