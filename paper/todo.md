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

## For the author

- [ ] **Two asserts fail on the full build** (2026-11-27, non-quick, no `fast_grid`:
  85 asserts, 83 pass, 2 fail; both `near` asserts pass, LoD at the fitted ρ 0.82% vs
  subsample LoD 1.0%). The render stops on these two, so each needs a wording or
  fact decision:
  - Table 1 caption, "≤0.5% on the same mixtures subsampled to the largest panel"
    (`results.claims.toml:257`): `subsample_lod_headline.mle_lod_2000x_400markers_pct`
    is 1.0%, censored at the lowest titration of the 7 mixtures that stop at 1%. The
    3 mixtures titrated to 0.5% give 0.51% at 2000x/400 markers (`to_0.5pct` row of
    `subsample_lod_summary.csv`), so even that subset is not ≤0.5%.
  - Figure S9, "only weakly correlated (r = …)" (`supplementary.claims.toml:142`,
    assert `|r| < 0.3`): `supp_synthetic.bias_stability_r` is 0.40. It comes from the
    empirical panel, not the simulation, so quick and full builds agree.

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
