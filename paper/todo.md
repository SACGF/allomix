# To do

Open work on the paper. Delete an item when it is done. How to work on the repo is in
`CLAUDE.md` and `docs/paper.md`.

## Pipeline

- [ ] **`panel_empirical` reads its one-row CSV with `csv.DictReader`**, not
  `vibepaper.facts.read_facts`, because vibepaper 0.13's `read_facts` upcasts every
  integer in a horizontal CSV to float (`1732` → `1732.0`). Switch to `read_facts` once
  vibepaper fixes it (SACGF/vibepaper#56), or make
  `paper/empirical_results/panel_empirical.csv` vertical.

## For the author

- [ ] **Word limits:** the target is JMD (the CSL in `paper.toml`), but `paper.toml` has no
  `[limits]` table yet. Add one from the JMD Regular Article limits (the author guide at
  jmdjournal.org / sciencedirect.com blocks automated fetches, so the limits need pasting
  in). `vibepaper wordcount` prints the current counts.
