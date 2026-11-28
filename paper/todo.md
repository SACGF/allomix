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

- [ ] **JMD conformance** (rules in `paper/jmd_guide.md`). Points where the draft may not
  follow the guide:
  - **Supplementary Methods (S1 to S11):** the guide says method detail in supplementary
    material is not acceptable. Decide whether to move them into Methods or keep them and
    risk an editor asking.
  - **First person:** the guide asks for minimal "we/our/us" and prefers the passive.
    There are about 45 uses (most in Methods, Results and Discussion).
  - **Abstract at the limit:** exactly 220 of 220 words, so any added word fails
    `vibepaper wordcount --check`.
  - **AI disclosure:** a "Declaration of AI and AI-assisted technologies in the writing
    process" section is needed if AI tools were used in writing.
  - **Article type:** Regular Article or Technical Advance (same limits and charges).
