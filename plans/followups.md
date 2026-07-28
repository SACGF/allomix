# Follow-ups

Repository and paper work arising from the donor-relatedness analysis. The
wetlab design itself is in `plans/wetlab_validation_design.md`.

## Paper

- **Add `"parent-child"` to the relatedness sweep.** `RELATEDNESS_LEVELS` in
  `paper/scripts/run_relatedness_validation.py:38` is
  `["unrelated", "cousin", "half-sibling", "sibling"]`, but `"parent-child"` is
  already implemented in `src/allomix/simulate.py:534`. One-line change plus a
  rerun of that single rule (not a full build). Gives Table 2 and Figure 5 a
  parent-child row, which is what the new Methods paragraph about types 0 and 1
  currently asserts without quantifying.

- **Quantify relatedness marker yield.** The Methods text now says full siblings
  retain "a small and variable number" of full-contrast markers. That is
  deliberately qualitative because the supporting counts are scratchpad
  simulations, not reproducible facts. If the paper is to give numbers, they need
  to come from a `paper/scripts/` fact like every other number in the text.

- **Figure 4 relatedness rows.** `plot_lod_grid.py:41` plots
  `["unrelated", "sibling"]` only. Worth deciding whether parent-child belongs
  there too once the sweep produces it, given the two relationships fail in
  different ways (siblings have fewer informative markers overall; parent-child
  deterministically has no full-contrast markers).

## Code

- **QC verdict does not flag a zero-marker presence test.** `qc.py:554` gates the
  host-presence REVIEW on `hp.n_markers > 0`, so a run where the donor set covers
  the recipient entirely can return an overall PASS with the presence readout
  silently absent and no message saying why. The report layer does the right
  thing already (`report.py:127-141` writes `NA` rather than the internal
  `lrt_pval=1.0`), so this is about the verdict and messaging, not the numbers.
  Candidate fix: REVIEW, or at least an informational note, when the presence
  detector selects zero markers from a non-empty informative set. Worth raising
  as a GitHub issue.

- **Pre-existing lint.** `ruff check` reports two `I001` import-order errors in
  `scripts/validate_contamination_lod_floor.py`. Present before this work and
  left alone; both are auto-fixable with `ruff check --fix`.

## Repository hygiene

- **Restore or retire the `claude/` design docs.** Nine references across
  `src/`, `tests/` and `scripts/` pointed at `claude/*.md` files that are not in
  the repository, not tracked, and not gitignored: `20_host_presence_detection_plan.md`,
  `further_improvements.md`, `step30_design.md`, `srp434573_figure_review_findings.md`.
  The dangling pointers have been removed from the source, but the rationale they
  held is genuinely gone. Some of it (the presence detector's acceptance gates,
  the contamination-table design decisions, the one-sided trim motivation) is
  worth reconstructing into `plans/` or `docs/` if any copies survive locally.

- **Track `plans/`.** The `claude/` case is the argument: docstrings that cite
  untracked design docs become dead references as soon as the working copy moves
  or is cleaned. Keep `plans/` in git rather than gitignored.
