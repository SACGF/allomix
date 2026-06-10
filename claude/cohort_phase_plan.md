# Cohort / multi-sample phase plan

Three improvements from the SRP434573 run (`further_improvements.md`, Observations
3, 4, and the pooled-bias-table part of 7) all block on the same missing piece: a
multi-sample entry point. `allomix monitor` takes one patient (host + donors +
serial admix samples); `allomix estimate-bias` takes one host/donor pair plus its
admix VCFs. Neither can aggregate across patients or across a sequencing batch.

Building that entry point once unlocks all three payloads, so they belong in one
phase rather than three separate efforts. This file is the plan; the evidence for
each payload lives in `further_improvements.md`.

## The shared dependency: a cohort entry point

Decide the interface first, because all three payloads consume it:

- **Input.** A manifest (TSV or small config) listing samples: for each, the role
  (host / donor / admix), the patient or pair it belongs to, the VCF path, and the
  run unit (already in the `##allomixRunUnit` header). A manifest is more honest than
  accepting many positional args and is reusable across the three payloads.
- **What it reads.** Per-sample informative markers and admix allele counts (same
  `genotype.classify_markers` path `monitor` uses), plus the per-sample contamination
  estimate (`contamination.py`) and run-unit metadata (`runmeta.py`).
- **Where it lives.** A new `allomix cohort` subcommand (or a thin module that the
  three payloads call), parallel to `monitor`. Keep `monitor` single-patient; do not
  overload it.

Confirm the current single-sample call paths before designing (they may have moved):
`cli.cmd_monitor`, `analysis.analyse_sample`, `cli.cmd_estimate_bias`.

## Payload 1: cohort-recurrence bad-site detection (Observation 4)

About 44 site-instances across the 10 two-person mixtures sat at ~42% minor allele
against a ~0.004% error floor and ~0.2% contamination, carrying ~92% of the pooled
minor reads. These are genotype miscalls or mapping artifacts at specific loci, and
because SRP434573 is the same 7 people across ~64 runs, a bad locus is bad in every
sample. **[data]**

- Use the repeated cohort to find loci that are systematically inconsistent (reference
  VAF contradicting the call, extreme strand/clip/position bias, or minor-allele
  behaviour that matches no contributor) and emit a panel-level blacklist.
- Both this panel (~1050 sites) and the deployment panel (76 rhAmpSeq SID SNPs) are
  small and fixed, so a one-time cohort QC pass is cheap and reusable.
- Output a blacklist file that `monitor` can consume (drop those loci panel-wide)
  rather than rediscovering them per sample. **[likely]**

This is the highest-value payload: the bad sites directly contaminate the
low-fraction signal that wave 1 and Observation 2 are trying to protect.

## Payload 2: batch / run-level contamination QC (Observation 3)

Contamination differed several-fold between the M3-F3 (v2) and F2-M1 (v1) libraries,
and the `##allomixRunUnit` header records each admix sample's flowcell + lane.
Nothing yet aggregates across that grouping. **[data]**

- Group admix samples by run unit, report the per-run contamination distribution, and
  flag a whole flowcell lane when its samples share an elevated floor (the signature
  of a bad multiplexing run rather than one bad sample).
- This closes the loop between the two issue-#12 features: same run **and** a
  cohort-wide elevated floor is much stronger evidence of index hopping than either
  the run-unit provenance flag or the in-data contamination estimate alone. **[likely]**

## Payload 3: pooled both-het bias table (Observation 7)

`estimate-bias --both-het` currently takes one host/donor pair plus its admix VCFs,
but a pair's both-het markers only inform *other* pairs, so the table is only useful
pooled across a cohort (it covered 0 of its own pair's 576 informative markers, but
66% of a held-out pair's informative set when pooled across 11 mixtures). **[data]**

- Accept a multi-pair manifest (the shared entry point above) and build the pooled
  table directly, instead of one pair per invocation.
- Note the wave-7 finding: bias correction is a small interior-fraction refinement and
  is neutral at low fractions, so this is a refinement payload, not a sensitivity lever.
  Build it because it makes the existing feature usable, not because it is expected to
  move the low-end result. **[likely]**

## Sequencing within the phase

1. The cohort entry point (manifest reader + `allomix cohort` skeleton).
2. Payload 1 (bad-site blacklist) first: highest value, and a clean consumer of the
   entry point.
3. Payload 2 (run-level contamination QC): reuses the entry point and the existing
   contamination estimator.
4. Payload 3 (pooled bias table): lowest priority per the Observation 7 evidence.

Validation for each payload should follow the project rule (in silico first, N>=5
replicates), and the bad-site blacklist in particular should be checked against the
clean-control work in Observation 6 so it does not blacklist good loci.
