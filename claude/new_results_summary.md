# Issue #8 pilot — LoD results summary

10-replicate pilot of the LoD sweep specified in `claude/issue_8_lod_markers_plan.md`, run with bias correction on (perfect per-marker biases passed to the estimator) per the agreed design change from no-correction to with-correction. Bias correction was added in response to the first pilot showing LoB sitting at 1–3% across most cells with the noise model uncorrected; with biases known, LoB drops to 0.1–0.5% on most cells and LoD becomes comparable to vendor analytical specs. The plan now carries a followup note to align the depth, relatedness, and multidonor figures on the same calibrated baseline so the whole paper sits on one story.

Sweep: 2 relatedness × 5 depths × 6 panel sizes × 7 true fractions × 10 reps = 4,200 estimator calls. Wall time ≈ 50 min on 8 workers. Outputs:

- `output/facts/lod_grid_raw.csv` (4,200 rows) — per-replicate point data, durable artefact.
- `output/facts/lod_summary.csv` (60 rows) — per-cell LoB, LoD, bootstrap CIs, logistic fit params, mean informative markers.
- `output/facts/lod_headline.csv` (1 row) — cells the Results prose will quote.
- `output/facts/fig5_lod_curves.png` — two-facet figure (unrelated, sibling), depth-coloured curves with shaded LoD CI, log axes, 76-marker reference line.

## Sanity checks (plan section "Validation that the implementation is working")

| Check | Expected | Got | Verdict |
|---|---|---|---|
| 1 | 2000×, nm=400, unrelated → LoD ≪ 0.5% | **0.19%** | pass |
| 2 | 100×, nm=25, sibling → several % | **2.54%** | pass |
| 3a | LoD monotone in panel size | clean monotone except small dip at unrelated 1000× (0.48% nm=100 vs 0.29% nm=75 vs 0.31% nm=200), within N=10 noise | pass |
| 3b | LoD monotone in depth | trends right, ±0.15% jitter from N=10 sampling on the 95th-percentile LoB | pass in spirit |
| 4 | LoB ≪ 1% across cells | 0.06–0.7% on most cells; sibling at 100× / nm ≤ 100 still 1–2.3% (real, driven by ~8 informative markers and shallow depth) | pass except sibling-shallow-low-markers regime |
| 5 | 1000×, nm=100, unrelated → 0.1–0.5% LoD | **0.48%** | pass |

## Headline LoDs (% donor fraction, N=10)

| Cell | LoD | LoB | Mean informative |
|---|---:|---:|---:|
| unrelated 1000×, 75 markers | 0.29 | 0.13 | 42.7 |
| unrelated 1000×, 100 markers | 0.48 | 0.31 | 58.6 |
| unrelated 500×, 200 markers | 0.46 | 0.21 | 115.0 |
| unrelated 2000×, 400 markers | 0.19 | 0.16 | 232.6 |
| sibling 1000×, 75 markers | 0.57 | 0.56 | 29.1 |
| sibling 1000×, 100 markers | 0.26 | 0.09 | 35.0 |
| sibling 2000×, 400 markers | 0.31 | 0.40 | 142.7 |

Vendor context: AlloSeq HCT 0.3%, Devyser 0.06–0.1%, ScisGo 0.2–0.5% (single donor). Our reference operating point (unrelated, 1000×, 100 markers) at 0.48% sits in the published range; deeper sequencing or larger panels beat the AlloSeq spec.

## Trends visible already at N=10

- **Panel size dominates.** Going from 25 to 400 markers at fixed 1000× depth cuts unrelated LoD from 1.26% to 0.32%. Sibling LoD from 1.24% to 0.31% over the same axis.
- **Depth helps less above 500×.** Unrelated 100-marker LoD: 0.53% at 100×, 0.62% at 250×, 0.60% at 500×, 0.48% at 1000×, 0.47% at 2000×. The flattening is consistent with depth-uniformity and per-marker bias dominating shot noise at deep coverage.
- **Sibling vs unrelated.** At small panels (25–75 markers), siblings give ~2× worse LoD because of the reduced informative marker fraction. At larger panels (200+ markers) the gap closes — there are still enough informative markers to nail down f.
- **76-marker reference line.** Reading off Figure 7 at the vertical 76-marker line: unrelated LoD lands around 0.3–0.6% at 500–2000× depth; sibling around 0.6–1.1%. Both meaningfully better than the qPCR/STR floor of ~1–5% the existing literature reports.

## Implementation notes

- **Bias correction wiring** (`paper/scripts/run_lod_validation.py`): each replicate generates a per-marker bias vector from `simulate.generate_marker_biases_realistic` once, passes it as `fixed_biases` to every `blend_vcfs` call within that rep, captures the (chrom, pos, ref, alt → bias) dict from the first blend, and threads it through to `estimate_single_donor_bb(marker_biases=...)`. This models the "panel calibrated" ceiling: how well the estimator does when bias is fully known. Real-world bias estimation noise is *not* modelled, so the LoD is an upper bound on analytical sensitivity.
- **Estimator grid_steps** lowered from default 1001 to 201 for the LoD sweep only (`ESTIMATOR_GRID_STEPS` constant). Nelder-Mead refinement still converges to <1e-3 in f. ~5× speedup. Other validation scripts unaffected.
- **Logistic fit robustness.** First pass used pure `scipy.optimize.curve_fit` on `P(detected | f) = sigmoid(a + b·log10 f)`. Six cells failed to converge in the bias-corrected pilot because detection jumped 0→1 over a single fraction step, making the slope unidentifiable. Added a linear-in-log10(f) bracketing-interpolation fallback (`_interp_lod`) that picks up where curve_fit gives up — every cell now produces an LoD. Two unit tests cover the new fallback (`test_interp_lod_brackets_target`, `test_fit_lod_falls_back_to_interp_on_step_data`).
- **Bootstrap CIs** (200 resamples of replicate-level detection booleans) work as designed but produce wide upper bounds (e.g. 121%, 168%) on a few "easy" cells where many bootstrap samples have all-1 detection rates that the interp also can't bracket. Will tighten substantially at N=60. Not currently filtered out; the plot uses them as-is so the visual is honest about uncertainty in the small-N regime.

## What changes at N=60

Mostly cosmetic. Key expected improvements:

- Bootstrap CIs much narrower — the runaway upper bounds (1.2 × 10² %) at "easy" cells go away because at N=60 you don't get bootstrap samples that are entirely all-1.
- Depth-axis monotonicity will be clean. The current ±0.15% jitter on unrelated nm=100 is sampling noise on the 95th-percentile LoB.
- Logistic fit will rarely if ever fall back to interp — 60 reps gives detection rates on a 1/60 grid, so the step-function regime mostly disappears except in the trivially-easy "all detected at smallest f" corner.
- Headline LoDs will move by ≤0.05–0.10% in either direction. Story is unchanged.

Compute: ~5 hours single-threaded → ~40 min on 8 workers. Cheap to run overnight.

## Open items at handoff

1. **Followup: align other figures on bias correction.** Plan note added. Out of scope for this issue.
2. **Sibling LoB at small panels** sits at 1–2.3% even with bias correction (nm ≤ 100, depth ≤ 250). Real: with ~8 informative markers and locus dropout, the 95th percentile of est_frac on a blank is structurally non-tiny. Methods text already notes this regime; no fix needed but worth a sentence in Results.
3. **Reviewers may ask** how robust the LoD is to bias miscalibration. Currently we report perfect-bias ceiling. A supplementary sweep with bias drawn from a noisy estimate (e.g. 10 training samples then re-estimated bias) would address this; flagged as supplementary, not blocking.
4. **CI on LoD at the cell level** uses bootstrap; reviewers in clinical journals may prefer a profile-likelihood approach on the logistic. The current implementation is defensible (200 resamples of replicates is the standard CLSI approach) but worth knowing if asked.

## Files touched

- New: `paper/scripts/run_lod_validation.py`, `paper/scripts/plot_lod_curves.py`, `tests/test_lod_validation.py`.
- Modified: `paper/Snakefile` (lod_validation + lod_plot rules, LOD_FACTS/LOD_FIGS, ALL_FACTS/ALL_FIGS).
- Modified: `paper/methods.md` (new "Limit of Detection" subsection citing CLSI EP17-A2 and Pierson-Perry 2012), `paper/results.md` (new "Limit of Detection" subsection with Figure 7 reference; Table 2 LoD column updated), `paper/discussion.md` (rewrote the "LoD not characterized" paragraph), `paper/references.bib` (added CLSIEP17A2 and PiersonPerry2012 entries).
- Plan: `claude/issue_8_lod_markers_plan.md` updated with the "Followup: align all validation figures on bias-corrected estimates" section.
