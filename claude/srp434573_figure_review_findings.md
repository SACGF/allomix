# SRP434573 figure-review findings (overnight run, 2026-06-24)

Experiments for Q1 and Q3 from `srp434573_figure_review_plan.md`. Scripts and data in
`output/figure_review/` (gitignored). All numbers are model-free or paper-faithful
(per-patient error tables loaded as in the paper run). Nothing committed to the paper.

## Q1: the 0.5% points read ~0.28% NOT because the mix is below label. Mixing is faithful; a donor-hom-specific contamination background is the driver.

Model-free realised host fraction = pooled host-allele read fraction at donor-homozygous
fully-informative markers (types 0/1), where in a fraction-f mix the host-allele fraction
equals f directly. Script: `output/figure_review/q1_realised_fraction.py`.

Across the full real ladder (10 mixtures, 0.5-10%):
- **raw donor-hom realised = 0.985 x nominal + 0.243%**. Slope ~1: the host DNA tracks
  the label. The mix is NOT below nominal.
- There is a **+0.24% additive background on the donor-hom host allele**, seen directly at
  the pure-donor endpoints (true host = 0%): median raw donor-hom = **0.29%**, per mixture
  0.14-0.48%.
- The **consensus-hom contamination floor reads ~0.00%** at these samples: it does NOT
  capture this background, because consensus-hom sites are chosen where neither contributor
  carries the minor allele, while the donor-hom host-allele markers are exactly where a
  co-pooled genome carrying the host allele inflates the count. This is the **Step 30**
  mechanism (per-marker contamination on the donor-absent/host allele).
- Per-mixture endpoint background correlates with the mixture: F-into-F pairs are highest
  (F3->F2 0.48%, M3->F2 0.45%, F2->F1 0.33%), F-into-M lowest (F2->M2 0.14%, M1->M2 0.17%,
  F2->M1 0.17%), consistent with co-pooled allele-sharing, not a uniform error floor.

The three 0.5% mixtures (F2->M1, F2->M2, M1->M2): raw donor-hom = 0.53-0.55%; their own
0% endpoint background = 0.14-0.17%; so host signal ~= 0.37-0.40%, roughly faithful to
the 0.5% label. The MLE/presence report ~0.28% because they subtract more background than
this; i.e. at 0.5% the estimators slightly OVER-attribute to background and under-report,
they do not reveal a half-strength mix.

Confirming cross-check: subtracting each mixture's OWN 0%-endpoint donor-hom background
from its titration reads gives realised-host-vs-nominal **slope 0.983, intercept -0.04%**
(near-perfect): once the donor-hom background is removed per mixture, true host tracks the
label. At 0.5% the background-subtracted host is ~0.37% (vs MLE/presence ~0.28%), so the
estimators slightly over-subtract at the lowest rung; they do not reveal a half mix.

**Bearing on the paper:** the committed `discussion.md` line ("most plausibly because the
realised mixing fraction departs slightly from the nominal one") is NOT supported by this
check. The realised mixing is ~faithful (slope 0.985). The real driver of the low-end
behaviour and the CI undercoverage is the donor-hom-specific co-pooled contamination
background that the consensus-hom floor misses. Recommend re-emphasising the discussion
toward (a) this contamination background (Step 30) and (b) real het overdispersion not
being an exact beta-binomial, and dropping the realised-vs-nominal framing. (Dave's call;
not edited.)

## Q2: where the real-reads / BAM figures are (no experiment)

- IN SILICO (simulator, binomial reads): Fig 1 (`fig5_lod_curves.png`) and S7
  (`fig_lod_saturation.png`), both from `run_lod_validation.py`.
- REAL READS / BAM-derived: Fig 4 (`fig_srp434573.png`, the dilution series + 3-person mix);
  **Fig 5 (`fig_subsample_lod_grid.png`) is the real-reads LoD** sub-sampled from the
  high-depth SRP434573 BAMs (`run_subsample_lod.py`), the direct counterpart of S7; and the
  semi-synthetic sub-0.5% ladder (`fig_srp434573_synthetic.png`).
- Suggested: a one-line cross-reference in the S7 caption to Fig 5 as the real-reads
  counterpart.

## Q3: the MLE's positive floor at zero host is mostly negligible; where it is significant it is real contamination, and the presence test already gates it correctly.

Profiled the host-fraction likelihood near 0 for each pure-donor endpoint (true host = 0%),
paper-faithful (per-patient error tables). Script: `output/figure_review/q3_mle_near_zero.py`.
"0 in 95% CI" = LL(host=0) within 1.92 of the profile peak.

Endpoints (true host = 0%):
- **6 of 10 are effectively zero** (MLE 0.000-0.015%, 0 inside the 95% CI): F1->F3, F2->M1,
  F2->M2, M3->F1, M3->F3, M3->M4. These are not false positives.
- **4 of 10 have a significant positive floor** (0 outside CI): F2->F1 0.037%, M1->M2 0.027%,
  F3->F2 0.101%, M3->F2 0.165%. The two largest (F3->F2, M3->F2) are exactly the two highest
  donor-hom contamination backgrounds from Q1 (0.48%, 0.45%). The MLE is correctly fitting a
  real excess of host-allele reads; the excess is co-pooled contamination, not host.

Contrast, real 0.5% host (F2->M1/M2, M1->M2): MLE 0.24-0.29%, 0 outside CI by 48-67 LL units,
far more significant than the endpoint floors (2-25 LL units). Real host at 0.5% is a strong
detection.

The host-presence test calls **not-detected (p=1.0) at ALL 10 endpoints**, including the two
contamination-heavy ones, because its null background is contamination-aware. The MLE point,
even with the per-marker error table, is not.

**Remedy recommendation (pick):**
1. **Gate the reported MLE host on the presence test** (cleanest, low-risk): when host-presence
   is not detected, report "host not detected / < LoD" instead of the MLE point. The presence
   call already separates all endpoints (not-detected) from real 0.5% host (detected). Surface
   the presence verdict next to the MLE in the report.
2. A pure CI/LRT gate on the MLE alone is INSUFFICIENT: F3->F2 and M3->F2 cross the 95%
   threshold on contamination, so the MLE would still false-positive.
3. Deeper fix = **Step 30** (apportion the contamination floor to the specific donor-hom
   markers a co-pooled genome inflates): this pulls the MLE point itself down, attacking the
   same background behind Q1. Strongest but more work.
4. If not fixing now, WRITE IT UP: the MLE is "how much given host present", the presence test
   is "is any present"; the raw MLE point must not be read as a detection near zero. Add the
   S12 caption sentence (the floating points are pure-donor endpoints; presence is the call).

## Follow-up (Dave's Q1 "why down not looser CI / is it a /2" and Q3 "was clamping a mistake")

Script: `output/figure_review/q1q3_followup.py`.

Q1 class decomposition at 0.5% nominal (host %):
| mixture | two-rho | shared-rho | donor-hom only | donor-het only |
|---|---|---|---|---|
| F2->M1 | 0.262 | 0.443 | 0.265 | 0.000 |
| F2->M2 | 0.287 | 0.834 | 0.287 | 0.292 |
| M1->M2 | 0.241 | 0.475 | 0.243 | 0.000 |

- **two-rho == donor-hom-only** (0.262 vs 0.265, etc). At low f two-rho down-weights the
  donor-het class to ~zero, so the estimate IS the donor-hom-only estimate.
- shared-rho is HIGHER (0.44-0.83) because it keeps the donor-het markers, whose symmetric
  near-0.5 overdispersion RECTIFIES upward into a positive host signal (the old #33 floor).
- So the "down" is a POINT/BIAS effect, not a variance effect: two-rho removes a DIRECTIONAL
  upward bias from the donor-het class, which moves the point (not just the CI). Contamination
  does not push the point down; it inflates the donor-hom reads (up). The point lands low
  because (a) two-rho drops the donor-het upward rectification and (b) the error model
  subtracts a per-marker background from the donor-hom reads (raw 0.54 -> MLE 0.26; model-free
  background-subtracted true ~0.38, so the error model over-subtracts ~0.1pp at this rung).

NOT a divide-by-2. two-rho estimate / nominal across levels (would be ~0.5 everywhere if /2):
0.5%: 0.48-0.57; 1%: 0.61-0.83; 2.5%: 0.78-0.88; 5%: 0.77-0.88; 10%: 0.85-0.99. The ratio
RISES toward 1 with fraction. The estimator fits ONE shared f; splitting rho does not divide
f. The ~half at 0.5% is a roughly FIXED background subtraction eating a level-dependent
fraction of signal (largest at the lowest rung) -> coincidence, not a /2.

Q3 clamp test (profile extended into negative host fraction):
- clean endpoint (F2->M1, true 0%): unconstrained peak @ **-0.020%** -> the clamp floors a
  ~0 (slightly negative) estimate at exactly 0. Negligible upward nudge.
- contamination endpoint (M3->F2, true 0%): unconstrained peak @ **+0.165%** -> a GENUINE
  positive interior maximum (real excess host-allele reads from contamination). The clamp is
  irrelevant here; unclamping would NOT remove it.
- So clamping was not the main mistake: clean endpoints sit ~0 with or without it; the visible
  S12 floors are contamination signal, not clamp pile-up. Valid kernel: a signed estimate
  would let true-0% samples center on 0 for diagnostics/presentation (removes the tiny
  boundary nudge), but for reporting the clamp is correct (negative % unphysical) and the real
  fix for the contamination cases is the presence gate / Step 30.

## Cross-cutting takeaway
Q1 and Q3 point to the SAME root cause: co-pooled contamination landing specifically on
donor-homozygous host-allele markers, which the consensus-hom contamination floor does not
see. It (a) puts a ~0.24% background under the low-fraction MLE/presence (Q1) and (b) produces
the significant MLE floors at the high-contamination zero-host endpoints (Q3). Step 30 (per-
marker contamination apportionment) is the principled fix for both; gating the MLE on the
presence test is the cheap interim fix for the Q3 cosmetics.
