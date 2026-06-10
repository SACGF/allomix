# allomix: interesting methods and results

A ranked inventory of what allomix actually does, built from a deep read of `src/allomix/`,
`scripts/`, and `paper/scripts/` (not from the paper, which is being rewritten separately).

The **Rank (1-10)** column is a methods-and-clinical-importance judgement: how novel, how
load-bearing for the tool's claims, and how much it would matter to a reviewer or a clinician.
10 = central and novel; 1 = minor implementation detail. Scores are my assessment, not the
authors'.

## Ranked table

| Rank | Feature | What we do | Why it is interesting / unexpected |
|---|---|---|---|
| 10 | **Host-presence detection (relapse early-warning)** | A presence/absence test for tiny residual host, separate from the magnitude MLE. Uses only donor-homozygous markers where the host carries the donor-absent allele, so that allele should appear only at the sequencing-error background in a pure-donor sample. Two statistics: a pooled one-sided Poisson test and a bounded-MLE likelihood ratio test on `q_i(f_h) = e_i + (h_i/2)·f_h`. Reports p-value, host fraction MLE, and a profile-likelihood CI. | This is the clinical headline: detecting host coming back *below* the resolution of the fraction estimate. The boundary `f_h >= 0` is handled correctly with a chi-bar-square (50:50 mixture of point mass at 0 and chi2_1) p-value, not naive Wilks. The CI is profile-likelihood clipped at 0. Getting boundary inference right here is the kind of thing routinely done wrong. |
| 10 | **Beta-binomial overdispersion as the real LoD limiter** | The whole likelihood is beta-binomial, not binomial, with a single concentration `rho` *jointly estimated* (profiled out at every grid point). Variance inflation factor `1 + (n-1)/(rho+1)` flows consistently into the GoF, the Fisher-information SE, and the CIs. | The simulator quantifies this as the single dominant control on achievable sensitivity at >1000x depth: at high depth binomial noise is negligible, so extra-binomial PCR/capture jitter sets the floor. LoD goes from 0.12% (binomial) to 0.41% at rho=100 (3.5x) to 0.75% at rho=30. This is the honest explanation for the gap between idealized in-silico LoD and real performance, and it ties directly to the SRP434573 REVIEW flags (GoF failures from overdispersion). |
| 9 | **SRP434573 real-data validation + forensic ~0.2% contamination floor** | Re-analyze a public titrated-mixture MIP panel dataset (10%→0.5% minor, 2- and 3-person mixes). Map the *minor* contributor to the HOST role so the dilution reads as a declining-chimerism / relapse trajectory. Establish a real co-pooled contamination floor (~0.1-0.2%) using only panel genotypes: at host+donor consensus-homozygous sites the minor allele can only be error or third-party DNA, and a carrier-vs-no-carrier dose-response separates index hopping from sequencing error and from mapping artifacts. | The role-flip reframes a generic dilution series as a clinically meaningful monitoring series. The floor is *established by dose-response*, not assumed: contamination scales with co-pooled allele dose, a site artifact does not, and a flowcell-wide index hop also touches the pure reference runs. Over the reliable range R2=0.95, MAE 0.63%; tracks to 1% host; collapses/scatters at 0.5% exactly where the contamination floor competes with truth. The minority contributor's low-fraction ceiling is bounded by a real property of co-pooled research data, not by the estimator. |
| 9 | **One-sided robust trim (preserves low-fraction host signal)** | The robust refit uses a median/MAD scale and iteratively drops outlier markers, but only trims residuals pointing *away* from host presence. Markers whose residual deviates in the host-present direction (`sign(host_alt/2 - expected_vaf)`) are protected and never trimmed. | A symmetric MAD cut destroys exactly the markers carrying the signal of interest: at low host fraction the host-carrying markers sit off a donor-dominated fit and read as outliers, collapsing the MLE to 0% when truth is ~1%. The asymmetric trim is a domain-aware robust estimator with the explicit trade "at the limit of detection we would rather keep a few artifacts than discard a real low-fraction host signal." Motivated by a real SRP434573 failure, then re-validated in silico against known truth. |
| 8 | **In-data contamination estimator (empirical p10 floor)** | Measure third-party signal at consensus-homozygous sites (neither host nor any donor can supply the minor allele). Headline estimate is the background-subtracted *median* per-site minor fraction; sequencing-error floor is the 10th percentile of per-site minor fractions (the no-carrier/error sites), so contamination is the heterogeneous *excess* over a uniform error floor. Cap sites above 10% as miscalls. | Three deliberate robustness choices: median not pooled mean (a few 40-100% miscall sites would dominate a mean), data-internal error floor (a uniform error elevation lifts the floor too and is correctly not called contamination), and the 10th percentile chosen so the floor still lands on no-carrier sites even on a densely pooled panel. Distinguishes contamination from real low chimerism by marker *geometry*, not magnitude: the three low-fraction signals (host presence, contamination, gross swap) read disjoint marker sets. |
| 8 | **Two-phase calling architecture (GATK GT + bcftools mpileup AD)** | Host/donor genotypes from GATK joint calling; admix allele depths from forced `bcftools mpileup` at panel sites. Documented and empirically verified rationale in `doc/joint_calling.md`. | The low-fraction signal lives in minority ALT reads that `HaplotypeCaller -ERC GVCF` strips at hom-ref blocks. Verified empirically: across ~9M reads at admix hom-ref calls, *zero* ALT reads retained in joint-called AD. No GATK flag recovers them. A somatic caller is also wrong because it is built to *reject* sub-1% events. The v1 paper still describes the old single-phase model; this is flagged as a correctness bug to fix before colleagues read it. |
| 8 | **CNV/LoH handling + robust-refit recovery** | Simulate sub-clonal host copy-number aberrations (CN-LoH, deletion CN1, gain CN3) with a copy-number-weighted VAF model (the local mixing fraction differs from genome-wide because the locus DNA mass changes). Validate donor-detection vs host-relapse-detection under aberration burden, with standard and robust estimators. | The recipient is a haematological-malignancy patient whose relapsing clone routinely carries CN changes, breaking the "divide by 2, everything is diploid" assumption. Result: relapse (host) detection is largely unaffected (aberration rides with the signal), but mixed-chimerism donor detection degrades badly (LoD pushed above the 20% ceiling), and the robust refit recovers most of it (gain-high donor LoD 2.72% → 0.20%). This is the quantitative justification for the robust trim. |
| 8 | **Profile-likelihood confidence intervals** | CIs from likelihood-ratio inversion (`chi2.ppf(0.95, df=1)`), with rho profiled out at each scan point, boundary pinning at 0/100%, and the max-LL reference re-derived from the *same* profiled optimizer to avoid brentq sign errors. | Not Wald/Fisher and not bootstrap. The boundary handling matters because chimerism near 0% or 100% is the norm. The self-consistent reference (profile_ll at f_mle, not the joint optimum) is a real numerical fix so the root-finder brackets a sign change. CI coverage validated at ~98% even at 100x depth; bias correction tightens CI width ~30%. |
| 7 | **Three-state QC with anti-gaming GoF and presence-vs-MLE cross-check** | PASS / REVIEW / FAIL. GoF is a beta-binomial-variance chi-squared computed twice (pre- and post-robust-trim), gated on the worse of the two. Cross-checks the host-presence LRT against the magnitude MLE and warns when low-level host is detected below the MLE's resolution. Marker-loss diagnosis names the dominant bottleneck. | REVIEW (computed but check a reliability flag) vs FAIL (unusable) is a deliberate clinical distinction. The pre/post-trim GoF prevents the robust trimmer from masking a genuinely bad fit by discarding its own outliers. The error-adjusted expected VAF fixes a real variance-floor blow-up at saturated homozygous markers. Honest about validation maturity: the presence-vs-MLE disagreement stays a soft warning because operating characteristics on real samples are still being mapped. |
| 7 | **Per-marker direction-specific empirical error null** | Estimate `e_refalt` and `e_altref` separately per site from a panel of normals (pooled reads, not averaged samples), with VAF guards against contamination and a 1e-5 floor. Feeds the host-presence H0 background per marker by direction. | REF→ALT and ALT→REF rates genuinely differ (oxidation/8-oxoG, strand bias, flanking context), and the dominant error *direction* sets the background for low-fraction detection. The planning doc calls the per-site empirical error null "the single highest-value algorithm change" remaining, expected to clear the overdispersion-driven REVIEW samples. Pooling reads is the correct sufficient-statistic MLE for a shared rate. |
| 7 | **Logit-space bias correction + both-het bias estimation** | Per-marker amplification bias measured at het sites as `median(VAF) - 0.5`, corrected multiplicatively in logit space `expit(logit(w) - logit(0.5 + bias))`, not additively. A second estimator reads bias from admix samples at host-and-all-donor-het markers (true VAF 0.5 regardless of mixing fraction). | Additive correction overcorrects at informative markers (expected VAF near 0/1, the norm at low chimerism); logit-space stays valid everywhere and reduces to `0.5 - bias` at a het. The both-het estimator calibrates *the same caller* being analyzed, sidestepping the caller-mismatch footgun of the two-phase pipeline. Honest finding from the plan: bias correction barely moves the point estimate but tightens CIs ~30% (it sharpens precision, not accuracy). |
| 6 | **EP17-A2 / Currie per-sample analytical LoB and LoD** | Fisher-information SE of the fraction (each informative marker contributes `(dp/df)^2 / var`, with the overdispersion inflation factor), then `LoB = z·SE(0)`, `LoD = LoB + z·SE(LoB)`, one-sided z=1.645. Reported per sample and floored up to the contamination estimate. | Blends a regulatory metrology framework (CLSI EP17-A2) with the model's own information matrix, substituting per-sample analytical SE for EP17's repeated blank measurements. Achieves 0.158% (unrelated) / 0.229% (sibling) LoD at 1000x/100 markers, 0.10% best corner. The code is explicit this is an analytical ceiling, not a validated assay LoD from replicates. |
| 6 | **Multi-donor estimation (2 donors + host)** | Simplex model with host as implicit complement, triangular grid over `(f1, f2)` with `f1+f2<=1`, renormalization on overshoot, per-donor informative-marker tracking, and per-donor profile-likelihood CIs. | Correctly ranks 16/16 asymmetric sibling mixes and resolves a 3-person mix cleanly. Honest weak spot recorded: per-donor CI coverage drops to ~50% because partitioning a mix between two *related* donors leaves few markers that separate each one individually (46/41 of 61 informative). Total and ranking stay accurate. |
| 6 | **Relatedness QC (allele-frequency-free kinship)** | Somalier-style robust coefficient `(shared_hets - 2·ibs0) / min(het_a, het_b)` on autosomal clean biallelic GTs. Asymmetric verdict logic: losing close relatedness → FAIL (swap signature); cousin→unrelated crossing only REVIEW; detected identical always FAIL (syngeneic donor is unmeasurable). | AF-free is the right choice for a panel-agnostic tool with unknown panel AFs. The asymmetry mirrors real failure modes rather than a symmetric "match declared" test. Relatedness matters because related donors share genotypes → fewer informative markers (sibling 36.8 vs unrelated 54.3) → worse LoD, while point-estimate accuracy stays flat (remaining markers are still unbiased). |
| 5 | **Read-level artifact filter (effect-size, panel-aware)** | Screen donor-homozygous markers for soft-clip-length bias (`|SCBZ|>3`), read-position bias (`|RPBZ|>6`), and strand bias (minor-strand fraction <0.10), using bcftools mpileup annotations. Auto-disable the strand test on single-strand amplicon panels when ≥90% of markers are one-strand. | Strand bias is judged by *effect size, not significance*: at >1000x a real allele's mild 55:45 skew is highly significant yet harmless, an artifact is extreme (~95:5) regardless of depth, so a significance-based filter would reject real alleles. The single-strand auto-detect (issue #18) stops the filter from silently destroying the test on the lab's actual amplicon panel. Caught a real TP53 intron-3 artifact. |
| 5 | **Self-consistent simulator (generative model matched to estimator)** | The synthetic-VCF generator uses the same 4-state error normalization (`p_alt/(1-2e/3)`), the same beta-binomial rho, the same logit-space bias injection, and the same Vynck typing as the likelihood. Layered noise: log-normal depth CV (0.43), heavy-tailed two-component bias mixture, overdispersion, allele/locus dropout, CNV/LoH, sibling-trio Mendelian segregation. | The credibility of the in-silico validation rests on this matching: the estimator is tested against data drawn from a deliberate *superset* of its own assumptions plus controlled mismatch terms (CNV/LoH, heavy tails, dropout, contamination) that the robust refit and flooring are designed to absorb. The sibling-trio path preserves true 3-way IBD correlation, which independent pairwise draws cannot. |
| 4 | **Index-hopping metadata flag, separate from in-data contamination** | An optional `##allomixRunUnit` VCF header records flowcell:lane. A pure-metadata flag warns when the admix sample shares a run unit with the host (index-hopping risk), kept separate from the in-data contamination estimate. | Clean separation of mechanism (metadata: sharing a run is a risk, not a defect) from measurement (in-data: whether it actually bit). On SRP434573 the metadata flag degrades to "cannot determine" because SRA stripped the read-group PU tag, while the in-data dose-response probe still recovers the floor: two independent routes to the same risk. |
| 4 | **Admix-consistency swap / third-genome test** | At consensus-homozygous sites the admixture must show that homozygote up to error regardless of mixing fraction. Per-site binomial tail at `error_rate` (conservative, not error/3), then a binomial-of-binomials `swap_pval` over discordant sites. | Catches a wrong-patient VCF or sample swap that the MLE goodness-of-fit cannot see, because the MLE only ever looks at *informative* markers and never at these consensus sites. Complements relatedness QC and contamination from a third marker geometry. |
| 3 | **Vynck marker-type informativeness taxonomy** | Classify each host/donor pair into 6 types by alt-dose (0,1 fully informative; 10,11,20,21 partial; equal = non-informative), tracked per donor in multi-donor sets, with a full marker-loss funnel (`MarkerCounts`) recording every drop reason. | Fully informative markers (host 0/100% of an allele) give a clean shift from 0 or 1; partial types start at 0.5 and give a half-size shift, which matters for per-marker information content and LoD. The drop-reason funnel exists purely for QC transparency; the estimator never reads it. |
| 3 | **Reproducible nested validation design** | LoD uses a two-level nested design (K donor/host pairs × M sequencing replicates, only the blend seed varies), strictly nested marker panels (50-marker = bit-identical prefix of 400), and SHA-256-derived process-stable seeds. N>=5 replicates everywhere. | The previous pooled design conflated IBD-sharing variation across sibling pairs with sequencing noise, leaking into the LoB and making LoD-vs-panel-size non-monotone. Separating the two sources fixes it. SHA-256 seeds fix a prior non-reproducibility bug (Python hash randomization). |
| 2 | **GT/AD consistency guard (panel side only)** | Drop reference-sample markers whose called GT contradicts the AD-derived VAF (het outside [0.35,0.65], etc.), requiring >=20 reads. Never applied to admix (a mixture's VAF is not expected at 0/0.5/1). | Stops a GATK het "rescued" from marginal evidence in a small 2-sample joint call from feeding systematic bias into the estimator. Bounds deliberately loose to tolerate genuine capture bias. |

## Top-level themes a reviewer would notice

A few cross-cutting ideas show up repeatedly and are worth stating once:

1. **Three low-fraction signals are separated by marker geometry, not by re-thresholding one
   statistic.** Host presence (donor-hom sites, host carries minor), contamination
   (consensus-hom sites, neither contributor carries minor), and gross swap (consensus-hom,
   minor individually significant) read disjoint marker sets, so they are genuinely orthogonal.

2. **Effect size over significance at high depth.** Both the strand-bias artifact filter and the
   contamination magnitude reporting deliberately prefer effect size to p-values, because at
   >1000x a significance test over-rejects harmless real signal.

3. **Boundary-aware inference throughout.** The `f >= 0` constraint is honored with chi-bar-square
   p-values, profile CIs clipped at 0, and boundary pinning, rather than Wald approximations that
   are wrong exactly where chimerism matters (near 0%).

4. **Real-data findings fed back into in-silico validation.** The two things SRP434573 revealed
   (a ~0.2% index-hopping contamination floor, and the symmetric-trim collapse at low host
   fraction) each became a dedicated validation script with known ground truth
   (`validate_contamination_lod_floor.py`, `validate_onesided_trim.py`).

5. **Honest about the gap between idealized and real performance.** Overdispersion is presented as
   the dominant real-data LoD limiter; the analytical LoB/LoD are labelled as ceilings, not
   validated assay limits; and several QC promotions are held as soft warnings pending real-sample
   operating characteristics.

## Caveats noted in the code worth not overstating in the paper

- The analytical LoB/LoD are Fisher-information idealizations on one sample, not empirical assay
  limits from replicate dilutions.
- Per-donor CIs under-cover (~50%) when splitting a mix between two related donors.
- The contamination floor is applied as a flat per-sample scalar (per-marker apportionment is a
  noted future refinement, Step 30, with a known risk of inverting the presence-vs-MLE QC gate).
- The LoD comparison against commercial assays is not head-to-head; allomix's number is an
  EP17-A2 analytical ceiling on near-binomial simulated data.

---

## Paper coverage audit (added 2026-06-10)

How many words the current paper (`paper/*.md`: abstract, introduction, methods, results,
discussion, supplementary) spends on each ranked feature above. Counts are *attributed
estimates*: I tagged the sentences/passages that describe each topic and summed their words,
rounding to the nearest ~10. Template placeholders (`{{ ... }}`) count as one word. Where a
passage serves two topics (e.g. the SRP contamination floor is both "real-data validation" and
"in-data contamination estimator") the words are split by emphasis, so column totals are
approximate, not exact. Section totals for reference: abstract 259, introduction 598, methods
3328, results 3908, discussion 1716, supplementary 1336.

### Per-topic word count

| Rank | Feature | Abs | Intro | Meth | Res | Disc | Supp | Total | Coverage |
|---:|---|--:|--:|--:|--:|--:|--:|--:|---|
| 10 | Host-presence detection (relapse early-warning) | 0 | 0 | ~15 | ~20 | 0 | 0 | **~35** | **Near-absent.** Only the Fig 9 "open squares" mention and the "host-presence strand-bias filter" aside. The Poisson test, bounded-MLE LRT, chi-bar-square boundary p-value, profile CI: **not described anywhere.** |
| 10 | Beta-binomial overdispersion as the LoD limiter | 0 | 0 | ~200 | ~90 | ~280 | ~290 | **~860** | Well covered. Dedicated Discussion section + Supp S7/S8 + Methods likelihood subsection. The best-served topic. |
| 9 | SRP434573 real-data validation + contamination floor | ~30 | 0 | ~290 | ~430 | ~110 | 0 | **~860** | Well covered (Methods dataset para + Results section + Discussion limitation). |
| 9 | One-sided robust trim (preserves low-fraction host) | 0 | 0 | ~50 | ~120 | 0 | 0 | **~170** | **Mis-described.** Methods calls it a generic "median/MAD outlier-resistant refit." The *one-sided / asymmetric / host-direction-protected* property (the actual novelty) is never stated. |
| 8 | In-data contamination estimator (p10 floor) | ~10 | 0 | ~10 | ~150 | ~20 | 0 | **~190** | Floor is *described as a dataset property* via dose-response in Results, but the estimator method (median per-site minor, p10 error floor, 10% miscall cap) is **not in Methods.** |
| 8 | Two-phase calling (GATK GT + bcftools AD) | ~5 | 0 | ~190 | 0 | ~120 | 0 | **~315** | Well covered. (Inventory flagged "v1 paper still single-phase" as a bug; **methods.md now describes two-phase correctly, so that bug is resolved.**) |
| 8 | CNV/LoH handling + robust-refit recovery | 0 | 0 | ~280 | ~430 | 0 | 0 | **~710** | Well covered (Methods + Results section + Fig 8). |
| 8 | Profile-likelihood confidence intervals | ~3 | 0 | ~180 | ~150 | ~120 | ~60 | **~510** | Well covered. |
| 7 | Three-state QC + anti-gaming GoF + presence-vs-MLE cross-check | ~6 | 0 | ~160 | ~60 | 0 | 0 | **~225** | **Partial.** Generic 5-point QC list is there, but PASS/**REVIEW**/FAIL three-state, pre/post-trim GoF, and the presence-vs-MLE cross-check are not described as such. |
| 7 | Per-marker direction-specific empirical error null | 0 | 0 | 0 | 0 | 0 | 0 | **0** | **Absent.** Paper uses only the flat ε=0.01 4-state model. The empirical per-site, per-direction null is not mentioned. |
| 7 | Logit-space bias correction + both-het estimation | ~5 | 0 | ~290 | ~180 | ~80 | ~40 | **~595** | Over-covered relative to impact (see commentary). |
| 6 | EP17-A2 analytical LoB/LoD | 0 | 0 | ~330 | ~450 | ~250 | ~270 | **~1300** | The single most-covered topic across the paper. |
| 6 | Multi-donor estimation (2 donors + host) | ~12 | 0 | ~180 | ~330 | ~120 | 0 | **~640** | Well covered. |
| 6 | Relatedness QC (AF-free kinship coefficient) | 0 | 0 | 0 | 0 | 0 | 0 | **0** | **Absent.** (NB: the "Effect of Donor-Host Relatedness" Results section, ~560 words, is a *different thing*: relatedness' effect on informative-marker count and accuracy, not the somalier-style kinship QC check.) |
| 5 | Read-level artifact filter (effect-size, panel-aware) | 0 | 0 | ~15 | 0 | 0 | 0 | **~15** | **Near-absent.** Only the "auto-disables on single-strand panel" aside in the SRP methods. SCBZ/RPBZ/strand-bias effect-size logic not described. |
| 5 | Self-consistent simulator | 0 | 0 | ~370 | ~30 | 0 | ~200 | **~600** | Well covered (Methods sim framework + Supp S1/S3/S4/S9). |
| 4 | Index-hopping metadata flag (separate from in-data) | 0 | 0 | 0 | ~20 | 0 | 0 | **~20** | **Near-absent.** "Index hopping" named once in Results as the likely floor mechanism; the `##allomixRunUnit` metadata flag and its separation from in-data contamination are not described. |
| 4 | Admix-consistency swap / third-genome test | 0 | 0 | 0 | 0 | 0 | 0 | **0** | **Absent.** |
| 3 | Vynck marker-type informativeness taxonomy | 0 | 0 | ~230 | ~30 | 0 | 0 | **~260** | Well covered (Methods marker classification). |
| 3 | Reproducible nested validation design | 0 | 0 | ~120 | ~40 | 0 | 0 | **~160** | Adequately covered inside the LoD methods. |
| 2 | GT/AD consistency guard (panel side) | 0 | 0 | 0 | 0 | 0 | 0 | **0** | **Absent.** |

### Cross-cutting themes (inventory's "themes a reviewer would notice")

| Theme | Words in paper | Coverage |
|---|--:|---|
| 1. Three low-fraction signals separated by marker *geometry* | ~0 | **Absent.** The unifying framing (host-presence vs contamination vs swap read disjoint marker sets) appears nowhere. This is the inventory's headline conceptual contribution. |
| 2. Effect size over significance at high depth | ~0 explicit | Only implicit in the bias-correction logit discussion; never stated as a principle. |
| 3. Boundary-aware inference (chi-bar-square, clipped CIs) | ~30 | Profile CI is described; the boundary handling (chi-bar-square p-value, pinning at 0) is barely surfaced. |
| 4. Real-data findings fed back into in-silico validation | ~60 | Present for the contamination floor; the symmetric-trim-collapse origin story is not told. |
| 5. Honest about idealized-vs-real gap | ~400 | Well covered (overdispersion framing + LoD-as-ceiling caveats in Discussion). |

### Large allocations that are NOT ranked features

- **Accuracy across sequencing depths** (Results "In Silico Validation Across Depths" + Table 1 +
  Figs 1-3): ~600 words. This is the paper's opening Results block and its single largest, yet it
  is not in the ranked inventory at all: it is generic "MAE is low" validation, the least novel
  thing the tool does.
- **Effect of donor-host relatedness on accuracy** (Results section + Table 3 + Fig 4): ~560
  words, also not a ranked feature (and easily confused with the absent relatedness *QC*).

## Commentary: what to fix and re-order

**1. Coverage is roughly *inverted* against the inventory ranking.** The two rank-10 clinical/
methods headlines split opposite ways: overdispersion gets ~860 words, but **host-presence
relapse detection (rank 10) gets ~35 and has no Methods description at all.** Meanwhile EP17 LoD
(rank 6) is the single most-covered topic (~1300 words) and bias correction (rank 7, which the
inventory and the paper itself both say "barely moves the point estimate") gets ~595 words across
four sections. The paper spends its words on the least differentiating, most idealized analyses
(depth-MAE, LoD ceiling, bias correction) and is nearly silent on the features that actually
distinguish allomix from a commercial kit.

**2. An entire class of features is missing from Methods.** Six ranked features are at or near
zero words: host-presence test (10), empirical error null (7), relatedness QC (6), artifact
filter (5), index-hopping flag (4), swap test (4), GT/AD guard (2). Five of these are the
"marker-geometry safety suite" that the inventory's theme #1 says is the conceptual core. Concrete
fix: add Methods subsections for **(a) host-presence detection**, **(b) contamination estimation**,
**(c) sample-swap / third-genome test**, **(d) relatedness QC**, and fold (e) the artifact filter,
(f) the GT/AD guard, and (g) the index-hopping metadata flag into a "QC and safety checks"
block. Today the Methods "Quality Control" section is a generic 5-point list that does not mention
any of them.

**3. Correct the robust-trim description (likely-wrong-as-written).** Methods describes the refit
as a symmetric "median/MAD outlier-resistant refit." The inventory's rank-9 point is that the trim
is deliberately **one-sided**: it protects markers whose residual points toward host presence. As
written, the paper describes the opposite of the actual (and more defensible) algorithm. This
should be fixed for correctness, not just emphasis.

**4. Surface the unifying frame.** Theme #1 ("three low-fraction signals separated by marker
geometry, not by re-thresholding one statistic") is the strongest single idea in the inventory and
appears nowhere. One paragraph stating that host-presence, contamination, and gross-swap read
disjoint marker sets would tie the three missing Methods subsections together and give the QC story
a spine.

**5. Re-order the Results to lead with the clinical and the real.** Current order: depth-MAE → bias
→ CI calibration → tool comparison → LoD → CNV → relatedness → multi-donor → timeline → real data.
This front-loads the least novel material and buries the two most compelling items (the real-read
SRP validation and the longitudinal/relapse use case) at the very end. Suggested order:
   1. **Longitudinal monitoring + host-presence relapse early-warning** (the actual clinical use
      case; promote the timeline figure and add the presence test).
   2. **Real-data validation on SRP434573** (the only real reads in the paper; currently last).
   3. **Accuracy and limit of detection** (depth, LoD, overdispersion as the honest limiter) as the
      analytical characterization that supports 1-2.
   4. **Robustness** (CNV/LoH, relatedness, multi-donor) as stress tests.
   5. **QC / safety suite** (swap, contamination, relatedness QC, artifact filter) as the
      deployment-readiness section.

**6. Trim what is over-weighted.** The depth-sweep accuracy block (~600 words) and the bias-
correction material (~595 words across four sections) can each be compressed once the QC and
host-presence material is added, keeping total length roughly flat. Bias correction in particular
is repeated in Methods, Results, Discussion, and Supp for an effect the text itself calls "modest
by design."

---

## Paper rewrite plan (decided 2026-06-10)

Decisions from a planning discussion, to be executed as a full in-place rewrite of `paper/*.md`
(the paper is too large for line-by-line edits; rewrite each section file in place on a branch and
let git history serve as the diff).

### Audience and voice

- Target: *Journal of Molecular Diagnostics* (per the repo CSL). Real readers are molecular
  pathologists and clinical lab directors, plus a few bioinformaticians.
- This audience is fluent in assay mechanics (capture bias, strand artifacts, sample swaps,
  contamination, LoD) but **formula-averse**: a page of log-likelihoods makes them stop reading.
- Rule: the **body carries no derivations.** Every formula currently in `methods.md` and
  `results.md` moves to **Supplementary Methods**. Keep exactly one anchoring *concept* in the body,
  stated as prose plus at most one simple line: the mixture model, "expected VAF = a weighted blend
  of the two known genotypes, ((1-f)*host_dose + f*donor_dose)/2." Everything heavier (beta-binomial
  Gamma-function likelihood, 4-state error algebra, profile-likelihood chi-square, logit bias
  formula, simplex grid search) goes to supp.
- Replace each removed formula with a one-sentence plain-language *mechanism* description (e.g.
  beta-binomial -> "allele counts are modeled so marker-to-marker scatter can exceed simple
  sampling noise; one overdispersion parameter is fit from the data and widens the confidence
  intervals on noisier panels").

### Worked-example figure (replaces the math for clinicians)

A single boxed worked example / figure walks through real A/C/G/T counts at two markers and teaches
the whole method by counting, no algebra:

- **Marker A, fully informative (host A/A, donor G/G).** Host never makes G, donor is all G, so the
  fraction of G reads *is* the donor fraction. `A=970,G=30` -> 3% donor. Flip to near-full-donor:
  `A=12,G=988` -> the 12 stray A reads can only be host -> 1.2% residual host still present. Same
  marker, two questions: "how much donor?" (majority shift) vs "is any host left?" (minority allele
  only one person could produce). This is where the **two complementary tests** are introduced.
- **Marker B, partially informative (host A/G, donor G/G).** Host already contributes 50% G alone,
  so the donor only nudges it: `A=485,G=515` -> 51.5% G -> 3% donor, but half the signal size.
  Motivates the informativeness idea without a Type-0/11 table in the body.
- Three one-sentence follow-ons: stray reads at 0% donor -> why there is a measured error floor;
  one marker is noisy -> why dozens are pooled into a confidence interval; markers disagree more
  than pure counting noise -> why the interval is widened (overdispersion).

### Headline framing and honesty guardrails

- One-sentence message: **sensitive chimerism monitoring from panels a lab already runs, no
  dedicated assay.** Supporting differentiators: sub-1% sensitivity, multi-donor, built-in QC/safety
  checks, and a dedicated residual-host presence test.
- **Two complementary tests** are a deliberate theme: a magnitude MLE ("how much donor?") and a
  residual-host presence/absence test ("is any host left, below the quantification limit?"). State
  they read different markers and answer different questions.
- Do **not** make "relapse early-warning" a demonstrated *result*. It is not clinically validated
  (presence-vs-MLE cross-check is still a soft warning; validation is in silico + one titrated
  research panel). Present the presence test as a *capability*, validate it on SRP down to 1% host,
  and keep the clinical relapse payoff in the Discussion as motivation only.
- Be openly bounded: "we started in silico because that is what we had," then bring the public
  real-data validation (SRP434573) up in prominence. Frame the EP17 LoD as an analytical ceiling,
  not a head-to-head assay LoD.

### Section spine (write the abstract first, then expand in this order)

1. **Tool overview** - VCFs in (already available), % chimerism out, single + multi-donor. Workflow
   schematic, no math.
2. **Accuracy + longitudinal monitoring** - lead with the timeline figure (clinical money shot);
   fold the depth-accuracy numbers in, compressed.
3. **Sensitivity vs commercial kits** - comparison table + LoD, framed as an analytical ceiling.
4. **Real titrated-mixture validation (SRP434573)** - promoted from last; the only real reads.
5. **Residual-host presence test** - the second of the two complementary tests; needs writing from
   scratch (currently ~0 words in the paper).
6. **Robustness** - relatedness, multi-donor, recipient CNV.
7. **Built-in QC and safety** - sample-swap / third-genome test, contamination estimate, relatedness
   QC. Sell this to wet-lab readers (mix-up paranoia is a feature, not a footnote).

### Must-ADD to the body (currently absent or near-absent, see audit above)

- Host-presence detection (the second complementary test).
- Contamination estimation (in-data, the dose-response idea in plain words).
- Sample-swap / third-genome consistency test.
- Relatedness QC (AF-free kinship check).
- The "three low-fraction signals separated by marker geometry" unifying frame (host-presence,
  contamination, swap read disjoint marker sets) - one short paragraph, the conceptual spine of the
  QC section.
- **Correctness fix:** describe the robust trim as **one-sided** (protects host-direction markers),
  not as a symmetric median/MAD refit, which is what `methods.md` currently (wrongly) says.

### Brutal cuts / compress (GitHub + supp absorb the detail)

- Depth-sweep accuracy block (~600 words, the current Results opener): compress hard; it is generic
  "MAE is low" and the least novel thing the tool does.
- Bias correction (~595 words across Methods/Results/Discussion/Supp): collapse to a few sentences;
  the text itself calls the effect "modest by design." Method detail -> supp.
- All equations -> Supplementary Methods.
- EP17 LoD (~1,300 words, the single most-covered topic) is over-weighted; keep the headline number
  and the analytical-ceiling caveat, push the sweep detail toward supp/figure.

---

## Paper coverage audit v2 (re-counted 2026-06-10, post-rewrite)

Re-count of the same per-topic attribution after the full in-place rewrite of `paper/*.md` was
executed. Same method as the v1 audit above (attributed estimates: tag the sentences/passages
describing each topic, sum, round to ~10; `{{ ... }}` placeholders count as one word; shared
passages split by emphasis, so column totals are approximate). The **v1 Total** column is copied
from the audit above so the two are directly comparable.

Current section totals (raw `wc -w`, including tables, captions, and figure legends; not directly
comparable to v1's prose-only section totals): abstract 330, introduction 727, methods 3518,
results 3388, discussion 1807, supplementary 2918.

### Per-topic word count (v1 -> v2)

| Rank | Feature | v1 Total | v2 Total | Change |
|---:|---|--:|--:|---|
| 10 | Host-presence detection (relapse early-warning) | ~35 | **~880** | **+845.** Now has its own Methods subsection (~180), Results section (~230), Discussion section (~180), Supp S7 (~180), plus abstract/intro framing. Was the single worst-covered headline; now one of the best. |
| 10 | Beta-binomial overdispersion as the LoD limiter | ~860 | **~1010** | +150. Still the best-served topic (dedicated Discussion section + Supp S3/S7/S8 + ablation S4). |
| 9 | SRP434573 real-data validation + contamination floor | ~860 | **~900** | +40. Methods dataset section + Results section (incl. dose-response floor) + Discussion limitation. |
| 9 | One-sided robust trim (preserves low-fraction host) | ~170 | **~250** | **Corrected.** Methods now describes it as explicitly *one-sided / host-direction-protected* (lines 68-70), fixing the v1 mis-description as a symmetric median/MAD refit. |
| 8 | In-data contamination estimator (p10 floor) | ~190 | **~350** | +160. Estimator method now in Methods QC bullet + Supp S8 (median per-site minor, p10 error floor, 10% miscall cap), not just as a Results dataset property. |
| 8 | Two-phase calling (GATK GT + bcftools AD) | ~315 | **~310** | ~flat. Methods section + Discussion workflow para. |
| 8 | CNV/LoH handling + robust-refit recovery | ~710 | **~490** | -220. Compressed; Methods + Results + Supp S9 retained, Fig 7. |
| 8 | Profile-likelihood confidence intervals | ~510 | **~310** | -200. Formula moved to Supp S4; body keeps the concept only. |
| 7 | Three-state QC + anti-gaming GoF + presence-vs-MLE cross-check | ~225 | **~330** | +105. PASS/REVIEW/FAIL, pre/post-trim GoF, and the presence-vs-MLE cross-check are now described as such (Methods QC intro + Results QC section). |
| 7 | Per-marker direction-specific empirical error null | 0 | **~37** | Still near-absent, but now *acknowledged*: Supp S7 states the symmetric global rate is a placeholder "pending a per-site, per-direction empirical error table," and Discussion lists it as a future priority. |
| 7 | Logit-space bias correction + both-het estimation | ~595 | **~400** | -195. Compressed per plan; method detail -> Supp S6, body calls the effect "modest by design." |
| 6 | EP17-A2 analytical LoB/LoD | ~1300 | **~640** | **-660.** No longer the most-covered topic; headline number + analytical-ceiling caveat kept in body, sweep detail pushed to Supp/figure. |
| 6 | Multi-donor estimation (2 donors + host) | ~640 | **~430** | -210. Methods + Results + Discussion + Supp S5, formula to supp. |
| 6 | Relatedness QC (AF-free kinship coefficient) | 0 | **~95** | **Now present.** Methods QC bullet describes the somalier-style AF-free coefficient and the asymmetric verdict logic; Results QC section names it. (Still distinct from the ~560-word relatedness-*effect* Results section.) |
| 5 | Read-level artifact filter (effect-size, panel-aware) | ~15 | **~105** | +90. SCBZ/RPBZ/strand-bias effect-size logic and single-strand auto-disable now in Methods QC + Results QC. |
| 5 | Self-consistent simulator | ~600 | **~450** | -150. Methods sim framework + Supp S9/S1-S6, trimmed. |
| 4 | Index-hopping metadata flag (separate from in-data) | ~20 | **~75** | +55. `##allomixRunUnit` header + separation from in-data contamination now in Methods QC bullet and Supp S8. |
| 4 | Admix-consistency swap / third-genome test | 0 | **~110** | **Now present.** Methods QC bullet + Supp S8 describe the consensus-site consistency test that catches a wrong-patient VCF the MLE cannot see. |
| 3 | Vynck marker-type informativeness taxonomy | ~260 | **~350** | +90. Methods "Which markers are informative" + Box 1 + Supp S1. |
| 3 | Reproducible nested validation design | ~160 | **~160** | ~flat. Inside the LoD methods + Supp S9. |
| 2 | GT/AD consistency guard (panel side) | 0 | **~60** | **Now present.** Methods QC closing sentence describes dropping reference-sample markers whose GT contradicts their AD, applied to references only. |

### New body content not in v1 (worked example)

| Item | v2 words | Note |
|---|--:|---|
| **Box 1: reading the donor fraction off the counts, no algebra** | ~350 | The plan's clinician-facing replacement for the math. Two worked markers (fully vs partially informative) teach the two complementary tests, the error floor, marker pooling, and overdispersion by counting. Serves rank-3 (informativeness), rank-10 (two tests), and the error-floor/overdispersion concepts at once. |

### What changed overall

The rewrite did what the plan at the top of this file set out to do, and the v1 "coverage is
roughly inverted against the ranking" finding is now largely resolved:

1. **The six absent/near-absent features all gained body coverage.** Host-presence test
   (0->880), relatedness QC (0->95), swap test (0->110), GT/AD guard (0->60), index-hopping
   flag (20->75), artifact filter (15->105). Five of these are the "marker-geometry safety
   suite" that theme #1 calls the conceptual core, and they now sit together in a QC and
   sample-integrity Methods block plus a Results QC section.
2. **The over-weighted, least-differentiating topics were compressed.** EP17 LoD (1300->640),
   CNV/LoH (710->490), profile CIs (510->310), bias correction (595->400), multi-donor
   (640->430), simulator (600->450). All equations moved to Supplementary Methods.
3. **The robust-trim mis-description is corrected** to one-sided/host-direction-protected.
4. **The unifying frame is now stated.** "Three low-fraction signals separated by marker
   geometry" appears explicitly in both Methods (QC intro) and Discussion (safety-suite
   section), where v1 had it at ~0 words.
5. **Honesty guardrails held.** Host-presence is presented as a validated *capability* (down to
   1% on SRP), not a demonstrated relapse result; EP17 LoD is repeatedly framed as an analytical
   ceiling, not a head-to-head assay limit; the presence-vs-MLE cross-check stays a soft warning.

Residual gap: the **per-site empirical error null** (rank 7) is still only future-work text
(~37 words), not an implemented-and-described method in the body. That is consistent with the
code (the presence test still uses a symmetric global rate), so this is an accurate omission
rather than a coverage failure.
