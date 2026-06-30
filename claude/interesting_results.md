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
| 8 | **Per-marker-type overdispersion (two-rho, het/hom split)** *(new since 2026-06-10, #33)* | Fit `f` jointly with two independent beta-binomial concentrations, one for donor-hom markers and one for donor-het markers, each profiled out at every grid `f`. Now the default; falls back to shared rho when either class has <30 markers (per-class rho not identifiable, recorded as a diagnostic field, not a QC warning). | The donor-het class (background VAF ~0.5) is more overdispersed, and at low fraction its symmetric amplification scatter rectifies into a small positive host fraction: a sub-0.5% floor that fakes residual host. Giving that class its own rho down-weights it (~0.22 pp offset traced directly to the het class) without touching the donor-hom markers that carry the real low-fraction signal. A targeted refinement of the rank-10 overdispersion model aimed exactly at the clinical operating point. |
| 7 | **Three-state QC with anti-gaming GoF and presence-vs-MLE cross-check** | PASS / REVIEW / FAIL. GoF is a beta-binomial-variance chi-squared computed twice (pre- and post-robust-trim), gated on the worse of the two. Cross-checks the host-presence LRT against the magnitude MLE and warns when low-level host is detected below the MLE's resolution. Marker-loss diagnosis names the dominant bottleneck. | REVIEW (computed but check a reliability flag) vs FAIL (unusable) is a deliberate clinical distinction. The pre/post-trim GoF prevents the robust trimmer from masking a genuinely bad fit by discarding its own outliers. The error-adjusted expected VAF fixes a real variance-floor blow-up at saturated homozygous markers. Honest about validation maturity: the presence-vs-MLE disagreement stays a soft warning because operating characteristics on real samples are still being mapped. |
| 7 | **Per-marker direction-specific empirical error null** | Estimate `e_refalt` and `e_altref` separately per site from a panel of normals (pooled reads, not averaged samples), with VAF guards against contamination and a 1e-5 floor. Feeds the host-presence H0 background per marker by direction. | REF→ALT and ALT→REF rates genuinely differ (oxidation/8-oxoG, strand bias, flanking context), and the dominant error *direction* sets the background for low-fraction detection. The planning doc calls the per-site empirical error null "the single highest-value algorithm change" remaining, expected to clear the overdispersion-driven REVIEW samples. Pooling reads is the correct sufficient-statistic MLE for a shared rate. |
| 7 | **Logit-space bias correction + both-het bias estimation** | Per-marker amplification bias measured at het sites as `median(VAF) - 0.5`, corrected multiplicatively in logit space `expit(logit(w) - logit(0.5 + bias))`, not additively. A second estimator reads bias from admix samples at host-and-all-donor-het markers (true VAF 0.5 regardless of mixing fraction). | Additive correction overcorrects at informative markers (expected VAF near 0/1, the norm at low chimerism); logit-space stays valid everywhere and reduces to `0.5 - bias` at a het. The both-het estimator calibrates *the same caller* being analyzed, sidestepping the caller-mismatch footgun of the two-phase pipeline. Honest finding from the plan: bias correction barely moves the point estimate but tightens CIs ~30% (it sharpens precision, not accuracy). |
| 7 | **Per-marker dose-response contamination correction (Step 30)** *(new since 2026-06-10, #30)* | Predict per-marker co-pooled contamination on donor-hom markers from the co-pooled carrier dose and subtract `slope·n_carriers·depth` from the host-allele count before the MLE. Per-flowcell gate on a significant positive consensus-hom dose-response slope; per-patient magnitude calibrated on the informative markers themselves, pooled across serial timepoints. Off by default. | Localizes contamination instead of taxing every donor-absent marker by the average: the host signal is identical at every donor-hom marker (independent of carrier count) while only the contamination scales with carrier dose, so the two separate cleanly. A clean flowcell has a flat slope and self-selects out (the correction becomes a byte-identical no-op). This is the implemented form of the floor the rank-9 SRP work characterized, turning a measured data property into an optional correction (resolves the old "flat per-sample scalar" caveat below). |
| 7 | **Real-data subsample LoD on SRP434573** *(new since 2026-06-10, #24)* | Rerun the panel-size / depth LoD sweep on real reads: throw away markers and reads from the high-depth titrated mixtures until LoD rises into the measurable 0.5-10% window, then report it for both the magnitude MLE (95% CI for host fraction excludes 0) and the presence test (`lrt_pval<0.05`), with blank-free per-sample detection rules. | Real reads keep what a simulator cannot fully reproduce: true per-marker capture bias, real between-marker overdispersion, and this dataset's co-pooled contamination floor. Complements the EP17-A2 *analytical* LoD (a best-case Fisher-information ceiling on near-binomial simulated data) with an empirical LoD curve on real data, which is the more defensible number for a reviewer. A faint information-theoretic overlay (per-sample analytical `lod_fraction`) is carried as a consistency check. |
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
- The contamination floor is applied as a flat per-sample scalar by default. Per-marker
  dose-response apportionment (Step 30, #30) is now implemented but ships **off by default**, so
  the default path is still the flat scalar. It carries a known risk of inverting the
  presence-vs-MLE QC gate, which is why it is gated and opt-in.
- The LoD comparison against commercial assays is not head-to-head; allomix's number is an
  EP17-A2 analytical ceiling on near-binomial simulated data.

---

## Code-change audit since the inventory was written (2026-06-10 -> 2026-06-30)

The ranked table above was written on 2026-06-10. Reviewing every `src/allomix` commit since
then (and the current uncommitted working tree):

- **Most of the working-tree diff is non-functional.** The uncommitted changes across 14 modules
  are comment/docstring trimming (verbosity, AI-tells, line length) with no logic change. The
  `use_sex_chroms` parameter, the one-sided trim, the relatedness verdict logic, and the
  contamination math are all unchanged. Nothing in the working tree warrants an inventory edit.
- **Three substantive features landed as commits** since 2026-06-10 and are now added to the table
  above: per-marker-type two-rho overdispersion (#33, rank 8), per-marker dose-response
  contamination correction (Step 30 / #30, rank 7), and real-data subsample LoD on SRP434573
  (#24, rank 7).
- **Non-inventory additions:** a standalone HTML report (#27) and "total markers alongside
  informative count" (#15) are reporting/UX, not methods-and-results, so they are not ranked here.
  An MLE speedup (hoisting f-invariant work out of the rho loops) and a forced hom-ref background
  call (#23) are plumbing already implied by existing rows (the latter feeds the rank-7 empirical
  error null).

## Paper coverage audit (rebuilt 2026-06-30)

How many words the current paper (`paper/*.md`) spends on each ranked feature. Counts are
*attributed estimates*: sentences/passages describing each topic were tagged and their words
summed, rounded to ~10. Template placeholders (`{{ ... }}`) count as one word. Where a passage
serves two topics the words are split by emphasis, so column totals are approximate. Per-section
audits were run independently, so cross-section splits are not coordinated. Section word totals
(raw `wc -w`, including markdown/figure scaffolding): abstract 379, introduction 820, methods
4209, results 4716, discussion 2271, supplementary 3613 (16,008 total). About 65% of words map to
a ranked feature; the rest is clinical framing, longitudinal-accuracy validation, and table/figure
scaffolding.

### Per-feature word count (sorted by total)

| Rank | Feature | Abs | Intro | Meth | Res | Disc | Supp | **Total** |
|---:|---|--:|--:|--:|--:|--:|--:|--:|
| 9 | SRP434573 real-data + contamination floor (data property) | 30 | 0 | 300 | 560 | 120 | 320 | **1330** |
| 10 | Host-presence detection (relapse early-warning) | 40 | 20 | 210 | 340 | 220 | 170 | **1000** |
| 10 | Beta-binomial overdispersion (shared-rho) | 0 | 0 | 190 | 40 | 240 | 410 | **880** |
| 5 | Self-consistent simulator | 0 | 0 | 340 | 70 | 10 | 420 | **840** |
| 6 | EP17-A2 **analytical** LoB/LoD | 30 | 0 | 110 | 420 | 90 | 90 | **740** |
| 7 | **Per-marker contamination correction (Step 30, new)** | 0 | 0 | 140 | 240 | 10 | 270 | **660** |
| 6 | Multi-donor estimation | 20 | 0 | 80 | 330 | 90 | 70 | **590** |
| 8 | CNV/LoH + robust-refit recovery | 0 | 0 | 170 | 210 | 0 | 100 | **480** |
| 8 | **Two-rho per-marker-type overdispersion (new)** | 0 | 0 | 190 | 180 | 70 | 0 | **440** |
| 7 | Logit-space bias correction | 0 | 0 | 120 | 60 | 60 | 190 | **430** |
| 8 | In-data contamination **estimator method** | 0 | 0 | 100 | 80 | 30 | 160 | **370** |
| 7 | **Real-data subsample LoD (new)** | 0 | 0 | 0 | 350 | 0 | 0 | **350** |
| 8 | Two-phase calling (GATK GT + bcftools AD) | 0 | 0 | 230 | 10 | 90 | 0 | **330** |
| 7 | Three-state QC + anti-gaming GoF + presence cross-check | 10 | 10 | 130 | 70 | 70 | 0 | **290** |
| 8 | Profile-likelihood CIs | 5 | 0 | 80 | 40 | 20 | 140 | **285** |
| 3 | Vynck marker-type taxonomy | 0 | 0 | 140 | 0 | 20 | 110 | **270** |
| 9 | One-sided robust trim | 0 | 0 | 150 | 70 | 0 | 0 | **220** |
| 3 | Reproducible nested validation design | 0 | 0 | 110 | 0 | 0 | 110 | **220** |
| 4 | Admix-consistency swap / third-genome test | 0 | 10 | 80 | 10 | 30 | 50 | **180** |
| 7 | Per-marker direction-specific empirical error null | 0 | 0 | 60 | 15 | 10 | 40 | **125** |
| 6 | Relatedness kinship QC (the *check*, not its LoD effect) | 0 | 10 | 90 | 15 | 10 | 0 | **125** |
| 5 | Read-level artifact filter | 0 | 0 | 90 | 20 | 0 | 0 | **110** |
| 4 | Index-hopping metadata flag | 0 | 0 | 40 | 10 | 0 | 20 | **70** |
| 2 | GT/AD consistency guard | 0 | 0 | 70 | 0 | 0 | 0 | **70** |
| | **Attributed** | 135 | 50 | 3220 | 3140 | 1190 | 2670 | **10,405** |

(The Results "effect of donor-host relatedness on marker count / LoD" subsection, ~200 words, is a
*different thing* from the rank-6 relatedness QC check and is left uncounted, as in the prior audit.)

### What changed vs the 2026-06-10 audit, and where to edit

**Good news, several prior gaps are now closed.** The previous audit flagged five features as
absent or mis-described; the rebuild shows them fixed:

- *One-sided robust trim* (rank 9): was "mis-described as a generic median/MAD refit." Methods now
  describes the asymmetric host-direction-protected property (~150 w). **Resolved.**
- *In-data contamination estimator* (rank 8): the estimator method (median per-site minor, p10
  floor, miscall cap) is now in Methods (~100 w), not just the Results dataset story. **Resolved.**
- *Relatedness kinship QC* (rank 6): was 0 / absent, now ~125 w (Methods ~90). **Resolved.**
- *Read-level artifact filter* (rank 5): ~15 -> 110 w. *Swap test* (rank 4): 0 -> 180 w.
  *GT/AD guard* (rank 2): 0 -> 70 w. **All now present.**
- *Two-phase calling* is correctly described as two-phase (the old "v1 still single-phase" bug
  stays fixed).

**The headline rebalanced in the right direction.** In the old audit the EP17 **analytical** LoD
was the single most-covered topic (~1300 w). It is now 740 w and demoted below the two genuinely
load-bearing, clinically novel topics: SRP real-data validation (1330) and host-presence detection
(1000). That is the correct ordering for a JMD adoption pitch.

**Spending too long on something boring (trim candidates):**

1. **Self-consistent simulator (840 w, rank 5)** is now the 4th-largest spend, behind only the
   three headline topics, and it is machinery, not a finding. The credibility argument (estimator
   tested against a superset of its own assumptions) is worth ~1-2 paragraphs; 340 w in Methods +
   420 w in Supp is more than that argument needs. Trim the Supp duplication.
2. **EP17 analytical LoD (740 w, rank 6)** is still heavy for a number the paper itself labels a
   best-case ceiling on near-binomial data. Now that **real-data subsample LoD (#24)** exists, the
   analytical version should *lean on* the real one and shrink: keep the EP17 framing, cut the
   per-corner analytical-number enumeration in Results (~420 w), and let Fig 5 (real LoD) carry it.
3. **Logit bias correction (430 w, rank 7)** was flagged "over-covered relative to impact" in the
   prior audit and still is: it is a precision-only refinement (tightens CIs ~30%, barely moves the
   point estimate). 430 w across four sections is too much for that. Consolidate to one place.

**Not bringing up something interesting enough (under-exposed):**

1. **Real-data subsample LoD (#24, 350 w, Results only)** has *no Methods description* and no
   Discussion tie-in. A real-data LoD is the most reviewer-defensible sensitivity number in the
   paper, yet the method for producing it (depth/marker subsampling of real reads, blank-free
   detection rule) is undocumented. **Add a short Methods paragraph and a Discussion sentence**,
   funded by trimming the analytical-LoD enumeration (point 2 above).
2. **Per-marker direction-specific empirical error null (rank 7, 125 w)** is still essentially
   *pending*: the paper mentions it only as a "future item" / "per-site error model" aside and
   still operates on the flat ε=0.01 4-state model. Either it is implemented (then describe the
   real e_refalt/e_altref null) or it is not (then stop implying it is close). This is the largest
   remaining "high-rank, low-words, and framed as vapor" gap.
3. **One-sided robust trim (rank 9, 220 w)** is now correctly *described* but still thin for a
   rank-9, genuinely novel domain-aware estimator. Given how much the low-fraction story leans on
   it (it is what recovers the CNV/LoH donor LoD), it could absorb ~50-80 w more, ideally a
   sentence in Discussion on *why* asymmetry is the right call at the LoD.
4. **New per-marker contamination correction (Step 30, 660 w)** is well covered but is **off by
   default** — make sure the paper frames it as optional/opt-in and does not let 660 w read as a
   headline default behaviour. The risk is over-claiming a correction the default pipeline does not
   apply.

**Net:** the paper is now spending its words on the right top-three (real data, host presence,
overdispersion). The remaining edits are second-order rebalancing: move ~300-400 w out of the
simulator + analytical-LoD + bias-correction machinery and into (a) a Methods paragraph for the
real-data subsample LoD and (b) resolving the per-marker error-null "is it in or not" ambiguity.
