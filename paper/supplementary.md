# Supplementary Data

## Supplementary Methods

This section gives the statistical detail summarised in plain language in the main
Methods. The notation is shared across subsections: for informative marker *i*,
$g_{r,i}$ and $g_{d,i}$ are the reference-allele doses (0, 1, or 2) of recipient and donor,
$n_i$ is the total read count, and $k_i$ is the alternative-allele count.

### S1. Mixture model and marker classification

For a proposed donor fraction *f*, the expected reference-allele weight at marker *i* is
a blend of the two known genotypes:

$$w_i(f) = (1 - f)\,\frac{g_{r,i}}{2} + f\,\frac{g_{d,i}}{2}$$

Markers are informative for the fraction estimate when recipient and donor genotypes differ.
Following Vynck et al.,[@Vynck2023bias] each informative marker is assigned one of six
types by recipient/donor alternative-allele dose (the six types and which part of allomix uses
each are tabulated in the main Methods, "Which markers are informative").

Markers where recipient and every donor share a genotype are non-informative for the fraction
but feed the QC checks (S8). At consensus-homozygous markers (all parties homozygous for
the same allele) the minority allele can only be a background artifact or foreign DNA,
which drives contamination estimation and the sample-swap test. At consensus-heterozygous
markers (all parties heterozygous) the admixture alternative-allele fraction should sit
near 0.5 whatever the mixing fraction, so a systematic skew is a separate signal for
contamination, allelic imbalance, or a sample mix-up. Default filters: recipient and donor GQ
$\geq$ 20, admixture DP $\geq$ 100, and at least three informative markers.

### S2. Sequencing-error model

To account for base substitutions and polymerase errors, the observed allele
probabilities use a 4-state (trinucleotide) model:

$$p_{alt,i} = (1 - w_i)(1 - \varepsilon) + w_i\frac{\varepsilon}{3}, \qquad p_{ref,i} = w_i(1 - \varepsilon) + (1 - w_i)\frac{\varepsilon}{3}$$

where $\varepsilon$ is the per-base error rate (default 0.01) and the factor of 3
distributes error among the three non-observed bases. Because VCF allele-depth fields
count only reference and alternative, the likelihood uses the conditional $\tilde{p}_i =
p_{alt,i} / (p_{ref,i} + p_{alt,i})$. This is the mixture genotype likelihood of Crysup
and Woerner[@CrysupWoerner2022] applied to the reverse problem (estimating the
fraction from known genotypes rather than genotyping at a known fraction). The single
rate $\varepsilon$ may be replaced per marker by a measured per-site, per-direction
empirical rate when an error table is supplied (S7, Methods); each marker then uses its
own REF-to-ALT or ALT-to-REF rate in place of $\varepsilon$, falling back to the global
value where a direction was not measured.

### S3. Beta-binomial likelihood and optimization

A binomial model assumes all variance comes from read sampling; in practice per-marker
amplification bias and depth variability produce overdispersion. allomix uses a
beta-binomial,[@HindeDemetrio1998] parameterised by *f* and a shared concentration $\rho
> 0$. The per-marker log-likelihood (up to a constant) is:

$$\ell_i(f, \rho) = \log\Gamma(k_i + \alpha_i) + \log\Gamma(n_i - k_i + \beta_i) - \log\Gamma(n_i + \rho) - \log\Gamma(\alpha_i) - \log\Gamma(\beta_i) + \log\Gamma(\rho)$$

with $\alpha_i = \tilde{p}_i\,\rho$ and $\beta_i = (1 - \tilde{p}_i)\,\rho$. As $\rho
\to \infty$ this converges to the binomial; smaller $\rho$ flattens the likelihood and
widens intervals. The total log-likelihood is $\mathcal{L}(f,\rho) =
\sum_{i=1}^{M}\ell_i(f,\rho)$. Both parameters are fit jointly: a grid search over 1,001
evenly spaced values of *f* in [0, 1] with $\rho$ profiled out (bounded Brent on
log-scale) at each point, followed by Nelder-Mead refinement over $(f, \log\rho)$ from
the grid maximum.

### S4. Profile-likelihood confidence intervals

The 95% interval for *f* inverts the profile log-likelihood $\mathcal{L}_P(f) =
\max_\rho \mathcal{L}(f,\rho)$, with bounds where

$$2\left[\mathcal{L}_P(\hat{f}) - \mathcal{L}_P(f)\right] = \chi^2_{1,\,0.95} \approx 3.84$$

found by Brent root-finding.[@Wilks1938] The reference maximum is re-derived from the
same profiled optimizer (the profile value at $\hat{f}$, not the joint optimum) so the
root-finder brackets a sign change. Bounds are pinned at 0 and 1, so an estimate near a
boundary does not produce an interval running past it.

### S5. Multi-donor extension

For two donors with fractions $f_1, f_2$ (recipient the remainder), the weight at marker *i*
is

$$w_i(f_1, f_2) = (1 - f_1 - f_2)\,\frac{g_{r,i}}{2} + f_1\,\frac{g_{d1,i}}{2} + f_2\,\frac{g_{d2,i}}{2}$$

A marker is informative if the recipient differs from any donor; per-donor informative counts
are tracked separately. Optimization is a triangular grid over the simplex $\{(f_1,
f_2): f_1, f_2 \geq 0,\; f_1 + f_2 \leq 1\}$ at 101 steps per dimension (~5,150
evaluations), then Nelder-Mead refinement. Per-donor 95% profile intervals use
$\chi^2_{1,0.95}$ (each interval profiles one donor while optimising the other).

### S6. Per-marker bias correction

Per-marker bias is estimated at heterozygous training observations as $b_i =
\text{median}(\text{VAF}_{het,i} - 0.5)$. A flat additive shift $w_i - b_i$ is valid
only near 0.5 and overcorrects at the extreme expected weights that dominate
low-fraction samples, so the correction is multiplicative in logit space:

$$w'_i = \text{expit}\!\left(\text{logit}(w_i) - \text{logit}(0.5 + b_i)\right)$$

clamped to $[10^{-6}, 1 - 10^{-6}]$. At a heterozygous site ($w_i = 0.5$) this reduces
to $0.5 - b_i$; at an extreme weight it is a small proportional shift. A marker is
correctable only where its bias was measured (where it was heterozygous), so the table
is built across other samples, either a reference cohort called the same way as the
admixture, or admixture samples at markers where recipient and every donor are heterozygous
(true VAF 0.5 regardless of mixing).

### S7. Residual-recipient presence test

The presence test uses only markers where the donor is homozygous and the recipient carries
the donor-absent allele (Vynck types 0, 1, 10, 11). Let $y_i$ be the donor-absent allele
count out of $n_i$ reads, $e_i$ the per-marker error background in that direction, and
$r_i$ the recipient dose of the donor-absent allele. Under a recipient fraction $f_r$ the expected
donor-absent allele probability is

$$q_i(f_r) = e_i + \frac{r_i}{2}\,f_r$$

Two statistics are reported. A pooled one-sided Poisson test uses $Y = \sum_i y_i$
against $\Lambda = \sum_i n_i e_i$, with $p = P(\text{Poisson}(\Lambda) \geq Y)$. A
bounded-MLE likelihood-ratio test maximises a per-marker binomial likelihood in
$q_i(f_r)$ over $f_r \geq 0$. Because the null sits on the boundary $f_r = 0$, the LRT
p-value uses a chi-bar-square reference (a 50:50 mixture of a point mass at 0 and
$\chi^2_1$) rather than naive Wilks, and the reported confidence interval is the
profile-likelihood interval clipped at 0. The test thus returns a p-value, a
recipient-fraction estimate $\hat{f}_r$, and a CI, and is calibrated against the per-marker
error background ($e_i$). When a per-site, per-direction empirical error table is
supplied (Methods, `estimate-errors`), $e_i$ is that marker's measured rate in the
donor-absent direction ($e_{refalt}$ where the donor is homozygous-reference,
$e_{altref}$ where homozygous-alternative), with a load-time floor so a zero observed
rate cannot make a stray read produce an infinite penalty; markers without a measured
rate fall back to the symmetric global $e_i = \varepsilon/3$. The in silico and
SRP434573 results use the global rate, because neither the simulator's uniform error
model nor this single public dataset provides a calibrated per-site table.

### S8. In-data contamination estimation

Contamination is measured at the consensus-homozygous markers (S1), where the minor
allele can only be a background artifact or foreign DNA. The headline estimate is the
background-subtracted median per-site
minor-allele fraction: the median is used rather than a pooled mean so a few gross
miscall sites do not dominate, sites above 10% minor fraction are capped as miscalls,
and the error floor is the 10th percentile of per-site minor fractions (the
no-carrier/error sites), so contamination is reported as the heterogeneous excess over a
uniform error floor (a uniform error elevation lifts the floor too and is correctly not
called contamination). Contamination is distinguished from real low-level chimerism by
marker geometry rather than magnitude: a dose-response in which the minor fraction rises
with the number of co-pooled panel individuals carrying that allele indicates foreign
reads, whereas a flat elevation indicates error. A separate sample-swap / foreign-genome
test runs at the same consensus sites: a per-site binomial tail at the error rate flags
sites where the minor allele is individually significant, combined into a swap p-value
over discordant sites, catching a wrong-patient VCF that the informative-marker
goodness-of-fit never sees. A complementary allele-balance check runs at
the consensus-heterozygous markers (S1) instead;
the fraction of markers skewed outside a 70:30 band is reported and promotes the sample
to manual review above a set count, flagging contamination, copy-number or allelic
imbalance, or a sample mix-up that the consensus-homozygous checks do not test. An
optional `##allomixRunUnit` VCF header (flowcell:lane) supports a pure-metadata
index-hopping flag, kept separate from the in-data estimate.

The optional per-marker contamination correction (off by default) acts on the magnitude
estimate rather than the contamination report. On a co-pooled run a donor-homozygous
informative marker carries extra reads on the recipient (donor-absent) allele from co-pooled
genomes that happen to carry it, scaling with the number of those co-pooled carriers;
the recipient signal is the same at every such marker while this contamination scales with
carrier dose. The correction subtracts a dose term, `slope * n_carriers * depth`, from
each donor-homozygous recipient-allele count before the fit, leaving the flat error floor to
the per-site error model so it is not double-counted. Two quantities are measured per
run, not assumed. The gate is the per-flowcell consensus-homozygous dose-response: the
minor-allele fraction at consensus-homozygous sites is regressed on the co-pooled
carrier count (weighted by depth, pooled across the run's admixtures with a
per-admixture intercept), and the correction is applied only when that slope is
significantly positive; a clean run has a flat slope and the correction is a no-op. The
magnitude is calibrated separately on the informative donor-homozygous markers
themselves (the same weighted dose regression), because the consensus-homozygous slope
predicts the informative-marker slope well enough to gate on but not to transfer
one-for-one. The carrier counts come from the cohort's joint-called genotypes, the same
input as the per-site error table. Only donor-homozygous markers (genotype-contrast
types where the donor is homozygous) are corrected; the magnitude slope is clamped
non-negative, so the correction can only lower the estimate, which is why it is gated
rather than applied unconditionally.

### S9. Simulation, limit of detection, and copy-number model

The simulator draws each marker's expected alternative-allele frequency as
$\text{VAF}_{expected} = [(1-f)a_h + f a_d]/2$ (with $a_h, a_d$ the alternative doses),
applies the four noise sources of the main Methods, and draws counts from a binomial (or
a beta-binomial at concentration $\rho$ for the overdispersion characterisation). The LoD
sweep uses a nested design: for each relatedness level, multiple donor/recipient pairs
({{ lod_headline.n_pairs_unrelated | int }} unrelated,
{{ lod_headline.n_pairs_sibling | int }} sibling) have fixed genotypes (MAF 0.2--0.5) and
per-marker bias reused across all depths and panel sizes, panels strictly nested (a
smaller panel is a bit-identical prefix of a larger one), and
{{ lod_headline.n_seq_reps | int }} sequencing replicates per cell varying only
read-sampling noise. Holding the pair fixed makes each pair's LoD curve monotonic in
panel size while the across-pair identity-by-descent spread is reported as a band rather
than leaking into the central estimate; seeds are SHA-256-derived for reproducibility.
For recipient copy-number variants (CNVs) the recipient is a mixture of normal diploid cells
and a CNV-bearing clone at a clonal fraction, with the expected allele fraction a
copy-number-weighted average over normal recipient, recipient clone, and donor; three
CNV types are produced by mutating one clone homolog (copy-neutral loss of heterozygosity
(cnLoH): retained homolog duplicated, heterozygous sites only; deletion: one copy; gain:
three copies), deletion and gain also changing the locus DNA contribution at homozygous
sites. The CNV is applied only to the admixture sample.

### S10. QC gating on effect size, not significance

The QC layer assigns each result PASS, REVIEW, or FAIL. Three of the REVIEW checks (the
model goodness-of-fit, the pre-trim guard against a robust trim hiding a poor fit, and
the consensus-homozygous swap test) are, at heart, significance tests, and a significance
test is the wrong instrument at panel depth. With hundreds of markers at more than 1000x,
a chi-squared statistic reaches p < 0.01 for a departure from the model far too small to
change the reported fraction, so a p-value threshold flags almost every real sample. On
the SRP434573 two-person mixtures the goodness-of-fit p-value is below 0.01 for
{{ qc_gating.gof_p_significant }} of {{ qc_gating.n_timepoints }} timepoints, yet the
beta-binomial model actually fits: the post-trim reduced chi-squared (chi-squared divided
by degrees of freedom, which is ~1 for a good fit and is depth-independent) has a median
of {{ qc_gating.gof_reduced_median }} and a maximum of {{ qc_gating.gof_reduced_max }}
across the series.

allomix therefore promotes a result to REVIEW only when the misfit is large as well as
significant. The goodness-of-fit and pre-trim checks require the reduced chi-squared to
exceed {{ qc_gating.gof_reduced_threshold }} (which sits just above the maximum observed
on this well-behaved dataset, so it does not fire on model-consistent scatter but will
catch a materially worse fit); the pre-trim guard additionally fires only when the recipient
fraction is within {{ qc_gating.pretrim_lod_multiple }}x the limit of detection, the
regime where the outlier-resistant trim could have discarded real low-fraction recipient
signal rather than artifacts; and the swap test requires the discordant fraction to
exceed {{ qc_gating.swap_review_fraction_pct }}%, well below the roughly one-half of
consensus-homozygous sites a genuine unrelated sample swap mismatches but far above the
handful of off-model sites (a genotyping miscall, a mapping artifact, or a site-specific
contaminant) that a co-pooled research dataset produces. These are the ``clinical_gating``
defaults; the legacy p-value-only rules remain available (``--no-clinical-gating``) for
reproducing the stricter behaviour.

The effect is to move the flags to where they change the interpretation. On the SRP434573
series clinical gating promotes {{ qc_gating.n_review_clinical }} of
{{ qc_gating.n_timepoints }} timepoints to REVIEW, all at the lowest recipient fractions
(0.5% to 1%), where the outlier-resistant refit set aside a large share of markers near
the detection floor. The legacy p-value-only rules would flag
{{ qc_gating.n_review_legacy }} of {{ qc_gating.n_timepoints }} on the identical fits,
including accurate estimates from 2.5% to 100% recipient that need no review.

### S11. Per-read and per-genotype weighting refinements considered and not adopted

Three per-marker weighting refinements were considered and left out, because each targets
a noise source that is not the low-fraction limiter at clinical depth.

Genotype-quality weighting would replace the hard `--min-gq` cutoff with a per-marker
weight from the Phred genotype posterior, keeping borderline calls but down-weighting
them. On a panel run above 1000x nearly every recipient and donor call is GQ 99, so the weight
is close to 1 at almost every marker and the method reduces to the current hard filter.
The only markers it recovers are the handful at GQ 10 to 20 on weak recipient genotyping, which
do not move the overdispersion and contamination floor that sets low-fraction accuracy.

Per-base-quality weighting (as in Conpair) would weight each read by its base quality
through a per-marker effective error rate. This is the most invasive change for the
smallest expected gain. The production VCF carries no per-base quality, so it would need
an extra `bcftools mpileup -a FORMAT/QS` pass at the panel sites, and on a Q30+ panel the
resulting per-marker error barely differs from the flat `--error-rate` already used. The
mean-Phred-to-error conversion also underestimates the true error under Jensen's inequality
when base qualities are heterogeneous. As with genotype-quality weighting, the dominant
low-fraction noise is overdispersion and contamination, neither of which this addresses.

An in-likelihood contamination term would add a per-marker contamination rate as a nuisance
parameter inside the MLE and joint-fit it with the recipient fraction, replacing the
pre-subtraction correction of Supplementary Methods S8. Recipient fraction and contamination
both add reads to the recipient allele and are partially confounded, so identifiability would
have to come from the co-pooled carrier-dose structure (contamination scales with the
number of co-pooled carriers, recipient signal does not), adding model complexity for an
interior refinement. The pre-subtraction together with the two independent readouts, the
residual-recipient presence test and the consensus-homozygous contamination report, already
separates these signals without that coupling.

## Supplementary Table S1. Empirical Panel Characterisation

Per-marker amplification bias, depth distribution, locus dropout, and allele dropout
were measured from {{ panel_empirical.n_vcfs | fmt('g') }} joint-called VCFs
({{ panel_empirical.n_samples | commas }} samples) generated from the
{{ panel_specs.n_markers_panel }} sample-identification SNPs of the IDT rhAmpSeq Sample
ID panel (Integrated DNA Technologies), incorporated into a custom hematology capture
panel and sequenced on Illumina instruments (for example NovaSeq) to a mean depth of
{{ panel_specs.typical_depth }} as part of routine clinical sequencing. All
{{ panel_empirical.n_bias_markers | fmt('g') }} biallelic
markers with heterozygous observations were included. Simulation parameters used
throughout this study were calibrated from these measurements.

| Parameter | Empirical Value | Simulation Default | Notes |
|:---|:---:|:---:|:---|
| **Amplification bias** | | | |
| Per-marker bias SD ($\sigma_{bias}$) | {{ panel_empirical.sd_bias }} | 0.02 | Close match |
| Mean \|bias\| | {{ panel_empirical.mean_abs_bias }} | | |
| Median \|bias\| | {{ panel_empirical.median_abs_bias }} | | |
| 95th percentile \|bias\| | {{ panel_empirical.p95_abs_bias }} | | Heavy tail |
| Max \|bias\| | {{ panel_empirical.max_abs_bias }} | | Single outlier marker |
| **Sequencing depth** | | | |
| Mean depth | {{ panel_empirical.mean_depth | commas }}x | 50–1,000x | Sims test lower depths |
| Median depth | {{ panel_empirical.median_depth | commas }}x | | |
| Min mean depth | {{ panel_empirical.min_depth | fmt('g') }}x | | Weakest marker |
| Max mean depth | {{ panel_empirical.max_depth | fmt('g') }}x | | |
| Per-sample depth CV | {{ panel_empirical.mean_sample_depth_cv }} | 0 (uniform) | Not yet modelled |
| **Locus dropout** | | | |
| Mean no-call rate | {{ panel_empirical.mean_nocall_pct }}% | 0% | Not yet modelled |
| Markers with >5% no-call | {{ panel_empirical.markers_gt5pct_nocall | fmt('g') }}/{{ panel_empirical.n_bias_markers | fmt('g') }} | | Single problematic marker |
| **Allele dropout** | | | |
| Mean het/HWE ratio | {{ panel_empirical.mean_het_ratio }} | | 1.0 = no ADO |
| Markers with ratio < 0.8 | {{ panel_empirical.markers_low_het | fmt('g') }}/{{ panel_empirical.n_bias_markers | fmt('g') }} | | Negligible at high depth |
| Estimated ADO rate | {{ panel_empirical.ado_estimate }} | 0 | Negligible |

## Supplementary Table S2. Per-Marker Detail

Per-marker statistics are available in the allomix repository at
`paper/empirical_results/panel_per_marker.tsv`. Fields include: number of observations,
call rate, genotype counts (hom-ref, het, hom-alt), observed-to-expected heterozygosity
ratio (HWE), mean depth, depth CV, and median amplification bias. Marker identities are
anonymised (sequential index only).

## Supplementary Table S3. Per-Sample Validation Results by Depth

Detailed per-sample validation results for each sequencing depth (50x, 100x, 200x, 500x,
1,000x) are available in the allomix repository at `output/depth_validation/`. For each
depth, the true donor fraction, estimated fraction, error, and 95% confidence interval
bounds are reported for all simulated mixture levels.

## Supplementary Table S4. In Silico Accuracy and CI Performance by Depth

We generated synthetic chimeric VCFs spanning 0% to 100% donor using 100 markers with
empirically calibrated per-marker bias, non-uniform depth (CV =
{{ sim_calibration.depth_cv }}), and {{ sim_calibration.locus_dropout_pct }}% locus
dropout, at five depths from 50x to 1,000x ({{ depth_50.n_replicates | dp(0) }}
replicates each). Mean absolute error stayed below 1% at every depth, improving from
{{ depth_50.mean_abs_error_pct | dp(2) }} ± {{ depth_50.mean_abs_error_sd_pct | dp(2) }}%
at 50x to {{ depth_1000.mean_abs_error_pct | dp(2) }} ±
{{ depth_1000.mean_abs_error_sd_pct | dp(2) }}% at 1,000x. The 95% profile-likelihood
intervals showed coverage of
{{ depth_1000.ci_coverage_pct }}--{{ depth_200.ci_coverage_pct }}%, close to the nominal
95% across the depth range, with the per-marker-type overdispersion model (separate
donor-homozygous and donor-heterozygous concentration) absorbing the extra-binomial
scatter that a single shared term left in the residuals. Per-marker bias correction is
available but modest by design: in a 2,000x experiment it left the point estimate
essentially unchanged (0% and 100% donor estimated at {{ bias_with_bias.est_0pct }}% and
{{ bias_with_bias.est_100pct }}%) while narrowing the mean interval from
{{ bias_no_bias.mean_ci_width_pct | dp(2) }}% to
{{ bias_with_bias.mean_ci_width_pct | dp(2) }}% (Supplementary Methods S6, Figure S4).

<!-- include-csv: output/facts/table_depth.csv
  align: center
-->


**Table S4.** allomix accuracy and confidence-interval performance across sequencing
depths (mean ± SD, N={{ depth_50.n_replicates | dp(0) }} replicates). MAE = mean absolute
error; RMSE = root mean square error. Error metrics are computed on interior fractions
(excluding 0% and 100%). The per-depth agreement scatter is Figure S14, and the
absolute-error boxplots and depth-summary panels are Figures S10 and S11.

## Supplementary Figures: Simulation Model Validation

### S1. Amplification Bias Distribution

![Figure S1]({{ facts_dir }}/fig_bias_distributions.png)

**Figure S1.** Per-marker amplification bias distribution across
{{ supp_synthetic.n_empirical_markers | dp(0) }} markers (median het VAF deviation from
0.5). The heavy-tailed mixture model used in simulation (95% N(0, 0.012), 5% N(0, 0.08))
tracks the empirical tail (95th percentile of |bias|
{{ supp_synthetic.empirical_p95_abs_bias }}), which a simple Gaussian underestimates. (A)
Histogram with kernel density estimates; (B) cumulative distribution of |bias|.

### S2. Depth Distribution

![Figure S2]({{ facts_dir }}/fig_depth_distributions.png)

**Figure S2.** Per-marker sequencing depth. (A) Empirical mean depth per marker vs
log-normal model draws at the same mean and CV. (B) Within-marker depth CV across
samples for each marker, showing the range of per-marker depth variability.

### S3. Heterozygous VAF Comparison

![Figure S3]({{ facts_dir }}/fig_het_vaf.png)

**Figure S3.** Violin plots of median heterozygous VAF per marker, empirical
({{ supp_synthetic.n_empirical_markers | dp(0) }} markers) versus simulated: both centred
on 0.5 with comparable spread, so the simulation reproduces the per-marker VAF
displacement seen in real data.

### S4. Noise Component Ablation

![Figure S4]({{ facts_dir }}/fig_ablation.png)

**Figure S4.** Effect of individual noise components on estimation accuracy (500x, 10
replicates, 7 conditions). (A) Overall RMSE by condition: ideal
{{ supp_synthetic.ablation_rmse_ideal_pct }}%; amplification bias alone
{{ supp_synthetic.ablation_rmse_bias_only_pct }}%, essentially unchanged by bias
correction ({{ supp_synthetic.ablation_rmse_bias_corrected_pct }}%) since injected biases
average near zero and their spread is absorbed by the overdispersion term. The full
binomial model (all noise sources, bias corrected) gives
{{ supp_synthetic.ablation_rmse_full_pct }}%, the no-overdispersion baseline; adding
per-marker overdispersion (beta-binomial at fitted rho = 100) raises it to
{{ supp_synthetic.ablation_rmse_overdispersion_pct }}%, a larger effect than any single
bias, depth, or sequencing-error component: the accuracy-side counterpart to
overdispersion controlling the LoD (Figures S7, S8). (B) Mean absolute error by true
donor fraction; dashed lines are bias-corrected conditions.

### S5. Confidence Interval Calibration

![Figure S5]({{ facts_dir }}/fig_ci_calibration.png)

**Figure S5.** CI calibration under the full noise model (100 replicates per fraction,
run as 10 parallel batches of 10 via Snakemake). (A) Observed 95% CI coverage rate by
true donor fraction; overall coverage is {{ supp_synthetic.cal_coverage_pct }}%. (B)
Mean CI width by true donor fraction, with standard deviation bars.

### S6. Per-Marker Residuals

![Figure S6]({{ facts_dir }}/fig_residuals.png)

**Figure S6.** Per-marker residuals (observed minus expected VAF) from a simulated 30%
donor mixture at 500x. (A) Residual histogram with normal fit. (B) Residuals plotted
against expected VAF, showing no systematic trend across the VAF range.

### S7. Limit of Detection vs Depth (Saturation)

![Figure S7]({{ facts_dir }}/fig_lod_saturation.png)

**Figure S7.** Limit of detection as a function of mean depth for the simulated
unrelated panels in the LoD sweep. Points are the in silico LoD per panel size; lines
are the LoD model $\mathrm{LoD} = (A/\sqrt{M})\sqrt{(n + \rho)/(n(\rho + 1))}$ fitted
across panels ($M$ = informative markers, $n$ = depth, $\rho$ = beta-binomial
overdispersion concentration). The simulator draws reads from a binomial (Methods), so
the fit returns a near-infinite $\rho$ and the LoD falls close to $1/\sqrt{n}$ with no
floor (dashed binomial reference). Under real, finite $\rho$ the per-marker variance
instead approaches $p(1-p)/(\rho + 1)$, so the LoD saturates at a floor and depth beyond
an effective cap of $\rho + 1$ reads yields diminishing returns. The in silico LoD
reported elsewhere is therefore an analytical best case under near-binomial sampling.

### S8. Effect of Overdispersion on the Limit of Detection

![Figure S8]({{ facts_dir }}/fig_overdispersion_lod.png)

**Figure S8.** In silico LoD as a function of the beta-binomial overdispersion
concentration $\rho$, at {{ overdispersion_lod_headline.depth }}x depth with
{{ overdispersion_lod_headline.n_markers }} informative markers (unrelated donor). Reads
were simulated beta-binomial across a grid of $\rho$ and the donor fraction estimated
with the standard pipeline; the analytic and simulated (tool) LoD agree closely. The LoD
rises from {{ overdispersion_lod_headline.lod_binomial_pct }}% under pure-binomial
sampling ($\rho \to \infty$) to {{ overdispersion_lod_headline.lod_rho100_pct }}% at
$\rho = 100$ (a {{ overdispersion_lod_headline.fold_rho100_vs_binomial }}-fold increase)
and {{ overdispersion_lod_headline.lod_rho30_pct }}% at $\rho = 30$. At clinical
coverage the overdispersion, not the depth, is the dominant control on the achievable
LoD. Fitting the estimator on the real SRP434573 mixtures
({{ overdispersion_lod_headline.n_real_mixed_samples }} two-component samples) gives a
median heterozygous-class concentration of $\rho \approx$
{{ overdispersion_lod_headline.real_rho_het_median }} (homozygous class $\rho \approx$
{{ overdispersion_lod_headline.real_rho_hom_median }}); the dashed line marks that fitted
value. At the fitted $\rho$ the in silico LoD is
{{ overdispersion_lod_headline.lod_at_real_rho_pct }}%, a
{{ overdispersion_lod_headline.fold_real_rho_vs_binomial }}-fold increase over the
pure-binomial value and close to the ~1% measured directly by subsampling the same
mixtures (Figure 3). Calibrating the simulation to the real per-sample fits, rather than
leaving it at the pure binomial, is what makes the in silico LoD a defensible performance
figure.

### S9. Fixed-Bias-Per-Marker Stability

![Figure S9]({{ facts_dir }}/fig_bias_stability.png)

**Figure S9.** Validation of the fixed-bias-per-marker assumption. Each of the
{{ supp_synthetic.n_empirical_markers | dp(0) }} markers is plotted by its absolute median
bias (systematic component) against its within-marker het-VAF SD across samples (random
component). The two are only weakly correlated (r =
{{ supp_synthetic.bias_stability_r | dp(2) }}), so bias behaves as a stable per-marker
offset rather than a quantity that grows with noise, supporting a fixed offset
(Supplementary Methods S6) with residual scatter absorbed by the overdispersion term.

### S10. Absolute Error by Depth (Boxplots)

![Figure S10]({{ facts_dir }}/fig_depth_boxplots.png)

**Figure S10.** Distribution of absolute estimation error by sequencing depth
(N={{ depth_50.n_replicates | dp(0) }} replicates per depth). Boxes show median and
interquartile range for interior fractions (excluding 0% and 100% donor). Whiskers
extend to 1.5× IQR. This is the per-fraction distribution behind the summary metrics in
Supplementary Table S4.

### S11. Depth-Performance Summary

![Figure S11]({{ facts_dir }}/fig_depth_summary.png)

**Figure S11.** allomix performance as a function of sequencing depth (mean ± SD,
N={{ depth_50.n_replicates | dp(0) }} replicates). Left: accuracy metrics (MAE, RMSE,
maximum error). Centre: 95% profile-likelihood CI coverage versus the nominal 95% level
(dashed). Right: mean CI width.

![Figure S12]({{ facts_dir }}/fig_srp434573_logy.png)

**Figure S12.** Confidence-interval view of the SRP434573 two-person dilution series
(main-text Figure 3A shows the same data as a log-log scatter). Each admixture is
plotted on a log recipient-fraction axis grouped by mixture, with the maximum-likelihood
estimate (filled circle, 100 minus donor%, with the per-marker contamination correction
applied) and the residual-recipient presence-test estimate (open square) each shown with its
95% confidence interval, against the known fraction (grey diamond). The dashed
horizontal line in each mixture is that mixture's independent in-data contamination
level, measured at consensus-homozygous markers (a marker class the magnitude estimate
never reads), so it is a floor estimated from different sites than the ones being
corrected: an estimate below its mixture's line is not separable from contamination. The
pure-donor (true-0%-recipient) endpoints, at the right of each group, sit at or near the 0
row after correction, at or below their contamination line. At the higher titration
levels the intervals are tight and bracket the known value; at the 0.5% level they widen
toward the contamination line, where the residual floor competes with the true recipient
signal (Results, Figure 3 caption).

![Figure S13]({{ facts_dir }}/fig_srp_contam.png)

**Figure S13.** Co-pooled contamination floor in SRP434573 as a dose-response, the
figure behind the median values in Results. Consensus-homozygous sites are those where
recipient and donor are both homozygous for the same allele, so the minor allele at that site
cannot come from either contributor. For each such site, reads were pooled across the
dilution samples (sites with pooled depth below 500 were dropped for a stable per-site
fraction) and the per-site minor-allele fraction was computed. Sites are binned on the
x-axis by the number of the five other co-pooled panel individuals that carry the minor
allele, het or homozygous counting equally as one carrier:
{{ srp_contam.n_nocarrier_sites | dp(0) }} sites have no carrier (the no-carrier bin)
and {{ srp_contam.n_carrier_sites | dp(0) }} have at least one. Boxes show the per-site
distribution on a log y-axis (median, interquartile range, 1.5x IQR whiskers, individual
outliers; exact-zero sites are drawn at a 0.001% floor so they render on the log axis);
the red line connects the per-bin medians and n below each box is the number of sites.
The median rises monotonically from the no-carrier floor (sequencing error,
{{ srp_contam.nocarrier_floor_pct | dp(3) }}%) through the carrier bins, the signature
of real reads from co-pooled material (most plausibly index hopping) rather than flat
sequencing error, which would not scale with co-pooled dose.

### S14. In Silico Accuracy Across Sequencing Depths

![Figure S14]({{ facts_dir }}/fig_depth_scatter.png)

**Figure S14.** In silico accuracy across sequencing depths, the per-depth agreement
behind Supplementary Table S4. Each panel shows true donor fraction (x-axis) versus
estimated donor fraction (y-axis) at the indicated depth, all replicates
(N={{ depth_50.n_replicates | dp(0) }}). Dashed line indicates perfect agreement. This is
the estimator's controlled-simulation accuracy across depth; the real-data accuracy is in
main-text Figures 2 and 3.

### S15. Error-Table Provenance and Low-Fraction Sensitivity

![Figure S15]({{ facts_dir }}/fig_error_table_arms.png)

**Figure S15.** Effect of the per-site error model on low-fraction performance,
on the SRP434573 public mixtures. Three background models are compared: the
built-in flat per-base `--error-rate` default; one `estimate-errors` table per
mixture, built from that mixture's two reference individuals; and one table
pooled across all seven reference individuals. (A) Presence-test detection rate
against nominal recipient fraction on the semi-synthetic ladder. Detection at
0.2% recipient rises from
{{ error_table_arms_headline.flat_detect_0p2 }} with the flat default to
{{ error_table_arms_headline.per_mixture_detect_0p2 }} with per-mixture tables
and {{ error_table_arms_headline.pooled_detect_0p2 }} when pooled, and full
detection is reached at 0.3% rather than 0.4%. The real titration moves the same
way: the median estimate at a known 0.5% recipient goes
{{ error_table_arms_headline.flat_real_0p5_mle_med_pct }}% to
{{ error_table_arms_headline.per_mixture_real_0p5_mle_med_pct }}% to
{{ error_table_arms_headline.pooled_real_0p5_mle_med_pct }}%. (B) The
corresponding false-signal floor, measured on the
{{ error_table_arms_headline.pooled_zero_n | dp(0) }} pure-donor endpoints, which are
genuine 0%-recipient samples carrying real reads and the real co-pooled
contamination. Bars are the per-arm maximum, dots the median; the grey band is
the contamination floor measured independently from the same data
({{ error_table_arms_headline.contam_floor_median_pct }}% to
{{ error_table_arms_headline.contam_floor_max_pct }}%).

A data-derived error table is worth building: the flat default detects nothing
at a 0.2% recipient fraction, where a pooled table detects two-thirds of
replicates. The gain does not cost specificity. No arm produced a false positive
on the pure-donor endpoints, and the maximum false signal stayed below the
independently measured contamination floor in all three
({{ error_table_arms_headline.flat_zero_mle_max_pct }}%,
{{ error_table_arms_headline.per_mixture_zero_mle_max_pct }}% and
{{ error_table_arms_headline.pooled_zero_mle_max_pct }}%).

No single arm dominates, and we report the comparison rather than a
recommendation. Pooling gives the best low-fraction sensitivity, but the
per-mixture tables give both the lowest false-signal floor and the closest
agreement on the semi-synthetic ladder (mean absolute deviation
{{ error_table_arms_headline.per_mixture_ladder_mean_abs_dev_pct }} against
{{ error_table_arms_headline.pooled_ladder_mean_abs_dev_pct }} percentage points
pooled and {{ error_table_arms_headline.flat_ladder_mean_abs_dev_pct }} flat).
The floor is also not ordered by how much data the table sees: pooling sits
close to the flat default, and only the per-mixture tables sit materially lower.
A laboratory monitoring near its limit of detection should therefore choose
deliberately, pooling when sensitivity governs and per-patient tables when a low
false-signal floor governs. Pooling is sound only across samples that share a
panel, run and chemistry; the substitution background is what pools, while
run-specific index hopping and cross-sample contamination are handled separately
(Methods).
