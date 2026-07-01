## Results

### Sensitivity compared with commercial kits

allomix is the only open-source chimerism tool, and to our knowledge the only tool of any
kind (open or commercial), that works with arbitrary marker panels from standard VCF
files (Table 1). To place its sensitivity alongside the limit-of-detection
figures cited by commercial vendors, we characterised the LoD across panel size,
sequencing depth, and donor-host relatedness under the EP17-A2 framework, for both
readouts allomix runs on the same sample (Methods). The LoD is the lowest fraction
recovered in at least 95% of replicates. The two readouts measure it in opposite
directions: the magnitude estimate quantifies the minor fraction (here the donor,
titrated against a pure-host background), and the residual-host presence test detects
residual host (the recipient titrated against a pure-donor background), the direction
that matters clinically when a patient is near full donor chimerism.

At an unrelated donor, 100 markers, and 1,000x mean depth, the in silico
magnitude-estimate LoD was {{ lod_headline.unrelated_lod_1000x_100markers_pct }}% and
the residual-host detection LoD was
{{ presence_lod_curve_headline.presence_unrelated_lod_1000x_100markers_pct }}%; with
full-sibling donors at the same panel and depth the magnitude LoD rose to
{{ lod_headline.sibling_lod_1000x_100markers_pct }}% (residual host
{{ presence_lod_curve_headline.presence_sibling_lod_1000x_100markers_pct }}%),
reflecting the smaller number of informative markers. Both readouts improve
monotonically with panel size and depth across the full sweep (Figure 1), letting a
laboratory read off the expected in silico LoD for its own assay.

These figures sit in the range reported for commercial NGS chimerism kits (0.06--0.5%),
but they are not a head-to-head comparison: the allomix numbers are best-case analytical
figures from the model's information on near-binomial simulated data, and a real assay's
LoD can only be higher, whereas the vendor numbers come from dilution series on real
DNA. The honest limiter on real panels is overdispersion, not depth (Discussion), and
the more defensible sensitivity number is the real-data LoD measured by subsampling the
SRP434573 mixtures (Figure 5, below) rather than this analytical ceiling.

| Tool | Markers | LoD | Open Source | Panel Agnostic | Input |
|:---|:---:|:---:|:---:|:---:|:---:|
| AlloSeq HCT | {{ tool_landscape.alloseq_n_markers }} SNPs | {{ tool_landscape.alloseq_lod }}%* | No | No | Proprietary |
| Devyser Chimerism | {{ tool_landscape.devyser_n_markers }} indels | {{ tool_landscape.devyser_lod }}% | No | No | Proprietary |
| NGStrack | 34 indels | 0.1% | No | No | Proprietary |
| ScisGo Chimerism MD | >200 SNPs + indels | 0.2% (single) / 0.5% (multi) | No | No | Proprietary |
| **allomix** | **Any biallelic** | **{{ presence_lod_curve_headline.presence_unrelated_lod_1000x_100markers_pct }}% (in silico, 100 markers, 1,000x)** | **Yes (MIT)** | **Yes** | **VCF** |

**Table 1.** NGS-based chimerism monitoring tools. LoD = limit of detection. Commercial
specifications are from published evaluations.[@Blouin2024comparison;
@Pedini2021devyser; @Kakodkar2023alloseq; @Qama2026devyser] *The AlloSeq HCT figure is
the vendor-stated LoD; the analytical sensitivity reported in independent evaluation
varies.[@Kakodkar2023alloseq] The allomix value shown is the residual-host detection LoD
(presence test) at the reference operating point (unrelated donor, 100 markers, 1,000x
mean depth), the readout comparable to the vendor detection limits; the
magnitude-estimate LoD at the same point is
{{ lod_headline.unrelated_lod_1000x_100markers_pct }}%. See Figure 1 for the full sweep
of both readouts. All LoD values refer to bulk extracted DNA; clinical sensitivity in a
given specimen further depends on the proportion of the lineage of interest and the
upstream cell-sorting workflow.

![**Figure 1.** Limit of detection as a function of panel size and sequencing depth, for both readouts allomix runs on the same sample. Columns: the MLE magnitude estimate (left) and the host-presence detection test (right). Rows: unrelated donor-host pairs (top) and full-sibling pairs (bottom). Each coloured curve is the median LoD across donor/host pairs at the indicated depth (100x to 2,000x); shaded bands are the 10th-90th percentile across pairs (the identity-by-descent spread, wide for siblings at small panels and narrowing as markers are added). The MLE panels also show the limit of blank (LoB) as a faint dashed line; the presence test has no blank because its null is the sequencing-error background (detection = host-presence likelihood-ratio test at p < 0.05). Dashed horizontal lines mark 0.5% and 1% minor (donor or residual-host) fraction. Both readouts come from the same simulated sweep (matched depths, panel sizes, donor/host pairs, and error model), so the columns are directly comparable. LoD was estimated per pair under the CLSI EP17-A2 workflow ({{ lod_headline.n_pairs_unrelated }} unrelated and {{ lod_headline.n_pairs_sibling }} sibling pairs, {{ lod_headline.n_seq_reps }} sequencing replicates each).]({{ facts_dir }}/fig_lod_curves.png)

### Longitudinal monitoring and accuracy

A primary clinical use of chimerism testing is serial monitoring. We simulated a
post-HSCT engraftment trajectory across {{ timeline.n_timepoints | dp(0) }} timepoints
from day +14 to day +365 (Figure 2). The scenario runs from early engraftment (15% donor
at day +14) through progressive engraftment to 95% donor by day +100, a clinically
relevant dip to 92% donor at day +180, and recovery to 97% donor at day +365. allomix
tracked the trajectory with mean absolute error {{ timeline.mae_pct | dp(2) }} ±
{{ timeline.mae_sd_pct | dp(2) }}% across all timepoints (maximum
{{ timeline.max_error_pct | dp(2) }} ± {{ timeline.max_error_sd_pct | dp(2) }}%), and
the day +180 dip was recovered in all replicates with
{{ timeline.dip_abs_error_pct | dp(2) }}% absolute error.

The accuracy behind that trajectory holds across sequencing depths. We generated
synthetic chimeric VCFs spanning 0% to 100% donor using 100 markers with empirically
calibrated per-marker bias, non-uniform depth (CV = {{ sim_calibration.depth_cv }}), and
{{ sim_calibration.locus_dropout_pct }}% locus dropout, at five depths from 50x to
1,000x ({{ depth_50.n_replicates | dp(0) }} replicates each). Mean absolute error stayed
below 1% at every depth, improving from {{ depth_50.mean_abs_error_pct | dp(2) }} ±
{{ depth_50.mean_abs_error_sd_pct | dp(2) }}% at 50x to
{{ depth_1000.mean_abs_error_pct | dp(2) }} ±
{{ depth_1000.mean_abs_error_sd_pct | dp(2) }}% at 1,000x (Table 2, Figure 3). The 95%
profile-likelihood intervals showed coverage of
{{ depth_1000.ci_coverage_pct }}--{{ depth_200.ci_coverage_pct }}%, close to the nominal
95% across the depth range, with the per-marker-type overdispersion model (separate
donor-homozygous and donor-heterozygous concentration) absorbing the extra-binomial
scatter that a single shared term left in the residuals. Per-marker bias correction is
available but its effect is modest by design: in a 2,000x experiment it left the point
estimate essentially unchanged (0% and 100% donor samples estimated at
{{ bias_with_bias.est_0pct }}% and {{ bias_with_bias.est_100pct }}%) while narrowing the
mean interval from {{ bias_no_bias.mean_ci_width_pct | dp(2) }}% to
{{ bias_with_bias.mean_ci_width_pct | dp(2) }}%, because in logit space the correction
is a small proportional adjustment at the extreme allele frequencies that dominate
low-fraction samples (Supplementary Figure S4).

![**Figure 2.** Simulated post-HSCT engraftment monitoring. True donor fraction trajectory (grey squares, dashed) versus allomix estimates (blue circles, solid) across six timepoints. Thin blue lines show individual replicates (N={{ timeline.n_replicates | dp(0) }}); the thick blue line is the mean estimate, and the shaded band is the mean 95% profile-likelihood CI. The trajectory includes a clinically relevant dip at day +180. Simulated with 100 markers, {{ timeline.depth | commas }}x depth, empirically calibrated per-marker bias, and {{ sim_calibration.locus_dropout_pct }}% locus dropout.]({{ facts_dir }}/fig_timeline.png)

| Depth | MAE (%) | RMSE (%) | Max Error (%) | CI Coverage (%) | Mean CI Width (%) |
|:---:|:---:|:---:|:---:|:---:|:---:|
| 50x | {{ depth_50.mean_abs_error_pct | dp(2) }} ± {{ depth_50.mean_abs_error_sd_pct | dp(2) }} | {{ depth_50.rmse_pct | dp(2) }} ± {{ depth_50.rmse_sd_pct | dp(2) }} | {{ depth_50.max_abs_error_pct | dp(2) }} ± {{ depth_50.max_abs_error_sd_pct | dp(2) }} | {{ depth_50.ci_coverage_pct }} ± {{ depth_50.ci_coverage_sd_pct }} | {{ depth_50.mean_ci_width_pct | dp(2) }} ± {{ depth_50.mean_ci_width_sd_pct | dp(2) }} |
| 100x | {{ depth_100.mean_abs_error_pct | dp(2) }} ± {{ depth_100.mean_abs_error_sd_pct | dp(2) }} | {{ depth_100.rmse_pct | dp(2) }} ± {{ depth_100.rmse_sd_pct | dp(2) }} | {{ depth_100.max_abs_error_pct | dp(2) }} ± {{ depth_100.max_abs_error_sd_pct | dp(2) }} | {{ depth_100.ci_coverage_pct }} ± {{ depth_100.ci_coverage_sd_pct }} | {{ depth_100.mean_ci_width_pct | dp(2) }} ± {{ depth_100.mean_ci_width_sd_pct | dp(2) }} |
| 200x | {{ depth_200.mean_abs_error_pct | dp(2) }} ± {{ depth_200.mean_abs_error_sd_pct | dp(2) }} | {{ depth_200.rmse_pct | dp(2) }} ± {{ depth_200.rmse_sd_pct | dp(2) }} | {{ depth_200.max_abs_error_pct | dp(2) }} ± {{ depth_200.max_abs_error_sd_pct | dp(2) }} | {{ depth_200.ci_coverage_pct }} ± {{ depth_200.ci_coverage_sd_pct }} | {{ depth_200.mean_ci_width_pct | dp(2) }} ± {{ depth_200.mean_ci_width_sd_pct | dp(2) }} |
| 500x | {{ depth_500.mean_abs_error_pct | dp(2) }} ± {{ depth_500.mean_abs_error_sd_pct | dp(2) }} | {{ depth_500.rmse_pct | dp(2) }} ± {{ depth_500.rmse_sd_pct | dp(2) }} | {{ depth_500.max_abs_error_pct | dp(2) }} ± {{ depth_500.max_abs_error_sd_pct | dp(2) }} | {{ depth_500.ci_coverage_pct }} ± {{ depth_500.ci_coverage_sd_pct }} | {{ depth_500.mean_ci_width_pct | dp(2) }} ± {{ depth_500.mean_ci_width_sd_pct | dp(2) }} |
| 1,000x | {{ depth_1000.mean_abs_error_pct | dp(2) }} ± {{ depth_1000.mean_abs_error_sd_pct | dp(2) }} | {{ depth_1000.rmse_pct | dp(2) }} ± {{ depth_1000.rmse_sd_pct | dp(2) }} | {{ depth_1000.max_abs_error_pct | dp(2) }} ± {{ depth_1000.max_abs_error_sd_pct | dp(2) }} | {{ depth_1000.ci_coverage_pct }} ± {{ depth_1000.ci_coverage_sd_pct }} | {{ depth_1000.mean_ci_width_pct | dp(2) }} ± {{ depth_1000.mean_ci_width_sd_pct | dp(2) }} |

**Table 2.** allomix accuracy and confidence-interval performance across sequencing
depths (mean ± SD, N={{ depth_50.n_replicates | dp(0) }} replicates). MAE = mean
absolute error; RMSE = root mean square error. Error metrics are computed on interior
fractions (excluding 0% and 100%). Absolute-error boxplots and the depth-summary panels
are in Supplementary Figures S10 and S11.

![**Figure 3.** In silico accuracy across sequencing depths. Each panel shows true donor fraction (x-axis) versus estimated donor fraction (y-axis) at the indicated depth, all replicates (N={{ depth_50.n_replicates | dp(0) }}). Dashed line indicates perfect agreement.]({{ facts_dir }}/fig_depth_scatter.png)

### Validation on real titrated mixtures

To test allomix on real reads rather than simulation alone, we ran it on a public
dataset of titrated DNA mixtures (SRA study SRP434573[@Chu2024mipsnp]), produced on a
panel and platform independent of those used to calibrate the simulator. The data are
{{ srp434573.n_individuals | dp(0) }} unrelated individuals on a
{{ srp434573.panel_n_snps | commas }}-SNP MIP sample-identification panel
({{ srp434573.platform }}, single-strand merged reads, raw on-panel depth
~{{ srp434573.raw_depth_min | commas }}-{{ srp434573.raw_depth_max | commas }}x;
{{ srp434573.n_intervals | commas }} amplicon intervals reconstructed from the reads),
combined into {{ srp434573.n_mixtures | dp(0) }} two-person mixtures titrated from
{{ srp434573.dilution_max_pct }}% down to {{ srp434573.dilution_min_pct }}% minor
contributor, plus one three-person mixture. We assigned the minor contributor to the
host role, so each dilution series reads as a declining-chimerism trajectory, the same
direction as relapse monitoring. Genotypes and admixture allele depths were produced
with the two-phase pipeline (Methods), and allomix was run with default parameters, with
the per-marker co-pooled contamination correction (Methods) additionally enabled for the
two-person series. That correction is off by default; it is turned on here because this
flowcell carries a significant contamination dose-response (the condition that gates the
correction on, characterised below), and it is a no-op on a run without one. After the
default genotype-quality and depth filters, a median of
{{ srp434573.markers_used_median | dp(0) }} informative markers per sample entered the
estimate (range {{ srp434573.markers_used_min | dp(0) }} to
{{ srp434573.markers_used_max | dp(0) }}). The dataset exercises the unrelated-donor case
only (no related or sibling donors).

Across the {{ srp434573.n_timepoints | dp(0) }} two-person admixtures, allomix
recovered the known host fraction with {{ srp434573.mae_reliable_pct }}% mean absolute
error over
the range at or above 2.5% (R^2 = {{ srp434573.r2_reliable }} for estimated versus
known), the range where the dataset's contamination floor is small relative to the true
fraction (Figure 4A). The maximum-likelihood host fraction tracked the dilution down to
1% host, with all {{ srp434573.n_onepct | dp(0) }} of the 1% admixtures estimated
between {{ srp434573.mle_onepct_min_pct }}% and {{ srp434573.mle_onepct_max_pct }}%. At
the lowest titration (0.5% host, {{ srp434573.n_lowest | dp(0) }} mixtures) the estimate
scattered between {{ srp434573.mle_lowest_min_pct }}% and
{{ srp434573.mle_lowest_max_pct }}% (mean {{ srp434573.mle_lowest_mean_pct }}%), where
the contamination floor competes directly with the true signal. For the three-person
mixture (1:3:5 of F2:M1:M2), allomix resolved all three contributors close to their
known fractions and ordered the two donors correctly: host F2 at
{{ srp434573.three_f2_est_pct }}% (known {{ srp434573.three_f2_known_pct }}%), donor M1
at {{ srp434573.three_m1_est_pct }}% (known {{ srp434573.three_m1_known_pct }}%), and
donor M2 at {{ srp434573.three_m2_est_pct }}% (known
{{ srp434573.three_m2_known_pct }}%) (Figure 4B).

This dataset carries a genuine low-level co-pooled contamination floor, most plausibly
index hopping on the patterned flowcell, which the source study independently identified
in these data as index misassignment.[@Chu2024mipsnp] allomix establishes it from the
data by a
dose-response argument rather than assuming it: at consensus-homozygous sites (host and
donor both homozygous for the same allele, so the minor allele cannot come from either
contributor), the median per-site minor-allele fraction rises with the number of other
co-pooled panel individuals carrying that allele, from
{{ srp_contam.nocarrier_floor_pct | dp(3) }}% at the
{{ srp_contam.n_nocarrier_sites | dp(0) }} sites with no co-pooled carrier (the
sequencing-error floor) to {{ srp_contam.dose_1carrier_pct | dp(2) }}% with one
co-pooled carrier and {{ srp_contam.dose_5carrier_pct | dp(2) }}% with five, across
{{ srp_contam.n_carrier_sites | dp(0) }} carrier sites (Supplementary Figure S13). That
monotonic rise is the signature of real reads from co-pooled material, not sequencing
error, which would not scale with co-pooled dose. This floor (a median
{{ srp_contam.carrier_median_pct | dp(2) }}% at carrier sites) sits on top of the sub-1%
host fractions, which is why the 0.5% estimates scatter. Given that significant
dose-response, the headline two-person estimates apply the per-marker contamination
correction (Methods): it subtracts the dose-predicted contamination from each
donor-homozygous host-allele count before the fit, calibrating the per-carrier rate on
each run's own consensus-homozygous and informative markers. All
{{ srp434573.correction_n_gated }} of {{ srp434573.correction_n_mixtures }} two-person
mixtures gate in (carrier-dose slopes up to {{ srp434573.correction_slope_max_pct }}%
per carrier), with the clean pairs self-selecting to a zero slope where the correction
is a no-op. It pulls the pure-donor (true-0%-host) endpoints down from a maximum of
{{ srp434573.endpoint_floor_max_baseline_pct }}% to
{{ srp434573.endpoint_floor_max_corrected_pct }}% (median
{{ srp434573.endpoint_floor_median_baseline_pct }}% to
{{ srp434573.endpoint_floor_median_corrected_pct }}%) and slightly tightens the
at-or-above-2.5% error ({{ srp434573.mae_reliable_baseline_pct }}% to
{{ srp434573.mae_reliable_pct }}%), while the 0.5% mean is essentially unchanged
({{ srp434573.mle_lowest_mean_baseline_pct }}% to {{ srp434573.mle_lowest_mean_pct }}%)
because the flat part of the floor is left to the per-site error model. The correction
only ever removes host-allele reads (its slope is clamped non-negative), so it is an
asymmetric downward adjustment that can under-report a real sub-0.1% host signal; that
is why it is gated on the dose-response and off by default. Supplementary Figure S12
draws the independent per-mixture contamination level (the consensus-homozygous floor,
{{ srp434573.contam_line_min_pct }}--{{ srp434573.contam_line_max_pct }}% across
mixtures), the level below which a corrected estimate is not separable from
contamination. Separately, {{ srp434573.n_review | dp(0) }} of
{{ srp434573.n_timepoints | dp(0) }} admixtures carry a REVIEW flag, predominantly from
the goodness-of-fit check rather than a biased point estimate: per-marker variance
exceeds the model expectation (overdispersion), so the fit is flagged even where the
fraction is recovered accurately, the same overdispersion gap quantified in the
limit-of-detection analysis (Supplementary Figures S7, S8). We inferred which named
individual is the titrated minor from the dataset's naming structure; the thesis does
not state it explicitly, and this affects only the host/donor labelling, not the
genotyping.

![**Figure 4.** allomix on the SRP434573 public titrated-mixture dataset (real reads). (A) Two-person dilution series: known host fraction versus allomix estimate (log-log) for the maximum-likelihood estimate (filled circles, 100 minus donor%, with the per-marker co-pooled contamination correction applied; Methods) and the residual-host presence test (open squares), across {{ srp434573.n_mixtures | dp(0) }} two-person mixtures ({{ srp434573.n_timepoints | dp(0) }} admixtures in total). Dashed line is perfect recovery. The 0.5% points still fall away from the line where the residual contamination floor competes with the true fraction. An alternative view of the same dilution series, with 95% confidence intervals on each estimate, is in Supplementary Figure S12. (B) The single three-person mixture (1:3:5 of F2:M1:M2): known versus estimated component fractions with 95% confidence intervals.]({{ facts_dir }}/fig_srp434573.png)

The real series floors at 0.5% host, where the co-pooled contamination competes with the
signal. To probe below that without the contamination confound, we built a
semi-synthetic sub-0.5% ladder from the same individuals: pure reference BAMs
sub-sampled and remixed at known host fractions from
{{ srp434573_synthetic.frac_min_pct }}% to {{ srp434573_synthetic.frac_max_pct }}%
({{ srp434573_synthetic.n_pairs }} donor-host pairs, {{ srp434573_synthetic.n_seeds }}
replicates each), depth-normalised on on-target reads so the realised minor fraction
matches the nominal one (real reads and genotyping path, artificial mixing ratio). On
this ladder the maximum-likelihood host estimate tracked the true fraction down to 0.1%
host, with the median estimate rising in step with the nominal one (median
{{ srp434573_synthetic.mle_med_0p1 }}% at 0.1% host,
{{ srp434573_synthetic.mle_med_0p5 }}% at 0.5%, near-unit slope) and no constant floor.
The residual-host presence test gives an independent low-fraction readout: it climbed
from a detection rate of {{ srp434573_synthetic.detect_rate_0p1 }} at 0.1% host to full
detection by 0.4% ({{ srp434573_synthetic.detect_rate_0p4 }}), reading only the
donor-homozygous markers. The same construction produced a host-plus-two-donor series
(host F2 titrated below 1% against donors M1 and M2), where allomix recovered the sub-1%
host and both donors with median absolute errors of
{{ srp434573_synthetic.three_host_med_abs_err_pct }} and
{{ srp434573_synthetic.three_donor_med_abs_err_pct }} percentage points. Reaching truth
at the low end depends on the overdispersion model, not on the data. On this dataset a
single overdispersion parameter shared across marker types leaves a near-constant positive
host offset of about 0.22 percentage points at these fractions, coming from the
donor-heterozygous markers, where the two contributors balance the alleles near 0.5 and
their symmetric extra-binomial scatter rectifies into a small positive host signal.
Fitting a separate overdispersion for the donor-heterozygous and donor-homozygous classes
(the default per-marker-type model) absorbs that scatter where it occurs and removes the
offset. The mixing itself is faithful (the realised minor fraction tracks the nominal one
with unit slope, and the donor-homozygous markers alone return the nominal fraction), and
the offset is neither a background-artifact floor nor an estimator boundary effect: a
plain binomial fit shows no offset, confirming it is a property of the shared-overdispersion
fit rather than of the reads.

### Residual-host presence test

The magnitude estimate above answers how much donor is present. The second test answers
whether any host remains, and reports a host-fraction estimate for it, reading only the
donor-homozygous markers where the host carries the donor-absent allele (Methods). On the
SRP434573 dilution series it returned a positive call (p < 0.05) at all
{{ srp434573.presence_n_detected | dp(0) }} of
{{ srp434573.presence_n_total | dp(0) }} two-person admixtures, reading a median of
{{ srp434573.presence_markers_median | dp(0) }} donor-homozygous markers per sample. The
same dataset supplies the specificity at true zero: each pair's pure-donor sample, the
donor's own DNA piled through the identical mpileup path, is a genuine 0%-host input, and
the test correctly withheld a host-present call at all
{{ srp434573.zero_host_presence_absent_n | dp(0) }} of {{ srp434573.zero_host_n | dp(0) }}
pure-donor controls (host-fraction estimate 0%), even with the real co-pooled
contamination floor present in those reads. The magnitude estimator floored there too,
with a maximum host estimate of {{ srp434573.zero_host_mle_max_pct }}% across those
{{ srp434573.zero_host_n | dp(0) }} controls. At
the 1% host level its host-fraction estimate ranged
{{ srp434573.presence_onepct_min_pct | dp(2) }} to
{{ srp434573.presence_onepct_max_pct | dp(2) }}%, tracking the dilution alongside the
magnitude estimate (Figure 4A, open squares); at the lowest 0.5% titration it cannot
separate residual host from the co-pooled contamination floor (below), so a positive
call there reflects both. We present this as a validated capability, not as a clinical
relapse-detection result: the test is calibrated against the background-artifact floor
and demonstrated in silico and on this titrated panel down to 1% host, but its operating
characteristics on real patient samples, and the clinical thresholds that would turn a
positive presence call into an action, remain to be established (Discussion). Its value
is that it is a separate readout from the magnitude estimate, designed to surface
residual host below the level the magnitude estimate can quantify, which is exactly the
regime where early relapse would first appear.

### Real-data limit of detection

The simulated LoD curves of Figure 1 have a real-data counterpart. We sub-sampled reads
and markers from the high-depth SRP434573 mixtures, reducing depth and panel size until
the LoD rose into the measurable window, and characterised the LoD across panel size and
sequencing depth on real reads for both readouts
({{ subsample_lod_headline.n_mixtures | dp(0) }} two-person mixtures,
{{ subsample_lod_headline.n_seeds | dp(0) }} sub-sampling replicates per cell; Figure
5). As in simulation, the LoD fell with both panel size and depth. In the seven mixtures
titrated only to 1% the median LoD reached at or below 1% (the lowest dilution they
carry); in the three titrated to {{ subsample_lod_headline.min_titration_pct | dp(1) }}%
host, both the magnitude and presence-test LoD reached at or below that level. The dilution series
stops at {{ subsample_lod_headline.min_titration_pct | dp(1) }}% (its lowest real
titration), and the dataset's co-pooled contamination floor (a median
{{ subsample_lod_headline.contamination_floor_pct | dp(2) }}% at the deepest, largest
panel) sits beneath it, so sub-1% cells are upper bounds rather than resolved values.
These are pseudo-replicates sub-sampled from one library, not independent low-depth
libraries, so the result confirms that the real-data LoD tracks the simulation within
the limits of the dilution grid rather than serving as an independent wet-lab LoD.

![**Figure 5.** Real-data limit of detection on the SRP434573 titrated mixtures, the real-data counterpart of Figure 1. Reads and markers were sub-sampled from the high-depth mixtures to bring the LoD into the measurable window. Columns: MLE magnitude estimate (left) and host-presence detection test (right). Rows: the seven mixtures titrated only to 1% (top) and the three titrated to {{ subsample_lod_headline.min_titration_pct | dp(1) }}% host (bottom), two disjoint sets. Only three of the {{ subsample_lod_headline.n_mixtures | dp(0) }} mixtures were diluted below 1%; keeping the rows disjoint stops those three from being buried in a top-row median that the 1%-floored mixtures would otherwise pin at 1%. Each coloured curve is the median LoD across mixtures at the indicated depth (100x to 2,000x); shaded bands are the 10th-90th percentile across mixtures. An X marks a cell where the LoD is at or below the lowest titration that mixture set carries (1% top row, 0.5% bottom row), not resolved lower. Per-mixture curves are constrained to be monotonic in panel size, which the nested marker panels justify. Points are jittered horizontally per depth. Sub-sampled pseudo-replicates, {{ subsample_lod_headline.n_seeds | dp(0) }} per cell, not independent libraries.]({{ facts_dir }}/fig_subsample_lod_grid.png)

### Stress tests (in silico): relatedness, multiple donors, and recipient copy-number changes

**Donor-host relatedness.** Related donors share genotypes with the host, reducing the
number of informative markers. Across four relatedness levels
({{ rel_unrelated.n_replicates }} replicate pairs each, {{ rel_unrelated.n_markers }}
markers, 500x), the mean number of informative markers fell from
{{ rel_unrelated.mean_informative }} (unrelated) to {{ rel_sibling.mean_informative }}
(full sibling), yet mean absolute error stayed below 2% throughout:
{{ rel_unrelated.mean_mae_pct }}% (unrelated), {{ rel_cousin.mean_mae_pct }}% (cousin),
{{ rel_half_sibling.mean_mae_pct }}% (half-sibling), and {{ rel_sibling.mean_mae_pct }}%
(sibling) (Table 3, Figure 6). Even with sibling donors the minimum informative-marker
count observed ({{ rel_sibling.min_informative }}) stayed well above the three required.
Relatedness costs sensitivity (fewer informative markers raise the LoD) but not
point-estimate accuracy, because the markers that remain are still unbiased.

| Relatedness | Mean Informative | Range | MAE (%) | RMSE (%) |
|:---|:---:|:---:|:---:|:---:|
| Unrelated | {{ rel_unrelated.mean_informative }} | {{ rel_unrelated.min_informative }}–{{ rel_unrelated.max_informative }} | {{ rel_unrelated.mean_mae_pct }} | {{ rel_unrelated.mean_rmse_pct }} |
| 1st cousin | {{ rel_cousin.mean_informative }} | {{ rel_cousin.min_informative }}–{{ rel_cousin.max_informative }} | {{ rel_cousin.mean_mae_pct }} | {{ rel_cousin.mean_rmse_pct }} |
| Half-sibling | {{ rel_half_sibling.mean_informative }} | {{ rel_half_sibling.min_informative }}–{{ rel_half_sibling.max_informative }} | {{ rel_half_sibling.mean_mae_pct }} | {{ rel_half_sibling.mean_rmse_pct }} |
| Full sibling | {{ rel_sibling.mean_informative }} | {{ rel_sibling.min_informative }}–{{ rel_sibling.max_informative }} | {{ rel_sibling.mean_mae_pct }} | {{ rel_sibling.mean_rmse_pct }} |

**Table 3.** Effect of donor-host relatedness on marker informativity and chimerism
accuracy. Each level: {{ rel_unrelated.n_replicates }} replicate donor-host pairs,
{{ rel_unrelated.n_markers }} markers, 500x depth.

![**Figure 6.** Effect of donor-host relatedness on allomix performance. Left: informative markers by relatedness level (dots = replicates, bars = means). Centre: mean absolute error. Right: truth versus estimated donor fraction across all replicates. Simulated with {{ rel_unrelated.n_markers }} markers, 500x mean depth (CV = {{ sim_calibration.depth_cv }}), {{ sim_calibration.seq_error_pct }}% sequencing error, empirically calibrated per-marker bias, and {{ sim_calibration.locus_dropout_pct }}% locus dropout.]({{ facts_dir }}/fig_relatedness.png)

**Multiple donors.** For the hardest multi-donor case, we generated a three-sibling
scenario (host and two donors sharing both parents) across
{{ multidonor.n_markers | dp(0) }} markers at {{ multidonor.depth | commas }}x, of which
{{ multidonor.n_informative_any | dp(0) }} were informative for at least one donor.
Across {{ multidonor.n_samples | dp(0) }} chimeric samples spanning the simplex,
per-donor mean absolute error was {{ multidonor.mae_d1_pct }}% (donor 1) and
{{ multidonor.mae_d2_pct }}% (donor 2), the total donor fraction was estimated with
{{ multidonor.mae_total_pct }}% MAE, and all {{ multidonor.n_asymmetric | dp(0) }}
asymmetric mixes were correctly ranked (Table 4, Figure 7). Per-donor CI coverage was
lower ({{ multidonor.ci_coverage_d1_pct }}% and {{ multidonor.ci_coverage_d2_pct }}%)
because partitioning a mix between two related donors leaves few markers that separate
one donor from the other; the total and the ranking stay accurate. Unrelated multi-donor
settings (such as cord-blood transplants) would yield more informative markers and
tighter per-donor intervals.

| Metric | Donor 1 | Donor 2 | Total |
|:---|:---:|:---:|:---:|
| MAE (%) | {{ multidonor.mae_d1_pct }} | {{ multidonor.mae_d2_pct }} | {{ multidonor.mae_total_pct }} |
| RMSE (%) | {{ multidonor.rmse_d1_pct }} | {{ multidonor.rmse_d2_pct }} | {{ multidonor.rmse_total_pct }} |
| Max error (%) | {{ multidonor.max_error_d1_pct }} | {{ multidonor.max_error_d2_pct }} | {{ multidonor.max_error_total_pct }} |
| CI coverage (%) | {{ multidonor.ci_coverage_d1_pct }} | {{ multidonor.ci_coverage_d2_pct }} | n/a |

**Table 4.** Multi-donor accuracy with sibling donors.
{{ multidonor.n_markers | dp(0) }} markers ({{ multidonor.n_informative_any | dp(0) }}
informative), {{ multidonor.depth | commas }}x depth. Error metrics on interior
fractions.

![**Figure 7.** Multi-donor chimerism estimation. (A) Per-donor accuracy: true versus estimated fraction for donor 1 (circles) and donor 2 (triangles), with 95% profile-likelihood CIs. (B) Two-dimensional log-likelihood surface for a representative mixture (60% host, 30% donor 1, 10% donor 2); contours show delta log-likelihood from the maximum, dashed line marks the 95% joint CI, grey region is infeasible (f1 + f2 > 1). Star = true value, circle = MLE.]({{ facts_dir }}/fig_multidonor.png)

**Recipient copy-number changes.** The relapsing recipient clone often carries
copy-number changes, which the diploid mixture model does not represent. We simulated
recipient copy-neutral loss of heterozygosity, deletion, and gain at controllable burden
and measured the LoD in two directions ({{ cnv_loh_headline.ref_markers }} markers,
{{ cnv_loh_headline.ref_depth }}x, {{ cnv_loh_headline.n_reps_per_frac }} replicates per
cell; Figure 8). For detecting the recipient clone as the minor component (the relapse
direction), aberrations had little effect: the LoD was
{{ cnv_loh_headline.relapse_lod_baseline_unrel_pct | dp(2) }}% (unrelated) at baseline
and stayed at or below {{ cnv_loh_headline.relapse_lod_max_pct | dp(2) }}% across every
aberration type and burden, because the aberration rides with the signal being detected.
The opposite direction, detecting a donor against an aberration-bearing recipient
background, degraded sharply: even a 10% loss-of-heterozygosity or deletion burden
pushed the unrelated donor LoD above the 20% probed ceiling. The one-sided
outlier-resistant refit recovered much of this at low to moderate burden (for example
returning the high-burden gain donor LoD from
{{ cnv_loh_headline.donor_lod_gain_high_unrel_std_pct | dp(2) }}% to
{{ cnv_loh_headline.donor_lod_gain_high_unrel_robust_pct | dp(2) }}%) while leaving
aberration-free samples unchanged, and the estimator flags the cases it cannot recover
for review (Methods).

![**Figure 8.** Effect of recipient copy-number aberrations on the limit of detection, by direction. Columns are aberration type (CN-LoH, deletion, gain); x-axis is burden (fraction of eligible recipient markers affected). Y-axis is the LoD of the minor component as a donor percentage, log-scaled. Solid lines: standard estimator; dashed lines: one-sided outlier-resistant refit; colours: relatedness. Top row: relapse detection (recipient clone is the minor component, insensitive to the aberration). Bottom row: donor detection against an aberration-bearing recipient background, where CN-LoH and deletion inflate the donor LoD past the 20% ceiling (dotted line) and the refit partly recovers it. Simulated at {{ cnv_loh_headline.ref_markers }} markers, {{ cnv_loh_headline.ref_depth }}x, same noise model as Figure 1.]({{ facts_dir }}/fig_cnv_loh.png)

### Built-in quality and sample-integrity checks

Because allomix is meant to run inside routine laboratory operations, it reports a
three-level verdict (PASS, REVIEW, FAIL) and a set of sample-integrity checks built from
the same marker data (Methods). Three of these target low-fraction signals the magnitude
estimate cannot see, kept separate by genotype geometry rather than by re-thresholding
one number: residual host (at donor-homozygous markers where the host carries the
donor-absent allele), contamination by a non-host, non-donor genome, and a gross sample
swap (both at the consensus-homozygous markers the magnitude estimate never reads). On SRP434573 it was
that contamination geometry, the consensus-homozygous markers, that exposed the
co-pooled floor by dose-response (above), reported as the excess minor signal over a
data-internal error floor. Two further checks guard identity at the sample level: a
relatedness check that flags unexpected kinship as a swap signature, and a read-level
artifact filter that judges artifacts by effect size and auto-disables its strand test
on single-strand amplicon panels (as on SRP434573). Together these turn the marker data
already used for the fraction estimate into a defence against the wrong-sample and
contamination errors that a quantitative chimerism result would otherwise carry
silently.
