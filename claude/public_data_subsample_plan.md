# Plan: real-data LoD curves by sub-sampling SRP434573 (issue #24)

Status: design only, for review. No implementation yet.

## 1. The ask (issue #24)

> Figure 1 changes the markers/depth to show how the LOD curve changes.
> The public data (SRP434573) has >1000 markers and ~2000x depth - we could
> subsample that down as a way to produce our LOD curve on real data rather
> than a sim.

So: reproduce the panel-size / depth LoD sweep that Figure 1 currently shows
**in silico**, but driven by the real SRP434573 reads, by throwing away markers
and reads from the high-depth real data and watching the LoD rise.

## 2. Why this is worth doing (where it fits the paper)

The current LoD figure (`fig5_lod_curves.png`, called "Figure 1" in
`results.md`) is entirely simulated. The results text already concedes the
weakness:

- "the allomix number is a best-case analytical figure from the model's
  information on **near-binomial simulated data**, and a real assay's LoD can
  only be higher"
- "The honest limiter on real panels is **overdispersion, not depth**."

Issue #24 is the experiment that lets us put those two sentences on real data.
Sub-sampling real reads keeps the things the simulator cannot fully reproduce:
the true per-marker capture bias, the real between-marker overdispersion, and
the known co-pooled contamination floor in this dataset (the thing that makes
the 0.5% points in Figure 4A fall off the line). The deliverable is a real-data
LoD-vs-(depth, panel size) curve we can lay next to the simulated one and say
"the sim is the floor; here is where real data lands."

This complements, not replaces, the existing SRP434573 work:
- `run_srp434573_allomix.py` -> Figure 4 (accuracy on the native full-depth data)
- `run_lod_validation.py` -> Figure 1/5 (simulated LoD sweep, CLSI EP17-A2)
- `run_presence_lod.py`, `run_overdispersion_lod.py` (simulated LoD variants)
- **new** `run_subsample_lod.py` -> real-data LoD sweep (this plan)

## 3. What the data gives us (and its limits)

Committed snapshot in `paper/public_data/SRP434573/genotypes/` (no `/tau`,
no FASTQ, no re-alignment needed):

- 11 mixtures (10 two-person + 1 three-person). Each has
  `<mix>.SRP434573.vcf.gz` (host+donor `GT`) and `<mix>.admix.vcf.gz`
  (per-timepoint `AD`/`DP` at panel sites), plus optional
  `<mix>.error_table.tsv`.
- Admix `FORMAT` carries `GT:PL:DP:AD` and INFO carries `DP4`/`*BZ` bias fields.
- Per-marker depth ~600-2100x (raw on-panel ~1000-1900x).
- ~1052 panel intervals; median informative markers per sample is reported in
  the paper facts (`srp434573.markers_used_median`).
- Titration fractions (minor = host): **10, 5, 2.5, 1.25, 1, 0.5 %**. Some
  mixtures have v1/v2 repeats per fraction.

Limits to design around:
- **Only 6 fractions**, geometric, bottoming at 0.5%. At full depth+markers the
  real LoD may already be **below 0.5%** (off the bottom of the titration), so
  there we can only report "< 0.5%". Sub-sampling is exactly the lever that
  pushes the LoD **up into** the 0.5-10% window we can actually measure. Good.
- **Few native replicates** per (mixture, fraction). Solved by the resampling
  design in section 5.
- **Unrelated donors only.** This produces the real analog of Figure 1's *left*
  panel (unrelated). The sibling panel stays simulation-only; the dataset has no
  related pairs (issue #16 notes CEPH-1463 as a future related-donor analog).
- **Contamination floor.** This dataset has a known co-pooled floor near ~0.5%.
  It will *cap* the achievable real LoD: more depth/markers cannot drive the LoD
  below the floor. That is a feature to report, not a bug, and it is the cleanest
  real-data illustration of the Discussion's "overdispersion/contamination, not
  depth, is the limiter" point.

## 4. Sub-sampling mechanics

Two independent knobs, both applied to the parsed informative-marker list so we
never touch the VCFs on disk and never call out to the aligner.

### 4a. Depth (binomial read thinning of AD at a single global rate)

Use **one global keep-rate per sample**, not a per-marker normalisation. To hit
a target *mean* depth `D` for a sample whose observed mean admix depth is
`D_obs`, set

```python
rate = min(1.0, D / D_obs)            # one rate for the whole sample
new_ref = rng.binomial(ad_ref, rate)  # applied identically to every marker
new_alt = rng.binomial(ad_alt, rate)
new_dp  = new_ref + new_alt
```

This is the faithful analog of `samtools view -s` / `seqtk sample` on the
FASTQ: a single keep-probability applied uniformly to every read. Because the
rate is the same everywhere, a 2000x amplicon lands near `2000*rate` and a 600x
one near `600*rate`, so the **real locus-to-locus depth CV is preserved** (deep
amplicons stay deeper, shallow ones drop out first). The figure axis is then
"mean depth", exactly how Figure 1 is labelled.

An earlier draft of this plan normalised each marker to the target depth
(`p = D / dp_i`), which forces every marker to the same depth and **flattens the
real depth distribution**, removing the locus-dropout variance. That is
optimistic and less faithful to FASTQ subsampling, so it is dropped in favour of
the single global rate above.

Why this is a fair stand-in for "subsample FASTQ then re-pileup" on *this*
dataset (the equivalence is exact, not approximate, under these conditions):

- Subsampling reads keeps each read independently with probability `rate`, so
  surviving ref/alt counts are `Binom(ad_ref, rate)` / `Binom(ad_alt, rate)` and
  independent. That **is** binomial thinning, not an approximation of it.
- The conditions that make it exact hold here: single-end pear-merged reads
  (read = fragment = independent unit), no UMI collapse and no dup-marking (so
  AD is a raw independent read count with no depth-nonlinear dedup), admix AD is
  a forced `bcftools mpileup` at fixed panel sites with no re-calling (lower
  depth changes counts only, never the site set or host/donor genotypes), and
  per-read-independent alignment (subsample-before vs -after align is the same).
- What thinning does *not* reproduce is independent-library overdispersion (PCR
  / capture jackpotting). But **subsampling the deposited FASTQ cannot reproduce
  that either** (same single library), so it is not a difference between the two
  methods, just a shared limit of any one-library subsample. The between-marker
  bias and overdispersion baked into the deposited reads are inherited
  identically by both.
- Residual second-order differences (mpileup BAQ/local-realignment depth
  dependence; indel/multiallelic context) are negligible for a biallelic-SNP
  panel.

**Re-apply the admix depth filter after thinning.** `classify_markers` applies
`min_dp` at parse time on the full-depth data; thinning happens after
classification, so the script must re-drop markers whose thinned `admix_dp`
falls below `min_dp`. This is the low-depth locus-dropout mechanism and it is
what makes the global-rate thinning behave like a real shallow pileup.

### 4b. Panel size (nested random marker subsets)

To get a monotone panel-size axis (as in `run_lod_validation.py`), draw one
random permutation of the informative markers per (mixture, seed) and take
prefixes:

```python
order = rng.permutation(len(markers))
for n in N_MARKERS_GRID:               # e.g. 25, 50, 75, 100, 200, 400
    if n > len(markers):
        continue                        # cap at this mixture's informative count
    subset = [markers[i] for i in order[:n]]
```

Nesting (prefix) guarantees that adding markers can only help, so each
(mixture, seed) curve is monotone and the across-mixture median is too. The max
panel size is capped per mixture by its informative-marker count (varies by
mixture); the grid auto-truncates.

### 4c. Replicate structure (how thin data becomes a sweep)

This is the key trick that makes 6 fractions usable. From one real admix sample
we generate many pseudo-replicates:

- `S_thin` thinning seeds (sequencing-noise replicates at fixed depth)
- `S_panel` permutation seeds (panel draws at fixed size)
- across the (up to) 10 two-person mixtures and their native v1/v2 repeats

These are **pseudo-replicates** (one real read draw underlies each thinning
family), so we will say so plainly in the figure caption and Methods. They are
valid for characterising the *sub-sampling response* (how LoD moves with
depth/markers); they are not 30 independent libraries. We mirror the
two-level design of `run_lod_validation.py`: treat each mixture as a "pair",
compute an LoD per mixture, report the **median across mixtures** as the curve
and the **10th-90th percentile across mixtures** as the band.

## 5. Two separate, complementary LoD analyses

allomix reports two independent low-fraction readouts, and they are not a
primary/secondary pair: they measure different things and each gets its own
sweep, facts, and figure (this mirrors Figure 4A, which already shows them side
by side, MLE as filled circles and presence as open squares):

1. **Magnitude (MLE) LoD** -- the beta-binomial donor/host-fraction estimate
   (`estimate_single_donor_bb`). Answers "how much minor is present".
2. **Presence-test LoD** -- the donor-homozygous residual-host detector
   (`host_presence_test`). Answers "is any minor present at all".

Both share the same per-cell pipeline (parse -> subset informative markers ->
global-rate thin -> re-apply `min_dp`) and run in the **same pass** over each
cell, so the only added cost over a single-test sweep is the second estimator
call. Each uses a blank-free, per-sample detection rule, so neither needs the
manufactured blanks discussed in section 9:

- **Magnitude (MLE) detection.** A cell is detected when the 95% profile-
  likelihood CI for the host (minor) fraction excludes 0 (equivalently the
  donor-fraction CI excludes 1): the estimator itself says it can separate the
  sample from "pure donor" at 95% confidence. (`ChimerismResult` also exposes a
  native per-sample analytical `lob_fraction` / `lod_fraction` from Fisher
  information; see the bonus curve below.)
- **Presence detection.** A cell is detected when
  `host_presence_test(...).lrt_pval < 0.05`, whose null is the sequencing-error
  background (optionally the per-site error table).

For **each** test independently, LoD = lowest titration fraction at which
`>= 95%` of replicate cells are detected, read from the same 2-parameter
logistic fit `P(detect | f) = sigmoid(a + b*log10 f)` over the 6 fractions that
`run_lod_validation.py` uses. Median across mixtures is the curve; 10/90 across
mixtures is the band. Two tests -> two curves, two figures.

**Bonus MLE curve (free).** Because `estimate_single_donor_bb` already returns a
per-sample analytical `lod_fraction` (the Fisher-information LoD given the
sample's markers and its *estimated* overdispersion `rho`), we can also plot the
median analytical `lod_fraction` vs depth/panel size with no titration fractions
involved. This is the information-theoretic LoD on the real markers' real bias
and real overdispersion, and it is conceptually the closest analog to how the
simulated Figure 1 LoD is defined. Worth emitting alongside the empirical MLE
curve as an internal consistency check.

**Caveat to watch on the MLE side.** This dataset has a co-pooled contamination
floor (~0.5%, the reason Figure 4A's 0.5% points fall off the line). The MLE
"CI excludes 0" rule will treat that floor as signal, so the empirical MLE LoD
may read artificially low at the bottom of the titration (it detects the
contamination, not the host). That is a real finding, not a bug: report it, and
let the analytical `lod_fraction` curve (which is floor-independent) sit next to
it to make the distinction visible. The presence test reads only donor-
homozygous markers and is affected differently, which is exactly why running the
two separately is informative.

Why not the EP17-A2 over-blanks route that Figure 1 uses internally: on real
data we have **no pristine zero-analyte sample** (titration bottoms at 0.5%, and
the pure single-individual runs are not in the committed snapshot). Manufacturing
blanks from cross-individual controls is biased (section 9, Q3). The per-sample
rules above are blank-free and apply uniformly to both tests.

## 6. Implementation

### 6a. New pure helper in `src/allomix/simulate.py` (unit-tested)

Thinning is small, reusable, and worth a test, so it goes in the library rather
than the script. It operates on `InformativeMarker` (from
`allomix.genotype`) and returns fresh copies:

```python
def thin_informative_markers(
    markers: list[InformativeMarker],
    rate: float,
    rng: np.random.Generator,
) -> list[InformativeMarker]:
    """Binomially down-sample admix AD by one global keep-rate (samtools -s).

    A single `rate` (0 < rate <= 1) is applied to every marker, so the real
    locus-to-locus depth CV is preserved. Preserves host/donor genotypes, marker
    type, and bias annotations; only the admix counts are resampled. rate >= 1
    is a no-op (cannot upsample).
    """
    if not 0.0 < rate <= 1.0:
        raise ValueError(f"rate must be in (0, 1], got {rate!r}")
    if rate == 1.0:
        return list(markers)
    out = []
    for m in markers:
        new_ref = int(rng.binomial(m.admix_ad_ref, rate))
        new_alt = int(rng.binomial(m.admix_ad_alt, rate))
        out.append(replace(  # dataclasses.replace
            m, admix_ad_ref=new_ref, admix_ad_alt=new_alt, admix_dp=new_ref + new_alt
        ))
    return out
```

The script turns a target mean depth into a rate per sample:

```python
def rate_for_mean_depth(markers_full, target_mean_depth):
    d_obs = float(np.mean([m.admix_dp for m in markers_full]))
    return min(1.0, target_mean_depth / d_obs)
```

Test (`tests/test_simulate.py`): thinned mean depth ~ rate * original; allele
ratio preserved in expectation over many seeds; the depth CV across markers is
preserved (not flattened); `rate=1.0` is a pass-through; `rng` makes it
deterministic.

### 6b. New script `paper/scripts/run_subsample_lod.py`

Structure mirrors `run_lod_validation.py` (argparse, `qval` quick-mode grid,
`ProcessPoolExecutor`, deterministic `derive_seed`, logistic LoD fit). The only
substantive change is that the per-cell marker list comes from **parsing +
thinning real data** instead of `blend_vcfs`.

Per mixture, parse once, then sweep:

```python
from allomix.genotype import parse_vcf, classify_markers, DEFAULT_MIN_DP
from allomix.chimerism import estimate_single_donor_bb
from allomix.detect import host_presence_test
from allomix.simulate import thin_informative_markers
from srp434573_common import resolve_srp434573_genotypes_dir

GEN = resolve_srp434573_genotypes_dir()

def load_mixture(name, host, donor):
    panel = GEN / f"{name}.SRP434573.vcf.gz"
    admix = GEN / f"{name}.admix.vcf.gz"
    host_md  = parse_vcf(panel, sample=host, min_gq=0, gt_ad_consistency=True)
    donor_md = parse_vcf(panel, sample=donor, min_gq=0, gt_ad_consistency=True)
    # one InformativeMarker list per admix timepoint (fraction)
    by_fraction = {}
    for s in VCF(str(admix)).samples:
        admix_md = parse_vcf(admix, sample=s, min_dp=0)
        mg = classify_markers(host_md, [donor_md], admix_md)
        by_fraction[known_host_pct(s)] = mg.informative   # 10,5,2.5,1.25,1,0.5
    return by_fraction

def cell_results(markers_full, frac, depth, n_markers, seed):
    """One subsampled cell -> both readouts (run in the same pass)."""
    rng = np.random.default_rng(derive_seed(frac, depth, n_markers, seed))
    rate = rate_for_mean_depth(markers_full, depth)   # one global rate per sample
    order = rng.permutation(len(markers_full))
    subset = [markers_full[i] for i in order[:n_markers]]
    thinned = thin_informative_markers(subset, rate, rng)
    thinned = [m for m in thinned if m.admix_dp >= DEFAULT_MIN_DP]  # low-depth dropout

    # Magnitude (MLE) readout: detection = host-fraction CI excludes 0.
    mle = estimate_single_donor_bb(thinned)
    host_ci_lo = 1.0 - mle.donor_fraction_ci[1]
    # Presence readout: detection = LRT p < 0.05.
    pres = host_presence_test(thinned)

    return {
        "n_used": len(thinned),
        "mle_host_frac": mle.host_fraction,
        "mle_detected": host_ci_lo > 0.0,
        "mle_analytical_lod": mle.lod_fraction,   # Fisher-info LoD (bonus curve)
        "presence_f": pres.f_host_mle,
        "presence_detected": pres.lrt_pval < 0.05,
    }
```

Marker subsetting is done on the **informative** marker set (the markers the
estimate actually uses), per Q5; the rate is computed from the full informative
set so it is a stable per-sample quantity, then applied to the subset. The
detection-rate -> logistic -> 95%-point reduction is then run **twice**, once on
`mle_detected` and once on `presence_detected`, producing two independent LoD
curves; `mle_analytical_lod` is aggregated directly (median across mixtures) for
the bonus information-theoretic MLE curve.

Then for each (mixture, depth, n_markers): compute detection rate per fraction
over `S_thin` x `S_panel` seeds, fit the logistic, read the per-mixture LoD;
aggregate to median + 10/90 band across mixtures.

Grids (quick-mode shrinks via `qval`, watermark applies):

```python
DEPTHS        = qval([100, 250, 500, 1000, 2000], [250, 1000])
N_MARKERS     = qval([25, 50, 75, 100, 200, 400], [50, 200])
FRACTIONS     = [10.0, 5.0, 2.5, 1.25, 1.0, 0.5]      # fixed by the titration
N_SEEDS       = qval(20, 3)                            # thin x panel seeds
```

**Quick-mode (required).** Like every other paper script, this must honour
`ALLOMIX_PAPER_QUICK` / the Snakefile `--config quick=1`. The grids above already
go through `qval(full, quick)`, which shrinks depths, panel sizes, and seed
counts so the rule drops from tens of minutes to ~a minute; both plot scripts
import `paper_quick` so their figures get the "QUICK BUILD" watermark and are not
used for publication. Add the new rules to the quick path the same way the
existing LoD rules are (the sweep runs both estimators, so its quick budget
should be sized like `run_lod_validation.py`'s, not double it).

Outputs (under `output/facts/`). A `test` column (`mle` | `presence`) keeps the
two analyses in the same files but cleanly separable:
- `subsample_lod_raw.csv` - one row per (mixture, fraction, depth, n_markers,
  seed): `mle_detected`, `presence_detected`, `mle_host_frac`, `presence_f`,
  `mle_analytical_lod`, `n_used`.
- `subsample_lod_per_mixture.csv` - per (test, mixture, depth, n_markers) LoD.
- `subsample_lod_summary.csv` - per (test, depth, n_markers): median + 10/90
  band across mixtures, mirroring `lod_summary.csv` columns so the existing plot
  helpers can be reused. The bonus analytical MLE LoD goes here too (e.g.
  `test=mle_analytical`).
- `subsample_lod_headline.csv` - named facts per test for later prose (e.g.
  `mle_lod_1000x_100markers_pct`, `presence_lod_1000x_100markers_pct`,
  `contamination_floor_pct`, `n_mixtures`, `n_seeds`).

### 6c. Figure `paper/scripts/plot_subsample_lod.py`

Reuse the look of `plot_lod_curves.py` so each reads as a Figure-1-style panel
(Q1): LoD (%) vs panel size, one curve per depth, 10/90 band across mixtures
shaded, dashed horizontal lines at 0.5% and 1%, log y-axis. Annotate the
contamination-floor plateau where the curve stops falling.

Produce the two tests as **separate figures** (they measure different things):
- `output/facts/fig_subsample_lod_mle.png` -- magnitude (MLE) LoD; optionally
  overlay the bonus analytical `lod_fraction` curve as a faint dashed line on
  the same axes (both are MLE-side, so this stays one figure).
- `output/facts/fig_subsample_lod_presence.png` -- presence-test LoD.

Per Q2, build these **standalone (real data only)** for now: no simulated
overlay. Whether to pair either with the simulated Figure 1, whether to combine
the two real-data panels, and what to do about the sibling panel (no real analog
in this unrelated-only dataset) are later decisions; separate self-contained
figures leave every option open.

### 6d. Snakefile wiring (`paper/Snakefile`)

Add a self-contained rule block alongside the other LoD rules:

```python
SUBSAMPLE_LOD_FACTS = [
    f"{FACTS_DIR}/subsample_lod_raw.csv",
    f"{FACTS_DIR}/subsample_lod_summary.csv",
    f"{FACTS_DIR}/subsample_lod_headline.csv",
]
SUBSAMPLE_LOD_FIGS = [
    f"{FACTS_DIR}/fig_subsample_lod_mle.png",
    f"{FACTS_DIR}/fig_subsample_lod_presence.png",
]

rule subsample_lod:
    input:
        # committed snapshot VCFs are the real inputs (declare for rerun-on-change)
        glob.glob("paper/public_data/SRP434573/genotypes/*.vcf.gz"),
    output: SUBSAMPLE_LOD_FACTS
    shell: "python paper/scripts/run_subsample_lod.py"

rule plot_subsample_lod:
    input: f"{FACTS_DIR}/subsample_lod_summary.csv"
    output: SUBSAMPLE_LOD_FIGS
    shell: "python paper/scripts/plot_subsample_lod.py"
```

Per Q2, wire these into the `ALL`/default target so `snakemake` builds them, but
**do not** add them to the `paper` rule's inputs yet: the facts and figure are
produced, but the rendered paper does not depend on them until you decide to pair
this with Figure 1. That keeps the build green while the figure is evaluated.

Two repo-specific gotchas: per the Snakemake-rerun rule, editing the script body
or grids forces a rerun of this (slow) rule, so flag it before committing; and
list every emitted fact key as a rule output so the build does not fail on a
stale facts CSV (the `supp_synthetic` gotcha) if/when it is later wired into the
paper.

### 6e. Prose + caption (deferred, per Q2)

No paper prose yet. The facts CSVs and figure are produced for evaluation only.
Once you decide whether/how to pair this with Figure 1, add a short paragraph to
the real-data results section, the honest framing from section 7 in the caption
and Methods, and the `{{ subsample_lod.* }}` fact placeholders.

## 7. Honest framing for the paper (must be explicit)

- "Sub-sampled real reads," not "an independent low-depth experiment." State the
  global-rate binomial-thinning method, that it is the exact analog of FASTQ /
  BAM read subsampling for this dataset (section 4a), and that it inherits the
  deposited library's bias and overdispersion (as would any one-library
  subsample) but does not add independent-library variance.
- "Pseudo-replicates from one real draw per mixture" (section 4c).
- Unrelated donors only; the real curve is the analog of Figure 1's left panel.
- The plateau at low fraction is the **contamination floor of this dataset**,
  not allomix's intrinsic limit; depth/markers cannot push below it. This is the
  concrete real-data version of the Discussion's "overdispersion/contamination,
  not depth, is the limiter."
- Expect the real LoD to sit **above** the simulated LoD at matched depth/markers
  (this is the predicted and desired result, confirming the "sim is best-case"
  statement). If it does not, that itself is a finding to chase, not paper over.

## 8. Suggested phasing

1. **Spike (half day):** `thin_informative_markers` (global-rate) + a throwaway
   script that, for one mixture, prints **both** the MLE host-fraction (with CI)
   and the presence detection rate per fraction at mean depths {2000, 500, 100}
   and panel sizes {full, 100, 25}. Confirms the LoDs move into the 0.5-10%
   window as the data is degraded, that the two tests separate, and that the
   contamination floor shows up (and how it hits the MLE empirical rule).
2. **Build the sweep:** full `run_subsample_lod.py` running both estimators per
   cell + facts CSVs (with the `test` column), quick-mode wired, unit test for
   the helper.
3. **Two standalone figures** (MLE and presence, real data only, no sim overlay).
4. **Snakefile wiring** into the default target but not the `paper` rule (Q2).

## 9. Decisions (resolved) and the one remaining question

Resolved in review:

- **Q1 (LoD definition / look).** Make it **look like Figure 1** (same axes,
  per-depth curves, bands, log y), on real data not sim. The MLE magnitude
  estimate and the presence test are **separate, complementary tests**, so each
  gets its own LoD sweep and figure (section 5), not a primary/secondary pair.
  Both use blank-free per-sample detection rules; the EP17/MLE-over-blanks
  construction is avoided because real data has no clean blank (Q3 below).
- **Q2 (figure scope).** **Standalone real-data figure** for now, no simulated
  overlay. Wire into Snakemake but **not** into the `paper` build rule. Pairing
  with Figure 1 and handling of the sibling panel are deferred to a later call.
- **Q4 (BAM down-sampling).** **Dropped.** The global-rate AD thinning is the
  exact statistical analog of FASTQ/BAM read subsampling for this dataset
  (section 4a), so a `samtools view -s` re-pileup is not worth the compute, and
  AD thinning keeps the whole experiment reproducible from a fresh checkout.
- **Q5 (marker scope).** **Sub-sample the informative markers** (what the
  estimate uses). The panel-intervals-first alternative is dropped; informative
  markers are the quantity everything else is a proxy for, and it matches
  `run_lod_validation.py`.

### Q3 expanded: why "blanks" are hard here (and why we avoid them)

EP17-A2 (the framework behind the simulated Figure 1) builds the LoD in two
steps: first a **Limit of Blank (LoB)** from samples that contain *zero* of the
thing being measured, then the LoD as the lowest true level reliably
distinguished from that LoB. The blank is the load-bearing part: it is the
estimator's output when the analyte is truly absent, and it sets the false-
positive threshold.

For SRP434573 the analyte is the minor (host) contributor, and **there is no
sample in the committed data with a true zero of a real contributor**: the
titration only goes down to 0.5%, and the pure single-individual runs are not in
the committed snapshot (only the mixture timepoints are). So to run the EP17
route we would have to *manufacture* blanks. The natural trick is a
**cross-individual negative control**: for a mixture of host `H` + donor `D`,
ask allomix to estimate the fraction of some third individual `G` who is not in
the mixture (all 7 individuals F1-F3, M1-M4 appear as host or donor in *some*
mixture, so each one's genotype is available somewhere). `G` contributed no DNA,
so in principle its estimated fraction is a blank.

Two problems make this a poor blank:

1. **Genotypes are scattered and batch-specific.** Each mixture's committed
   genotype VCF holds only its own two individuals, joint-called in a separate
   GATK run. Pulling `G` from a different mixture's VCF means the site sets and
   filtering differ between batches; `classify_markers` intersects by
   (chrom,pos,ref,alt) so it still runs, but on fewer markers. A uniform genotype
   set for all 7 would need a one-off joint-call of the pure individuals together
   (a `/tau` job), which Q4's "no extra compute" steer also argues against.
2. **The blank is not actually clean (the real killer).** `G` is unrelated to
   `H` and `D`, but unrelated humans still share roughly half their common-SNP
   alleles by chance. At the markers where `G` carries an allele `D` lacks, the
   real contributor `H` will, at a chunk of those markers, happen to carry the
   same allele, so the admix shows signal there. The "blank" therefore measures
   sequencing error **plus** `H`/`G` allele-sharing, biasing the LoB upward by an
   amount that depends on which `G` you pick. That is not the quantity EP17 wants.

The presence test sidesteps all of this: its null hypothesis *is* "host absent,
reads are sequencing error", evaluated per sample from the error-rate background
(optionally the per-site error table). That is precisely what an empirical LoB
tries to estimate, but computed analytically and per sample, with no genotype
sourcing and no allele-sharing contamination. Hence the presence-test detection
call is both the cleaner and the more defensible LoD definition on real data,
and it still plots as a Figure-1-style curve.

**Resolved:** use the blank-free per-sample detection rules for both tests
(MLE: host-fraction CI excludes 0; presence: LRT p < 0.05). The EP17-from-blanks
construction is not pursued (no clean real blank, and the manufactured
cross-individual blanks are biased). Nothing here is gated on `/tau`. All
decisions are now settled; this plan is ready to implement.
</content>
</invoke>
