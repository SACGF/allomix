# Wetlab validation: sample choice and mixture design

Design notes for the titrated-mixture wetlab arm. The analysis behind these
choices is kept because it is what justifies them, but the parts that have
already been written into the paper have been dropped, and repository and paper
follow-ups live in `plans/followups.md`.

Status key: **[verified]** = computed from allomix's own code, **[derived]** =
closed-form, checked by hand.

---

## Recommendation

Order the **GIAB Ashkenazim trio**: NA24143 (mother), NA24149 (father),
NA24385 (son). Preferably as **NIST RM 8392**, which is those three individuals
as a reference material: three vials, ~10 ug each, in TE, extracted from three
large LCL growths (confirmed from the RM 8392 report of investigation). Buy 3
units in one lot, up front.

Three mixture arms:

| Arm | Series | Purpose |
|--:|:--|:--|
| 1 | `mum <- dad` | unrelated baseline, on our own panel and run |
| 2 | `kid <- mum` | parent-child, sex-mismatched |
| 3 | `mum <- dad + kid` | multi-donor, mixed relatedness |

Notation is `recipient <- donor(s)`.

Arm 1 is not optional. Without it the only unrelated comparator is SRP434573, a
different panel, platform and laboratory, so any related-versus-unrelated
difference is confounded with all of that. Arm 1 is what makes arms 2 and 3
interpretable.

Drop `kid <- dad`: statistically near-identical to arm 2, differing only in being
sex-matched.

Optionally a small `kid <- mum + dad` arm as a **negative control**, to show the
tool reports no presence result rather than a p-value when the donor set covers
the recipient (see C).

Why this trio specifically:

- Sex structure suits the sex-mismatch requirement: mother F, father M, son M, so
  `mum <- dad` and `kid <- mum` are both mismatched.
- GIAB benchmark VCFs give an independent check on the genotyping phase, which
  matters because `docs/joint_calling.md` is load-bearing for the whole method.
  Check first that the panel markers fall inside the GIAB high-confidence regions;
  outside them there is no truth genotype.
- Open-consent PGP material, so the resulting FASTQs and VCFs can be released.
  That makes the series a reusable community resource, which fits the adoption
  goal. Note this covers the *data*, not the physical DNA, which stays under
  Coriell/NIST transfer terms.
- Buying one lot up front matters because LCL DNA acquires culture CNVs, and a
  lot change mid-study would be indistinguishable from a real effect given the
  CNV stress test in the paper.

Arm 2 has no full-contrast markers, so the co-pool contamination correction is
unavailable there (see B). Use unique dual indexes and keep the low-fraction
libraries off a heavily co-pooled flowcell.

---

## A. Minimum input material

Two separate constraints govern low-fraction measurement:

1. **Per-locus unique molecular depth** governs *detection*. Pooling across
   markers helps here.
2. **Input genome equivalents** governs *how accurately the nominal fraction is
   realised*. This error is **coherent across every marker**: the same minor
   genomes carry all loci, so if the tube realised 0.083% instead of 0.1%, every
   marker is biased identically. Pooling across markers cannot reduce it, and
   neither can sequencing deeper.

### The number **[derived]**

At 6.6 pg per diploid genome, input mass `M` (ng) supplies `151.5 x M` genome
equivalents, of which the minor component contributes `n = f x 151.5 x M`.
Poisson relative CV on that count is `1/sqrt(n)`, so for a target CV:

```
M_min (ng)  =  0.0066 / (f x CV^2)
```

Equivalently: **at least `1/CV^2` genome equivalents of the minor component**, so
100 for 10% relative precision.

| Minor fraction | Genome equivalents (200 ng) | Poisson CV | M for 10% CV | M for 5% CV |
|--:|--:|--:|--:|--:|
| 1% | 303 | 5.7% | 66 ng | 264 ng |
| 0.5% | 152 | 8.1% | 132 ng | 528 ng |
| 0.25% | 76 | 11.5% | 264 ng | 1.1 ug |
| 0.1% | 30 | 18.3% | 660 ng | 2.6 ug |

This is the ideal limit at 100% library conversion. Only a fraction `e` of input
molecules become sequenceable, so the real requirement is `M_min / e`. `e` is
measurable from our own data as deduplicated per-locus coverage over input genome
equivalents, which is why unique (deduplicated) depth is the number to report,
not raw depth.

### Design consequence

Pool size, not library input, pins the true fraction. Make mixture pools at ug
scale and aliquot from them: a 5 ug pool at 0.1% holds ~758 minor genome
equivalents (3.6% CV), where mixing at 200 ng library scale holds 30 (18%).
Replicate libraries from one pool then separate library-plus-sequencing variance
from mixture variance, which is the decomposition the overdispersion argument in
the Discussion currently makes by inference rather than measurement.

Per-library input: 200 ng is adequate at >=1%; 500 ng at 0.25-0.5%; 1 ug at 0.1%,
treated as exploratory. Confirm the capture kit's validated input range first,
since it may cap this regardless.

---

## B. Why parent-child rather than siblings

A parent and child share exactly one allele IBD at every autosomal locus, so
opposite homozygotes cannot occur: marker types 0 and 1 (the full-contrast
classes) have probability **zero**, not merely low. The co-pool contamination
correction acts only on those types, so it is unavailable with a parent-child
donor. The presence test still runs on types 10/11, all half-contrast.

### Marker yield, 76-marker panel, MAF ~U(0.2,0.5), 20,000 pairs **[verified]**

Computed through `MarkerType.classify` and `simulate.RELATEDNESS_IBD`.

| Relationship | Informative [5-95%] | Presence-capable | Full-contrast (0/1) | P(zero full-contrast) |
|:--|:--|--:|:--|--:|
| unrelated | 44.5 [37-51] | 26.0 | 7.4 [3-12] | 0% |
| parent-child | 33.5 [26-41] | 16.7 | **0.0** | **100%** |
| half-sibling | 39.0 [32-46] | 21.4 | 3.7 [1-7] | 2% |
| sibling | 27.8 [21-35] | 14.8 | 1.9 [0-4] | 15% |

### The N=1 argument

This is the deciding reason for a parent-child arm over a sibling arm.

Parent-child sharing is deterministic (IBD = 1 at every locus), so **one real
parent-child pair generalises to all parent-child pairs**. Sibling sharing is
stochastic (IBD 0/1/2 with probability 1/4, 1/2, 1/4), so a single sibling pair
is one draw from a wide distribution and characterises nothing. Averaging it out
would need several independent sibling pairs.

Note the two relationships are hard in different ways: siblings have fewer
informative markers overall (27.8 vs 33.5), so they remain the harder case for
magnitude LoD, while parent-child deterministically loses the full-contrast
classes. Haploidentical donors (parent to child, child to parent) are a growing
share of related-donor HSCT, so the parent-child arm is a real clinical
configuration rather than a contrived worst case.

### Available sibling material, if a sibling arm is ever wanted

- **NIGMS00018**: NA28756 (proband, F, 18y) and NA28757 (brother, M, 10y), GAMT
  deficiency. The only unambiguous sex-mismatched sibling pair with DNA. No
  parental DNA, so the sibship cannot be verified independently.
- **Family 3410**: NA27241 (M, 2y) and NA27461 (M, 3y), both affected, with both
  parents (NA27242 M 37y, NA27243 F 41y) also available as DNA. FOXG1 congenital
  Rett variant. Probably a genuine sibling pair plus parents, i.e. a quad, but
  sex-matched.

Some catalogue families do carry duplicate records for one individual: family
3508 lists two "fathers" both aged 38 YR. The duplicates share an identical age
at sampling (3222: 18/18, 3525: 16/16, 3508: 38/38), which is how 3410 (ages 2
and 3) can be told apart from them. Worth confirming with Coriell before relying
on either.

---

## C. Donor sets that cover the recipient

The constraint that decides the arm-3 configuration.

### The principle **[derived]**

What matters is how much of the recipient's genotype the **donor set** covers
between them:

- **Parent-child pair**: the donor covers *one* of the recipient's alleles at
  every locus. Kills opposite homozygotes. Presence test survives on
  half-contrast markers.
- **Child recipient, both parents as donors**: the donor set covers *both* of the
  recipient's alleles at every locus. A child has no private allele relative to
  its parents. Kills the presence test entirely.

### The presence-test result **[verified]**

`host_presence.py` selects markers where every donor is homozygous for the same
allele and the recipient carries the donor-absent allele. If both parents are
donors and the child is the recipient, that set is empty by construction: for
both parents to be B/B the child must have received B from each, so the child is
B/B and carries no donor-absent allele.

Confirmed through `select_donor_hom_markers` over 5,000 simulated Mendelian
trios: **zero presence-test markers in 100% of trios**.

| Configuration | Informative | Presence-test | P(zero) |
|:--|--:|--:|--:|
| kid <- mum + dad | 51.9 | **0.0** | **100%** |
| mum <- dad + kid | 51.9 | 9.3 | 0% |
| dad <- mum + kid | 51.9 | 9.2 | 0% |

The parent-as-recipient configurations recover markers because a heterozygous
parent has an *untransmitted* allele that neither the other parent nor the child
need carry.

### The magnitude result **[verified]**

The magnitude estimate does still work in the kid-as-recipient configuration.
Where both parents are homozygous (~25% of markers at MAF 0.5) the child's dosage
is exactly `(mum + dad)/2`, so adding child at fraction `f` is algebraically
identical to adding `f/2` to each parent, and those markers carry no information
about the child's fraction. Heterozygous-parent markers break the degeneracy. Via
`estimate_multi_donor`, 200 markers, 1000x, 40 reps:

| True recipient | kid <- mum + dad | mum <- dad + kid |
|--:|:--|:--|
| 50% | 49.70 +/- 0.32 | 50.01 +/- 0.26 |
| 10% | 9.94 +/- 0.51 | 10.27 +/- 0.22 |
| 2% | 1.90 +/- 0.37 | 2.43 +/- 0.12 |

So **2.3x the scatter at 10%, 3.1x at 2%**, and worse on a 76-marker panel where
the informative set is smaller and the degenerate quarter is a fixed proportion.

A 3x precision penalty would be tolerable; losing a headline readout entirely is
not. Hence arm 3 is `mum <- dad + kid`.

### Behaviour on the empty case **[verified]**

Relevant to the optional negative-control arm. `host_presence_test`
short-circuits on an empty marker set, and the report layer writes `NA` rather
than a p-value into the presence cells, so "could not test" cannot be misread as
"confidently no residual recipient". The QC verdict does not flag it, which is
tracked in `plans/followups.md`.

---

## D. Open items

- Confirm whether the Haem panel covers any chrY. If it does, a sex-mismatched
  arm gives an orthogonal estimate of the male fraction independent of the
  autosomal MLE. If not, the sex mismatch only exercises the default X/Y
  exclusion and `--use-sex-chroms`.
- Confirm the capture kit's validated input range (see A).
- Check the panel marker positions against the GIAB high-confidence BED before
  ordering.
- Confirm with Coriell whether family 3410's two probands are genuinely two
  children (see B).
