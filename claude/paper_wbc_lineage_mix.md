# Whole-blood vs lineage chimerism: what allomix should say

Tracking issue: SACGF/allomix#6
Source papers (in `claude/papers/`):
- Clark et al., BJH 2025 (UK BSBMTCT/UKCGG consensus guidelines)
- Kakodkar et al., Front. Genet. 2023 (AlloSeq HCT validation, including CD3/CD66 vs WB comparison)

## The gap in the current paper

The paper presents `% donor chimerism` as a single scalar estimated from a VCF. Every result and every claim about sensitivity (MAE <1% across depths, 0.6% MAE benchmark, CI coverage, etc.) is framed against bulk DNA mixtures. The methods and simulation never mention that the input DNA usually comes from a heterogeneous cell population, and the discussion never addresses the difference between *analytical* sensitivity (assay LoD on bulk DNA) and *clinical* sensitivity (smallest shift in a clinically informative cell subset that the assay can resolve in the analysed specimen).

The issue thread frames the gap clearly:

- Whole blood is a blend of lineages (T, B, NK, myeloid, erythroid); each has its own donor fraction.
- A measured WB donor fraction is a lineage-abundance-weighted average. Two patients with identical WB chimerism can have very different lineage compositions, and the clinically meaningful shift can be hidden in a small lineage.
- Worked example from Dave: if the lineage of interest is 5% of WBCs and a clinically meaningful shift in that lineage is 5 percentage points, the corresponding WB shift is 0.05 × 0.05 = 0.25%. To resolve that in WB you need sub-0.25% sensitivity, comfortably below STR's 1–5% LoD and tight against modern NGS LoDs (0.06–0.3%).
- David Ross's reply: there is no international standardisation, no agreed clinically meaningful threshold for sub-STR detection. UK consensus (Clark 2025) is to test WB and store sorted fractions for later if needed.
- Wendy's reply: STR has high precision variability (cited as ~20% CV); the WB STR/NGS discrepancy in our example sample is likely within CI.

The paper currently makes none of this explicit. A reader could reasonably conclude that <1% MAE in WB is sufficient for clinical chimerism monitoring, when in fact the *adequacy* of any analytical sensitivity depends on the upstream specimen choice (WB vs sorted lineage) and on which clinical question is being asked.

## Empirical anchor from the Kakodkar paper

Kakodkar 2023 directly quantifies the lineage effect using AlloSeq HCT on matched WB and sorted CD3+/CD66+ fractions:

- 36 patients with mixed chimerism: CD3-enriched showed mean **7.1 ± 7.0-fold higher host % detection** than WB (range 1× to 38.9×). This is a lineage-amplification factor, exactly Dave's `1 / lineage_fraction` scaling.
- CD66+ (myeloid) was a near-useless surrogate (Pearson 0.73 with STR, R² = 0.535) because neutrophils reconstitute earliest and were ~100% donor in nearly all patients. Lineage choice matters: not every lineage is informative at every timepoint.
- The same NGS assay achieved 0.3% LoD on bulk DNA. The 7-fold lineage gain on top of that is what brings the *clinical* sensitivity down to the regime Dave's calculation requires.

These numbers should appear in our discussion as concrete justification for why WB-only validation is incomplete.

## What allomix actually does (and doesn't) do

Worth being clear-eyed about this:

- allomix is panel-agnostic and **lineage-agnostic**. It estimates the donor fraction in whatever DNA is in the VCF. Run it on a CD3-sorted VCF and you get CD3 chimerism; run it on WB and you get WB chimerism. The tool itself does not know or care.
- allomix cannot recover lineage composition from bulk DNA. The lab/clinical workflow choice (sort or not, which markers, which fractions) determines what the number means.
- Therefore the question "what sensitivity is required" is not a tool question; it is a workflow question. The tool's job is to be at least as good as the lineage-specific assay sensitivity demands. Our <1% WB MAE is already below STR LoD; whether it is below the *required* sensitivity depends entirely on the cell-mix and the lineage of interest.

## Recommendations: changes to make to the paper

### 1. Introduction — add one sentence acknowledging cellular composition

After the paragraph describing STR limitations, add a brief framing sentence that chimerism is intrinsically a property of a cell population, not of "blood", and that clinical practice often involves lineage-specific testing. One sentence, not a paragraph. Cite Clark 2025 (UK consensus) and Kakodkar 2023.

### 2. Methods — caveat the simulation framework

The simulation section currently says "synthetic chimeric VCFs by blending two genotype VCFs at a specified donor fraction." That is true at the DNA level but glosses over what the donor fraction physically represents. Add a brief note that the simulated fraction corresponds to the donor proportion in the analysed DNA, regardless of whether that DNA was extracted from unfractionated whole blood or from a sorted lineage. The same statistical framework applies; only the interpretation changes.

### 3. Results — no changes needed to existing results

The depth/relatedness/multi-donor/timeline results are valid as bulk-DNA characterisations of the estimator. Don't restate them in lineage terms; that would overclaim. Just frame them carefully in the discussion (next item).

### 4. Discussion — add a new subsection: "Cellular composition and clinical sensitivity"

Insert before "Limitations and Future Directions". Suggested content:

> Chimerism is a property of a cell population, not of whole blood per se. A donor fraction measured in unfractionated peripheral blood is a lineage-abundance-weighted average of the donor fractions across T cells, B cells, NK cells, granulocytes, monocytes, and other constituent lineages. Lineage proportions vary substantially between patients and across time post-transplant, and a clinically meaningful change in a small lineage may produce only a fractional change in the whole-blood signal. As a worked example, a 5-percentage-point shift in a CD3+ T-cell population that constitutes 5% of post-transplant lymphopenic whole blood corresponds to a whole-blood donor-fraction change of approximately 0.25%, comfortably below the 1–5% sensitivity of STR-based chimerism monitoring and at the edge of current NGS analytical sensitivities.
>
> Lineage-specific testing addresses this directly. Kakodkar et al. reported a mean 7.1 ± 7.0-fold (range 1× to 38.9×) higher host detection in CD3+-enriched samples compared to matched whole-blood inputs from 36 mixed-chimerism cases, with CD3 the most informative subset for relapse-relevant signal. The same study found CD66+ (myeloid) chimerism a poor surrogate (Pearson 0.73 versus STR), since granulocytes reconstitute earliest and were near-fully donor in most patients regardless of the underlying clinical state.[@Kakodkar2023alloseq]
>
> The 2025 UK BSBMTCT/UKCGG consensus recommends whole-blood STR chimerism at day +30 with lineage-specific testing (CD3, CD15) for mixed cases or specific clinical contexts (haemoglobinopathies, inborn errors), and explicitly notes that intervention thresholds derived from STR may not translate cleanly to higher-sensitivity methods. Clinical utility of sub-STR microchimerism detection remains under active investigation.[@Clark2025bjh]
>
> allomix is lineage-agnostic by design: it estimates the donor fraction in the DNA represented by the input VCF, whether that DNA was extracted from whole blood or from a sorted cell subset. The tool itself imposes no constraint on the upstream specimen, and the same validation results apply to lineage-sorted inputs. However, the choice of specimen, and any decision to sort cells before sequencing, remains a clinical and laboratory workflow decision that determines what the resulting chimerism number means clinically. allomix's analytical sensitivity (sub-1% MAE in our in silico bulk-DNA validation) is a necessary but not sufficient condition for adequate clinical sensitivity in any given monitoring scenario; the latter depends on the lineage composition of the analysed specimen and on the clinical question being asked.

### 5. Limitations — add an explicit caveat

Add a bullet or sentence to "Limitations and Future Directions":

> Our in silico validation characterises the estimator on bulk DNA mixtures and does not address the cellular composition of the source specimen. Clinical validation should include matched whole-blood and lineage-sorted (CD3, CD15/CD66) samples to characterise concordance and to confirm that allomix's analytical sensitivity translates into the clinical sensitivity required by specific monitoring contexts (e.g., early lymphoid mixed chimerism in T-depleted reduced-intensity transplants).

### 6. Comparison table — qualify the LoD column

The Table 2 LoD figures (0.06%, 0.3%, etc. for commercial tools) and the allomix entry (~0.6% MAE) are all bulk-DNA numbers. Add a footnote to Table 2 along the lines of: "All LoD values refer to detection in bulk extracted DNA. Clinical sensitivity in a given specimen further depends on the proportion of the lineage of interest and on the upstream cell-sorting workflow."

### 7. References to add to `references.bib`

- Clark A et al., BJH 2025 (UK BSBMTCT consensus). DOI 10.1111/bjh.70061. Citation key suggestion: `Clark2025bjh`.
- The Kakodkar reference (`Kakodkar2023alloseq`) is already cited; we will be re-using it for the lineage-fold finding.

### 8. What we are deliberately NOT changing

- We are not adding any new simulation. The lineage point is conceptual and clinical, not a computational gap; running a "simulated lineage" experiment would be circular (we would just be re-blending DNA at a different fraction).
- We are not redefining or proposing new clinical thresholds. The UK consensus is explicit that no agreed sub-STR thresholds exist; we should match that humility.
- We are not weakening the existing analytical claims. <1% MAE in WB is real and reportable; the new framing simply contextualises what that number does and doesn't guarantee.
- We are not making a recommendation on whether labs should sort. That is a deployment decision (David Ross's GMP point in our internal review), not a tool decision.

## Quick wording check: avoid AI-tells per CLAUDE.md

The draft text above uses "comfortably below", "directly", "explicitly", "actively". I've avoided em-dashes, "robust", "comprehensive", "leveraging", "underscores", "landscape", "crucial", "pivotal", "nuanced". When porting to `discussion.md`, re-check for "furthermore", "moreover", "notably", "importantly".

## Summary of action items

1. One-sentence framing in `introduction.md` after the STR-limitations paragraph.
2. One-sentence caveat in `methods.md` simulation section about what the simulated fraction represents.
3. New subsection "Cellular composition and clinical sensitivity" in `discussion.md` (before Limitations).
4. Explicit caveat bullet in `discussion.md` Limitations.
5. Footnote on Table 2 LoD column.
6. Add `Clark2025bjh` to `references.bib`.
