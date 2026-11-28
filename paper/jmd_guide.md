# JMD author guide: what applies to this paper

The target journal is The Journal of Molecular Diagnostics (JMD). This file summarises the
parts of the JMD Guide for Authors that affect this paper. It is a summary, not a copy.
Source: the guide on ScienceDirect
(https://www.sciencedirect.com/journal/the-journal-of-molecular-diagnostics/publish/guide-for-authors),
revised April 2026 and saved on 2026-11-28. Re-check the guide before submission; it changes.

The paper fits the **Regular Article** type. **Technical Advance** (a methodology or
analytical tool, with a detailed description of the method's design and of how it improves
practice) is the alternative, and it has the same limits and charges.

## Length and limits

| Item | Regular Article / Technical Advance |
|---|---|
| Abstract | **220 words or fewer** (the only hard limit) |
| Main text | No official word limit |
| References | No official limit |
| Figures and tables | No official limit |
| Page charges | US $75 per typeset page for AMP/ASIP regular members, $100 otherwise. Opting into Open Access replaces page charges with the APC (articles accepted from 2026-01-01) |

The editors can ask for the text to be shortened or for data to be removed. Since the
charges are per page, length costs money even without a limit.

`paper.toml` `[limits]` enforces the abstract limit (`vibepaper wordcount --check`).

## Abstract

- 220 words or fewer, on its own page, readable without the text.
- Cover the background, methods, results and significance.
- **Not structured**: no section headings.
- Avoid abbreviations, website links and citations.

## Structure and style

- Section order: Introduction, Materials and Methods (with any ethics statement), Results
  (with any data-sharing information), Discussion. These are optional and come after the
  Discussion: Conclusions, Acknowledgements, Author Contributions.
- Key Points are optional: 3 or 4 bullets at the start.
- A graphical abstract is optional (TIF, JPG or PDF; generative AI is not allowed for it).
- **Keep first-person pronouns (I, we, our, us) to a minimum.** Prefer the passive voice.
- **Tense**: past tense in Results, present tense in Results headings and in figure
  legends (except when describing what was done).
- Spell out an abbreviation at first use. Gene and locus symbols are italic; protein
  symbols are not.
- Statistics: name each test, the sample size and the measures reported. Sensitivity and
  specificity need confidence intervals. Do not use "discrepant analysis". Say whether the
  sample size was set statistically, and give the replicates per condition.
- **Methods must be complete in the main text.** The guide says method and reagent detail
  in supplementary material is not acceptable, and online protocol links are not
  permitted. Analytical software counts as a novel material to be described well enough to
  reproduce.
- Code and data: publication implies sharing, on reasonable request, of whatever is needed
  to verify the claims, including analysis code.
- STARD: optional; submit the checklist if the paper should carry the STARD footnote.
- Sex and gender: define the terms used, or discuss the gap as a limitation (SAGER).
- **AI disclosure**: if AI tools were used in writing (beyond grammar and spelling
  checks), add a section titled "Declaration of AI and AI-assisted technologies in the
  writing process" at the end, using the guide's set statement. AI must not be used to
  create or alter figures, unless it is part of the research method.

## Title page

Concise title, not a question. Full author names and affiliations. The number of text
pages, tables and figures. A running head of 40 characters or fewer. Funding and grant
numbers. The corresponding author's contact details (at most two co-corresponding
authors). Conflict-of-interest disclosures.

## References

- Numbered in order of citation, including citations in tables and legends. Citations
  that first appear in supplementary data are numbered by where the item is first called
  out in the text. There is no separate supplementary reference list.
- List every author; use "et al" only for more than 35 authors (then the first 6 names).
- Websites, product inserts and unpublished data are cited in the text only, not in the
  reference list.
- The CSL in `paper.toml` (`paper/the-journal-of-molecular-diagnostics.csl`) handles the
  format.

## Figures

- White background. Choose grayscale or colour so that readers can tell the data points
  apart. No pattern fills and no 3D bars for single data sets.
- Width: one column (8 cm) or two columns (17 cm), with a maximum height of 22.5 cm. Each
  figure and all its panels (A, B, C) on one page.
- A plain font such as Arial. Labels readable at print size, with one capitalisation style
  across all figures.
- **Numbers of 10,000 or more take commas.** Put a space before units (2 kDa) and around
  math operators (2 + 2 = 4).
- Axis labels in one consistent form, for example "Depth (x)" or "Depth, x", not a bare
  unit.
- Legends go in a list in the manuscript. They must explain every symbol, line style and
  significance marker.
- Files: PDF, TIFF or JPG in RGB. No PowerPoint.

## Tables

Word tables placed after the main text, one per page. Black and white text only, with no
shading or images. Emphasis may be shown with bold, italic or underline, and must be
explained in the legend. Every symbol and abbreviation is defined in the legend.

## Supplementary data

Peer reviewed and published online; the editors may move items into the main article.
Supplementary figures are PDF/TIFF/JPG, tables are Word or Excel, and legends are Word.
Keep it small. Methods and reagent details must not be supplementary (see above).

## Submission

- The first submission can be in any format; style rules apply from the first revision.
  Text for revisions must be in Word, double-spaced, with page and line numbers.
- Cover letter: conflicts of interest, no prior or concurrent publication, and all authors
  and acknowledged contributors have approved. Disclose a preprint if there is one; JMD
  does not prefer preprints but may accept them.
- Review is blinded. Minor revisions are due in 30 days and major revisions in 90.
