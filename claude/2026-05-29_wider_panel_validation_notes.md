# 2026-05-29 — Wider Panel Validation + Pipeline Fixes

End-of-day write-up of the work that delivered Step 22 (pileup-only commit) and
Step 23 (wider force-output panel). Captures what changed, what the numbers
look like, and the open items to follow up on.

## What landed (in dependency order)

1. **Pipeline additive BED/panel semantics.** `pipeline/Snakefile`: the
   `intervals` BED is required for `-L` discovery scope, and the
   `panel_alleles_vcf` is purely additive via `--force-output-intervals`.
   Previously the panel VCF overrode the BED, restricting GATK to only the
   force-output sites. Config comments in `pipeline/config.yaml` and
   `doc/joint_calling.md` updated to match.

2. **Step 22 — pileup-only CLI.** `--vcf` removed from `monitor` and
   `timeline`; `--panel-vcf` + `--admix-vcf` are now required.
   `_resolve_vcf_inputs` deleted. `estimate-bias` / `estimate-errors`
   unchanged. Integration + multidonor tests rewritten as panel/admix
   pairs (joint VCF passed twice for synthetic fixtures). README, CLAUDE.md,
   doc/joint_calling.md updated.

3. **gnomAD-derived force-output panel build.** New script
   `scripts/build_force_output_panel.sh` + reference data
   `pipeline/gnomad_refseq_to_hg38_chrs.tsv`. Takes a BED + gnomAD v4.1
   sites VCF + AF threshold and emits a filtered, chr-renamed, indexed
   panel VCF. Filters: PASS / EXOMES_FILTERED / GENOMES_FILTERED (drops
   BOTH_FILTERED only), biallelic SNPs, INFO/AF >= threshold.

   Source: `/data/annotation/VEP/annotation_data/GRCh38/gnomad4.1_GRCh38_contigs.vcf.gz`
   (RefSeq-accessioned, 58 GB, tabix-indexed). The chr-rename file maps
   `NC_000001.11` → `chr1` etc.

   Built locally (haem BED at `output/union_sid_haem_vendor_probes.bed`,
   432 regions):

   - `output/union_sid_haem_gnomad_af05.vcf.gz` — 258 sites at AF >= 0.05
   - `output/union_sid_haem_gnomad_af01.vcf.gz` — 329 sites at AF >= 0.01
     (strict superset, contains all 71 existing SID panel sites)

   Verified all 71 autosomal SID SNPs are in gnomAD at AF >= 0.01; the
   5 sex-chromosome forensic markers (AMEL/SRY/ZFXY) are not in gnomAD
   and not chimerism panel members anyway.

4. **Pipeline join bug fix: `bcftools call -A`.** Without `-A`,
   `bcftools call -C alleles` strips the panel ALT to `.` at admix sites
   with zero ALT reads, breaking the `(chr,pos,ref,alt)` join with the
   panel VCF. Affected ~48 informative SNPs per patient. Fix is one
   character in `pipeline/Snakefile:401`.

5. **Pipeline panel_tsv fix: drop `ALT="."` rows.**
   `pipeline/Snakefile:354-355` — adds `-e 'ALT="."'` to `bcftools query`
   so the targets TSV excludes force-output sites that GATK genotyped as
   hom-ref in both host and donor. Two reasons: `bcftools call -A`
   produces malformed PL records at those rows, and those sites are
   uninformative anyway (host == donor == hom-ref by GATK's evidence).

6. **allomix: drop indels at parse_vcf.**
   `src/allomix/genotype.py:118-126` — skip records where
   `len(REF) != 1 or len(ALT) != 1`. Indels can't be reliably counted by
   straight pileup the way local-reassembly GATK does, producing
   systematic admix=0-ALT at panel-het indel sites.

7. **allomix: GT/AD consistency check (panel samples only).**
   `src/allomix/genotype.py:165-181` — new `gt_ad_consistency=True` flag
   on `parse_vcf`, wired through `cmd_monitor` and `cmd_timeline` in
   `cli.py` for host/donor (not admix). Drops markers where the called
   GT contradicts the AD VAF: het outside [0.35, 0.65], hom-ref VAF > 0.05,
   hom-alt VAF < 0.95. Catches GATK's miscalls (e.g. 21% VAF called het
   from marginal evidence in a 2-sample joint call).

## Why each filter was needed

run8 (just `-A` + the wider panel) recovered markers but biased the
estimate hard against flow truth (BHOA donor% dropped 99.86 → 82.15,
flow says 100%). GoF was a misleading 1.0 because the MLE found an
internally consistent fraction for the data, just one that disagreed
with biology. Investigation showed the newly-recovered markers were:

- **Indels** that pileup couldn't count (donor het by GATK, admix
  AD=2043/0 because mpileup doesn't track the deletion reads)
- **Marginal hets called by GATK from ~20% VAF reads** (donor AD=432/114
  at chr21:34792306, called 0/1 but VAF says it's not a real het)

Filtering these out at parse time recovers the marker-count win without
the bias.

## Numbers: run6 (pre-Step 23) vs run9 (post-Step 23)

Flow truth = 100% donor on CD45/CD3/CD13 for all samples below except
RCAR (mixed) and BCOL (CD13 88.30%).

| Sample | run6 n / donor% | run9 n / donor% | Δ markers |
|---|---|---|---|
| QUDO TP1 | 47 / 100.00 | 57 / 100.00 | +10 |
| QUDO TP2 | 47 / 98.74 | 57 / 99.03 | +10 |
| RCAR | 131 / 41.12 | 131 / 41.12 | 0 |
| PNOL | 29 / 100.00 | 49 / 100.00 | +20 |
| NDAD | 91 / 99.68 | 146 / 99.75 | +55 |
| BHOA | 103 / 99.86 | 144 / 99.91 | +41 |
| GBRI | 85 / 100.00 | 128 / 100.00 | +43 |
| PCAH | 95 / 99.04 | 132 / 99.17 | +37 |
| LNAN | 60 / 99.94 | 76 / 99.80 | +16 |
| BCOL | 98 / 97.47 | 150 / 97.59 | +52 |

run9 batch: `output/validation_run9/batch.tsv`.

Host-presence ("any detection") run9:

| Sample | host_f_est | p | call |
|---|---|---|---|
| QUDO TP1 (pure donor) | 0% | 1 | not detected ✓ |
| QUDO TP2 | 0.87% | 4.4e-44 | **detected** |
| PNOL (pure donor) | 0% | 1 | not detected ✓ |
| NDAD | 0.22% | 3.1e-19 | **detected** |
| BHOA | 0.08% | 1.3e-4 | **detected** |
| GBRI (pure donor) | 0% | 1 | not detected ✓ |
| PCAH | 0.76% | 3.6e-124 | **detected** |
| LNAN | 0.09% | 0.16 | not detected (borderline) |
| BCOL | 2.38% | 0 | **detected** |

Compared with run6, p-values are stronger across all "detected" samples
because of the higher marker count; magnitudes are preserved.

## Files added / changed

Pipeline:
- `pipeline/Snakefile` — additive BED/panel semantics, `-e 'ALT="."'` in
  panel_tsv, `-A` on bcftools call.
- `pipeline/config.yaml` — `intervals:` documented as required; panel
  documented as additive.
- `pipeline/gnomad_refseq_to_hg38_chrs.tsv` — new, RefSeq → UCSC map.

Scripts:
- `scripts/build_force_output_panel.sh` — new, reproducible panel build.

Source:
- `src/allomix/genotype.py` — indel skip; `gt_ad_consistency` option.
- `src/allomix/cli.py` — `--vcf` removed; `gt_ad_consistency=True` for
  panel samples in monitor and timeline.

Tests:
- `tests/test_integration.py` — 21 CLI tests rewritten as panel/admix
  pairs; dual-mode tests deleted.
- `tests/test_multidonor.py` — 4 CLI tests rewritten the same way.

Docs:
- `README.md`, `CLAUDE.md`, `doc/joint_calling.md` — two-VCF input now
  the only mode.

Built artefacts (local only, gitignored):
- `output/union_sid_haem_gnomad_af05.vcf.gz` (+ tbi) — recommended panel.
- `output/union_sid_haem_gnomad_af01.vcf.gz` (+ tbi) — wider alternative.

Pipeline config used on frgeneseq04:
- `pipeline/frgeneseq04_haem_big_bed.yaml` (user-side) currently points
  `panel_alleles_vcf` at the af05 panel.

## Open items

### To follow up with the post-doc (user's note)

- How did the expanded BED + gnomAD force-output panel perform on real
  samples? run9 numbers (above) are what to share. Specifically:
  - Marker count gains per patient (run6 → run9 in the table above).
  - That n_informative no longer correlates with cohort size — same per
    patient regardless of who's on the run.
- How did the "any detection" host-presence test perform? It detected
  sub-1% host in QUDO TP2, NDAD, BHOA, PCAH; not detected on pure-donor
  samples (QUDO TP1, PNOL, GBRI). LNAN at p=0.16 is borderline — worth
  discussing whether to treat as positive or negative.

### Possible Step 24 — calibrate beta-binomial overdispersion to wider panel

NDAD, BHOA, PCAH have `gof_pval = 0.0000` in run9 → flagged QC=REVIEW.
The model fits the chimerism fraction fine but residuals exceed
beta-binomial expectation. The wider panel has more diverse per-marker
behaviour than the overdispersion `rho` is currently calibrated for.
Two ways forward:

- Refit `rho` on the wider panel using `scripts/diagnose_sample.py`.
- Move toward per-marker (or per-marker-type) overdispersion as noted
  in Step 21 TODO.

Not urgent — the QC=REVIEW flag still surfaces these for human review,
and the chimerism numbers themselves match flow truth.

### Plan-doc tidying

The Step 22 section in `claude/allomix_overall_plan.md` has the original
"if we commit to pileup-only, the work to do:" subsection still embedded
under the ✅ COMPLETE header. Worth cleaning up next time someone is in
the file. Not urgent.

### Comparison plot needs run-replacement

`scripts/plot_chimerism_comparison.py` was used today to refresh the
run1-vs-run2-vs-run3-vs-run6 plot. Once Step 23 lands properly with the
post-doc's review, regenerate as run1-vs-run2-vs-run3-vs-run9 (or
similar) and replace `output/run1_vs_run2_vs_run3_vs_run6.png` in the
canonical comparison.

### Pipeline-rerun checklist for frgeneseq04 (one-liner)

If reruns are needed after a Snakefile change, the dependencies are:
1. Edit Snakefile.
2. Delete the rule's output files (Snakemake caches by mtime + content).
3. Rerun snakemake with the same config.
4. Copy back to `output/joint_called/` locally (note: directory has the
   trailing 'd' — `joint_called` not `joint_call`).
