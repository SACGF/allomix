# Step 2: BAM vs VCF — Primary Input Format Decision

## What allomix needs at each marker

For every marker position in an admixture sample, allomix requires:

| Field | Purpose | Essential? |
|-------|---------|------------|
| **Ref allele count** | Numerator/denominator for VAF calculation | Yes |
| **Alt allele count** | Numerator/denominator for VAF calculation | Yes |
| **Total depth (DP)** | Confidence weighting in MLE, QC filtering | Yes |
| **Genotype (GT)** | Classify donor/host as hom-ref, het, hom-alt | Yes (genotyping samples) |
| **Genotype quality (GQ)** | Filter low-confidence genotype calls in donor/host | Yes (genotyping samples) |
| **Phred-scaled likelihoods (PL)** | Genotype uncertainty modelling | Nice to have |
| **Per-sample VAF (AF)** | Convenience — can derive from AD | Nice to have |
| **Strand bias (SB/FS)** | QC: flag markers with strand bias artifacts | Nice to have |

The critical inputs are the **allele depths (AD)** — the integer ref and alt read counts. Everything else is either derivable from AD+DP or used only for QC/filtering.

For genotyping samples (pure donor, pure host), GT and GQ are also essential: we need confident genotype calls to classify which markers are informative. For admixture samples, GT is less important (GATK will call hom-ref at low donor fractions), but the AD counts are everything.

---

## Are VCF AD fields sufficient?

**Yes.** The existing VCFs contain exactly what we need:

```
# Het variant site — has ref and alt counts
GT:AD:DP:GQ:PL:AF    0/1:1182,1280:2462:99:24424,0,22256:0.5199

# Hom-alt site — has ref and alt counts
GT:AD:DP:GQ:PL:AF    1/1:2,1999:2001:99:51143,5937,0:0.999

# Hom-ref site — has ref count ONLY
GT:AD:DP:GQ:PL        0/0:2926:2926:99:0
```

At variant sites, AD provides both ref and alt counts as comma-separated integers. At hom-ref sites, AD has only the ref count (single value, no alt). This is standard GATK GenotypeGVCFs behaviour — when there is no ALT allele, there is no alt count to report.

**For genotyping samples (pure donor/host):** This is fine. Hom-ref means "no alt allele present." The AD ref-only value plus the GT call gives us everything we need to classify the marker.

**For admixture samples:** This is the critical concern — see next section.

---

## The GATK low-fraction problem and its solution

### The problem

When a post-HSCT admixture sample has a low donor fraction (say 1%), the donor's alleles will be present at very low VAF:

- Host = A/A, Donor = A/B, donor fraction = 1% => expected B VAF = 0.5%
- Host = A/A, Donor = B/B, donor fraction = 1% => expected B VAF = 1.0%

At 0.5-1% VAF with 2000x depth, we expect 10-20 alt reads. GATK HaplotypeCaller, running on a single sample, will likely:

1. Call the position as **hom-ref (0/0)** because the alt reads look like sequencing error
2. Report AD with **only the ref count** (no alt count)
3. The alt read information is **lost from the VCF**

This is the scenario where a naive "just use the VCFs" approach breaks down.

### The solution: pipeline adjustment

The user controls the GATK pipeline and can modify it. Two approaches fix this, and both are standard GATK practice:

#### Option A: Joint calling with genotyping samples (recommended)

Include the admixture samples in the same GenomicsDBImport/GenotypeGVCFs joint call as the donor and host samples. When GATK joint-calls across a cohort, it evaluates ALT alleles discovered in **any** sample at **every** sample's position. If the donor is het A/B, the B allele is discovered from the donor's data, and GATK will report AD ref,alt counts for every sample at that position — including the admixture sample where B is at 0.5%.

This is already how the existing pipeline works (the GenomicsDBImport command in the example VCF header imports 100+ samples). The key requirement is that donor and host genotyping samples must be in the same joint call as the admixture samples. The pipeline already does this.

**Confirmed with real data:** The joint-called multi-sample VCF at `data/joint_called_example.vcf` (114 samples, 9 variant sites) demonstrates this directly. At variant sites, hom-ref samples have two-element AD fields with explicit zero alt counts:

```
# Hom-ref sample at a variant site in joint-called VCF — AD has ref,alt (3608 ref, 0 alt)
0/0:3608,0:3608:99:0,120,1800:0
```

This is exactly what allomix needs: every sample has ref AND alt counts at every variant position, regardless of genotype. The joint calling approach works.

#### Option B: Force-calling at known sites

Use GATK GenotypeGVCFs with `--force-call-filtered-alleles` or `--include-non-variant-sites`, or use `bcftools mpileup` at a BED file of the 76 marker positions. This guarantees allele counts at every position regardless of whether an ALT allele was independently discovered.

#### Verification

Either approach ensures that the output VCF contains AD with both ref and alt counts at every marker, even when the alt fraction is very low. The pipeline is already nearly there — the only requirement is that admixture samples are joint-called alongside genotyping samples, which is the existing workflow.

---

## What would BAM-level access add?

Working directly from BAMs (via pileup) would give:

| Feature | Available from BAM | Available from VCF | Verdict |
|---------|-------------------|-------------------|---------|
| Allele counts at any position | Yes (pileup) | Yes (if joint-called or force-called) | VCF sufficient |
| Strand-specific counts | Yes (per-read) | Partially (SB field if present) | BAM adds value, but SB in VCF is adequate for QC |
| Base quality per read | Yes | No (aggregated into QD/BaseQRankSum) | BAM adds minor value for error modelling |
| Mapping quality per read | Yes | No (aggregated into MQ/MQRankSum) | BAM adds minor value |
| UMI/molecular barcode dedup | Yes (if UMIs used) | No | BAM required **only if** UMI dedup is needed |
| Read-level phasing | Yes | Partial (PGT/PID) | Not relevant for biallelic SNP chimerism |
| Custom allele counting (e.g., min BQ filter) | Yes | No — counts are pre-computed | BAM adds flexibility |
| Indel realignment control | Yes | No — GATK already handles this | Marginal |

**Summary:** BAMs provide more flexibility and granularity, but for our use case (biallelic SNPs at >1000x depth with GATK-called VCFs from a controlled pipeline), the VCF AD fields capture the essential information. The scenarios where BAM access genuinely helps are:

1. **UMI deduplication** — not currently used in the rhAmpSeq pipeline, and would need to be handled upstream of allomix anyway
2. **Custom base-quality filtering** — marginal benefit at >1000x depth where random errors average out
3. **Troubleshooting outlier markers** — useful in development, not needed in production

None of these justify requiring BAMs as primary input.

---

## How reference tools handle input

| Tool | Primary input | What it extracts | Notes |
|------|--------------|-----------------|-------|
| **Demixtify** | BAM + BCF panel | Pileup allele counts at known sites | BAM-first because it operates on forensic samples where VCFs may not exist or may be called differently |
| **EuroForMix** | CSV (allele, height/count) | Pre-extracted allele counts | Input-agnostic — consumes counts, not raw data |
| **All-FIT** | TSV (VAF, depth, ploidy) | Pre-extracted VAFs | VCF-derived in practice |
| **chimerism_smmip** | BAM | Pileup counts at smMIP targets | BAM-first because smMIPs require UMI-aware dedup |
| **Conpair** | BAM or VCF | Pileup or VCF genotypes | Supports both |

The pattern: tools that work from BAMs do so because they either (a) need UMI dedup, (b) need to operate in environments where VCFs don't exist, or (c) want complete control over allele counting. Tools that work from VCFs or pre-extracted counts do so because the upstream pipeline is trusted to produce correct counts.

In our case, the upstream pipeline is controlled and trusted — we are building both the pipeline and the analysis tool.

---

## Minimum VCF FORMAT field requirements

allomix should require the following FORMAT fields in input VCFs:

### Required

| Field | VCF spec | Used for |
|-------|----------|----------|
| **GT** | Standard | Genotype classification (donor/host samples); ignored for admixture |
| **AD** | Standard (GATK) | Allele depths — the primary data for chimerism calculation |
| **DP** | Standard | Total depth — QC filtering and confidence weighting |

### Recommended (used if present, not required)

| Field | VCF spec | Used for |
|-------|----------|----------|
| **GQ** | Standard | Genotype quality filtering for donor/host samples |
| **PL** | Standard (GATK) | Genotype likelihood — could inform uncertain genotype handling |
| **AF** | bcftools | Per-sample VAF — convenience cross-check against AD-derived VAF |

### Validation rules allomix should enforce

1. **AD must have ref AND alt counts** at every marker in admixture samples (i.e., two comma-separated values). A single-value AD at a marker means the ALT allele was not evaluated — this is the pipeline problem described above.
2. **DP must be above a configurable minimum** (default: 100x) for the marker to be used.
3. **GT must be callable (not ./.)** for donor/host genotyping samples.
4. **GQ >= 20** (if present) for donor/host genotypes. Low-confidence genotypes risk misclassifying informative markers.

---

## Decision: VCF as primary input

**allomix will take VCFs as primary input.** BAM support is not planned for v1.

### Reasoning

1. **The data is already there.** The existing GATK pipeline produces VCFs with GT, AD, DP, GQ, PL, and AF at all 76 markers. No additional variant calling or data extraction is needed.

2. **The pipeline is controlled and adjustable.** The user owns the GATK/bcftools pipeline. If a VCF is missing alt counts at a position (the low-fraction concern), the fix is in the pipeline (joint calling or force-calling), not in allomix. Requiring allomix to parse BAMs to work around a pipeline deficiency is solving the wrong problem.

3. **VCFs are the natural interface boundary.** The GATK pipeline produces VCFs. VariantGrid consumes VCFs. allomix sits between them. VCF-in, results-out keeps the architecture clean and the tool pipeline-agnostic.

4. **Panel agnosticism.** A VCF-first tool works with any panel (not just rhAmpSeq) and any variant caller (not just GATK). A BAM-first tool would need to know the marker positions and handle pileup — coupling it to the panel definition. VCFs already encode the marker positions implicitly.

5. **Simpler implementation and testing.** VCF parsing (via cyvcf2) is ~50 lines. BAM pileup at known sites is ~200 lines and introduces edge cases (MAPQ filtering, base quality thresholds, duplicate handling, indel realignment). The BAM approach adds complexity without adding information we don't already have.

6. **GMP auditability.** A tool that reads a standard file format (VCF) and performs a statistical calculation is easier to validate than one that re-extracts data from BAMs. The variant calling is already validated as part of the existing pipeline.

7. **Reference tool precedent.** EuroForMix (the most statistically rigorous mixture tool) works on pre-extracted allele counts, not BAMs. All-FIT works on VAFs from VCFs. The BAM-first tools (Demixtify, chimerism_smmip) have specific reasons for BAM access (forensic context, UMI dedup) that don't apply here.

### What allomix will NOT do

- allomix will not call variants. It consumes pre-called VCFs.
- allomix will not perform pileup. If allele counts are wrong, the pipeline should be fixed.
- allomix will not require BAMs. If a future use case demands BAM-level data (e.g., UMI-aware counting), that would be a separate preprocessing step that produces a VCF or counts file, not a change to allomix's core input.

### Fallback: if VCFs prove insufficient

If during validation we discover that GATK's allele counts are systematically biased in a way that BAM-level access would fix (e.g., GATK's local reassembly in HaplotypeCaller alters counts at certain markers), we can:

1. Add a `bcftools mpileup` preprocessing step that generates raw pileup counts at the 76 positions — output is still a VCF
2. Provide a utility script that extracts allele counts from BAMs into a simple TSV that allomix can also consume
3. In extremis, add a `--bam` input mode — but this is unlikely to be needed

This fallback costs nothing to keep in reserve and doesn't affect the v1 architecture.

---

## Pipeline requirements for admixture VCFs

To ensure allomix gets what it needs, the GATK pipeline for admixture samples should:

1. **Joint-call admixture samples with donor/host genotyping samples** in the same GenomicsDBImport + GenotypeGVCFs run. This ensures ALT alleles discovered in donors are evaluated in admixture samples.

2. **Do not apply hard filters that remove low-VAF variants.** The current hard filters (QD<2, QUAL<30, etc.) are applied at the site level across the joint call. A site that is clearly variant in the donor will pass filters even though the admixture sample has low VAF at that site. This is already correct.

3. **Retain AD with both ref and alt counts.** The `bcftools view --trim-alt-alleles` step in the current pipeline trims unused ALT alleles, which is fine — but verify it doesn't collapse a biallelic site to ALT="." if the sample's ALT count is very low. If it does, use `bcftools view` without `--trim-alt-alleles`, or use `--min-ac 0` to retain all alleles from the joint call.

4. **Do not use bcftools setGT to convert low-depth to no-call at the genotype level** for admixture samples. The current `setGT -t q -n ./. -i FMT/DP=0` converts zero-depth sites to no-call, which is correct (no data = no call). But ensure this threshold isn't raised to something that would discard low-count sites.

---

## VCF audit script

The following script can be run on /tau to verify that existing VCFs have the fields allomix requires. It outputs only summary statistics — never raw patient data.

```python
#!/usr/bin/env python3
"""
Audit VCFs in a directory for allomix compatibility.

Checks that VCFs contain the required FORMAT fields (GT, AD, DP) and reports
summary statistics about field completeness, depth distribution, and AD format.

Usage:
    python audit_vcfs.py /tau/data/clinical_hg38/idt_rhampseq_sid/

Output: summary statistics only — no sample identifiers, no genotypes, no
genomic coordinates, no individual-level data.
"""

import sys
import os
from pathlib import Path
from collections import Counter

try:
    from cyvcf2 import VCF
except ImportError:
    sys.exit("ERROR: cyvcf2 not installed. Run: pip install cyvcf2")


REQUIRED_FORMAT = {"GT", "AD", "DP"}
RECOMMENDED_FORMAT = {"GQ", "PL", "AF"}
MIN_DEPTH = 100


def audit_single_vcf(vcf_path: str) -> dict:
    """Audit a single VCF and return summary statistics."""
    stats = {
        "n_variants": 0,
        "has_required_fields": True,
        "missing_fields": set(),
        "ad_single_value_count": 0,   # AD with only ref (no alt) — hom-ref issue
        "ad_two_value_count": 0,      # AD with ref,alt — what we need
        "ad_multi_value_count": 0,    # AD with >2 values — multiallelic
        "depths": [],
        "low_depth_count": 0,
        "nocall_count": 0,
        "format_fields_seen": set(),
        "parse_error": None,
    }

    try:
        vcf = VCF(str(vcf_path))
    except Exception as e:
        stats["parse_error"] = str(e)
        return stats

    # Check header for FORMAT fields
    for fmt in vcf.header_iter():
        if fmt["HeaderType"] == "FORMAT":
            stats["format_fields_seen"].add(fmt["ID"])

    for field in REQUIRED_FORMAT:
        if field not in stats["format_fields_seen"]:
            stats["has_required_fields"] = False
            stats["missing_fields"].add(field)

    # Iterate over records
    for variant in vcf:
        stats["n_variants"] += 1

        # Check GT
        gts = variant.genotypes
        if gts and gts[0][0] == -1:  # no-call
            stats["nocall_count"] += 1

        # Check AD format (single sample expected)
        ad = variant.format("AD")
        if ad is not None:
            ad_values = ad[0]  # first (only) sample
            n_alleles = sum(1 for v in ad_values if v >= 0)
            if n_alleles <= 1:
                stats["ad_single_value_count"] += 1
            elif n_alleles == 2:
                stats["ad_two_value_count"] += 1
            else:
                stats["ad_multi_value_count"] += 1

        # Check DP
        dp = variant.format("DP")
        if dp is not None:
            depth = int(dp[0][0])
            stats["depths"].append(depth)
            if depth < MIN_DEPTH:
                stats["low_depth_count"] += 1

    vcf.close()
    return stats


def summarise_depths(all_depths: list[int]) -> dict:
    """Compute depth summary statistics."""
    if not all_depths:
        return {"count": 0}
    import statistics
    sorted_d = sorted(all_depths)
    n = len(sorted_d)
    return {
        "count": n,
        "min": sorted_d[0],
        "q25": sorted_d[n // 4],
        "median": sorted_d[n // 2],
        "q75": sorted_d[3 * n // 4],
        "max": sorted_d[-1],
        "mean": round(statistics.mean(sorted_d), 1),
    }


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <vcf_directory_or_file> [--glob '*.vcf.gz']")
        sys.exit(1)

    target = Path(sys.argv[1])
    glob_pattern = "*.vcf.gz"

    # Parse optional --glob argument
    if "--glob" in sys.argv:
        idx = sys.argv.index("--glob")
        if idx + 1 < len(sys.argv):
            glob_pattern = sys.argv[idx + 1]

    # Collect VCF files
    if target.is_file():
        vcf_files = [target]
    elif target.is_dir():
        vcf_files = sorted(target.glob(glob_pattern))
        if not vcf_files:
            # Try uncompressed
            vcf_files = sorted(target.glob("*.vcf"))
    else:
        sys.exit(f"ERROR: {target} is not a file or directory")

    if not vcf_files:
        sys.exit(f"ERROR: No VCF files found in {target} matching {glob_pattern}")

    print(f"Auditing {len(vcf_files)} VCF file(s) in {target}")
    print(f"Glob pattern: {glob_pattern}")
    print("=" * 70)

    # Aggregate statistics
    total_files = len(vcf_files)
    files_ok = 0
    files_missing_fields = 0
    files_with_errors = 0
    total_variants = 0
    total_ad_single = 0
    total_ad_two = 0
    total_ad_multi = 0
    total_nocall = 0
    total_low_depth = 0
    all_depths = []
    missing_field_counts = Counter()
    format_field_counts = Counter()

    for vcf_path in vcf_files:
        stats = audit_single_vcf(str(vcf_path))

        if stats["parse_error"]:
            files_with_errors += 1
            continue

        total_variants += stats["n_variants"]
        total_ad_single += stats["ad_single_value_count"]
        total_ad_two += stats["ad_two_value_count"]
        total_ad_multi += stats["ad_multi_value_count"]
        total_nocall += stats["nocall_count"]
        total_low_depth += stats["low_depth_count"]
        all_depths.extend(stats["depths"])

        for field in stats["format_fields_seen"]:
            format_field_counts[field] += 1

        if stats["has_required_fields"]:
            files_ok += 1
        else:
            files_missing_fields += 1
            for field in stats["missing_fields"]:
                missing_field_counts[field] += 1

    # Report
    print()
    print("FILES")
    print(f"  Total VCFs scanned:      {total_files}")
    print(f"  Pass (all required):     {files_ok}")
    print(f"  Missing required fields: {files_missing_fields}")
    print(f"  Parse errors:            {files_with_errors}")

    if missing_field_counts:
        print()
        print("MISSING REQUIRED FORMAT FIELDS")
        for field, count in missing_field_counts.most_common():
            print(f"  {field}: missing in {count} file(s)")

    print()
    print("FORMAT FIELDS PRESENT (across all files)")
    for field, count in sorted(format_field_counts.items()):
        marker = ""
        if field in REQUIRED_FORMAT:
            marker = " [REQUIRED]"
        elif field in RECOMMENDED_FORMAT:
            marker = " [RECOMMENDED]"
        print(f"  {field}: {count}/{total_files} files{marker}")

    print()
    print("VARIANT RECORDS")
    print(f"  Total records across all files: {total_variants}")
    print(f"  No-call (./.) genotypes:        {total_nocall}")
    avg = total_variants / max(total_files - files_with_errors, 1)
    print(f"  Avg markers per file:           {avg:.1f}")

    print()
    print("AD FIELD FORMAT (allele depth)")
    total_ad = total_ad_single + total_ad_two + total_ad_multi
    if total_ad > 0:
        print(f"  Single value (ref only):   {total_ad_single} "
              f"({100*total_ad_single/total_ad:.1f}%) — hom-ref, no alt count")
        print(f"  Two values (ref,alt):      {total_ad_two} "
              f"({100*total_ad_two/total_ad:.1f}%) — biallelic, ideal for chimerism")
        print(f"  >2 values (multiallelic):  {total_ad_multi} "
              f"({100*total_ad_multi/total_ad:.1f}%)")
    else:
        print("  No AD data found")

    print()
    print("  NOTE: Single-value AD at hom-ref sites is expected for genotyping")
    print("  samples (pure donor/host). For ADMIXTURE samples, all markers must")
    print("  have two-value AD (ref,alt). If admixture samples show single-value")
    print("  AD, the pipeline needs adjustment (joint-call with donor/host samples).")

    print()
    print("DEPTH (DP)")
    depth_summary = summarise_depths(all_depths)
    if depth_summary["count"] > 0:
        print(f"  Min:    {depth_summary['min']}")
        print(f"  Q25:    {depth_summary['q25']}")
        print(f"  Median: {depth_summary['median']}")
        print(f"  Q75:    {depth_summary['q75']}")
        print(f"  Max:    {depth_summary['max']}")
        print(f"  Mean:   {depth_summary['mean']}")
        print(f"  Below {MIN_DEPTH}x: {total_low_depth} records "
              f"({100*total_low_depth/depth_summary['count']:.1f}%)")
    else:
        print("  No DP data found")

    print()
    print("=" * 70)
    print("ALLOMIX COMPATIBILITY VERDICT")
    print()

    issues = []
    if files_missing_fields > 0:
        issues.append(f"{files_missing_fields} file(s) missing required FORMAT fields")
    if total_ad == 0:
        issues.append("No AD data found in any file")
    if depth_summary.get("count", 0) > 0 and depth_summary["median"] < MIN_DEPTH:
        issues.append(f"Median depth ({depth_summary['median']}x) below {MIN_DEPTH}x threshold")

    if not issues:
        print("  PASS — VCFs appear compatible with allomix.")
        print()
        print("  Next step: verify that ADMIXTURE sample VCFs (not just genotyping")
        print("  samples) have two-value AD at all marker positions. Run this script")
        print("  on a known admixture VCF and check for single-value AD counts.")
    else:
        print("  ISSUES FOUND:")
        for issue in issues:
            print(f"    - {issue}")
        print()
        print("  Review pipeline configuration before using these VCFs with allomix.")

    print()


if __name__ == "__main__":
    main()
```

Save this script as `src/allomix/scripts/audit_vcfs.py` or run it standalone. Example:

```bash
# Audit all per-sample VCFs in the rhAmpSeq SID directory
python audit_vcfs.py /tau/data/clinical_hg38/idt_rhampseq_sid/

# Audit a specific file
python audit_vcfs.py /tau/data/clinical_hg38/idt_rhampseq_sid/HID_SAMPLE.gatk.vcf.gz

# Use a custom glob pattern
python audit_vcfs.py /tau/data/clinical_hg38/idt_rhampseq_sid/ --glob '*.gatk.vcf.gz'
```

---

## Summary

| Question | Answer |
|----------|--------|
| Primary input format? | **VCF** |
| Required FORMAT fields? | GT, AD, DP |
| Recommended FORMAT fields? | GQ, PL, AF |
| BAM support needed? | No (v1) |
| Pipeline changes needed? | Verify joint calling includes admixture + genotyping samples together; verify `--trim-alt-alleles` doesn't drop low-count ALT alleles |
| Risk if VCFs are insufficient? | Low — fallback is `bcftools mpileup` preprocessing, still outputs VCF |
| Blocking concern? | None. Existing VCFs have all required fields. Pipeline already joint-calls across samples. |
