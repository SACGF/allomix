"""Separate genotype/mapping artifact sites from genuine low-level contamination (issue #12).

probe_contam_mechanism showed the pooled "2%" third-party signal is dominated by
~120 sites sitting at 10-100% minor fraction (read-weighted mean inflated by a few
sites), while the median carrier site is ~0.2%. A site at ~100% minor allele cannot
be low-level contamination; it is a wrong genotype call or a mapping artifact at
that locus. This script flags those sites using the PURE host/donor reference runs
(which must be near-zero at a correctly-called hom site) and re-measures the
contamination signal on the clean remainder.

Definitions:
  carrier site   = host & donor both hom for the same allele, and >=1 other panel
                   individual carries the minor allele.
  no-carrier site= same but NO panel individual carries the minor allele (pure
                   sequencing-error control).
  BAD site       = a carrier-or-no-carrier site whose pooled PURE-reference minor
                   fraction exceeds REF_BAD_FRAC at adequate depth: the genotype
                   is mis-called or the locus mismaps, independent of any pool.
  CLEAN site     = reference minor fraction at/below the error floor.

On CLEAN sites only, the admix carrier-vs-no-carrier gap is the genuine low-level
third-party (pooled-contamination / index-hopping) signal, with the genotype and
mapping artifacts removed.

Reads only output/genotypes/SRP434573. Pooled statistics only.
"""

import sys
from collections import defaultdict
from pathlib import Path

from cyvcf2 import VCF

GEN = Path("output/genotypes/SRP434573")
REF_BAD_FRAC = 0.02      # reference minor frac above this => mis-called / artifact site
REF_MIN_DP = 50          # need this much reference depth to judge a site

MIXES = {
    "mix_F1_into_F3": ("F1", "F3"), "mix_F2_into_F1": ("F2", "F1"),
    "mix_F2_into_M1": ("F2", "M1"), "mix_F2_into_M2": ("F2", "M2"),
    "mix_F3_into_F2": ("F3", "F2"), "mix_M1_into_M2": ("M1", "M2"),
    "mix_M3_into_F1": ("M3", "F1"), "mix_M3_into_F2": ("M3", "F2"),
    "mix_M3_into_F3": ("M3", "F3"), "mix_M3_into_M4": ("M3", "M4"),
}
RATIO_PCT = {9: 10.0, 19: 5.0, 39: 2.5, 79: 1.25, 99: 1.0, 199: 0.5}
HOM_REF, HET, UNKNOWN, HOM_ALT = 0, 1, 2, 3


def known_host_pct(sample: str):
    parts = sample.split("_")
    try:
        return RATIO_PCT.get(int(parts[1]))
    except (IndexError, ValueError):
        return None


def build_genotype_matrix():
    gt = defaultdict(dict)
    conflicts = defaultdict(set)
    for panel in sorted(GEN.glob("*.SRP434573.vcf.gz")):
        vcf = VCF(str(panel))
        samples = vcf.samples
        for v in vcf:
            if len(v.ALT) != 1:
                continue
            key = (v.CHROM, v.POS, v.REF, v.ALT[0])
            for i, s in enumerate(samples):
                t = int(v.gt_types[i])
                if t == UNKNOWN:
                    continue
                if key in gt[s] and gt[s][key] != t:
                    conflicts[s].add(key)
                    continue
                gt[s][key] = t
    for s, keys in conflicts.items():
        for k in keys:
            gt[s].pop(k, None)
    return gt


def minor_reads_of(ad_row, minor_is_alt):
    ref_n, alt_n = int(ad_row[0]), int(ad_row[1])
    dp = ref_n + alt_n
    return (alt_n if minor_is_alt else ref_n), dp


def main() -> int:
    gt = build_genotype_matrix()
    all_indiv = sorted(gt)

    # Pooled accumulators, split by site cleanliness.
    acc = {
        ("clean", "carrier"): [0, 0], ("clean", "nocarrier"): [0, 0],
        ("bad", "carrier"): [0, 0], ("bad", "nocarrier"): [0, 0],
        ("unknown", "carrier"): [0, 0], ("unknown", "nocarrier"): [0, 0],
    }
    n_sites = defaultdict(int)  # (clean/bad/unknown, carrier/nocarrier) -> site count
    # clean carrier signal by known host %, and dose-response on clean sites
    clean_by_pct = defaultdict(lambda: [0, 0])      # pct -> [minor, dp] (carrier, clean)
    clean_nocarrier_by_pct = defaultdict(lambda: [0, 0])
    clean_dose = defaultdict(lambda: [0, 0])        # n_carriers -> [minor, dp]

    for name, (host, donor) in MIXES.items():
        others = [s for s in all_indiv if s not in (host, donor)]
        marker_info = {}
        for key, hg in gt[host].items():
            dg = gt[donor].get(key)
            if dg is None:
                continue
            if hg == HOM_REF and dg == HOM_REF:
                minor_is_alt = True
            elif hg == HOM_ALT and dg == HOM_ALT:
                minor_is_alt = False
            else:
                continue
            n_carriers = sum(
                1 for o in others
                if key in gt[o] and (
                    (minor_is_alt and gt[o][key] in (HET, HOM_ALT))
                    or (not minor_is_alt and gt[o][key] in (HET, HOM_REF))
                )
            )
            marker_info[key] = (minor_is_alt, n_carriers)

        # Reference minor fraction per site -> clean/bad/unknown label.
        pvcf = VCF(str(GEN / f"{name}.SRP434573.vcf.gz"))
        ps = pvcf.samples
        hi, di = ps.index(host), ps.index(donor)
        site_label = {}
        for v in pvcf:
            key = (v.CHROM, v.POS, v.REF, v.ALT[0] if v.ALT else None)
            info = marker_info.get(key)
            if info is None:
                continue
            minor_is_alt, _ = info
            ad = v.format("AD")
            if ad is None:
                site_label[key] = "unknown"
                continue
            m = d = 0
            for idx in (hi, di):
                mr, dp = minor_reads_of(ad[idx], minor_is_alt)
                m += mr
                d += dp
            if d < REF_MIN_DP:
                site_label[key] = "unknown"
            elif m / d > REF_BAD_FRAC:
                site_label[key] = "bad"
            else:
                site_label[key] = "clean"

        for key, lab in site_label.items():
            _, n_carriers = marker_info[key]
            n_sites[(lab, "carrier" if n_carriers > 0 else "nocarrier")] += 1

        # Admix AD bucketed by site label.
        avcf = VCF(str(GEN / f"{name}.admix.vcf.gz"))
        asamples = avcf.samples
        for v in avcf:
            key = (v.CHROM, v.POS, v.REF, v.ALT[0] if v.ALT else None)
            info = marker_info.get(key)
            if info is None:
                continue
            minor_is_alt, n_carriers = info
            lab = site_label.get(key, "unknown")
            cbucket = "carrier" if n_carriers > 0 else "nocarrier"
            ad = v.format("AD")
            if ad is None:
                continue
            for i, s in enumerate(asamples):
                mr, dp = minor_reads_of(ad[i], minor_is_alt)
                if dp <= 0:
                    continue
                acc[(lab, cbucket)][0] += mr
                acc[(lab, cbucket)][1] += dp
                if lab == "clean":
                    pctk = known_host_pct(s)
                    if cbucket == "carrier":
                        clean_by_pct[pctk][0] += mr
                        clean_by_pct[pctk][1] += dp
                        b = min(n_carriers, 5)
                        clean_dose[b][0] += mr
                        clean_dose[b][1] += dp
                    else:
                        clean_nocarrier_by_pct[pctk][0] += mr
                        clean_nocarrier_by_pct[pctk][1] += dp

    def pct(a, b):
        return a / b * 100 if b else float("nan")

    print("Site classification (per mixture-site, pooled over 10 mixtures):")
    for lab in ("clean", "bad", "unknown"):
        c = n_sites[(lab, "carrier")]
        nc = n_sites[(lab, "nocarrier")]
        print(f"  {lab:>8}: carrier={c:>5}  no-carrier={nc:>5}")
    print()

    print("Admix minor-allele fraction by site class:")
    print(f"  {'class':>8} {'carrier%':>10} {'no-carrier%':>12} {'car_dp':>11}")
    for lab in ("clean", "bad", "unknown"):
        cr = acc[(lab, "carrier")]
        nc = acc[(lab, "nocarrier")]
        print(f"  {lab:>8} {pct(*cr):>10.4f} {pct(*nc):>12.4f} {cr[1]:>11}")
    print()

    print("CLEAN sites only — genuine low-level signal by known host %:")
    print(f"  {'host%':>6} {'carrier%':>10} {'no-carrier%':>12} {'excess%':>9}")
    for pk in sorted([p for p in clean_by_pct if p is not None]):
        cr = clean_by_pct[pk]
        nc = clean_nocarrier_by_pct[pk]
        print(f"  {pk:>6} {pct(*cr):>10.4f} {pct(*nc):>12.4f} "
              f"{pct(*cr) - pct(*nc):>9.4f}")
    print()

    print("CLEAN sites only — dose-response (minor frac vs # other carriers):")
    print(f"  {'#carriers':>10} {'minor%':>10} {'depth':>12}")
    for b in sorted(clean_dose):
        print(f"  {b:>10} {pct(*clean_dose[b]):>10.4f} {clean_dose[b][1]:>12}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
