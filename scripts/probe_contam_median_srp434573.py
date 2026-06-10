"""Robust (median) estimate of the low-level third-party signal in SRP434573 (issue #12).

The pooled read-weighted means are dominated by a few genotype-miscall / mapping
artifact sites (one site at 40-100% outweighs hundreds at <0.5%). To estimate the
TYPICAL low-level signal, aggregate reads per site across the dilution samples,
take that site's minor fraction, then report the MEDIAN over sites. A genuine
co-pooled-contamination signal should: (a) raise the carrier median above the
no-carrier (pure sequencing error) median, and (b) rise with the number of
co-pooled genomes carrying the allele (dose-response). A site artifact would do
neither.

Reads only output/genotypes/SRP434573. Pooled statistics only.
"""

from collections import defaultdict
from pathlib import Path
from statistics import median

from cyvcf2 import VCF

GEN = Path("output/genotypes/SRP434573")
HOM_REF, HET, UNKNOWN, HOM_ALT = 0, 1, 2, 3
MIXES = {
    "mix_F1_into_F3": ("F1", "F3"), "mix_F2_into_F1": ("F2", "F1"),
    "mix_F2_into_M1": ("F2", "M1"), "mix_F2_into_M2": ("F2", "M2"),
    "mix_F3_into_F2": ("F3", "F2"), "mix_M1_into_M2": ("M1", "M2"),
    "mix_M3_into_F1": ("M3", "F1"), "mix_M3_into_F2": ("M3", "F2"),
    "mix_M3_into_F3": ("M3", "F3"), "mix_M3_into_M4": ("M3", "M4"),
}


def build_gt():
    gt = defaultdict(dict)
    conflicts = defaultdict(set)
    for panel in sorted(GEN.glob("*.SRP434573.vcf.gz")):
        vcf = VCF(str(panel))
        for v in vcf:
            if len(v.ALT) != 1:
                continue
            key = (v.CHROM, v.POS, v.REF, v.ALT[0])
            for i, s in enumerate(vcf.samples):
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


def main() -> int:
    gt = build_gt()
    all_indiv = sorted(gt)

    # Per site, aggregate admix minor reads/depth across all dilution samples.
    nocarrier_fracs = []
    carrier_fracs_by_dose = defaultdict(list)
    carrier_fracs_all = []

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

        site_acc = defaultdict(lambda: [0, 0])  # key -> [minor, dp] summed over admix samples
        avcf = VCF(str(GEN / f"{name}.admix.vcf.gz"))
        n_samples = len(avcf.samples)
        for v in avcf:
            key = (v.CHROM, v.POS, v.REF, v.ALT[0] if v.ALT else None)
            info = marker_info.get(key)
            if info is None:
                continue
            minor_is_alt, _ = info
            ad = v.format("AD")
            if ad is None:
                continue
            for i in range(n_samples):
                ref_n, alt_n = int(ad[i][0]), int(ad[i][1])
                dp = ref_n + alt_n
                if dp <= 0:
                    continue
                mr = alt_n if minor_is_alt else ref_n
                site_acc[key][0] += mr
                site_acc[key][1] += dp

        for key, (mr, dp) in site_acc.items():
            if dp < 500:  # need enough pooled depth for a stable per-site fraction
                continue
            frac = mr / dp
            n_carriers = marker_info[key][1]
            if n_carriers == 0:
                nocarrier_fracs.append(frac)
            else:
                carrier_fracs_all.append(frac)
                carrier_fracs_by_dose[min(n_carriers, 5)].append(frac)

    def summ(xs):
        xs = sorted(xs)
        n = len(xs)
        if n == 0:
            return "n=0"
        return (f"n={n:>4}  median={median(xs)*100:7.4f}%  "
                f"p75={xs[int(n*0.75)]*100:7.4f}%  p95={xs[int(n*0.95)]*100:7.4f}%")

    print("Per-site minor fraction (reads pooled across dilution samples, dp>=500):")
    print(f"  no-carrier (pure seq error): {summ(nocarrier_fracs)}")
    print(f"  carrier (>=1 co-pooled)    : {summ(carrier_fracs_all)}")
    print()
    print("Dose-response (MEDIAN per-site minor fraction by # co-pooled carriers):")
    print(f"  {'#carriers':>10} {'median%':>10} {'p75%':>9} {'n_sites':>8}")
    for b in sorted(carrier_fracs_by_dose):
        xs = sorted(carrier_fracs_by_dose[b])
        n = len(xs)
        print(f"  {b:>10} {median(xs)*100:>10.4f} {xs[int(n*0.75)]*100:>9.4f} {n:>8}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
