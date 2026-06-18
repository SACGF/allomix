"""Robust (median) estimate of the low-level third-party signal in SRP434573 (issue #12).

The pooled read-weighted means are dominated by a few genotype-miscall / mapping
artifact sites (one site at 40-100% outweighs hundreds at <0.5%). To estimate the
TYPICAL low-level signal, aggregate reads per site across the dilution samples,
take that site's minor fraction, then report the MEDIAN over sites. A genuine
co-pooled-contamination signal should: (a) raise the carrier median above the
no-carrier (pure sequencing error) median, and (b) rise with the number of
co-pooled genomes carrying the allele (dose-response). A site artifact would do
neither.

Reads the committed genotype snapshot in paper/public_data/SRP434573/genotypes
(or a freshly joint-called output/genotypes/SRP434573 if present). Pooled
statistics only.
"""

import csv
from collections import defaultdict
from pathlib import Path
from statistics import median

from cyvcf2 import VCF
from srp434573_common import resolve_srp434573_genotypes_dir

GEN = resolve_srp434573_genotypes_dir()
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
    persite_rows = []  # (n_carriers, n_alleles, minor_frac) per pooled site, for the figure

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
            # n_carriers: each other individual carrying the minor allele counts 1
            # (het or hom). n_alleles: dose-weighted (het=1, hom=2), the number of
            # co-pooled minor alleles, which is the more mechanistically honest
            # x-axis for an index-hopping floor.
            n_carriers = 0
            n_alleles = 0
            for o in others:
                g = gt[o].get(key)
                if g is None:
                    continue
                hom_minor = HOM_ALT if minor_is_alt else HOM_REF
                if g == HET:
                    n_carriers += 1
                    n_alleles += 1
                elif g == hom_minor:
                    n_carriers += 1
                    n_alleles += 2
            marker_info[key] = (minor_is_alt, n_carriers, n_alleles)

        site_acc = defaultdict(lambda: [0, 0])  # key -> [minor, dp] summed over admix samples
        avcf = VCF(str(GEN / f"{name}.admix.vcf.gz"))
        n_samples = len(avcf.samples)
        # Endpoint columns (pure host/donor, 0%/100% host) are single-source
        # samples, not co-pooled dilution mixtures; exclude them so the floor is
        # estimated only from the titration samples.
        skip_idx = {i for i, s in enumerate(avcf.samples) if s in (host, donor)}
        for v in avcf:
            key = (v.CHROM, v.POS, v.REF, v.ALT[0] if v.ALT else None)
            info = marker_info.get(key)
            if info is None:
                continue
            minor_is_alt = info[0]
            ad = v.format("AD")
            if ad is None:
                continue
            for i in range(n_samples):
                if i in skip_idx:
                    continue
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
            _, n_carriers, n_alleles = marker_info[key]
            persite_rows.append((n_carriers, n_alleles, frac))
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

    # Write headline facts for the paper (single-row CSV -> vibepaper namespace
    # ``srp_contam``). These back the contamination dose-response described in
    # the Results section.
    facts_dir = Path("output/facts")
    facts_dir.mkdir(parents=True, exist_ok=True)
    facts = {
        "n_nocarrier_sites": len(nocarrier_fracs),
        "n_carrier_sites": len(carrier_fracs_all),
        "nocarrier_floor_pct": f"{median(nocarrier_fracs) * 100:.3f}",
        "carrier_median_pct": f"{median(carrier_fracs_all) * 100:.2f}",
        "dose_1carrier_pct": f"{median(carrier_fracs_by_dose[1]) * 100:.2f}",
        "dose_5carrier_pct": f"{median(carrier_fracs_by_dose[5]) * 100:.2f}",
    }
    with open(facts_dir / "srp_contam.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(facts.keys())
        w.writerow(facts.values())
    print(f"\nWrote {facts_dir / 'srp_contam.csv'}")

    # Long-format per-site export for the contamination-floor figure (#19). One row
    # per pooled site (dp>=500, same site set as the medians above), so the figure
    # and the cited medians describe the same data. No genomic coordinates needed.
    persite_path = facts_dir / "srp_contam_persite.csv"
    with open(persite_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["n_carriers", "n_alleles", "minor_frac"])
        for n_carriers, n_alleles, frac in persite_rows:
            w.writerow([n_carriers, n_alleles, f"{frac:.8f}"])
    print(f"Wrote {persite_path}  ({len(persite_rows)} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
