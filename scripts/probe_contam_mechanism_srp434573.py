"""Work out the mechanism behind the ~2% third-party signal in SRP434573 (issue #12).

Builds on probe_index_hopping_srp434573.py. That probe showed that at markers where
both host and donor are homozygous for the same allele, the minor allele sits at
~2% wherever one of the other panel individuals carries it, vs ~0.02% (the
sequencing-error floor) where none do. This script tries to distinguish the
candidate explanations using only genotypes + AD we already have:

  H1  Genotyping miscall: host/donor are actually het at some "hom" sites, so a
      few sites carry ~45% minor reads and the pooled 2% is just those.
      -> Test A: per-site distribution. If H1, the excess concentrates in a few
         high-fraction sites; if real low-level contamination, it is spread
         uniformly at ~2% across many sites.

  H2  Reference-mapping artifact at polymorphic loci (paralog cross-mapping):
      a fixed property of the site, present even without a pool.
      -> Test C: the pure single-source HOST and DONOR reference runs are hom at
         these sites. If H2, they show the same ~2% minor fraction. If they are
         clean (~error floor), the 2% is not a site/mapping property.

  H3  Pooled cross-contamination (index hopping at the sequencer, or physical
      library/pooling carryover): reads from co-loaded genomes that carry the
      allele.
      -> Test B: dose-response. Bin carrier sites by how many of the other panel
         individuals carry the minor allele. Contamination scales with co-pooled
         allele dose; a site artifact does not.
      -> Test C also separates index hopping from admix-only physical contam: a
         flowcell-wide hop hits the pure references too (they show ~2%); contam
         introduced only during admix dilution prep leaves references clean.

Reads only output/genotypes/SRP434573. Emits pooled statistics, no coordinates.
"""

import sys
from collections import defaultdict
from pathlib import Path

from cyvcf2 import VCF

GEN = Path("output/genotypes/SRP434573")

MIXES = {
    "mix_F1_into_F3": ("F1", "F3"),
    "mix_F2_into_F1": ("F2", "F1"),
    "mix_F2_into_M1": ("F2", "M1"),
    "mix_F2_into_M2": ("F2", "M2"),
    "mix_F3_into_F2": ("F3", "F2"),
    "mix_M1_into_M2": ("M1", "M2"),
    "mix_M3_into_F1": ("M3", "F1"),
    "mix_M3_into_F2": ("M3", "F2"),
    "mix_M3_into_F3": ("M3", "F3"),
    "mix_M3_into_M4": ("M3", "M4"),
}
HOM_REF, HET, UNKNOWN, HOM_ALT = 0, 1, 2, 3


def build_genotype_matrix() -> dict[str, dict[tuple, int]]:
    gt: dict[str, dict[tuple, int]] = defaultdict(dict)
    conflicts: dict[str, set] = defaultdict(set)
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
                prev = gt[s].get(key)
                if prev is not None and prev != t:
                    conflicts[s].add(key)
                    continue
                gt[s][key] = t
    for s, keys in conflicts.items():
        for k in keys:
            gt[s].pop(k, None)
    return gt


def minor_dose(gt_type: int, minor_is_alt: bool) -> int:
    """Copies of the minor allele (0/1/2) this genotype carries."""
    if minor_is_alt:
        return {HOM_REF: 0, HET: 1, HOM_ALT: 2}.get(gt_type, 0)
    return {HOM_REF: 2, HET: 1, HOM_ALT: 0}.get(gt_type, 0)


def main() -> int:
    gt = build_genotype_matrix()
    all_indiv = sorted(gt)
    sys.stderr.write(f"Individuals: {', '.join(all_indiv)}\n\n")

    # Accumulators pooled across all 10 two-person mixtures.
    # Test A: per-site admix minor fractions (one entry per (mixture, site, sample)).
    persite_admix: list[float] = []
    excess_reads_lowfrac = 0   # reads from sites with per-site frac <= 10%
    excess_reads_highfrac = 0  # reads from sites with per-site frac > 10% (miscall-like)
    # Test B: dose-response, binned by number of OTHER panel carriers (1..5).
    dose_bin_reads: dict[int, int] = defaultdict(int)
    dose_bin_depth: dict[int, int] = defaultdict(int)
    # Test C: pure reference (host & donor own runs) at carrier vs no-carrier sites.
    ref_carrier_reads = ref_carrier_depth = 0
    ref_nocarrier_reads = ref_nocarrier_depth = 0
    admix_carrier_reads = admix_carrier_depth = 0
    admix_nocarrier_reads = admix_nocarrier_depth = 0

    for name, (host, donor) in MIXES.items():
        others = [s for s in all_indiv if s not in (host, donor)]

        # Resolve per-marker status from genotypes.
        marker_info: dict[tuple, tuple[bool, int]] = {}  # key -> (minor_is_alt, n_other_carriers)
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
                if key in gt[o] and minor_dose(gt[o][key], minor_is_alt) > 0
            )
            marker_info[key] = (minor_is_alt, n_carriers)

        # --- Reference AD (Test C): host & donor own runs from the panel VCF ---
        pvcf = VCF(str(GEN / f"{name}.SRP434573.vcf.gz"))
        psamples = pvcf.samples
        hi = psamples.index(host)
        di = psamples.index(donor)
        for v in pvcf:
            key = (v.CHROM, v.POS, v.REF, v.ALT[0] if v.ALT else None)
            info = marker_info.get(key)
            if info is None:
                continue
            minor_is_alt, n_carriers = info
            ad = v.format("AD")
            if ad is None:
                continue
            for idx in (hi, di):
                ref_n, alt_n = int(ad[idx][0]), int(ad[idx][1])
                dp = ref_n + alt_n
                if dp <= 0:
                    continue
                minor_reads = alt_n if minor_is_alt else ref_n
                if n_carriers > 0:
                    ref_carrier_reads += minor_reads
                    ref_carrier_depth += dp
                else:
                    ref_nocarrier_reads += minor_reads
                    ref_nocarrier_depth += dp

        # --- Admix AD (Tests A, B, C) ---
        avcf = VCF(str(GEN / f"{name}.admix.vcf.gz"))
        asamples = avcf.samples
        for v in avcf:
            key = (v.CHROM, v.POS, v.REF, v.ALT[0] if v.ALT else None)
            info = marker_info.get(key)
            if info is None:
                continue
            minor_is_alt, n_carriers = info
            ad = v.format("AD")
            if ad is None:
                continue
            for i in range(len(asamples)):
                ref_n, alt_n = int(ad[i][0]), int(ad[i][1])
                dp = ref_n + alt_n
                if dp <= 0:
                    continue
                minor_reads = alt_n if minor_is_alt else ref_n
                if n_carriers > 0:
                    admix_carrier_reads += minor_reads
                    admix_carrier_depth += dp
                    frac = minor_reads / dp
                    persite_admix.append(frac)
                    if frac > 0.10:
                        excess_reads_highfrac += minor_reads
                    else:
                        excess_reads_lowfrac += minor_reads
                    b = min(n_carriers, 5)
                    dose_bin_reads[b] += minor_reads
                    dose_bin_depth[b] += dp
                else:
                    admix_nocarrier_reads += minor_reads
                    admix_nocarrier_depth += dp

    def pct(a: int, b: int) -> float:
        return a / b * 100 if b else float("nan")

    print("=" * 70)
    print("TEST A — is the 2% a few miscalled sites or spread uniformly?")
    print("=" * 70)
    n = len(persite_admix)
    persite_admix.sort()
    hi_sites = sum(1 for f in persite_admix if f > 0.10)
    print(f"  carrier site-observations: {n}")
    print(f"  median per-site minor frac: {persite_admix[n // 2] * 100:.3f}%")
    print(f"  90th pct: {persite_admix[int(n * 0.9)] * 100:.3f}%   "
          f"99th pct: {persite_admix[int(n * 0.99)] * 100:.3f}%")
    print(f"  site-observations with frac > 10% (miscall-like): "
          f"{hi_sites} ({pct(hi_sites, n):.2f}% of obs)")
    print(f"  share of carrier minor reads from >10% sites: "
          f"{pct(excess_reads_highfrac, excess_reads_highfrac + excess_reads_lowfrac):.1f}%")
    print()

    print("=" * 70)
    print("TEST B — dose-response: minor frac vs # other panel carriers")
    print("=" * 70)
    print(f"  {'#carriers':>10} {'minor_frac%':>12} {'depth':>12}")
    for b in sorted(dose_bin_reads):
        print(f"  {b:>10} {pct(dose_bin_reads[b], dose_bin_depth[b]):>12.4f} "
              f"{dose_bin_depth[b]:>12}")
    print()

    print("=" * 70)
    print("TEST C — pure HOST/DONOR reference runs vs admix, at the same sites")
    print("=" * 70)
    print(f"  {'sample set':>18} {'carrier%':>10} {'no-carrier%':>12}")
    print(f"  {'pure references':>18} {pct(ref_carrier_reads, ref_carrier_depth):>10.4f} "
          f"{pct(ref_nocarrier_reads, ref_nocarrier_depth):>12.4f}")
    print(f"  {'admix (pooled)':>18} {pct(admix_carrier_reads, admix_carrier_depth):>10.4f} "
          f"{pct(admix_nocarrier_reads, admix_nocarrier_depth):>12.4f}")
    print()
    print(f"  reference carrier depth: {ref_carrier_depth}, "
          f"no-carrier depth: {ref_nocarrier_depth}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
