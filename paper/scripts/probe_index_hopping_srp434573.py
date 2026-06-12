"""Empirical index-hopping / third-party-contamination probe for SRP434573 (issue #12).

SRP434573 was sequenced on an Illumina HiSeq 3000 (patterned flowcell + ExAmp),
the platform class most prone to index hopping, with many barcoded samples pooled
together. If index hopping leaks reads between co-pooled samples, a two-person
admix sample should carry trace alleles from individuals who are neither its host
nor its donor.

Discriminating test (uses only the genotype panels we already have):
  At markers where BOTH host and donor are homozygous for the SAME allele, the
  minor allele cannot come from either contributor. Its reads are sequencing
  error OR a co-pooled third party. Split those markers by whether any of the
  OTHER five panel individuals carries the minor allele:

    - "third-party carrier" sites: some other individual is HET/HOM for the minor
      allele, so cross-sample leakage would deposit it here.
    - "no-carrier" sites: NO panel individual carries the minor allele, so the
      only source is sequencing error.

  Pure sequencing error is allele- and individual-blind, so it gives the SAME
  minor-allele fraction in both buckets. Index hopping (or any pooled-library
  cross-contamination) deposits real reads only where a co-pooled genome carries
  the allele, so it lifts the carrier bucket above the no-carrier baseline. The
  carrier-minus-baseline excess is an estimate of the third-party leakage rate.

Reads only output/genotypes/SRP434573; writes a TSV to output/. Touches no /tau
data and emits no genomic coordinates, only pooled summary statistics.
"""

import sys
from collections import defaultdict
from pathlib import Path

from cyvcf2 import VCF
from srp434573_common import resolve_srp434573_genotypes_dir

GEN = resolve_srp434573_genotypes_dir()
OUT = Path("output/index_hopping_probe.tsv")

# name -> (host = minor, [donors = major(s)]); same mapping as run_srp434573_allomix.
MIXES = {
    "mix_F1_into_F3": ("F1", ["F3"]),
    "mix_F2_into_F1": ("F2", ["F1"]),
    "mix_F2_into_M1": ("F2", ["M1"]),
    "mix_F2_into_M2": ("F2", ["M2"]),
    "mix_F3_into_F2": ("F3", ["F2"]),
    "mix_M1_into_M2": ("M1", ["M2"]),
    "mix_M3_into_F1": ("M3", ["F1"]),
    "mix_M3_into_F2": ("M3", ["F2"]),
    "mix_M3_into_F3": ("M3", ["F3"]),
    "mix_M3_into_M4": ("M3", ["M4"]),
}
RATIO_PCT = {9: 10.0, 19: 5.0, 39: 2.5, 79: 1.25, 99: 1.0, 199: 0.5}

# cyvcf2 gt_types: 0 HOM_REF, 1 HET, 2 UNKNOWN, 3 HOM_ALT.
HOM_REF, HET, UNKNOWN, HOM_ALT = 0, 1, 2, 3


def known_host_pct(sample: str) -> float | None:
    parts = sample.split("_")
    try:
        return RATIO_PCT.get(int(parts[1]))
    except (IndexError, ValueError):
        return None


def build_genotype_matrix() -> dict[tuple, dict[str, int]]:
    """individual -> {marker_key: gt_type}, unioned over all panel VCFs.

    Only biallelic, confidently-genotyped (non-UNKNOWN) sites are stored. Markers
    where an individual appears in several panels with disagreeing calls are
    dropped for that individual (rare; conservative).
    """
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


def carries_minor(gt_type: int, minor_is_alt: bool) -> bool:
    """True if this genotype carries the minor allele (ALT if minor_is_alt else REF)."""
    if minor_is_alt:
        return gt_type in (HET, HOM_ALT)
    return gt_type in (HET, HOM_REF)


def main() -> int:
    gt = build_genotype_matrix()
    all_indiv = sorted(gt)
    sys.stderr.write(
        f"Genotype matrix: {len(all_indiv)} individuals "
        f"({', '.join(all_indiv)})\n"
    )

    rows = []
    for name, (host, donors) in MIXES.items():
        if len(donors) != 1:
            continue
        donor = donors[0]
        others = [s for s in all_indiv if s not in (host, donor)]
        admix_vcf = VCF(str(GEN / f"{name}.admix.vcf.gz"))
        asamples = admix_vcf.samples

        # Pre-resolve, per marker, the "both-hom-same-allele" status and whether a
        # third party carries the minor allele. Done once per mixture.
        marker_info: dict[tuple, tuple[bool, bool]] = {}  # key -> (minor_is_alt, third_carrier)
        for key in gt[host]:
            hg = gt[host][key]
            dg = gt[donor].get(key)
            if dg is None:
                continue
            if hg == HOM_REF and dg == HOM_REF:
                minor_is_alt = True
            elif hg == HOM_ALT and dg == HOM_ALT:
                minor_is_alt = False
            else:
                continue
            third = any(
                carries_minor(gt[o][key], minor_is_alt)
                for o in others
                if key in gt[o]
            )
            marker_info[key] = (minor_is_alt, third)

        # Accumulate pooled minor-allele reads / depth per admix sample per bucket.
        for i, s in enumerate(asamples):
            acc = {
                "carrier": [0, 0],     # [minor_reads, depth]
                "nocarrier": [0, 0],
            }
            n_carrier = n_nocarrier = 0
            admix_vcf2 = VCF(str(GEN / f"{name}.admix.vcf.gz"))
            for v in admix_vcf2:
                key = (v.CHROM, v.POS, v.REF, v.ALT[0] if v.ALT else None)
                info = marker_info.get(key)
                if info is None:
                    continue
                minor_is_alt, third = info
                ad = v.format("AD")
                if ad is None:
                    continue
                ref_n, alt_n = int(ad[i][0]), int(ad[i][1])
                dp = ref_n + alt_n
                if dp <= 0:
                    continue
                minor_reads = alt_n if minor_is_alt else ref_n
                bucket = "carrier" if third else "nocarrier"
                acc[bucket][0] += minor_reads
                acc[bucket][1] += dp
                if third:
                    n_carrier += 1
                else:
                    n_nocarrier += 1

            cr_reads, cr_dp = acc["carrier"]
            nc_reads, nc_dp = acc["nocarrier"]
            cr_frac = cr_reads / cr_dp if cr_dp else None
            nc_frac = nc_reads / nc_dp if nc_dp else None
            excess = (cr_frac - nc_frac) if (cr_frac is not None and nc_frac is not None) else None
            rows.append({
                "mixture": name,
                "sample": s,
                "host": host,
                "donor": donor,
                "known_host_pct": known_host_pct(s),
                "n_carrier_sites": n_carrier,
                "n_nocarrier_sites": n_nocarrier,
                "carrier_minor_frac_pct": cr_frac * 100 if cr_frac is not None else None,
                "nocarrier_minor_frac_pct": nc_frac * 100 if nc_frac is not None else None,
                "excess_pct": excess * 100 if excess is not None else None,
                "carrier_reads": cr_reads,
                "carrier_depth": cr_dp,
                "nocarrier_reads": nc_reads,
                "nocarrier_depth": nc_dp,
            })

    cols = ["mixture", "sample", "host", "donor", "known_host_pct",
            "n_carrier_sites", "n_nocarrier_sites",
            "carrier_minor_frac_pct", "nocarrier_minor_frac_pct", "excess_pct",
            "carrier_reads", "carrier_depth", "nocarrier_reads", "nocarrier_depth"]
    lines = ["\t".join(cols)]
    for r in rows:
        vals = []
        for c in cols:
            v = r.get(c)
            if isinstance(v, float):
                v = f"{v:.4f}"
            vals.append("" if v is None else str(v))
        lines.append("\t".join(vals))
    OUT.write_text("\n".join(lines) + "\n")
    sys.stderr.write(f"Wrote {len(rows)} rows to {OUT}\n")

    # Print a compact summary to stderr: pooled across mixtures by known host %.
    by_pct: dict[float, list[dict]] = defaultdict(list)
    for r in rows:
        if r["known_host_pct"] is not None:
            by_pct[r["known_host_pct"]].append(r)
    sys.stderr.write(
        "\nPooled by known host %% (carrier vs no-carrier minor-allele fraction):\n"
    )
    sys.stderr.write(
        f"{'host%':>6} {'carrier%':>10} {'nocarrier%':>11} {'excess%':>9} "
        f"{'car_dp':>9} {'noc_dp':>9}\n"
    )
    for pct in sorted(by_pct):
        rs = by_pct[pct]
        cr_reads = sum(r["carrier_reads"] for r in rs)
        cr_dp = sum(r["carrier_depth"] for r in rs)
        nc_reads = sum(r["nocarrier_reads"] for r in rs)
        nc_dp = sum(r["nocarrier_depth"] for r in rs)
        cr = cr_reads / cr_dp * 100 if cr_dp else float("nan")
        nc = nc_reads / nc_dp * 100 if nc_dp else float("nan")
        sys.stderr.write(
            f"{pct:>6} {cr:>10.4f} {nc:>11.4f} {cr - nc:>9.4f} "
            f"{cr_dp:>9} {nc_dp:>9}\n"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
