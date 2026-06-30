#!/usr/bin/env python3
"""Write a VCF of donor-homozygous markers carrying host-presence signal.

Diagnostic script (run from `scripts/`, not part of the installed allomix package).

For one admixture sample, this lists the donor-homozygous markers used by the
presence test, annotated with the per-marker host signal, and flags the ones
that are "upregulated": markers carrying significantly more donor-absent reads
than the sample-wide pooled host fraction predicts. Those are the candidates
for host CNV/LOH (a residual disease clone over-represented at some loci) versus
a uniform low-level host fraction.

The output is a genuine VCF (CHROM/POS/REF/ALT), sorted by genomic position so
clustering is visible, and carries INFO fields so it can be annotated by VEP /
intersected with driver-gene panels. It contains coordinates, so it is written
to a LOCAL file only (see CLAUDE.md); nothing here is printed to stdout beyond
de-identified counts.

Per-marker model: the host carries the donor-absent allele at dose h (1 if host
het, 2 if host hom), so under a single host fraction f the expected donor-absent
VAF is e + f*(h/2), with e the per-direction background (error_rate/3 under the
global fallback). A marker is UPREG if the observed donor-absent count y exceeds
that pooled expectation at a Bonferroni-corrected one-sided binomial test.

Usage:
    python scripts/host_presence_markers_vcf.py \
        --vcf-dir output/joint_called \
        --samples-csv-dir pipeline/sample_csvs \
        --panel-suffix .union_sid_haem_vendor_probes.vcf.gz \
        --admix SAMPLE_AAAA \
        --out-dir output/host_presence_markers
"""

import argparse
import csv
import glob
import os
from pathlib import Path

from scipy.stats import binom

from allomix.detect import host_presence_test, select_donor_hom_markers
from allomix.genotype import classify_markers, parse_vcf

ALPHA = 0.05


def index_csvs(csv_dir: str) -> dict[str, dict]:
    """Map each ADMIX sample id to its patient stem, host and donor ids."""
    out: dict[str, dict] = {}
    for path in sorted(glob.glob(os.path.join(csv_dir, "*.csv"))):
        stem = Path(path).stem
        host = donor = None
        admixes: list[str] = []
        with open(path) as f:
            for r in csv.DictReader(f):
                stype = r["sample_type"].strip().upper()
                sid = r["sample_id"].strip()
                if stype == "HOST":
                    host = sid
                elif stype == "DONOR":
                    donor = sid
                elif stype == "ADMIX":
                    admixes.append(sid)
        for a in admixes:
            out[a] = {"stem": stem, "host": host, "donor": donor}
    return out


def _chrom_key(chrom: str) -> tuple[int, int | str]:
    """Sort key giving chr1..chr22, chrX, chrY, chrM natural order."""
    c = chrom[3:] if chrom.lower().startswith("chr") else chrom
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if c.isdigit():
        return (0, int(c))
    return (1, order.get(c.upper(), 99))


def build_markers(
    meta: dict,
    admix_sample: str,
    vcf_dir: str,
    panel_suffix: str,
    min_gq: int,
    min_dp: int,
    error_rate: float,
) -> tuple[list[dict], dict]:
    """Return (per-marker dicts, pooled dict) for one admix sample."""
    panel = os.path.join(vcf_dir, meta["stem"] + panel_suffix)
    admix_vcf = os.path.join(vcf_dir, meta["stem"] + ".admix.vcf.gz")
    host = parse_vcf(panel, sample=meta["host"], min_gq=min_gq, gt_ad_consistency=True)
    donor = parse_vcf(panel, sample=meta["donor"], min_gq=min_gq, gt_ad_consistency=True)
    admix = parse_vcf(admix_vcf, sample=admix_sample, min_dp=0)
    genotypes = classify_markers(
        host, [donor], admix, min_dp=min_dp, min_gq=min_gq, use_sex_chroms=True
    )

    result = host_presence_test(genotypes.informative, error_rate=error_rate)
    f = result.f_host_mle
    e = error_rate / 3.0
    rows = select_donor_hom_markers(genotypes.informative)
    n_tot = len(rows)
    bonf = ALPHA / n_tot if n_tot else ALPHA

    out: list[dict] = []
    for m in rows:
        chrom, pos, ref, alt = m.key
        coef = m.h / 2.0
        vaf = m.y / m.n if m.n else 0.0
        implied = max(0.0, vaf - e) / coef
        # Expected donor-absent VAF under the pooled host fraction, and a
        # one-sided binomial p for observing at least y given that expectation.
        p0 = min(max(e + f * coef, 1e-9), 0.5)
        p_up = float(binom.sf(m.y - 1, m.n, p0)) if m.y > 0 else 1.0
        # Donor-absent allele = REF for alt->ref directions, else ALT.
        host_allele = ref if m.direction == "alt->ref" else alt
        out.append(
            {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "host_allele": host_allele,
                "y": m.y,
                "n": m.n,
                "h": m.h,
                "vaf": vaf,
                "implied": implied,
                "pooled": f,
                "fold": (implied / f) if f > 0 else float("nan"),
                "p_up": p_up,
                "upreg": p_up < bonf,
            }
        )
    out.sort(key=lambda d: (_chrom_key(d["chrom"]), d["pos"]))
    pooled = {"f": f, "p": result.lrt_pval, "n": n_tot, "bonf": bonf}
    return out, pooled


_HEADER = """\
##fileformat=VCFv4.2
##source=allomix host_presence_markers_vcf.py
##INFO=<ID=HOSTY,Number=1,Type=Integer,Description="Donor-absent allele read count">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Admixture depth at marker">
##INFO=<ID=DOSE,Number=1,Type=Integer,Description="Host dose of donor-absent allele (1 het, 2 hom)">
##INFO=<ID=HOSTALLELE,Number=1,Type=String,Description="The donor-absent (host) allele">
##INFO=<ID=RAWVAF,Number=1,Type=Float,Description="Observed donor-absent VAF">
##INFO=<ID=IMPLIEDF,Number=1,Type=Float,Description="Dose-normalised implied host fraction">
##INFO=<ID=POOLEDF,Number=1,Type=Float,Description="Sample-wide pooled MLE host fraction">
##INFO=<ID=FOLD,Number=1,Type=Float,Description="Implied host fraction / pooled host fraction">
##INFO=<ID=PUP,Number=1,Type=Float,Description="One-sided binomial p for excess over pooled">
##INFO=<ID=UPREG,Number=0,Type=Flag,Description="Above pooled host fraction (Bonferroni)">
"""


def write_vcf(rows: list[dict], pooled: dict, sample: str, path: Path, upreg_only: bool) -> int:
    """Write the marker VCF; return the number of records written."""
    selected = [r for r in rows if r["upreg"]] if upreg_only else rows
    with open(path, "w") as fh:
        fh.write(_HEADER)
        fh.write(f"##sample={sample}\n")
        fh.write(
            f"##pooled_host_fraction={pooled['f']:.6f}\t"
            f"presence_p={pooled['p']:.3e}\tn_donor_hom_markers={pooled['n']}\n"
        )
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for r in selected:
            info = (
                f"HOSTY={r['y']};DP={r['n']};DOSE={r['h']};HOSTALLELE={r['host_allele']};"
                f"RAWVAF={r['vaf']:.5f};IMPLIEDF={r['implied']:.5f};"
                f"POOLEDF={r['pooled']:.5f};FOLD={r['fold']:.2f};PUP={r['p_up']:.2e}"
            )
            if r["upreg"]:
                info += ";UPREG"
            filt = "UPREG" if r["upreg"] else "PASS"
            fh.write(f"{r['chrom']}\t{r['pos']}\t.\t{r['ref']}\t{r['alt']}\t.\t{filt}\t{info}\n")
    return len(selected)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf-dir", default="output/joint_called")
    ap.add_argument("--samples-csv-dir", default="pipeline/sample_csvs")
    ap.add_argument("--panel-suffix", default=".union_sid_haem_vendor_probes.vcf.gz")
    ap.add_argument("--admix", nargs="+", required=True, help="ADMIX sample ids")
    ap.add_argument("--min-gq", type=int, default=20)
    ap.add_argument("--min-dp", type=int, default=20)
    ap.add_argument("--error-rate", type=float, default=0.01)
    ap.add_argument("--out-dir", type=Path, default=Path("output/host_presence_markers"))
    ap.add_argument(
        "--all",
        action="store_true",
        help="Write every donor-homozygous marker (flagged), not only UPREG ones.",
    )
    args = ap.parse_args()

    csv_index = index_csvs(args.samples_csv_dir)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    for admix_sample in args.admix:
        meta = csv_index.get(admix_sample)
        if meta is None:
            raise SystemExit(f"ADMIX sample {admix_sample!r} not found")
        rows, pooled = build_markers(
            meta,
            admix_sample,
            args.vcf_dir,
            args.panel_suffix,
            args.min_gq,
            args.min_dp,
            args.error_rate,
        )
        n_up = sum(1 for r in rows if r["upreg"])
        out = args.out_dir / f"{admix_sample}.host_presence.vcf"
        written = write_vcf(rows, pooled, admix_sample, out, upreg_only=not args.all)
        # De-identified summary only (no coordinates).
        print(
            f"{admix_sample}: pooled f={pooled['f'] * 100:.3f}%  "
            f"{pooled['n']} donor-hom markers, {n_up} UPREG "
            f"(Bonferroni p<{pooled['bonf']:.1e})  ->  {written} records  {out}"
        )


if __name__ == "__main__":
    main()
