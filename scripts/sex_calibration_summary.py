#!/usr/bin/env python3
"""De-identified sex-inference calibration summary from joint-called VCFs.

Calibration data for the sex-inference rules in ``plans/sex_markers.md``. Runs
on the internal data (where the VCFs and BAMs live) and prints only aggregate,
de-identified numbers, so the output can leave the share:

- samples are renamed ``s001``, ``s002``, ... (the mapping is written locally
  with ``--mapping-out`` and never printed);
- no genomic coordinates are printed; chrX sites are referred to by rank.

Two tables go to stdout, separated by a blank line:

1. ``per_sample``: one row per sample with declared sex (if a ``--sex-tsv`` is
   given), counts of called / heterozygous non-PAR chrX sites, chrX and
   autosomal het rates, median chrX and autosomal depth from the VCF, and, when
   ``--bam-list`` and ``--sex-bed`` are given, mean read depth at each named
   sex-typing region (SRY, ZFY, AMELY, chrX_AMELX, ...) plus the relative chrY
   depth and the imputed sex using the thresholds of the lab's amplicon
   pipeline (``combine_target_region_mean_depth_for_amplicon_and_impute_sex``).
2. ``per_x_site``: one row per non-PAR chrX site (by rank, no position) with how
   many samples were called and how many were het, split by declared sex. A
   site that is het in several declared males is a paralog or mapping artifact
   and should be excluded from inference.

Depth at the sex-typing regions comes from ``samtools depth`` over the BED
rows (any BED; chrX/chrY rows are the sex regions, the rest give the
autosomal reference depth, as in the lab pipeline). Needs ``samtools`` on PATH
for that part only.

Usage:
    python scripts/sex_calibration_summary.py \\
        --vcf joint_called_a.vcf.gz [--vcf joint_called_b.vcf.gz ...] \\
        [--sex-tsv sample_sex.tsv] \\
        [--bam-list sample_bam.tsv --sex-bed idt_rhampseq_sid_SNPsQC.bed] \\
        [--mapping-out sample_mapping.local.tsv] \\
        > sex_calibration_summary.tsv

``sample_sex.tsv`` and ``sample_bam.tsv`` are two-column, tab-separated:
``sample_name<TAB>F|M`` and ``sample_name<TAB>/path/to.bam``. Sample names must
match the VCF header.
"""

import argparse
import csv
import statistics
import subprocess
import sys
from pathlib import Path

from cyvcf2 import VCF

# Pseudoautosomal regions, 1-based inclusive, exact union of GRCh37 and GRCh38
# so the mask is conservative whichever build the VCF is on (PAR1 merges, the
# two builds' PAR2 intervals do not overlap and stay separate). Same values as
# allomix.contigs. Source: NCBI assembly_regions.txt via the bioutils fork (see
# plans/sex_markers.md, Phase 0).
PAR_X = [(10001, 2781479), (154931044, 155260560), (155701383, 156030895)]
PAR_Y = [(10001, 2781479), (56887903, 57217415), (59034050, 59363566)]

MIN_DP_DEFAULT = 20
MIN_GQ_DEFAULT = 20
# Lab amplicon-pipeline thresholds on max(chrY depth) / median(non-sex depth).
CHRY_MALE = 0.5
CHRY_BORDERLINE_MALE = 0.15
CHRY_BORDERLINE_FEMALE = 0.05


def _norm_chrom(chrom: str) -> str:
    c = chrom[3:] if chrom.lower().startswith("chr") else chrom
    return c.upper()


def _in_par(chrom: str, pos: int) -> bool:
    c = _norm_chrom(chrom)
    regions = PAR_X if c == "X" else PAR_Y if c == "Y" else []
    return any(s <= pos <= e for s, e in regions)


def _read_two_col(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    with open(path, encoding="utf-8") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) >= 2 and row[0] and not row[0].startswith("#"):
                out[row[0]] = row[1]
    return out


def _norm_sex(value: str | None) -> str:
    if value is None:
        return "NA"
    v = value.strip().lower()
    if v in ("f", "female"):
        return "F"
    if v in ("m", "male"):
        return "M"
    return "NA"


def _median(values: list[float]) -> str:
    return f"{statistics.median(values):.0f}" if values else "NA"


def _rate(num: int, den: int) -> str:
    return f"{num / den:.3f}" if den else "NA"


def collect_vcf_stats(vcf_paths: list[Path], min_dp: int, min_gq: int) -> tuple[dict, dict]:
    """Per-sample chrX/autosome genotype tallies and per-chrX-site tallies."""
    per_sample: dict[str, dict] = {}
    per_site: dict[tuple[str, int], dict[str, dict[str, int]]] = {}
    seen: set[tuple[str, str, int]] = set()

    for path in vcf_paths:
        vcf = VCF(str(path))
        samples = vcf.samples
        for s in samples:
            per_sample.setdefault(
                s,
                {
                    "x_called": 0,
                    "x_het": 0,
                    "x_par": 0,
                    "x_dp": [],
                    "auto_called": 0,
                    "auto_het": 0,
                    "auto_dp": [],
                    "y_called": 0,
                    "y_dp": [],
                },
            )
        for var in vcf:
            if len(var.ALT) != 1:
                continue
            chrom = _norm_chrom(var.CHROM)
            is_x = chrom == "X"
            is_y = chrom == "Y"
            is_mt = chrom in ("M", "MT")
            if is_mt:
                continue
            is_par = (is_x or is_y) and _in_par(var.CHROM, var.POS)
            gts = var.genotypes
            dps = var.format("DP")
            gqs = var.gt_quals
            for i, s in enumerate(samples):
                # A reference sample can recur across per-patient joint VCFs;
                # count each (sample, site) once.
                if (s, var.CHROM, var.POS) in seen:
                    continue
                seen.add((s, var.CHROM, var.POS))
                st = per_sample[s]
                alleles = gts[i][:-1]
                if len(alleles) != 2 or any(a < 0 for a in alleles):
                    continue
                dp = int(dps[i][0]) if dps is not None and dps[i][0] >= 0 else -1
                gq = float(gqs[i]) if gqs is not None else -1.0
                if dp < min_dp or (gq >= 0 and gq < min_gq):
                    continue
                het = alleles[0] != alleles[1]
                if is_x:
                    if is_par:
                        st["x_par"] += 1
                        continue
                    st["x_called"] += 1
                    st["x_het"] += int(het)
                    st["x_dp"].append(dp)
                    site = per_site.setdefault((var.CHROM, var.POS), {})
                    site.setdefault(s, {"called": 1, "het": int(het)})
                elif is_y:
                    if is_par:
                        continue
                    st["y_called"] += 1
                    st["y_dp"].append(dp)
                else:
                    st["auto_called"] += 1
                    st["auto_het"] += int(het)
                    st["auto_dp"].append(dp)
    return per_sample, per_site


def region_depths(bam: str, bed: Path) -> dict[str, float]:
    """Mean depth per BED region (keyed by the BED name column, or chrom:rank)."""
    regions: list[tuple[str, int, int, str]] = []
    with open(bed, encoding="utf-8") as fh:
        for n, line in enumerate(fh):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            name = f[3].split(";")[0] if len(f) > 3 and f[3] else f"region{n}"
            regions.append((f[0], int(f[1]), int(f[2]), name))
    proc = subprocess.run(
        ["samtools", "depth", "-a", "-b", str(bed), bam],
        check=True,
        capture_output=True,
        text=True,
    )
    depth: dict[tuple[str, int], int] = {}
    for line in proc.stdout.splitlines():
        chrom, pos, d = line.split("\t")
        depth[(chrom, int(pos))] = int(d)
    out: dict[str, float] = {}
    for chrom, start, end, name in regions:
        vals = [depth.get((chrom, p), 0) for p in range(start + 1, end + 1)]
        out[name] = sum(vals) / len(vals) if vals else 0.0
    return out


def impute_sex_from_depth(
    reg: dict[str, float], sex_names: set[str], y_names: set[str]
) -> tuple[str, str]:
    """Lab-pipeline rule: max(chrY) / median(non-sex regions)."""
    non_sex = [v for k, v in reg.items() if k not in sex_names]
    y_vals = [reg[k] for k in y_names if k in reg]
    if not non_sex or not y_vals or statistics.median(non_sex) <= 0:
        return "NA", "NA"
    rel = max(y_vals) / statistics.median(non_sex)
    if rel > CHRY_MALE:
        label = "M"
    elif rel > CHRY_BORDERLINE_MALE:
        label = "M*"
    elif rel > CHRY_BORDERLINE_FEMALE:
        label = "F*"
    else:
        label = "F"
    return f"{rel:.3f}", label


def main() -> None:
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--vcf", action="append", required=True, type=Path, help="Joint-called VCF(s)")
    ap.add_argument("--sex-tsv", type=Path, help="sample<TAB>F|M declared sex")
    ap.add_argument("--bam-list", type=Path, help="sample<TAB>bam path, for region depth")
    ap.add_argument("--sex-bed", type=Path, help="BED whose chrX/chrY rows are the sex regions")
    ap.add_argument("--mapping-out", type=Path, help="Write index<TAB>sample mapping here (local)")
    ap.add_argument("--min-dp", type=int, default=MIN_DP_DEFAULT)
    ap.add_argument("--min-gq", type=int, default=MIN_GQ_DEFAULT)
    args = ap.parse_args()
    if bool(args.bam_list) != bool(args.sex_bed):
        ap.error("--bam-list and --sex-bed go together")

    declared = (
        {k: _norm_sex(v) for k, v in _read_two_col(args.sex_tsv).items()} if args.sex_tsv else {}
    )
    bams = _read_two_col(args.bam_list) if args.bam_list else {}

    per_sample, per_site = collect_vcf_stats(args.vcf, args.min_dp, args.min_gq)
    names = sorted(per_sample)
    index = {s: f"s{i + 1:03d}" for i, s in enumerate(names)}
    if args.mapping_out:
        with open(args.mapping_out, "w", encoding="utf-8") as fh:
            for s in names:
                fh.write(f"{index[s]}\t{s}\n")

    # Region depth (optional). Sex-region names are whatever the BED calls them.
    depth_rows: dict[str, dict[str, float]] = {}
    sex_names: set[str] = set()
    y_names: set[str] = set()
    if args.sex_bed:
        with open(args.sex_bed, encoding="utf-8") as fh:
            for n, line in enumerate(fh):
                if not line.strip() or line.startswith(("#", "track", "browser")):
                    continue
                f = line.rstrip("\n").split("\t")
                c = _norm_chrom(f[0])
                if c in ("X", "Y"):
                    name = f[3].split(";")[0] if len(f) > 3 and f[3] else f"region{n}"
                    sex_names.add(name)
                    if c == "Y":
                        y_names.add(name)
        for s in names:
            if s in bams:
                depth_rows[s] = region_depths(bams[s], args.sex_bed)
    sex_cols = sorted(sex_names)

    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    header = [
        "sample",
        "declared_sex",
        "x_sites_called",
        "x_het",
        "x_het_rate",
        "x_par_sites_skipped",
        "auto_sites_called",
        "auto_het_rate",
        "median_dp_x",
        "median_dp_auto",
        "y_sites_called_in_vcf",
        "median_dp_y_vcf",
    ]
    if args.sex_bed:
        header += [f"depth_{c}" for c in sex_cols] + ["chry_rel_depth", "imputed_sex_depth"]
    w.writerow(["#per_sample"])
    w.writerow(header)
    for s in names:
        st = per_sample[s]
        row = [
            index[s],
            declared.get(s, "NA"),
            st["x_called"],
            st["x_het"],
            _rate(st["x_het"], st["x_called"]),
            st["x_par"],
            st["auto_called"],
            _rate(st["auto_het"], st["auto_called"]),
            _median(st["x_dp"]),
            _median(st["auto_dp"]),
            st["y_called"],
            _median(st["y_dp"]),
        ]
        if args.sex_bed:
            reg = depth_rows.get(s)
            if reg is None:
                row += ["NA"] * (len(sex_cols) + 2)
            else:
                row += [f"{reg.get(c, 0.0):.0f}" for c in sex_cols]
                row += list(impute_sex_from_depth(reg, sex_names, y_names))
        w.writerow(row)

    w.writerow([])
    w.writerow(["#per_x_site"])
    w.writerow(
        ["x_site_rank", "n_called", "n_het", "n_called_M", "n_het_M", "n_called_F", "n_het_F"]
    )
    for rank, key in enumerate(sorted(per_site, key=lambda k: (k[0], k[1])), 1):
        calls = per_site[key]
        n_called = len(calls)
        n_het = sum(v["het"] for v in calls.values())
        m = [v for s, v in calls.items() if declared.get(s) == "M"]
        f = [v for s, v in calls.items() if declared.get(s) == "F"]
        w.writerow(
            [
                rank,
                n_called,
                n_het,
                len(m),
                sum(v["het"] for v in m),
                len(f),
                sum(v["het"] for v in f),
            ]
        )


if __name__ == "__main__":
    main()
