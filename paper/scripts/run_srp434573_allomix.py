"""Run allomix monitor on all SRP434573 mixtures and compare to ground truth.

Role mapping matches the regenerated CSVs (see paper/public_data/SRP434573):
the minor (titrated) contributor is the HOST (the residual / recurring patient
we monitor), the major (background) contributor is the DONOR. So the quantity we
validate is the HOST fraction, and:

  - MLE estimate of the monitored fraction = 100 - donor_pct (two-component) or
    the reported host_pct (three-person).
  - Presence-test estimate = host_f_est (native orientation now; the strand-bias
    artifact filter auto-skips on this single-strand panel, issue #18, so no
    --no-artifact-filter is needed).

Ground truth: admix alias 1_<N>_<X>-<Y> => minor (= host) fraction 1/(1+N).
Three-person 1_3_5_F2-M1-M2 (1:3:5 of F2:M1:M2) => host F2 (1/9), donors M1
(3/9) and M2 (5/9).

Writes (the three-person sample is kept but split out of the two-person accuracy
series; the semi-synthetic outputs appear only when that snapshot is present):

  output/srp434573_two_person.tsv    one row per two-person dilution timepoint
  output/srp434573_three_person.tsv  one row per component of the real 3-person mix
  output/srp434573_synthetic.tsv     semi-synthetic two-person sub-0.5% series (#5)
  output/srp434573_synthetic_three_person.tsv  semi-synthetic host + 2 donor mix (#5)

Genotype VCFs are read from the committed snapshot in
``paper/public_data/SRP434573/genotypes`` so the paper builds out of the box,
unless a freshly joint-called ``output/genotypes/SRP434573`` is present (full
from-scratch reproduction), which takes precedence. Writes nothing to /tau.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path

from cyvcf2 import VCF
from srp434573_common import (
    resolve_srp434573_genotypes_dir,
    resolve_srp434573_synthetic_dir,
)

GEN = resolve_srp434573_genotypes_dir()
OUT = Path("output")
ALLOMIX = [shutil.which("allomix") or ".venv/bin/allomix", "monitor"]

# Per-marker-type overdispersion (issue #33) is the estimator default, so every
# monitor invocation uses it. Set ALLOMIX_NO_MARKER_TYPE_OVERDISPERSION=1 to pass
# --no-marker-type-overdispersion instead, recovering the legacy shared-rho
# baseline (used to regenerate the pre-#33 ladder numbers for comparison).
NO_MARKER_TYPE_OVERDISPERSION = os.environ.get("ALLOMIX_NO_MARKER_TYPE_OVERDISPERSION") == "1"

# name -> (host = minor, [donors = major(s)]); two-person first, three-person last
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
    "mix_F2_M1_into_M2": ("F2", ["M1", "M2"]),
}

# ratio token N -> minor (= host) percent
RATIO_PCT = {9: 10.0, 19: 5.0, 39: 2.5, 79: 1.25, 99: 1.0, 199: 0.5}

# three-person 1:3:5 of F2:M1:M2 -> known component percents
THREE_PERSON_KNOWN = {"F2": 100.0 / 9, "M1": 300.0 / 9, "M2": 500.0 / 9}


def known_host_pct(sample: str) -> float | None:
    """Expected minor (= host) percent from a two-person admix alias."""
    parts = sample.split("_")
    try:
        return RATIO_PCT.get(int(parts[1]))
    except (IndexError, ValueError):
        return None


def parse_synthetic_sample(sample: str) -> tuple[float, int] | None:
    """Parse ``syn_<minor>-<major>_f<pct>_rep<n>`` into ``(known_pct, seed_rep)``.

    The semi-synthetic admix samples (issue #5) encode their known host (minor)
    fraction and replicate directly in the name, mirroring how the real
    ``1_N_X-Y`` aliases encode theirs. Returns ``None`` for any name that does
    not match (so non-synthetic samples are ignored).
    """
    if not sample.startswith("syn_"):
        return None
    try:
        frac_tok = next(p for p in sample.split("_") if p.startswith("f"))
        rep_tok = next(p for p in sample.split("_") if p.startswith("rep"))
        return float(frac_tok[1:]), int(rep_tok[3:])
    except (StopIteration, ValueError):
        return None


def parse_synthetic3_sample(sample: str) -> tuple[float, list[float], int] | None:
    """Parse ``syn3_<h>-<d1>-<d2>_h<hpct>_d<d1pct>-<d2pct>_rep<n>``.

    The three-person semi-synthetic mixtures (issue #5, host + 2 donors) encode
    the known host percent and both donor percents directly in the name (written
    by ``make_semisynthetic_srp434573.sample_name3``), so the ground truth needs
    no separate split table. Returns ``(host_pct, [donor1_pct, donor2_pct], rep)``
    or ``None`` for any name that does not match.
    """
    if not sample.startswith("syn3_"):
        return None
    parts = sample.split("_")
    try:
        h_tok = next(p for p in parts if p.startswith("h") and p[1:2].isdigit())
        d_tok = next(p for p in parts if p.startswith("d") and "-" in p[1:])
        rep_tok = next(p for p in parts if p.startswith("rep"))
        hpct = float(h_tok[1:])
        d1pct, d2pct = (float(x) for x in d_tok[1:].split("-"))
        return hpct, [d1pct, d2pct], int(rep_tok[3:])
    except (StopIteration, ValueError):
        return None


def admix_samples(vcf: Path) -> list[str]:
    return list(VCF(str(vcf)).samples)


def run_mix(
    name: str,
    host: str,
    donors: list[str],
    panel: Path | None = None,
    admix: Path | None = None,
    error_table: Path | None = None,
) -> list[dict]:
    """Run allomix monitor (TSV) and return one parsed dict per admix sample.

    Args:
        name: Mixture id (e.g. ``mix_F2_into_M1``).
        host: Host (minor) sample id.
        donors: Donor (majority) sample id(s).
        panel: Genotype VCF (defaults to the real ``GEN/<name>.SRP434573.vcf.gz``).
        admix: Admix VCF (defaults to the real ``GEN/<name>.admix.vcf.gz``).
        error_table: Per-patient error table (defaults to ``GEN/<name>.error_table.tsv``
            when present). Pass an explicit path to reuse the real table for the
            semi-synthetic run (the host/donor individuals are unchanged).
    """
    if panel is None:
        panel = GEN / f"{name}.SRP434573.vcf.gz"
    if admix is None:
        admix = GEN / f"{name}.admix.vcf.gz"
    if error_table is None:
        error_table = GEN / f"{name}.error_table.tsv"
    cmd = [
        *ALLOMIX,
        "--panel-vcf", str(panel),
        "--admix-vcf", str(admix),
        "--host-sample", host,
        "--format", "tsv",
    ]
    if NO_MARKER_TYPE_OVERDISPERSION:
        cmd += ["--no-marker-type-overdispersion"]
    # Per-patient empirical error table (issue #23). When the committed snapshot
    # carries one (built TAU-side by the pipeline's phase-1b reference pileup),
    # pass it so the host-presence background is data-derived per site instead of
    # the flat --error-rate default, which over-attributes signal to error at the
    # lowest dilutions. Absent (fresh checkout before the table is generated) the
    # run falls back to the default, matching the previous behaviour.
    if error_table.exists():
        cmd += ["--error-table", str(error_table)]
    for d in donors:
        cmd += ["--donor-sample", d, "--expected-relatedness", "unrelated"]
    for s in admix_samples(admix):
        cmd += ["--sample", s]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        sys.stderr.write(f"[{name}] FAILED rc={res.returncode}\n{res.stderr}\n")
        return []
    rows, header = [], None
    for line in res.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if fields[0] == "sample":
            header = fields
            continue
        if header is not None:
            rows.append(dict(zip(header, fields)))
    return rows


def fnum(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def write_tsv(path: Path, cols: list[str], rows: list[dict]) -> None:
    lines = ["\t".join(cols)]
    for r in rows:
        vals = []
        for c in cols:
            v = r.get(c)
            if isinstance(v, float):
                v = f"{v:.4f}"
            vals.append("" if v is None else str(v))
        lines.append("\t".join(vals))
    path.write_text("\n".join(lines) + "\n")


def two_person_row(name: str, host: str, donor: str, known: float | None, rec: dict) -> dict:
    """Build one two-person result row dict from an allomix monitor record."""
    sample = rec.get("sample")
    # MLE host fraction = 100 - donor_pct (two-component).
    dpct = fnum(rec.get("donor_pct"))
    dlo, dhi = fnum(rec.get("ci_lo")), fnum(rec.get("ci_hi"))
    mle = (100.0 - dpct) if dpct is not None else None
    mle_lo = (100.0 - dhi) if dhi is not None else None
    mle_hi = (100.0 - dlo) if dlo is not None else None
    # Presence-test host fraction (native orientation, fixed filter).
    pf = fnum(rec.get("host_f_est"))
    plo, phi = fnum(rec.get("host_f_ci_lo")), fnum(rec.get("host_f_ci_hi"))
    return {
        "mixture": name,
        "sample": sample,
        "host": host,
        "donor": donor,
        "known_pct": known,
        "mle_pct": mle,
        "mle_ci_lo": mle_lo,
        "mle_ci_hi": mle_hi,
        "presence_pct": pf * 100 if pf is not None else None,
        "presence_ci_lo": plo * 100 if plo is not None else None,
        "presence_ci_hi": phi * 100 if phi is not None else None,
        "presence_p": fnum(rec.get("host_present_p")),
        "presence_markers": rec.get("host_detect_markers"),
        "n_used": rec.get("n_used"),
        "mean_depth": rec.get("mean_depth"),
        "gof_pval": rec.get("gof_pval"),
        "qc": rec.get("qc_status"),
    }


def run_synthetic(syn_dir: Path) -> list[dict]:
    """Run allomix on the committed semi-synthetic admix VCFs (issue #5).

    Reuses each pair's real per-patient error table (the host/donor individuals
    are unchanged). Returns one row per synthetic mixture sample, with the known
    host fraction and replicate parsed from the ``syn_..._f<pct>_rep<n>`` name.
    """
    rows: list[dict] = []
    for name, (host, donors) in MIXES.items():
        if len(donors) != 1:
            continue  # synthetic series is two-person only
        donor = donors[0]
        panel = syn_dir / f"{name}.synthetic.SRP434573.vcf.gz"
        admix = syn_dir / f"{name}.synthetic.admix.vcf.gz"
        if not admix.exists():
            continue
        # Reuse the real per-patient error table (same host/donor individuals).
        error_table = GEN / f"{name}.error_table.tsv"
        recs = run_mix(name, host, donors, panel=panel, admix=admix,
                       error_table=error_table)
        for rec in recs:
            parsed = parse_synthetic_sample(rec.get("sample") or "")
            if parsed is None:
                continue
            known, rep = parsed
            row = two_person_row(name, host, donor, known, rec)
            row["frac_pct"] = known
            row["seed"] = rep
            rows.append(row)
    return rows


def run_synthetic_three(syn_dir: Path) -> list[dict]:
    """Run allomix on the semi-synthetic host + 2 donor mixtures (issue #5).

    The three-person trio (host F2, donors M1 + M2) is titrated with the host at
    the low fraction ladder while the two donors split the background. Returns one
    row per (sample, component) with the known and estimated percent, mirroring the
    real three-person rows but with the ground truth decoded from the ``syn3_``
    name. Reuses the trio's real per-patient error table (same individuals).
    """
    rows: list[dict] = []
    for name, (host, donors) in MIXES.items():
        if len(donors) != 2:
            continue
        panel = syn_dir / f"{name}.synthetic.SRP434573.vcf.gz"
        admix = syn_dir / f"{name}.synthetic.admix.vcf.gz"
        if not admix.exists():
            continue
        error_table = GEN / f"{name}.error_table.tsv"
        recs = run_mix(name, host, donors, panel=panel, admix=admix,
                       error_table=error_table)
        for rec in recs:
            parsed = parse_synthetic3_sample(rec.get("sample") or "")
            if parsed is None:
                continue
            hpct, (d1pct, d2pct), rep = parsed
            comps = [
                (host, "host", hpct, fnum(rec.get("host_pct")), None, None),
                (donors[0], "donor", d1pct, fnum(rec.get("donor1_pct")),
                 fnum(rec.get("donor1_ci_lo")), fnum(rec.get("donor1_ci_hi"))),
                (donors[1], "donor", d2pct, fnum(rec.get("donor2_pct")),
                 fnum(rec.get("donor2_ci_lo")), fnum(rec.get("donor2_ci_hi"))),
            ]
            for indiv, role, known, est, lo, hi in comps:
                rows.append({
                    "mixture": name,
                    "sample": rec.get("sample"),
                    "component": indiv,
                    "role": role,
                    "host_known_pct": hpct,
                    "known_pct": known,
                    "est_pct": est,
                    "ci_lo": lo,
                    "ci_hi": hi,
                    "seed": rep,
                    "qc": rec.get("qc_status"),
                })
    return rows


def main() -> int:
    two, three = [], []
    for name, (host, donors) in MIXES.items():
        recs = run_mix(name, host, donors)
        if len(donors) == 1:
            donor = donors[0]
            for rec in recs:
                sample = rec.get("sample")
                known = known_host_pct(sample) if sample else None
                two.append(two_person_row(name, host, donor, known, rec))
        else:  # three-person: host = F2, donor1 = M1, donor2 = M2
            rec = recs[0]
            comps = [
                (host, "host", fnum(rec.get("host_pct")), None, None),
                (donors[0], "donor", fnum(rec.get("donor1_pct")),
                 fnum(rec.get("donor1_ci_lo")), fnum(rec.get("donor1_ci_hi"))),
                (donors[1], "donor", fnum(rec.get("donor2_pct")),
                 fnum(rec.get("donor2_ci_lo")), fnum(rec.get("donor2_ci_hi"))),
            ]
            for indiv, role, est, lo, hi in comps:
                known = THREE_PERSON_KNOWN.get(indiv)
                three.append(
                    {
                        "mixture": name,
                        "sample": rec.get("sample"),
                        "component": indiv,
                        "role": role,
                        "known_pct": known,
                        "est_pct": est,
                        "ci_lo": lo,
                        "ci_hi": hi,
                        "qc": rec.get("qc_status"),
                    }
                )

    two_cols = ["mixture", "sample", "host", "donor", "known_pct", "mle_pct",
                "mle_ci_lo", "mle_ci_hi", "presence_pct", "presence_ci_lo",
                "presence_ci_hi", "presence_p", "presence_markers", "n_used",
                "mean_depth", "gof_pval", "qc"]
    three_cols = ["mixture", "sample", "component", "role", "known_pct",
                  "est_pct", "ci_lo", "ci_hi", "qc"]
    write_tsv(OUT / "srp434573_two_person.tsv", two_cols, two)
    write_tsv(OUT / "srp434573_three_person.tsv", three_cols, three)
    sys.stderr.write(
        f"Wrote {len(two)} two-person rows and {len(three)} three-person rows.\n"
    )

    # Semi-synthetic sub-0.5% mixtures (issue #5), written to a separate TSV so
    # the real-data accuracy metrics above are untouched. Only present once the
    # committed snapshot exists (TAU-side generation step); otherwise skipped.
    syn_dir = resolve_srp434573_synthetic_dir()
    if syn_dir is not None:
        syn = run_synthetic(syn_dir)
        syn_cols = two_cols + ["frac_pct", "seed"]
        write_tsv(OUT / "srp434573_synthetic.tsv", syn_cols, syn)
        sys.stderr.write(f"Wrote {len(syn)} semi-synthetic two-person rows from {syn_dir}.\n")

        # Three-person host + 2 donor semi-synthetic mixtures (issue #5). Separate
        # TSV so the two-person low-fraction series stays clean; absent when only
        # two-person synthetic VCFs were generated.
        syn3 = run_synthetic_three(syn_dir)
        syn3_cols = ["mixture", "sample", "component", "role", "host_known_pct",
                     "known_pct", "est_pct", "ci_lo", "ci_hi", "seed", "qc"]
        write_tsv(OUT / "srp434573_synthetic_three_person.tsv", syn3_cols, syn3)
        sys.stderr.write(f"Wrote {len(syn3)} semi-synthetic three-person rows.\n")
    else:
        sys.stderr.write("No semi-synthetic snapshot found; skipping synthetic run.\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
