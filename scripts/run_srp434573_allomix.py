"""Run allomix monitor on all SRP434573 mixtures and compare to ground truth.

Ground truth: admix alias 1_<N>_<X>-<Y> => minor (first, = donor) fraction
1/(1+N). File mix_<DONOR>_into_<HOST> sets host/donor. All 7 individuals are
unrelated. Three-person mix_F2_M1_into_M2 (1:3:5 of F2:M1:M2) => host M2,
donors F2 (1/9) and M1 (3/9).

Writes two files (the three-person sample is kept but split out of the
two-person accuracy series):

  output/srp434573_two_person.tsv    one row per two-person dilution timepoint
  output/srp434573_three_person.tsv  one row per component of the 3-person mix

Reads only output/genotypes/SRP434573; writes nothing to /tau.
"""

import subprocess
import sys
from pathlib import Path

from cyvcf2 import VCF

GEN = Path("output/genotypes/SRP434573")
OUT = Path("output")
ALLOMIX = [".venv/bin/allomix", "monitor"]

# name -> (host, [donors in order]); two-person first, three-person last
MIXES = {
    "mix_F1_into_F3": ("F3", ["F1"]),
    "mix_F2_into_F1": ("F1", ["F2"]),
    "mix_F2_into_M1": ("M1", ["F2"]),
    "mix_F2_into_M2": ("M2", ["F2"]),
    "mix_F3_into_F2": ("F2", ["F3"]),
    "mix_M1_into_M2": ("M2", ["M1"]),
    "mix_M3_into_F1": ("F1", ["M3"]),
    "mix_M3_into_F2": ("F2", ["M3"]),
    "mix_M3_into_F3": ("F3", ["M3"]),
    "mix_M3_into_M4": ("M4", ["M3"]),
    "mix_F2_M1_into_M2": ("M2", ["F2", "M1"]),
}

# ratio token N -> minor (single-donor) percent
RATIO_PCT = {9: 10.0, 19: 5.0, 39: 2.5, 79: 1.25, 99: 1.0, 199: 0.5}

# three-person 1:3:5 of F2:M1:M2 -> known component percents
THREE_PERSON_KNOWN = {"F2": 100.0 / 9, "M1": 300.0 / 9, "M2": 500.0 / 9}


def known_donor_pct(sample: str) -> float | None:
    """Expected minor/donor percent from a two-person admix alias."""
    parts = sample.split("_")
    try:
        return RATIO_PCT.get(int(parts[1]))
    except (IndexError, ValueError):
        return None


def admix_samples(vcf: Path) -> list[str]:
    return list(VCF(str(vcf)).samples)


def run_mix(name: str, host: str, donors: list[str]) -> list[dict]:
    """Run allomix monitor; return one parsed TSV dict per admix sample."""
    panel = GEN / f"{name}.SRP434573.vcf.gz"
    admix = GEN / f"{name}.admix.vcf.gz"
    cmd = [
        *ALLOMIX,
        "--panel-vcf", str(panel),
        "--admix-vcf", str(admix),
        "--host-sample", host,
        "--format", "tsv",
    ]
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


def main() -> int:
    two, three = [], []
    for name, (host, donors) in MIXES.items():
        recs = run_mix(name, host, donors)
        if len(donors) == 1:
            donor = donors[0]
            for rec in recs:
                sample = rec.get("sample")
                known = known_donor_pct(sample) if sample else None
                est = fnum(rec.get("donor_pct"))
                two.append(
                    {
                        "mixture": name,
                        "sample": sample,
                        "donor": donor,
                        "host": host,
                        "known_pct": known,
                        "est_pct": est,
                        "error": (est - known) if (est is not None and known is not None) else None,
                        "ci_lo": fnum(rec.get("ci_lo")),
                        "ci_hi": fnum(rec.get("ci_hi")),
                        "n_used": rec.get("n_used"),
                        "mean_depth": rec.get("mean_depth"),
                        "gof_pval": rec.get("gof_pval"),
                        "qc": rec.get("qc_status"),
                    }
                )
        else:  # three-person: donor1 = donors[0], donor2 = donors[1], host
            rec = recs[0]
            comps = [
                (donors[0], "donor", fnum(rec.get("donor1_pct")),
                 fnum(rec.get("donor1_ci_lo")), fnum(rec.get("donor1_ci_hi"))),
                (donors[1], "donor", fnum(rec.get("donor2_pct")),
                 fnum(rec.get("donor2_ci_lo")), fnum(rec.get("donor2_ci_hi"))),
                (host, "host", fnum(rec.get("host_pct")), None, None),
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
                        "error": (est - known) if (est is not None and known is not None) else None,
                        "ci_lo": lo,
                        "ci_hi": hi,
                        "qc": rec.get("qc_status"),
                    }
                )

    two_cols = ["mixture", "sample", "donor", "host", "known_pct", "est_pct",
                "error", "ci_lo", "ci_hi", "n_used", "mean_depth", "gof_pval", "qc"]
    three_cols = ["mixture", "sample", "component", "role", "known_pct",
                  "est_pct", "error", "ci_lo", "ci_hi", "qc"]
    write_tsv(OUT / "srp434573_two_person.tsv", two_cols, two)
    write_tsv(OUT / "srp434573_three_person.tsv", three_cols, three)
    sys.stderr.write(
        f"Wrote {len(two)} two-person rows and {len(three)} three-person rows.\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
