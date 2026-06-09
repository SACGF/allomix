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


def admix_samples(vcf: Path) -> list[str]:
    return list(VCF(str(vcf)).samples)


def run_mix(name: str, host: str, donors: list[str]) -> list[dict]:
    """Run allomix monitor (TSV) and return one parsed dict per admix sample."""
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
                known = known_host_pct(sample) if sample else None
                # MLE host fraction = 100 - donor_pct (two-component).
                dpct = fnum(rec.get("donor_pct"))
                dlo, dhi = fnum(rec.get("ci_lo")), fnum(rec.get("ci_hi"))
                mle = (100.0 - dpct) if dpct is not None else None
                mle_lo = (100.0 - dhi) if dhi is not None else None
                mle_hi = (100.0 - dlo) if dlo is not None else None
                # Presence-test host fraction (native orientation, fixed filter).
                pf = fnum(rec.get("host_f_est"))
                plo, phi = fnum(rec.get("host_f_ci_lo")), fnum(rec.get("host_f_ci_hi"))
                two.append(
                    {
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
                )
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
