"""Sex-mismatch cross-check arm of the SRP434573 real-data validation (#48).

Runs ``allomix detect`` on the five sex-mismatched two-person SRP434573
titrations (M3 into F1, M3 into F2, M3 into F3, F2 into M1, F2 into M2) with
the default settings and the same error table as the headline run, and
collects, per titration, the known minor fraction, the autosomal MLE with its
CI, and the independent chrX copy-number estimate the sex-mismatch cross-check
produces (``sexchrom_*`` columns), with the number of chrX markers it used and
the concordance flag. The chrY depth readout cannot run here: the public panel
has no non-PAR chrY content and no midpoint pileup, so the basis is ``chrX-cn``
throughout.

Role mapping follows ``run_srp434573_allomix.py``: the titrated minor
contributor is the HOST, so every fraction is reported as the host percent
(``100 - donor``) against the known dilution. The pure reference samples in
each admix VCF are kept with ``known_pct`` empty, as in the chrX arm.

Writes:

  output/srp434573_sexmismatch.tsv          one row per mixture x titration
  output/facts/srp434573_sexmismatch.csv    the same table as a facts CSV
"""

import csv
import sys
from pathlib import Path

from run_srp434573_allomix import MIXES, OUT, fnum, known_host_pct, run_mix, write_tsv

FACTS_DIR = Path("output/facts")

# Sex-mismatched two-person mixtures, male minor first, then female minor.
MISMATCHED = [
    "mix_M3_into_F1",
    "mix_M3_into_F2",
    "mix_M3_into_F3",
    "mix_F2_into_M1",
    "mix_F2_into_M2",
]

COLS = [
    "mixture",
    "host",
    "donor",
    "sample",
    "known_pct",
    "sex_pair",
    "host_sex",
    "donor_sex",
    "mle_pct",
    "mle_ci_lo",
    "mle_ci_hi",
    "sexchrom_pct",
    "sexchrom_ci_lo",
    "sexchrom_ci_hi",
    "sexchrom_basis",
    "sexchrom_n",
    "sexchrom_concordant",
    "sexchrom_minus_mle_pp",
    "n_informative",
    "qc_status",
]


def _host_pct_from_ci(est: float | None, lo: float | None, hi: float | None):
    """Donor-fraction (0-1 or percent) triple -> host-percent triple (estimate, lo, hi)."""
    h_est = 100.0 - est if est is not None else None
    h_lo = 100.0 - hi if hi is not None else None
    h_hi = 100.0 - lo if lo is not None else None
    return h_est, h_lo, h_hi


def _sexchrom(rec: dict) -> tuple[float | None, float | None, float | None]:
    """Host percent (estimate, lo, hi) from the ``sexchrom_frac`` / ``sexchrom_ci`` cells."""
    frac = fnum(rec.get("sexchrom_frac"))
    ci = rec.get("sexchrom_ci") or ""
    lo = hi = None
    if "," in ci:
        a, b = ci.split(",", 1)
        lo, hi = fnum(a), fnum(b)
    return _host_pct_from_ci(
        100.0 * frac if frac is not None else None,
        100.0 * lo if lo is not None else None,
        100.0 * hi if hi is not None else None,
    )


def collect() -> list[dict]:
    rows: list[dict] = []
    for name in MISMATCHED:
        host, (donor,) = MIXES[name]
        for rec in run_mix(name, host, [donor]):
            mle, mle_lo, mle_hi = _host_pct_from_ci(
                fnum(rec.get("donor_pct")), fnum(rec.get("ci_lo")), fnum(rec.get("ci_hi"))
            )
            sx, sx_lo, sx_hi = _sexchrom(rec)
            diff = sx - mle if sx is not None and mle is not None else None
            rows.append(
                {
                    "mixture": name,
                    "host": host,
                    "donor": donor,
                    "sample": rec["sample"],
                    "known_pct": known_host_pct(rec["sample"]),
                    "sex_pair": rec.get("sex_pair"),
                    "host_sex": rec.get("host_sex"),
                    "donor_sex": rec.get("donor_sex"),
                    "mle_pct": mle,
                    "mle_ci_lo": mle_lo,
                    "mle_ci_hi": mle_hi,
                    "sexchrom_pct": sx,
                    "sexchrom_ci_lo": sx_lo,
                    "sexchrom_ci_hi": sx_hi,
                    "sexchrom_basis": rec.get("sexchrom_basis"),
                    "sexchrom_n": rec.get("sexchrom_n"),
                    "sexchrom_concordant": rec.get("sexchrom_concordant"),
                    "sexchrom_minus_mle_pp": diff,
                    "n_informative": rec.get("n_informative"),
                    "qc_status": rec.get("qc_status"),
                }
            )
    return rows


def write_facts(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS)
        w.writeheader()
        for r in rows:
            w.writerow({c: (f"{v:.4f}" if isinstance(v, float) else v) for c, v in r.items()})


def main() -> int:
    rows = collect()
    write_tsv(OUT / "srp434573_sexmismatch.tsv", COLS, rows)
    write_facts(FACTS_DIR / "srp434573_sexmismatch.csv", rows)
    sys.stderr.write(
        f"Wrote {len(rows)} sex-mismatched titration rows -> "
        f"{OUT / 'srp434573_sexmismatch.tsv'} and {FACTS_DIR / 'srp434573_sexmismatch.csv'}\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
