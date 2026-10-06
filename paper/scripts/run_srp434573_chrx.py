"""chrX routing arm of the SRP434573 real-data validation (#46).

Runs ``allomix detect`` on the five same-sex two-person SRP434573 titrations
(F1 into F3, F3 into F2, F2 into F1, M1 into M2, M3 into M4) twice: with the
default sex-aware contig policy (non-PAR chrX routed in for a sex-matched
pair; hom/hom sites only for the male pairs) and with ``--contig-policy
autosomes_only``. Both arms use the same error table as the headline run and
no Step 30 contamination table, so the two estimates differ only in the chrX
routing. The sex-mismatched titrations are not run here: the gate excludes
chrX for them, so both arms would be identical.

Role mapping follows ``run_srp434573_allomix.py``: the titrated minor
contributor is the HOST, so the reported quantity is the host percent
(``100 - donor_pct``) against the known dilution.

Writes:

  output/srp434573_chrx.tsv          one row per mixture x titration, both arms
  output/facts/srp434573_chrx.csv    the same table as a facts CSV for the paper
"""

import csv
import sys
from pathlib import Path

from run_srp434573_allomix import MIXES, OUT, fnum, known_host_pct, run_mix, write_tsv

FACTS_DIR = Path("output/facts")

# Same-sex two-person mixtures, in the order the paper lists them.
SAME_SEX = [
    "mix_F1_into_F3",
    "mix_F3_into_F2",
    "mix_F2_into_F1",
    "mix_M1_into_M2",
    "mix_M3_into_M4",
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
    "mle_pct_chrx",
    "mle_ci_lo_chrx",
    "mle_ci_hi_chrx",
    "mle_pct_auto",
    "mle_ci_lo_auto",
    "mle_ci_hi_auto",
    "n_informative_chrx",
    "n_informative_auto",
    "n_used_chrx",
    "n_used_auto",
    "n_chrx_used",
    "n_chrx_male_het_dropped",
    "n_informative_sex_chrom_excluded_auto",
    "qc_chrx",
    "qc_auto",
]


def _host_pct(rec: dict) -> tuple[float | None, float | None, float | None]:
    """(estimate, lo, hi) of the host percent from a detect TSV record."""
    dpct, dlo, dhi = (fnum(rec.get(k)) for k in ("donor_pct", "ci_lo", "ci_hi"))
    est = 100.0 - dpct if dpct is not None else None
    lo = 100.0 - dhi if dhi is not None else None
    hi = 100.0 - dlo if dlo is not None else None
    return est, lo, hi


def collect() -> list[dict]:
    rows: list[dict] = []
    for name in SAME_SEX:
        host, (donor,) = MIXES[name]
        on = {
            r["sample"]: r
            for r in run_mix(name, host, [donor], extra_args=["--contig-policy", "sex_aware"])
        }
        off = {
            r["sample"]: r
            for r in run_mix(name, host, [donor], extra_args=["--contig-policy", "autosomes_only"])
        }
        for sample in on:
            if sample not in off:
                sys.stderr.write(f"[{name}] {sample}: missing from the autosomes-only arm\n")
                continue
            a, b = on[sample], off[sample]
            est_a, lo_a, hi_a = _host_pct(a)
            est_b, lo_b, hi_b = _host_pct(b)
            rows.append(
                {
                    "mixture": name,
                    "host": host,
                    "donor": donor,
                    "sample": sample,
                    "known_pct": known_host_pct(sample),
                    "sex_pair": a.get("sex_pair"),
                    "host_sex": a.get("host_sex"),
                    "donor_sex": a.get("donor_sex"),
                    "mle_pct_chrx": est_a,
                    "mle_ci_lo_chrx": lo_a,
                    "mle_ci_hi_chrx": hi_a,
                    "mle_pct_auto": est_b,
                    "mle_ci_lo_auto": lo_b,
                    "mle_ci_hi_auto": hi_b,
                    "n_informative_chrx": a.get("n_informative"),
                    "n_informative_auto": b.get("n_informative"),
                    "n_used_chrx": a.get("n_used"),
                    "n_used_auto": b.get("n_used"),
                    "n_chrx_used": a.get("n_chrx_used"),
                    "n_chrx_male_het_dropped": a.get("n_chrx_male_het_dropped"),
                    "n_informative_sex_chrom_excluded_auto": b.get(
                        "n_informative_sex_chrom_excluded"
                    ),
                    "qc_chrx": a.get("qc_status"),
                    "qc_auto": b.get("qc_status"),
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
    write_tsv(OUT / "srp434573_chrx.tsv", COLS, rows)
    write_facts(FACTS_DIR / "srp434573_chrx.csv", rows)
    sys.stderr.write(
        f"Wrote {len(rows)} same-sex titration rows (chrX on vs off) -> "
        f"{OUT / 'srp434573_chrx.tsv'} and {FACTS_DIR / 'srp434573_chrx.csv'}\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
