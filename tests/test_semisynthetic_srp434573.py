"""Tests for the semi-synthetic SRP434573 mixture helpers (issue #5).

Covers the two pieces that have to stay consistent for the synthetic fractions to
mean what they say:

  - run_srp434573_allomix.parse_synthetic_sample: the ``syn_..._f<pct>_rep<n>``
    name -> (known host %, replicate) parser the paper runner uses to label points.
  - make_semisynthetic_srp434573.read_pair_csv: reading a pair's HOST/DONOR BAMs,
    and the host(minor)=mix_bams-DONOR convention mapping (the easy thing to get
    backwards, which would silently invert every fraction).
"""

import csv
import sys
from pathlib import Path

import pytest

# Make the paper/scripts and scripts modules importable.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "paper" / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))

import make_semisynthetic_srp434573 as gen  # noqa: E402
import run_srp434573_allomix as runner  # noqa: E402


@pytest.mark.parametrize(
    "name,expected",
    [
        ("syn_F2-M1_f0.1_rep3", (0.1, 3)),
        ("syn_F1-F3_f0.5_rep5", (0.5, 5)),
        ("syn_M3-M4_f0.05_rep1", (0.05, 1)),
    ],
)
def test_parse_synthetic_sample_ok(name, expected):
    assert runner.parse_synthetic_sample(name) == expected


@pytest.mark.parametrize(
    "name",
    [
        "1_199_F2-M1_v1",   # real titration alias, not synthetic
        "100_0_F1-F3",      # real endpoint control
        "syn_bad",          # missing fraction/rep tokens
        "syn_F2-M1_fx_rep1",  # non-numeric fraction
        "",
    ],
)
def test_parse_synthetic_sample_rejects(name):
    assert runner.parse_synthetic_sample(name) is None


def test_sample_name_round_trips_through_parser():
    sid = gen.sample_name("F2", "M1", 0.3, 4)
    assert sid == "syn_F2-M1_f0.3_rep4"
    assert runner.parse_synthetic_sample(sid) == (0.3, 4)


def _write_pair_csv(path: Path, rows: list[tuple[str, str, str]]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample_id", "bam_filename", "sample_type"])
        w.writerows(rows)


def test_read_pair_csv_returns_host_minor_and_donor_major(tmp_path):
    csv_path = tmp_path / "mix_F2_into_M1.csv"
    _write_pair_csv(
        csv_path,
        [
            ("F2", "/tau/F2.bam", "HOST"),   # minor (titrated) contributor
            ("M1", "/tau/M1.bam", "DONOR"),  # majority background
            ("1_199_F2-M1_v1", "/tau/admix.bam", "ADMIX"),
        ],
    )
    pair = gen.read_pair_csv(csv_path)
    assert pair is not None
    host, donor = pair
    assert host == ("F2", "/tau/F2.bam")
    assert donor == ("M1", "/tau/M1.bam")


def test_read_pair_csv_skips_three_person(tmp_path):
    csv_path = tmp_path / "mix_F2_M1_into_M2.csv"
    _write_pair_csv(
        csv_path,
        [
            ("F2", "/tau/F2.bam", "HOST"),
            ("M1", "/tau/M1.bam", "DONOR"),
            ("M2", "/tau/M2.bam", "DONOR"),
            ("1_3_5_F2-M1-M2", "/tau/admix.bam", "ADMIX"),
        ],
    )
    assert gen.read_pair_csv(csv_path) is None


def test_build_pair_maps_host_minor_to_mix_bams_donor_arg(tmp_path, capsys):
    """The allomix host (minor) BAM must be passed as mix_bams' DONOR_BAM.

    mix_bams.sh subsamples its DONOR_BAM to the small fraction, so getting this
    mapping wrong would put the *majority* individual at the titrated fraction
    and silently invert the whole series.
    """
    host = ("F2", "/tau/F2.bam")    # minor
    donor = ("M1", "/tau/M1.bam")   # majority
    n = gen.build_pair(
        "mix_F2_into_M1", host, donor,
        fractions_pct=[0.1], reps=1,
        bam_dir=tmp_path / "bam", out_csv_dir=tmp_path / "csv",
        mix_script=Path("scripts/mix_bams.sh"), dry_run=True,
    )
    assert n == 1
    out = capsys.readouterr().out
    # The single dry-run mix line: HOST_BAM arg (majority) then DONOR_BAM (minor).
    mix_line = next(line for line in out.splitlines() if line.strip().startswith("mix:"))
    parts = mix_line.split()
    # parts: ['mix:', <script>, HOST_BAM, DONOR_BAM, frac, out, name, seed]
    assert parts[2] == "/tau/M1.bam"   # mix_bams HOST_BAM = majority (donor)
    assert parts[3] == "/tau/F2.bam"   # mix_bams DONOR_BAM = minor (allomix host)
    assert parts[4] == "0.001"         # 0.1% as a fraction
    assert parts[6] == "syn_F2-M1_f0.1_rep1"
