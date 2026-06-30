"""Tests for the semi-synthetic SRP434573 mixture helpers (issue #5).

Covers the pieces that have to stay consistent for the synthetic fractions to mean
what they say:

  - run_srp434573_allomix.parse_synthetic_sample / parse_synthetic3_sample: the
    ``syn_..._f<pct>_rep<n>`` and ``syn3_..._h<pct>_d<d1>-<d2>_rep<n>`` name parsers
    the paper runner uses to label points with their known fractions.
  - make_semisynthetic_srp434573.read_mix_csv: reading a mixture's HOST/DONOR BAMs
    (one or two donors).
  - make_semisynthetic_srp434573.build_mix: the mix_bams.sh command shape, where
    each component carries its explicit target fraction (so there is no minor/major
    argument-order trap to get backwards).
"""

import csv
import sys
from pathlib import Path

import pytest

# Make the paper/scripts modules importable. This file lives in paper/tests/, so
# paper/scripts is one level up.
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
        "1_199_F2-M1_v1",  # real titration alias, not synthetic
        "100_0_F1-F3",  # real endpoint control
        "syn_bad",  # missing fraction/rep tokens
        "syn_F2-M1_fx_rep1",  # non-numeric fraction
        "syn3_F2-M1-M2_h0.5_d49.75-49.75_rep1",  # three-person, not two-person
        "",
    ],
)
def test_parse_synthetic_sample_rejects(name):
    assert runner.parse_synthetic_sample(name) is None


def test_sample_name_round_trips_through_parser():
    sid = gen.sample_name("F2", "M1", 0.3, 4)
    assert sid == "syn_F2-M1_f0.3_rep4"
    assert runner.parse_synthetic_sample(sid) == (0.3, 4)


def test_sample_name3_round_trips_through_parser():
    sid = gen.sample_name3("F2", "M1", "M2", 0.5, 49.75, 49.75, 2)
    assert sid == "syn3_F2-M1-M2_h0.5_d49.75-49.75_rep2"
    assert runner.parse_synthetic3_sample(sid) == (0.5, [49.75, 49.75], 2)


@pytest.mark.parametrize(
    "name",
    [
        "syn_F2-M1_f0.1_rep3",  # two-person, not three
        "syn3_bad",
        "1_3_5_F2-M1-M2",  # real three-person alias
        "",
    ],
)
def test_parse_synthetic3_sample_rejects(name):
    assert runner.parse_synthetic3_sample(name) is None


def _write_csv(path: Path, rows: list[tuple[str, str, str]]) -> None:
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["sample_id", "bam_filename", "sample_type"])
        w.writerows(rows)


def test_read_mix_csv_returns_host_minor_and_one_donor(tmp_path):
    csv_path = tmp_path / "mix_F2_into_M1.csv"
    _write_csv(
        csv_path,
        [
            ("F2", "/tau/F2.bam", "HOST"),  # minor (titrated) contributor
            ("M1", "/tau/M1.bam", "DONOR"),  # majority background
            ("1_199_F2-M1_v1", "/tau/admix.bam", "ADMIX"),
        ],
    )
    mix = gen.read_mix_csv(csv_path)
    assert mix is not None
    host, donors = mix
    assert host == ("F2", "/tau/F2.bam")
    assert donors == [("M1", "/tau/M1.bam")]


def test_read_mix_csv_returns_host_and_two_donors(tmp_path):
    csv_path = tmp_path / "mix_F2_M1_into_M2.csv"
    _write_csv(
        csv_path,
        [
            ("F2", "/tau/F2.bam", "HOST"),
            ("M1", "/tau/M1.bam", "DONOR"),
            ("M2", "/tau/M2.bam", "DONOR"),
            ("1_3_5_F2-M1-M2", "/tau/admix.bam", "ADMIX"),
        ],
    )
    mix = gen.read_mix_csv(csv_path)
    assert mix is not None
    host, donors = mix
    assert host == ("F2", "/tau/F2.bam")
    assert donors == [("M1", "/tau/M1.bam"), ("M2", "/tau/M2.bam")]


def _mix_lines(out: str) -> list[list[str]]:
    return [line.split() for line in out.splitlines() if line.strip().startswith("mix:")]


def test_build_mix_two_person_passes_explicit_component_fractions(tmp_path, capsys):
    """Each component is passed with its explicit target fraction.

    The mix_bams.sh args are OUTPUT_BAM PANEL_BED SAMPLE SEED then (BAM FRAC)
    pairs: donor (background, 1 - t) then host (minor, t). Getting the fraction on
    the wrong BAM would invert the series, so pin the mapping.
    """
    host = ("F2", "/tau/F2.bam")  # minor
    donor = ("M1", "/tau/M1.bam")  # majority
    n = gen.build_mix(
        "mix_F2_into_M1",
        host,
        [donor],
        fractions_pct=[0.1],
        reps=1,
        donor_splits=gen.DEFAULT_DONOR_SPLITS,
        bam_dir=tmp_path / "bam",
        out_csv_dir=tmp_path / "csv",
        mix_script=Path("scripts/mix_bams.sh"),
        panel_bed=Path("panel.bed"),
        dry_run=True,
    )
    assert n == 1
    parts = _mix_lines(capsys.readouterr().out)[0]
    # parts: ['mix:', <script>, OUT, BED, SAMPLE, SEED, donor_bam, 0.999, host_bam, 0.001]
    assert parts[4] == "syn_F2-M1_f0.1_rep1"
    assert parts[6] == "/tau/M1.bam"  # donor (background)
    assert parts[7] == "0.999"  # 1 - 0.001
    assert parts[8] == "/tau/F2.bam"  # host (minor)
    assert parts[9] == "0.001"  # 0.1% as a fraction


def test_build_mix_three_person_emits_one_bam_per_split_with_donor_fractions(tmp_path, capsys):
    host = ("F2", "/tau/F2.bam")
    donors = [("M1", "/tau/M1.bam"), ("M2", "/tau/M2.bam")]
    n = gen.build_mix(
        "mix_F2_M1_into_M2",
        host,
        donors,
        fractions_pct=[0.5],
        reps=1,
        donor_splits={"eq": (1.0, 1.0), "2to1": (2.0, 1.0)},
        bam_dir=tmp_path / "bam",
        out_csv_dir=tmp_path / "csv",
        mix_script=Path("scripts/mix_bams.sh"),
        panel_bed=Path("panel.bed"),
        dry_run=True,
    )
    assert n == 2  # one mix per split at this single fraction
    lines = _mix_lines(capsys.readouterr().out)
    assert len(lines) == 2

    # Equal split: host 0.5%, donors 49.75% each, in component order d1 d2 host.
    # The split label is in the output BAM path (parts[2]); the sample name
    # (parts[4]) carries the resolved donor percents.
    eq = next(p for p in lines if "_eq_" in p[2])
    assert eq[4] == "syn3_F2-M1-M2_h0.5_d49.75-49.75_rep1"
    assert eq[6:] == ["/tau/M1.bam", "0.4975", "/tau/M2.bam", "0.4975", "/tau/F2.bam", "0.005"]
    # The encoded donor percents must parse back to the same component fractions.
    hpct, dpcts, rep = runner.parse_synthetic3_sample(eq[4])
    assert hpct == 0.5 and dpcts == [49.75, 49.75] and rep == 1

    # 2:1 split: donors at 2/3 and 1/3 of the 99.5% background.
    two = next(p for p in lines if "_2to1_" in p[2])
    assert two[6] == "/tau/M1.bam" and two[8] == "/tau/M2.bam"
    assert abs(float(two[7]) - 0.995 * 2 / 3) < 1e-6
    assert abs(float(two[9]) - 0.995 * 1 / 3) < 1e-6
