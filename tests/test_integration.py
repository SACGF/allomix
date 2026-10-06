"""Integration tests -- full pipeline from joint VCFs through to results.

Uses the joint-called test VCFs (produced by build_joint_vcf) and runs
genotype -> chimerism -> qc -> report, verifying the output.
"""

import json
from pathlib import Path

import pytest

from allomix.cli import main
from allomix.estimate.chimerism import estimate_single_donor_bb
from allomix.genotype import classify_markers, parse_vcf
from allomix.qc.qc import assess_quality
from allomix.simulate import build_joint_vcf, write_joint_vcf

TEST_DATA_DIR = Path(__file__).resolve().parent / "test_data"
JOINT_VCF = TEST_DATA_DIR / "joint_single_donor.vcf"
JOINT_MULTI_VCF = TEST_DATA_DIR / "joint_multi_donor.vcf"
SEXCHROM_DIR = TEST_DATA_DIR / "sexchrom"


def _run_pipeline(
    vcf_path,
    host_sample="HOST",
    donor_sample="DONOR",
    admix_sample="ADMIX_F0.10",
    min_dp=0,
    min_gq=0,
):
    """Run the full genotype -> chimerism -> qc pipeline from a joint VCF."""
    host = parse_vcf(vcf_path, sample=host_sample, min_dp=0, min_gq=0)
    donor = parse_vcf(vcf_path, sample=donor_sample, min_dp=0, min_gq=0)
    admix = parse_vcf(vcf_path, sample=admix_sample, min_dp=0, min_gq=0)

    genotypes = classify_markers(
        host,
        [donor],
        admix,
        min_dp=min_dp,
        min_gq=min_gq,
        pass_only=False,
    )
    genotypes.sample_name = admix_sample

    result = estimate_single_donor_bb(genotypes.informative, error_rate=0.01)
    qc = assess_quality(result, genotypes)
    return result, qc, genotypes


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestFullPipeline:
    """End-to-end: joint VCF -> genotype -> chimerism -> QC.

    Estimation accuracy across the fraction range is covered cheaply in
    test_chimerism.py; this keeps one case as the parse -> classify -> estimate
    wiring check.
    """

    def test_ten_percent(self):
        """f=0.10: estimate should be near 10%."""
        result, _, _ = _run_pipeline(JOINT_VCF, admix_sample="ADMIX_F0.10")
        assert 0.05 < result.donor_fraction < 0.20


class TestCLIIntegration:
    """Test CLI wiring runs without error."""

    def test_monitor_tsv_and_json(self, tmp_path):
        """A single detect run emits both a valid TSV and the JSON envelope."""
        out_tsv = tmp_path / "cli_out.tsv"
        out_json = tmp_path / "cli_out.json"
        rc = main(
            [
                "detect",
                "--genotype-vcf",
                str(JOINT_VCF),
                "--admix-vcf",
                str(JOINT_VCF),
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR",
                "--sample",
                "ADMIX_F0.10",
                "--tsv",
                str(out_tsv),
                "--json",
                str(out_json),
                "--min-dp",
                "0",
                "--min-gq",
                "0",
            ]
        )
        assert rc == 0
        assert "donor_pct" in out_tsv.read_text()
        # monitor --json writes the report envelope; the per-sample analysis is
        # nested under "analysis".
        data = json.loads(out_json.read_text())
        assert "donor_pct" in data["analysis"]

    def test_detect_exclude_sites(self, tmp_path):
        """--exclude-sites drops the BED marker positions before analysis."""
        markers = parse_vcf(JOINT_VCF, sample="HOST", min_dp=0, min_gq=0)
        n_total = len(markers)
        bed = tmp_path / "exclude.bed"
        # Exclude the first two marker positions (BED is 0-based half-open).
        bed.write_text(
            "".join(f"{m.chrom}\t{m.pos - 1}\t{m.pos}\n" for m in markers[:2]),
            encoding="utf-8",
        )
        out = tmp_path / "excl.tsv"
        rc = main(
            [
                "detect",
                "--genotype-vcf",
                str(JOINT_VCF),
                "--admix-vcf",
                str(JOINT_VCF),
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR",
                "--sample",
                "ADMIX_F0.10",
                "--tsv",
                str(out),
                "--min-dp",
                "0",
                "--min-gq",
                "0",
                "--exclude-sites",
                str(bed),
            ]
        )
        assert rc == 0
        rows = out.read_text().splitlines()
        header = rows[0].split("\t")
        values = rows[1].split("\t")
        n_total_col = header.index("n_total_markers")
        assert int(values[n_total_col]) == n_total - 2

    def test_detect_exclude_and_include_mutually_exclusive(self, tmp_path):
        bed = tmp_path / "sites.bed"
        bed.write_text("chr1\t0\t1\n", encoding="utf-8")
        with pytest.raises(SystemExit):
            main(
                [
                    "detect",
                    "--genotype-vcf",
                    str(JOINT_VCF),
                    "--admix-vcf",
                    str(JOINT_VCF),
                    "--host-sample",
                    "HOST",
                    "--donor-sample",
                    "DONOR",
                    "--sample",
                    "ADMIX_F0.10",
                    "--exclude-sites",
                    str(bed),
                    "--include-sites",
                    str(bed),
                ]
            )

    def test_monitor_estimate_bias_inline(self, tmp_path):
        """--estimate-bias runs bias estimation inline (issue #11), no table file.

        The JSON report records that inline bias estimation was active, so assert
        that provenance rather than just a non-zero exit.
        """
        out = tmp_path / "cli_bias.json"
        rc = main(
            [
                "detect",
                "--genotype-vcf",
                str(JOINT_VCF),
                "--admix-vcf",
                str(JOINT_VCF),
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR",
                "--sample",
                "ADMIX_F0.10",
                "--estimate-bias",
                "--json",
                str(out),
                "--min-dp",
                "0",
                "--min-gq",
                "0",
            ]
        )
        assert rc == 0
        report = json.loads(out.read_text())
        assert report["params"]["estimate_bias"] is True
        assert "--estimate-bias" in report["params"]["command"]
        assert "donor_pct" in report["analysis"]

    def test_estimate_bias_both_het_table_builder(self, tmp_path):
        """estimate-bias --both-het builds a pooled table from admix VCFs (issue #11)."""
        out = tmp_path / "both_het_bias.tsv"
        rc = main(
            [
                "estimate-bias",
                str(JOINT_VCF),
                "--both-het",
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR",
                "--admix-vcfs",
                str(JOINT_VCF),
                "--output",
                str(out),
            ]
        )
        assert rc == 0
        # Header plus at least one both-het marker row (a leading # allomixCaller
        # comment may precede the column header, issue #42).
        lines = [ln for ln in out.read_text().splitlines() if not ln.startswith("#")]
        assert lines[0].split("\t") == ["chrom", "pos", "ref", "alt", "bias", "n_het"]

    def test_estimate_bias_both_het_requires_inputs(self):
        """--both-het without the required genotype/admix inputs is an error."""
        with pytest.raises(SystemExit):
            main(["estimate-bias", str(JOINT_VCF), "--both-het"])

    def test_estimate_bias_samples_file(self, tmp_path):
        """--samples-file reads names (skipping blanks/# comments), merged with --sample."""
        pass_list = tmp_path / "pass.txt"
        pass_list.write_text("# cohort pass-list\nHOST\n\n")
        out = tmp_path / "bias.tsv"
        rc = main(
            [
                "estimate-bias",
                str(JOINT_VCF),
                "--samples-file",
                str(pass_list),
                "--sample",
                "DONOR",
                "--output",
                str(out),
            ]
        )
        assert rc == 0
        # File contributed HOST, --sample added DONOR: a joint-mode run over both.
        # A leading # allomixCaller comment may precede the header (issue #42).
        lines = [ln for ln in out.read_text().splitlines() if not ln.startswith("#")]
        assert lines[0].split("\t") == [
            "chrom",
            "pos",
            "ref",
            "alt",
            "bias",
            "n_het",
        ]

    def test_estimate_bias_samples_file_missing(self):
        """A missing --samples-file is a clean error, not a traceback."""
        with pytest.raises(SystemExit):
            main(["estimate-bias", str(JOINT_VCF), "--samples-file", "/no/such/file.txt"])

    def test_monitor_estimate_bias_conflicts_with_table(self, tmp_path):
        """--estimate-bias and --bias-table together is an error."""
        table = tmp_path / "bias.tsv"
        table.write_text("chrom\tpos\tref\talt\tbias\tn_het\n")
        with pytest.raises(SystemExit):
            main(
                [
                    "detect",
                    "--genotype-vcf",
                    str(JOINT_VCF),
                    "--admix-vcf",
                    str(JOINT_VCF),
                    "--host-sample",
                    "HOST",
                    "--donor-sample",
                    "DONOR",
                    "--sample",
                    "ADMIX_F0.10",
                    "--estimate-bias",
                    "--bias-table",
                    str(table),
                    "--min-dp",
                    "0",
                    "--min-gq",
                    "0",
                ]
            )

    def test_monitor_json_records_command(self, tmp_path):
        """The report JSON records the analysis invocation, minus output flags."""
        out = tmp_path / "cmd.json"
        rc = main(
            [
                "detect",
                "--genotype-vcf",
                str(JOINT_VCF),
                "--admix-vcf",
                str(JOINT_VCF),
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR",
                "--sample",
                "ADMIX_F0.10",
                "--no-host-presence",
                "--json",
                str(out),
            ]
        )
        assert rc == 0
        command = json.loads(out.read_text())["params"]["command"]
        assert command.startswith("allomix detect ")
        assert "--no-host-presence" in command  # analysis flag kept
        assert "--host-sample HOST" in command
        assert "--json" not in command  # output flag stripped (keeps it reproducible)

    def test_timeline(self, tmp_path):
        out = tmp_path / "timeline.json"
        rc = main(
            [
                "timeline",
                "--genotype-vcf",
                str(JOINT_VCF),
                "--admix-vcf",
                str(JOINT_VCF),
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR",
                "--sample",
                "ADMIX_F0.00",
                "--sample",
                "ADMIX_F0.10",
                "--json",
                str(out),
                "--min-dp",
                "0",
                "--min-gq",
                "0",
            ]
        )
        assert rc == 0
        data = json.loads(out.read_text())
        assert len(data["timepoints"]) == 2

    def test_monitor_multiple_samples(self, tmp_path):
        """Monitor with multiple admixture samples."""
        out = tmp_path / "multi.tsv"
        rc = main(
            [
                "detect",
                "--genotype-vcf",
                str(JOINT_VCF),
                "--admix-vcf",
                str(JOINT_VCF),
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR",
                "--sample",
                "ADMIX_F0.00",
                "--sample",
                "ADMIX_F0.10",
                "--sample",
                "ADMIX_F0.50",
                "--tsv",
                str(out),
                "--min-dp",
                "0",
                "--min-gq",
                "0",
            ]
        )
        assert rc == 0
        content = out.read_text()
        for sample in ("ADMIX_F0.00", "ADMIX_F0.10", "ADMIX_F0.50"):
            assert sample in content, f"{sample} missing from multi-sample TSV"

    def test_detect_admix_samples_split_across_vcfs(self, tmp_path):
        """--admix-vcf is repeatable; each --sample resolves to its containing VCF."""
        from cyvcf2 import VCF, Writer

        def subset(sample: str, dest):
            src = VCF(str(JOINT_VCF))
            src.set_samples([sample])
            w = Writer(str(dest), src)
            for rec in src:
                w.write_record(rec)
            w.close()
            src.close()
            return dest

        vcf_a = subset("ADMIX_F0.10", tmp_path / "a.vcf")
        vcf_b = subset("ADMIX_F0.50", tmp_path / "b.vcf")
        out = tmp_path / "split.tsv"
        rc = main(
            [
                "detect",
                "--genotype-vcf",
                str(JOINT_VCF),
                "--admix-vcf",
                str(vcf_a),
                "--admix-vcf",
                str(vcf_b),
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR",
                "--sample",
                "ADMIX_F0.10",
                "--sample",
                "ADMIX_F0.50",
                "--tsv",
                str(out),
                "--min-dp",
                "0",
                "--min-gq",
                "0",
            ]
        )
        assert rc == 0
        # Both samples resolved (one from each VCF) and were estimated.
        content = out.read_text()
        assert "ADMIX_F0.10" in content
        assert "ADMIX_F0.50" in content

    def test_detect_admix_sample_not_in_any_vcf(self, tmp_path):
        """A --sample absent from every --admix-vcf is a clear error."""
        with pytest.raises(SystemExit, match="not found in any --admix-vcf"):
            main(
                [
                    "detect",
                    "--genotype-vcf",
                    str(JOINT_VCF),
                    "--admix-vcf",
                    str(JOINT_VCF),
                    "--host-sample",
                    "HOST",
                    "--donor-sample",
                    "DONOR",
                    "--sample",
                    "NO_SUCH_TIMEPOINT",
                    "--min-dp",
                    "0",
                    "--min-gq",
                    "0",
                ]
            )

    def test_invalid_sample_name(self):
        """CLI should fail with clear error for bad sample name."""
        with pytest.raises(SystemExit):
            main(
                [
                    "detect",
                    "--genotype-vcf",
                    str(JOINT_VCF),
                    "--admix-vcf",
                    str(JOINT_VCF),
                    "--host-sample",
                    "NONEXISTENT",
                    "--donor-sample",
                    "DONOR",
                    "--sample",
                    "ADMIX_F0.10",
                ]
            )


class TestHostPresenceCli:
    """Smoke tests for the host-presence detector wiring in the CLI.

    Sits next to TestCLIIntegration rather than in a dedicated tests/test_cli.py
    because the existing integration suite already owns the CLI smoke surface
    (JOINT_VCF fixture, sample naming conventions).
    """

    _HP_COLS = (
        "host_present_p",
        "host_f_est",
        "host_f_ci_lo",
        "host_f_ci_hi",
        "host_detect_markers",
        "host_err_source",
    )

    def _run_monitor(self, tmp_path, extra_args=()):
        out = tmp_path / "cli_hp.tsv"
        rc = main(
            [
                "detect",
                "--genotype-vcf",
                str(JOINT_VCF),
                "--admix-vcf",
                str(JOINT_VCF),
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR",
                "--sample",
                "ADMIX_F0.10",
                "--tsv",
                str(out),
                "--min-dp",
                "0",
                "--min-gq",
                "0",
                *extra_args,
            ]
        )
        assert rc == 0
        return out.read_text()

    def test_monitor_emits_host_presence_columns(self, tmp_path):
        content = self._run_monitor(tmp_path)
        header = content.splitlines()[0].split("\t")
        for col in self._HP_COLS:
            assert col in header, f"{col} missing from TSV header"
        data = content.splitlines()[1].split("\t")
        # host_err_source cell is one of the documented sentinels.
        idx = header.index("host_err_source")
        assert data[idx] in {"per-site", "global-fallback", "mixed", "none", "NA"}

    def test_no_host_presence_suppresses_values(self, tmp_path):
        """--no-host-presence keeps headers (TSV stays rectangular) but the
        cells should all be the NA sentinel."""
        content = self._run_monitor(tmp_path, ["--no-host-presence"])
        header = content.splitlines()[0].split("\t")
        data = content.splitlines()[1].split("\t")
        for col in self._HP_COLS:
            assert col in header
            assert data[header.index(col)] == "NA", f"expected NA for {col}"


class TestDynamicJointVcf:
    """Test building a joint VCF on the fly and running pipeline through it."""

    def test_round_trip(self, tmp_path):
        """Build joint VCF, write it, run pipeline, verify estimate."""
        host_vcf = TEST_DATA_DIR / "host.vcf"
        donor_vcf = TEST_DATA_DIR / "donor.vcf"

        result = build_joint_vcf(
            host_path=str(host_vcf),
            donor_paths=[str(donor_vcf)],
            admix_fractions=[0.20],
            admix_sample_names=["TP1"],
            host_sample_name="H",
            donor_sample_names=["D"],
            target_depth=2000,
            seed=99,
        )
        joint_path = tmp_path / "joint.vcf"
        write_joint_vcf(result, joint_path)

        # Run through the CLI
        out = tmp_path / "out.tsv"
        rc = main(
            [
                "detect",
                "--genotype-vcf",
                str(joint_path),
                "--admix-vcf",
                str(joint_path),
                "--host-sample",
                "H",
                "--donor-sample",
                "D",
                "--sample",
                "TP1",
                "--tsv",
                str(out),
                "--min-dp",
                "0",
                "--min-gq",
                "0",
            ]
        )
        assert rc == 0
        rows = out.read_text().splitlines()
        header = rows[0].split("\t")
        values = rows[1].split("\t")
        donor_pct = float(values[header.index("donor_pct")])
        # Truth is f=0.20; deep (2000x), 100-marker clean data recovers it tightly.
        assert donor_pct == pytest.approx(20.0, abs=3.0), (
            f"round-trip recovered donor_pct={donor_pct}, expected ~20"
        )


class TestSexCli:
    """Sex flags (#50): declared sex parsing, the retired --use-sex-chroms, columns."""

    @staticmethod
    def _argv(*extra: str) -> list[str]:
        return [
            "detect",
            "--genotype-vcf",
            str(JOINT_VCF),
            "--admix-vcf",
            str(JOINT_VCF),
            "--host-sample",
            "HOST",
            "--donor-sample",
            "DONOR",
            "--sample",
            "ADMIX_F0.10",
            "--min-dp",
            "0",
            "--min-gq",
            "0",
            *extra,
        ]

    def test_use_sex_chroms_is_retired(self, capsys):
        with pytest.raises(SystemExit) as exc:
            main(self._argv("--use-sex-chroms"))
        assert exc.value.code == 2
        err = capsys.readouterr().err
        assert "--use-sex-chroms has been retired" in err
        assert "--recipient-sex" in err and "--donor-sex" in err
        assert "removed in a later release" in err

    def test_use_sex_chroms_hidden_from_help(self, capsys):
        with pytest.raises(SystemExit) as exc:
            main(["detect", "--help"])
        assert exc.value.code == 0
        out = capsys.readouterr().out
        assert "--use-sex-chroms" not in out
        assert "--donor-sex" in out
        assert "--recipient-sex" in out

    def test_donor_sex_count_mismatch_errors(self):
        with pytest.raises(SystemExit, match="--donor-sex given 2 value"):
            main(self._argv("--donor-sex", "F", "--donor-sex", "M", "--tsv", "-"))

    def test_declared_recipient_sex_on_autosome_only_panel(self, tmp_path):
        """No chrX content: inference is unavailable, the declaration resolves it."""
        out = tmp_path / "sex.tsv"
        rc = main(self._argv("--recipient-sex", "F", "--tsv", str(out)))
        assert rc == 0
        header, values = (ln.split("\t") for ln in out.read_text().splitlines()[:2])
        row = dict(zip(header, values, strict=True))
        assert row["host_sex"] == "F"
        assert row["sex_source"] == "declared"
        assert row["donor_sex"] == "unavailable"
        assert row["donor_sex_source"] == "unavailable"
        assert row["sex_pair"] == "unknown"
        assert row["n_chrx_used"] == "0"
        assert row["n_par_excluded"] == "0"
        assert row["n_other_contig_excluded"] == "0"
        assert row["qc_status"] != "FAIL"

    def test_donor_sex_na_and_json_sex_object(self, tmp_path):
        out = tmp_path / "sex.json"
        rc = main(self._argv("--recipient-sex", "male", "--donor-sex", "NA", "--json", str(out)))
        assert rc == 0
        data = json.loads(out.read_text())
        sex = data["analysis"]["sex"]
        assert sex["host"]["declared"] == "M"
        assert sex["host"]["effective"] == "M"
        assert sex["host"]["source"] == "declared"
        assert sex["donors"][0]["declared"] is None
        assert sex["pair"] == "unknown"
        # NA is "nothing declared", so no donor sex is shown in the header meta.
        assert data["meta"]["donors"][0]["sex"] is None
        assert data["meta"]["sex"] == "male"
        assert data["params"]["contig_policy"] == "sex_aware"
        assert "use_sex_chroms" not in data["params"]

    def test_unparseable_recipient_sex_warns_and_runs(self, tmp_path, capsys):
        out = tmp_path / "sex.tsv"
        rc = main(self._argv("--recipient-sex", "intersex", "--tsv", str(out)))
        assert rc == 0
        assert "not F/female or M/male" in capsys.readouterr().err
        header, values = (ln.split("\t") for ln in out.read_text().splitlines()[:2])
        row = dict(zip(header, values, strict=True))
        assert row["host_sex"] == "unavailable"
        assert row["sex_source"] == "unavailable"


# ---------------------------------------------------------------------------
# Sex-aware chrX routing on the sexchrom fixtures (#46)
# ---------------------------------------------------------------------------


def _sexchrom_truth() -> dict[tuple[str, str], dict[str, str]]:
    rows = (SEXCHROM_DIR / "truth_table.tsv").read_text().splitlines()
    header = rows[0].split("\t")
    out = {}
    for line in rows[1:]:
        rec = dict(zip(header, line.split("\t"), strict=True))
        out[(rec["vcf"], rec["sample_name"])] = rec
    return out


def _detect_rows(tmp_path, vcf_name: str, *extra: str) -> dict[str, dict[str, str]]:
    """Run ``detect`` on one sexchrom fixture; return TSV rows keyed by sample."""
    vcf = SEXCHROM_DIR / vcf_name
    out = tmp_path / f"{vcf_name}.tsv"
    rc = main(
        [
            "detect",
            "--genotype-vcf",
            str(vcf),
            "--admix-vcf",
            str(vcf),
            "--host-sample",
            "HOST",
            "--donor-sample",
            "DONOR",
            "--sample",
            "ADMIX_F0.95",
            "--sample",
            "ADMIX_F0.80",
            "--min-dp",
            "0",
            "--min-gq",
            "0",
            "--tsv",
            str(out),
            *extra,
        ]
    )
    assert rc == 0
    lines = out.read_text().splitlines()
    header = lines[0].split("\t")
    # The multi-sample TSV repeats the header before every row.
    data = [ln for ln in lines if ln and not ln.startswith("sample\t")]
    return {
        rec["sample"]: rec
        for rec in (dict(zip(header, ln.split("\t"), strict=True)) for ln in data)
    }


class TestSexChromFixtures:
    """``detect`` on the FF / MM / MF fixtures: routing, counts and accuracy."""

    # Estimate tolerance against the blended truth, in percentage points. The
    # fixtures have 29-45 autosomal plus up to 10 chrX informative markers at
    # depth 1000, so the estimate sits well inside 1 pp of truth.
    TOL_PP = 1.0

    def _check_common(self, rows, truth, vcf_name):
        for sample, rec in rows.items():
            t = truth[(vcf_name, sample)]
            assert rec["n_par_excluded"] == t["n_chrx_par"] == "2"
            assert rec["n_other_contig_excluded"] == t["n_other_contig"] == "1"
            assert rec["host_sex"] == t["host_sex"]
            assert rec["qc_status"] != "FAIL"
            est = float(rec["donor_pct"])
            assert abs(est - 100 * float(t["true_donor_fraction"])) < self.TOL_PP, (sample, est)

    def test_ff_uses_chrx(self, tmp_path):
        truth = _sexchrom_truth()
        rows = _detect_rows(tmp_path, "joint_FF.vcf")
        self._check_common(rows, truth, "joint_FF.vcf")
        for sample, rec in rows.items():
            t = truth[("joint_FF.vcf", sample)]
            assert rec["donor_sex"] == "F"
            assert rec["sex_pair"] == "matched_female"
            assert rec["sex_source"] == rec["donor_sex_source"] == "inferred"
            n_chrx = int(rec["n_chrx_used"])
            assert n_chrx > 0
            assert n_chrx == int(t["n_chrx_nonpar_informative"])
            assert int(rec["n_informative"]) == int(t["n_autosomal_informative"]) + n_chrx
            assert rec["n_informative_sex_chrom_excluded"] == "0"
            assert rec["n_chrx_male_het_dropped"] == "0"

    def test_ff_autosomes_only_policy_matches_truth_without_chrx(self, tmp_path):
        truth = _sexchrom_truth()
        rows = _detect_rows(tmp_path, "joint_FF.vcf", "--contig-policy", "autosomes_only")
        self._check_common(rows, truth, "joint_FF.vcf")
        for sample, rec in rows.items():
            t = truth[("joint_FF.vcf", sample)]
            assert rec["sex_pair"] == "matched_female"  # inference still runs
            assert rec["n_chrx_used"] == "0"
            assert rec["n_informative"] == t["n_autosomal_informative"]
            assert int(rec["n_informative_sex_chrom_excluded"]) > 0

    def test_mm_hom_hom_used_and_spurious_hets_dropped(self, tmp_path):
        """Male/male pair: hom/hom chrX contrasts used, het calls dropped.

        The MM donor carries 2 spurious hets in 20 non-PAR chrX sites, which the
        likelihood ratio (male model het rate 0.02) leaves ambiguous, so the
        donor sex is declared; the declaration resolves the pair for routing.
        """
        truth = _sexchrom_truth()
        rows = _detect_rows(tmp_path, "joint_MM.vcf", "--recipient-sex", "M", "--donor-sex", "M")
        self._check_common(rows, truth, "joint_MM.vcf")
        for sample, rec in rows.items():
            t = truth[("joint_MM.vcf", sample)]
            assert rec["donor_sex"] == "M"
            assert rec["sex_pair"] == "matched_male"
            assert rec["sex_source"] == "inferred+declared"
            assert rec["donor_sex_source"] == "declared"
            n_drop = int(rec["n_chrx_male_het_dropped"])
            assert n_drop >= 1
            # Every spurious het site is dropped (host and donor hets are at
            # distinct sites in this fixture).
            n_spurious = int(t["n_chrx_host_spurious_het"]) + int(t["n_chrx_donor_spurious_het"])
            assert n_drop == n_spurious
            assert int(rec["n_chrx_used"]) > 0
            assert rec["n_informative_sex_chrom_excluded"] == "0"

    def test_mm_without_declaration_is_unknown_and_excludes_chrx(self, tmp_path):
        truth = _sexchrom_truth()
        rows = _detect_rows(tmp_path, "joint_MM.vcf")
        self._check_common(rows, truth, "joint_MM.vcf")
        for rec in rows.values():
            assert rec["donor_sex"] == "ambiguous"
            assert rec["sex_pair"] == "unknown"
            assert rec["n_chrx_used"] == "0"
            assert rec["n_chrx_male_het_dropped"] == "0"
            assert int(rec["n_informative_sex_chrom_excluded"]) > 0

    def test_mf_excludes_chrx(self, tmp_path):
        truth = _sexchrom_truth()
        rows = _detect_rows(tmp_path, "joint_MF.vcf")
        self._check_common(rows, truth, "joint_MF.vcf")
        for sample, rec in rows.items():
            t = truth[("joint_MF.vcf", sample)]
            assert (rec["host_sex"], rec["donor_sex"]) == ("M", "F")
            assert rec["sex_pair"] == "mismatched"
            assert rec["n_chrx_used"] == "0"
            assert rec["n_chrx_male_het_dropped"] == "0"
            assert int(rec["n_informative_sex_chrom_excluded"]) > 0
            assert rec["n_informative"] == t["n_autosomal_informative"]

    def test_mf_stderr_explains_exclusion(self, tmp_path, capsys):
        _detect_rows(tmp_path, "joint_MF.vcf")
        err = capsys.readouterr().err
        assert "host and donor sexes differ (mismatched: host M, donor F)" in err
        assert "pseudoautosomal" in err and "non-primary contig" in err

    def test_ff_stderr_reports_chrx_use(self, tmp_path, capsys):
        _detect_rows(tmp_path, "joint_FF.vcf")
        err = capsys.readouterr().err
        assert "informative non-PAR chrX marker(s) (matched_female: host F, donor F)" in err


class TestSexMismatchCrossCheck:
    """The sex-mismatch cross-check (#48) on the MF fixture through ``detect``."""

    def test_mf_fixture_chrx_estimate(self, tmp_path):
        truth = _sexchrom_truth()
        rows = _detect_rows(tmp_path, "joint_MF.vcf")
        for sample, rec in rows.items():
            t = truth[("joint_MF.vcf", sample)]
            assert rec["sex_pair"] == "mismatched"
            assert rec["sexchrom_basis"] == "chrX-cn"
            assert int(rec["sexchrom_n"]) >= 5
            est = float(rec["sexchrom_frac"])
            assert abs(est - float(t["true_donor_fraction"])) < 0.02, (sample, est)
            lo, hi = (float(x) for x in rec["sexchrom_ci"].split(","))
            assert lo <= est <= hi
            assert rec["sexchrom_concordant"] == "true"
            assert rec["qc_status"] != "FAIL"

    def test_mf_stderr_reports_cross_check(self, tmp_path, capsys):
        _detect_rows(tmp_path, "joint_MF.vcf")
        err = capsys.readouterr().err
        assert "sex-chromosome cross-check (chrX-cn, n=" in err
        assert "(concordant)" in err

    def test_matched_fixture_columns_na(self, tmp_path):
        rows = _detect_rows(tmp_path, "joint_FF.vcf")
        for rec in rows.values():
            assert rec["sexchrom_basis"] == "NA"
            assert rec["sexchrom_frac"] == "NA"
            assert rec["sexchrom_concordant"] == "NA"

    def test_mf_json_and_html(self, tmp_path):
        vcf = SEXCHROM_DIR / "joint_MF.vcf"
        out_json = tmp_path / "mf.json"
        out_html = tmp_path / "mf.html"
        rc = main(
            [
                "detect",
                "--genotype-vcf",
                str(vcf),
                "--admix-vcf",
                str(vcf),
                "--host-sample",
                "HOST",
                "--donor-sample",
                "DONOR",
                "--sample",
                "ADMIX_F0.80",
                "--min-dp",
                "0",
                "--min-gq",
                "0",
                "--json",
                str(out_json),
                "--html",
                str(out_html),
            ]
        )
        assert rc == 0
        sm = json.loads(out_json.read_text())["analysis"]["sex_mismatch"]
        assert sm["basis"] == "chrX-cn" and sm["concordant"] is True
        assert sm["male_party"] == "host"
        assert "Sex-chromosome cross-check" in out_html.read_text()
