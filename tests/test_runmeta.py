"""Tests for allomix.runmeta — reading ##allomixRunUnit header lines."""

from allomix.runmeta import RunUnitInfo, read_run_units

_VCF_BODY = (
    "##fileformat=VCFv4.2\n"
    "##contig=<ID=chr1,length=1000>\n"
    "{runmeta}"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\n"
    "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/0\t0/1\t1/1\n"
)


def _write_vcf(path, runmeta_lines):
    body = "".join(line + "\n" for line in runmeta_lines)
    path.write_text(_VCF_BODY.format(runmeta=body))
    return str(path)


def test_reads_all_fields(tmp_path):
    vcf = _write_vcf(
        tmp_path / "a.vcf",
        ["##allomixRunUnit=<ID=s1,RunUnit=FC1:1,Source=RG:PU,SharesRunWithHost=true>"],
    )
    m = read_run_units(vcf)
    assert m == {"s1": RunUnitInfo(run_unit="FC1:1", source="RG:PU", shares_run_with_host=True)}


def test_bool_variants(tmp_path):
    vcf = _write_vcf(
        tmp_path / "a.vcf",
        [
            "##allomixRunUnit=<ID=s1,RunUnit=FC1:1,Source=RG:PU,SharesRunWithHost=true>",
            "##allomixRunUnit=<ID=s2,RunUnit=FC2:2,Source=RG:PU,SharesRunWithHost=false>",
            "##allomixRunUnit=<ID=s3,RunUnit=FC3:3,Source=readname,SharesRunWithHost=unknown>",
        ],
    )
    m = read_run_units(vcf)
    assert m["s1"].shares_run_with_host is True
    assert m["s2"].shares_run_with_host is False
    assert m["s3"].shares_run_with_host is None  # "unknown" -> None
    assert m["s3"].source == "readname"


def test_multi_unit_value(tmp_path):
    vcf = _write_vcf(
        tmp_path / "a.vcf",
        ["##allomixRunUnit=<ID=s1,RunUnit=FC1:1;FC1:2,Source=RG:PU,SharesRunWithHost=true>"],
    )
    assert read_run_units(vcf)["s1"].run_unit == "FC1:1;FC1:2"


def test_no_metadata_returns_empty(tmp_path):
    # The optional case: a VCF without run metadata yields no entries.
    vcf = _write_vcf(tmp_path / "a.vcf", [])
    assert read_run_units(vcf) == {}


def test_ignores_other_header_lines(tmp_path):
    vcf = _write_vcf(
        tmp_path / "a.vcf",
        [
            "##SAMPLE=<ID=s1,Description=something>",
            "##allomixRunUnit=<ID=s2,RunUnit=FC2:2,Source=RG:PU,SharesRunWithHost=false>",
        ],
    )
    m = read_run_units(vcf)
    assert set(m) == {"s2"}


def test_strips_quoted_values(tmp_path):
    # Defensive: a writer that quotes a value should still parse.
    vcf = _write_vcf(
        tmp_path / "a.vcf",
        ['##allomixRunUnit=<ID=s1,RunUnit="FC1:1",Source=RG:PU,SharesRunWithHost=true>'],
    )
    assert read_run_units(vcf)["s1"].run_unit == "FC1:1"
