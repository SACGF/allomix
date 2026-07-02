"""Tests for allomix.qc.caller — detecting the variant caller from a VCF header."""

from allomix.calibration.bias import (
    MarkerBias,
    load_bias_table,
    read_bias_table_caller,
    save_bias_table,
)
from allomix.qc.caller import (
    Caller,
    caller_from_token,
    detect_caller,
    detect_caller_from_header,
)

# Header snippets modelled on the real SRP434573 VCFs (paper/public_data).
_GATK_HEADER = (
    "##fileformat=VCFv4.2\n"
    '##GATKCommandLine=<ID=GenotypeGVCFs,CommandLine="GenotypeGVCFs ...",Version="4.6.2.0">\n'
    "##source=HaplotypeCaller\n"
    '##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="...">\n'
    '##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="...">\n'
    '##INFO=<ID=FS,Number=1,Type=Float,Description="...">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHOST\n"
)

_MPILEUP_HEADER = (
    "##fileformat=VCFv4.2\n"
    "##bcftoolsVersion=1.23.1\n"
    "##bcftoolsCommand=mpileup -f ref.fasta -a FORMAT/AD,FORMAT/DP sample.bam\n"
    '##INFO=<ID=DP4,Number=4,Type=Integer,Description="...">\n'
    "##bcftools_callCommand=call -m -A -C alleles\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tADMIX\n"
)


def test_gatk_from_command_line():
    info = detect_caller_from_header(_GATK_HEADER)
    assert info.caller is Caller.GATK


def test_mpileup_from_command_line():
    info = detect_caller_from_header(_MPILEUP_HEADER)
    assert info.caller is Caller.MPILEUP


def test_gatk_from_info_fingerprint_only():
    """Command lines stripped: GATK-only INFO IDs, no DP4."""
    header = (
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=SOR,Number=1,Type=Float,Description="...">\n'
        '##INFO=<ID=MQRankSum,Number=1,Type=Float,Description="...">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tHOST\n"
    )
    assert detect_caller_from_header(header).caller is Caller.GATK


def test_mpileup_from_info_fingerprint_only():
    header = (
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=DP4,Number=4,Type=Integer,Description="...">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tADMIX\n"
    )
    assert detect_caller_from_header(header).caller is Caller.MPILEUP


def test_explicit_stamp_overrides_sniffing():
    """An explicit ##allomixCaller wins even against contradictory signatures."""
    header = "##allomixCaller=mpileup\n" + _GATK_HEADER
    info = detect_caller_from_header(header)
    assert info.caller is Caller.MPILEUP
    assert "allomixCaller" in info.evidence


def test_unrecognised_explicit_value_falls_through():
    header = "##allomixCaller=freebayes\n" + _MPILEUP_HEADER
    assert detect_caller_from_header(header).caller is Caller.MPILEUP


def test_unknown_when_no_signature():
    header = (
        "##fileformat=VCFv4.2\n"
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="...">\n'
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
    )
    assert detect_caller_from_header(header).caller is Caller.UNKNOWN


def test_caller_from_token():
    assert caller_from_token("gatk") is Caller.GATK
    assert caller_from_token("MPILEUP") is Caller.MPILEUP
    assert caller_from_token(None) is Caller.UNKNOWN
    assert caller_from_token("nonsense") is Caller.UNKNOWN


def test_detect_caller_from_file(tmp_path):
    p = tmp_path / "admix.vcf"
    p.write_text(_MPILEUP_HEADER)
    assert detect_caller(str(p)).caller is Caller.MPILEUP


def _bias():
    return {("chr1", 100, "A", "G"): MarkerBias("chr1", 100, "A", "G", bias=0.02, n_het=5)}


def test_bias_table_records_and_reads_caller(tmp_path):
    p = tmp_path / "bias.tsv"
    save_bias_table(_bias(), p, caller="gatk")
    assert read_bias_table_caller(p) == "gatk"
    # The comment line does not disturb the data rows.
    loaded = load_bias_table(p)
    assert loaded == {("chr1", 100, "A", "G"): 0.02}


def test_bias_table_without_caller_reads_none(tmp_path):
    p = tmp_path / "bias.tsv"
    save_bias_table(_bias(), p)
    assert read_bias_table_caller(p) is None
    assert load_bias_table(p) == {("chr1", 100, "A", "G"): 0.02}
