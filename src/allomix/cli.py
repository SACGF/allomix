"""Command-line interface for allomix."""

import argparse
import json
import sys

from cyvcf2 import VCF

from allomix import __version__
from allomix.analysis import analyse_sample
from allomix.bias import (
    biases_to_simple_dict,
    estimate_biases,
    load_bias_table,
    save_bias_table,
)
from allomix.chimerism import PanelCalibration
from allomix.constants import (
    DEFAULT_ERROR_RATE,
    DEFAULT_MIN_DP,
    DEFAULT_MIN_GQ,
    ROBUST_K_DEFAULT,
)
from allomix.error_rates import (
    estimate_error_rates,
    load_error_table,
    save_error_table,
)
from allomix.genotype import parse_vcf
from allomix.relatedness import VALID_DECLARATIONS
from allomix.report import timeline_json, to_json, to_tsv


def _expected_relatedness_value(value: str) -> str:
    """Validate one ``--expected-relatedness`` value (argparse ``type``).

    Accepts the relationship declarations plus NA (case-insensitive), returning
    the lowercased form. Rejects "identical" with an explanation: host and donor
    are only identical for a monozygotic-twin (syngeneic) donor, which has no
    host/donor genetic differences to measure, so genotype-based chimerism does
    not apply and there is nothing to declare.
    """
    v = value.strip().lower()
    if v == "identical":
        raise argparse.ArgumentTypeError(
            "'identical' is not a valid expected relatedness. Host and donor are "
            "only identical for a monozygotic-twin (syngeneic) donor, which has "
            "no host/donor genetic differences to measure, so genotype-based "
            "chimerism does not apply. (If samples do come back identical, "
            "allomix fails QC and says so.)"
        )
    allowed = {*VALID_DECLARATIONS, "na"}
    if v not in allowed:
        valid = ", ".join([*VALID_DECLARATIONS, "NA"])
        raise argparse.ArgumentTypeError(
            f"invalid expected relatedness {value!r}; choose from {valid}"
        )
    return v


def _add_common_args(parser: argparse.ArgumentParser) -> None:
    """Add arguments shared between monitor and timeline."""
    parser.add_argument(
        "--panel-vcf",
        required=True,
        help="Panel VCF with host/donor genotypes (typically GATK joint-called; "
             "see doc/joint_calling.md)",
    )
    parser.add_argument(
        "--admix-vcf",
        required=True,
        help="Admix VCF with raw pileup AD (typically bcftools mpileup output)",
    )
    parser.add_argument("--host-sample", required=True, help="Host sample name in VCF")
    parser.add_argument(
        "--donor-sample",
        required=True,
        action="append",
        metavar="SAMPLE_NAME",
        help="Donor sample name in VCF (repeat for multi-donor)",
    )
    parser.add_argument(
        "--expected-relatedness",
        action="append",
        metavar="RELATIONSHIP",
        type=_expected_relatedness_value,
        help="Declared host-vs-donor relationship for the QC relatedness check, "
             "one per --donor-sample in the same order (repeat to match). One of "
             f"{', '.join(VALID_DECLARATIONS)} or NA (no expectation). A declared "
             "relationship that crosses the related/unrelated boundary fails QC. "
             "'identical' is rejected: an identical-twin (syngeneic) donor cannot "
             "be monitored by genotype.",
    )
    parser.add_argument(
        "--relatedness-tolerance",
        type=int,
        default=1,
        help="Allowed degree distance before a declared-vs-detected relatedness "
             "mismatch is flagged for review (default: 1)",
    )
    parser.add_argument(
        "--sample",
        required=True,
        action="append",
        metavar="SAMPLE_NAME",
        help="Admixture sample name in VCF (repeat for multiple timepoints)",
    )
    parser.add_argument("--output", "-o", default="-", help="Output file (default: stdout)")
    parser.add_argument(
        "--min-dp",
        type=int,
        default=DEFAULT_MIN_DP,
        help=f"Minimum depth (default: {DEFAULT_MIN_DP})",
    )
    parser.add_argument(
        "--min-gq",
        type=int,
        default=DEFAULT_MIN_GQ,
        help=f"Minimum GQ (default: {DEFAULT_MIN_GQ})",
    )
    parser.add_argument(
        "--use-sex-chroms",
        action="store_true",
        help="Include sex / mitochondrial contigs (X/Y/M). Off by default: in "
             "sex-mismatched transplants the host/donor dosage on chrX/chrY is "
             "wrong. Enable per run only once host and donor sex are known to "
             "match. The informative sex-chrom markers dropped are reported.",
    )
    parser.add_argument(
        "--error-rate",
        type=float,
        default=DEFAULT_ERROR_RATE,
        help=f"Sequencing error rate (default: {DEFAULT_ERROR_RATE})",
    )
    parser.add_argument(
        "--robust",
        choices=["off", "auto", "force"],
        default="auto",
        help="Robust refit: iteratively drop residual-outlier markers "
             "(host copy-number / LoH, genotyping errors) and refit. 'auto' "
             "(default) keeps a marker floor and is a no-op on clean data; "
             "'force' trims further; 'off' disables. A large exclusion is "
             "flagged for REVIEW.",
    )
    parser.add_argument(
        "--robust-k",
        type=float,
        default=ROBUST_K_DEFAULT,
        help="Robust residual cut in robust SDs (median/MAD) for --robust "
             f"(default: {ROBUST_K_DEFAULT})",
    )
    parser.add_argument("--verbose", action="store_true", help="Include per-marker detail")
    parser.add_argument(
        "--bias-table",
        default=None,
        help="Per-marker bias table TSV (from allomix estimate-bias or simulation)",
    )
    parser.add_argument(
        "--no-bias-correction",
        action="store_true",
        help="Disable bias correction even when a bias table is provided",
    )
    parser.add_argument(
        "--estimate-bias",
        action="store_true",
        help="Estimate per-marker bias inline from all samples in --panel-vcf, "
             "held in memory (no separate `estimate-bias` step or table file). "
             "Mutually exclusive with --bias-table. Estimate from data called "
             "the same way as the admix; works best when the panel VCF holds "
             "many samples.",
    )
    parser.add_argument(
        "--estimate-bias-min-het",
        type=int,
        default=1,
        help="Minimum het observations per marker for inline --estimate-bias "
             "(default: 1).",
    )
    parser.add_argument(
        "--error-table",
        default=None,
        help="Per-site empirical error-rate table TSV (from "
             "`allomix estimate-errors`). Sites with per-direction rates "
             "override --error-rate; missing sites or missing directions "
             "fall back to --error-rate.",
    )
    parser.add_argument(
        "--no-error-correction",
        action="store_true",
        help="Disable empirical error-rate correction even when an error "
             "table is provided",
    )
    parser.add_argument(
        "--no-host-presence",
        action="store_true",
        help="Disable the host-presence detection test (see "
             "`allomix.detect`). On by default; cheap to run.",
    )
    parser.add_argument(
        "--no-artifact-filter",
        action="store_true",
        help="Disable the read-level artifact filter in the host-presence "
             "test (strand/soft-clip/read-position bias). On by default; "
             "drops alignment-artifact markers (see `allomix.detect`).",
    )


def _validate_expected_relatedness(args: argparse.Namespace) -> None:
    """Check --expected-relatedness count matches the number of donors.

    Raises SystemExit with a clear message rather than letting a count mismatch
    surface later as a strict-zip error in the QC step.
    """
    er = args.expected_relatedness
    if er is not None and len(er) != len(args.donor_sample):
        raise SystemExit(
            f"--expected-relatedness given {len(er)} value(s) but there are "
            f"{len(args.donor_sample)} donor(s); provide exactly one per "
            "--donor-sample, in the same order (use NA for no expectation)"
        )


def _validate_sample_names(vcf_path: str, required: list[str]) -> None:
    """Check that all required sample names exist in the VCF header.

    Raises SystemExit with a clear error listing available names if any
    are missing.
    """
    vcf = VCF(vcf_path)
    available = list(vcf.samples)
    vcf.close()
    missing = [s for s in required if s not in available]
    if missing:
        raise SystemExit(f"Sample(s) not found in {vcf_path}: {missing}\nAvailable: {available}")


def _run_single_sample(
    host: list,
    donors: list[list],
    vcf_path: str,
    admix_sample: str,
    min_dp: int,
    min_gq: int,
    error_rate: float,
    calibration: PanelCalibration | None = None,
    run_host_presence: bool = True,
    use_sex_chroms: bool = False,
    artifact_filter: bool = True,
    robust: str = "off",
    robust_k: float = ROBUST_K_DEFAULT,
    expected_relatedness: list[str] | None = None,
    relatedness_tolerance: int = 1,
) -> tuple:
    """Run the chimerism pipeline for one admixture sample.

    Takes pre-parsed host and donor markers to avoid redundant VCF reads, then
    delegates to ``allomix.analysis.analyse_sample`` (shared with the
    diagnostic scripts). Automatically uses multi-donor estimation when more
    than one donor is provided.

    Returns (ChimerismResult | MultiDonorResult, QCReport, MarkerGenotypes).
    """
    admix = parse_vcf(vcf_path, sample=admix_sample, min_dp=0)

    analysis = analyse_sample(
        host,
        donors,
        admix,
        min_dp=min_dp,
        min_gq=min_gq,
        error_rate=error_rate,
        calibration=calibration,
        run_host_presence=run_host_presence,
        use_sex_chroms=use_sex_chroms,
        artifact_filter=artifact_filter,
        sample_name=admix_sample,
        robust=robust,
        robust_k=robust_k,
        expected_relatedness=expected_relatedness,
        relatedness_tolerance=relatedness_tolerance,
    )

    if not use_sex_chroms and analysis.genotypes.n_sex_chrom_excluded:
        print(
            f"{admix_sample}: excluded {analysis.genotypes.n_sex_chrom_excluded} "
            "informative sex-chromosome marker(s) (use --use-sex-chroms to keep them)",
            file=sys.stderr,
        )

    return analysis.result, analysis.qc, analysis.genotypes


def _open_output(path: str):
    """Open output file or return stdout."""
    if path == "-":
        return sys.stdout
    return open(path, "w", encoding="utf-8")


def _load_calibration(args: argparse.Namespace) -> PanelCalibration:
    """Build the per-marker calibration from the CLI table options.

    Bias comes from --bias-table, or is estimated inline from the panel VCF
    samples when --estimate-bias is set (issue #11), or is empty. Per-site error
    comes from --error-table or is empty. The --no-*-correction flags force the
    respective correction off.
    """
    if getattr(args, "estimate_bias", False) and not args.no_bias_correction:
        if args.bias_table:
            raise SystemExit("Use either --bias-table or --estimate-bias, not both")
        samples = list(VCF(args.panel_vcf).samples)
        marker_lists = [
            parse_vcf(args.panel_vcf, sample=s, min_dp=0, min_gq=0) for s in samples
        ]
        biases = biases_to_simple_dict(
            estimate_biases(marker_lists, min_het=args.estimate_bias_min_het)
        )
        sys.stderr.write(
            f"Estimated per-marker bias for {len(biases)} marker(s) from "
            f"{len(samples)} panel sample(s)\n"
        )
    elif args.bias_table and not args.no_bias_correction:
        biases = load_bias_table(args.bias_table)
    else:
        biases = {}
    errors = (
        load_error_table(args.error_table)
        if args.error_table and not args.no_error_correction
        else {}
    )
    return PanelCalibration(biases=biases, errors=errors)


def cmd_monitor(args: argparse.Namespace) -> int:
    """Run the monitor subcommand."""
    _validate_expected_relatedness(args)
    _validate_sample_names(args.panel_vcf, [args.host_sample] + args.donor_sample)
    _validate_sample_names(args.admix_vcf, args.sample)

    calibration = _load_calibration(args)

    # Parse host and donors once — they're the same for every timepoint.
    # gt_ad_consistency=True is the panel-side miscall guard: drops
    # markers where the called GT contradicts the AD VAF (e.g. GATK
    # called het from 20% VAF reads in a 2-sample joint call). Without
    # it, the wider gnomAD-derived panel recovers markers that bias the
    # estimator toward false host signal — see Step 23 verification.
    host = parse_vcf(
        args.panel_vcf, sample=args.host_sample, min_gq=args.min_gq, gt_ad_consistency=True
    )
    donors = [
        parse_vcf(args.panel_vcf, sample=d, min_gq=args.min_gq, gt_ad_consistency=True)
        for d in args.donor_sample
    ]

    out = _open_output(args.output)
    try:
        for sample_name in args.sample:
            result, qc, genotypes = _run_single_sample(
                host,
                donors,
                args.admix_vcf,
                sample_name,
                args.min_dp,
                args.min_gq,
                args.error_rate,
                calibration=calibration,
                run_host_presence=not args.no_host_presence,
                use_sex_chroms=args.use_sex_chroms,
                artifact_filter=not args.no_artifact_filter,
                robust=args.robust,
                robust_k=args.robust_k,
                expected_relatedness=args.expected_relatedness,
                relatedness_tolerance=args.relatedness_tolerance,
            )

            if args.format == "json":
                data = to_json(result, qc, sample_name=genotypes.sample_name)
                out.write(json.dumps(data, indent=2) + "\n")
            else:
                to_tsv(result, qc, out, verbose=args.verbose, sample_name=genotypes.sample_name)
    finally:
        if out is not sys.stdout:
            out.close()

    return 0


def cmd_timeline(args: argparse.Namespace) -> int:
    """Run the timeline subcommand."""
    _validate_expected_relatedness(args)
    _validate_sample_names(args.panel_vcf, [args.host_sample] + args.donor_sample)
    _validate_sample_names(args.admix_vcf, args.sample)

    calibration = _load_calibration(args)

    # Parse host and donors once. See cmd_monitor for gt_ad_consistency.
    host = parse_vcf(
        args.panel_vcf, sample=args.host_sample, min_gq=args.min_gq, gt_ad_consistency=True
    )
    donors = [
        parse_vcf(args.panel_vcf, sample=d, min_gq=args.min_gq, gt_ad_consistency=True)
        for d in args.donor_sample
    ]

    results = []
    for sample_name in args.sample:
        result, qc, genotypes = _run_single_sample(
            host,
            donors,
            args.admix_vcf,
            sample_name,
            args.min_dp,
            args.min_gq,
            args.error_rate,
            calibration=calibration,
            run_host_presence=not args.no_host_presence,
            use_sex_chroms=args.use_sex_chroms,
            artifact_filter=not args.no_artifact_filter,
            robust=args.robust,
            robust_k=args.robust_k,
            expected_relatedness=args.expected_relatedness,
            relatedness_tolerance=args.relatedness_tolerance,
        )
        results.append((genotypes.sample_name, result, qc))

    data = timeline_json(results)

    out = _open_output(args.output)
    try:
        out.write(json.dumps(data, indent=2) + "\n")
    finally:
        if out is not sys.stdout:
            out.close()

    return 0


def cmd_estimate_bias(args: argparse.Namespace) -> int:
    """Run the estimate-bias subcommand."""
    if args.vcfs and args.vcf:
        raise SystemExit("Use either --vcfs or --vcf/--samples, not both")
    if not args.vcfs and not args.vcf:
        raise SystemExit("One of --vcfs or --vcf is required")
    if args.vcf and not args.samples:
        raise SystemExit("--samples is required when using --vcf")

    marker_lists = []
    if args.vcfs:
        for vcf_path in args.vcfs:
            markers = parse_vcf(vcf_path, min_dp=0, min_gq=0)
            marker_lists.append(markers)
        n_source = f"{len(args.vcfs)} VCFs"
    else:
        _validate_sample_names(args.vcf, args.samples)
        for sample in args.samples:
            markers = parse_vcf(args.vcf, sample=sample, min_dp=0, min_gq=0)
            marker_lists.append(markers)
        n_source = f"{len(args.samples)} samples from {args.vcf}"

    biases = estimate_biases(marker_lists, min_het=args.min_het)
    save_bias_table(biases, args.output)
    print(
        f"Estimated bias for {len(biases)} markers from {n_source} -> {args.output}",
        file=sys.stderr,
    )
    return 0


def cmd_estimate_errors(args: argparse.Namespace) -> int:
    """Run the estimate-errors subcommand."""
    if args.vcfs and args.vcf:
        raise SystemExit("Use either --vcfs or --vcf/--samples, not both")
    if not args.vcfs and not args.vcf:
        raise SystemExit("One of --vcfs or --vcf is required")
    if args.vcf and not args.samples:
        raise SystemExit("--samples is required when using --vcf")

    marker_lists = []
    if args.vcfs:
        for vcf_path in args.vcfs:
            markers = parse_vcf(vcf_path, min_dp=0, min_gq=args.min_gq)
            marker_lists.append(markers)
        n_source = f"{len(args.vcfs)} VCFs"
    else:
        _validate_sample_names(args.vcf, args.samples)
        for sample in args.samples:
            markers = parse_vcf(args.vcf, sample=sample, min_dp=0, min_gq=args.min_gq)
            marker_lists.append(markers)
        n_source = f"{len(args.samples)} samples from {args.vcf}"

    errors = estimate_error_rates(
        marker_lists,
        min_reads=args.min_reads,
        max_vaf_homref=args.max_vaf_homref,
        min_vaf_homalt=args.min_vaf_homalt,
    )
    save_error_table(errors, args.output)
    print(
        f"Estimated error rates for {len(errors)} sites from {n_source} "
        f"-> {args.output}",
        file=sys.stderr,
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    """Entry point for the allomix CLI."""
    parser = argparse.ArgumentParser(
        prog="allomix",
        description="NGS-based donor chimerism monitoring for HSCT",
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    subparsers = parser.add_subparsers(dest="command")

    monitor_parser = subparsers.add_parser(
        "monitor",
        help="Calculate chimerism for one or more samples",
    )
    _add_common_args(monitor_parser)
    monitor_parser.add_argument(
        "--format",
        choices=["tsv", "json"],
        default="tsv",
        help="Output format (default: tsv)",
    )

    timeline_parser = subparsers.add_parser(
        "timeline",
        help="Generate chimerism timeline across timepoints",
    )
    _add_common_args(timeline_parser)

    bias_parser = subparsers.add_parser(
        "estimate-bias",
        help="Estimate per-marker amplification bias from VCFs",
    )
    bias_input = bias_parser.add_mutually_exclusive_group()
    bias_input.add_argument(
        "--vcfs",
        nargs="+",
        metavar="VCF",
        help="Per-sample VCFs, one per file (reads first sample from each)",
    )
    bias_input.add_argument(
        "--vcf",
        metavar="VCF",
        help="Joint-called multi-sample VCF (use with --samples)",
    )
    bias_parser.add_argument(
        "--samples",
        nargs="+",
        metavar="SAMPLE_NAME",
        help="Sample names to extract from --vcf",
    )
    bias_parser.add_argument(
        "--output",
        "-o",
        default="bias_table.tsv",
        help="Output bias table TSV (default: bias_table.tsv)",
    )
    bias_parser.add_argument(
        "--min-het",
        type=int,
        default=1,
        help="Minimum het observations per marker (default: 1)",
    )

    err_parser = subparsers.add_parser(
        "estimate-errors",
        help="Estimate per-site empirical error rates from VCFs",
    )
    err_input = err_parser.add_mutually_exclusive_group()
    err_input.add_argument(
        "--vcfs",
        nargs="+",
        metavar="VCF",
        help="Per-sample VCFs, one per file (reads first sample from each)",
    )
    err_input.add_argument(
        "--vcf",
        metavar="VCF",
        help="Joint-called multi-sample VCF (use with --samples)",
    )
    err_parser.add_argument(
        "--samples",
        nargs="+",
        metavar="SAMPLE_NAME",
        help="Sample names to extract from --vcf",
    )
    err_parser.add_argument(
        "--output",
        "-o",
        default="error_table.tsv",
        help="Output error table TSV (default: error_table.tsv)",
    )
    err_parser.add_argument(
        "--min-reads",
        type=int,
        default=1000,
        help="Minimum total reads per direction to retain a site's estimate "
             "(default: 1000)",
    )
    err_parser.add_argument(
        "--max-vaf-homref",
        type=float,
        default=0.10,
        help="Drop hom-ref training observations with vaf > this "
             "(default: 0.10)",
    )
    err_parser.add_argument(
        "--min-vaf-homalt",
        type=float,
        default=0.90,
        help="Drop hom-alt training observations with vaf < this "
             "(default: 0.90)",
    )
    err_parser.add_argument(
        "--min-gq",
        type=int,
        default=DEFAULT_MIN_GQ,
        help=f"Minimum GQ for training calls (default: {DEFAULT_MIN_GQ})",
    )

    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 1

    if args.command == "monitor":
        return cmd_monitor(args)
    if args.command == "timeline":
        return cmd_timeline(args)
    if args.command == "estimate-bias":
        return cmd_estimate_bias(args)
    if args.command == "estimate-errors":
        return cmd_estimate_errors(args)
    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
