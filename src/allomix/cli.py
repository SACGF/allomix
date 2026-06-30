"""Command-line interface for allomix."""

import argparse
import datetime
import importlib.util
import json
import sys
from pathlib import Path

from cyvcf2 import VCF

from allomix import __version__
from allomix.analysis import analyse_sample
from allomix.bias import (
    biases_to_simple_dict,
    estimate_biases,
    estimate_biases_both_het,
    load_bias_table,
    save_bias_table,
)
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
from allomix.html.render import render_single
from allomix.likelihood import PanelCalibration
from allomix.marker_contamination import (
    DEFAULT_DOSE_CAP,
    DEFAULT_GATE_ALPHA,
    DEFAULT_GATE_MIN_SLOPE,
    estimate_contamination_table,
    load_contamination_table,
    save_contamination_table,
)
from allomix.relatedness import VALID_DECLARATIONS
from allomix.report import (
    DonorMeta,
    ReportMeta,
    report_data,
    timeline_report_data,
    to_marker_csv,
    to_tsv,
)
from allomix.runmeta import RunUnitInfo, read_run_units


def _expected_relatedness_value(value: str) -> str:
    """Validate one ``--expected-relatedness`` value (argparse ``type``).

    Accepts the relationship declarations plus NA (case-insensitive), returns the
    lowercased form. Rejects "identical" (see the raised error for why).
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
    parser.add_argument(
        "--marker-csv",
        type=Path,
        default=None,
        metavar="PATH",
        help="Also write per-marker detail (allele depths, observed vs expected "
        "VAF, residual, include flag) to this CSV. This is the "
        "bioinformatician-facing detail that the clinician HTML report omits; "
        "one row per marker per sample, with a leading sample column.",
    )
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
    parser.add_argument(
        "--marker-type-overdispersion",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Fit a separate beta-binomial concentration (rho) for the "
        "donor-homozygous and donor-heterozygous marker classes instead of "
        "one shared rho (single-donor only). Removes the sub-0.5%% MLE "
        "floor (issue #33). On by default; use --no-marker-type-overdispersion "
        "for the legacy shared-rho estimator. Falls back to shared rho for a "
        "sample when a class has too few markers.",
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
        help="Minimum het observations per marker for inline --estimate-bias (default: 1).",
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
        help="Disable empirical error-rate correction even when an error table is provided",
    )
    parser.add_argument(
        "--contamination-table",
        default=None,
        help="Per-marker co-pooled contamination table TSV (from "
        "`allomix build-contamination-table`). Used only when "
        "--contamination-correction is set (Step 30, issue #30).",
    )
    parser.add_argument(
        "--contamination-correction",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Subtract dose-predicted co-pooled contamination from "
        "donor-homozygous host-allele counts before the MLE, using "
        "--contamination-table (Step 30, issue #30). OFF by default. A "
        "table built on a clean flowcell gates itself out, so this is a "
        "no-op there even when enabled.",
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


def _add_report_meta_args(parser: argparse.ArgumentParser) -> None:
    """Add the optional HTML-report metadata flags.

    All optional: they populate the HTML header band and are ignored by the
    tsv/json formats.
    """
    group = parser.add_argument_group("HTML report metadata (optional)")
    group.add_argument("--recipient-id", help="Recipient identifier for the report header")
    group.add_argument("--recipient-name", help="Recipient display name")
    group.add_argument("--recipient-sex", help="Recipient sex (shown verbatim)")
    group.add_argument("--recipient-dob", help="Recipient date of birth (shown verbatim)")
    group.add_argument(
        "--transplant-type",
        default="HSCT",
        help="Transplant type label (default: HSCT)",
    )
    group.add_argument(
        "--transplant-date",
        help="Transplant date, ISO (YYYY-MM-DD) for the days-post-transplant derivation",
    )
    group.add_argument(
        "--donor-relationship",
        action="append",
        metavar="RELATIONSHIP",
        help="Declared donor relationship shown in the header, one per "
        "--donor-sample in the same order (free text, e.g. 'unrelated', "
        "'sibling'). Separate from --expected-relatedness, which drives QC.",
    )
    group.add_argument(
        "--sample-date",
        action="append",
        metavar="DATE",
        help="Sample collection date, one per --sample in the same order, ISO "
        "(YYYY-MM-DD) for the days-post-transplant derivation and timeline "
        "x-axis.",
    )
    group.add_argument(
        "--report-timestamp",
        help="Override the report generation time shown in the header / footer "
        "(any string). Defaults to the current local time. Pin it for "
        "byte-reproducible report output.",
    )


def _build_report_meta(args: argparse.Namespace) -> ReportMeta:
    """Assemble a ReportMeta from the optional CLI metadata flags."""
    rels = args.donor_relationship or []
    donors = [
        DonorMeta(donor_id=d, relationship=rels[i] if i < len(rels) else None)
        for i, d in enumerate(args.donor_sample)
    ]
    dates = args.sample_date or []
    sample_dates = {s: dates[i] for i, s in enumerate(args.sample) if i < len(dates)}
    return ReportMeta(
        recipient_id=args.recipient_id,
        recipient_name=args.recipient_name,
        sex=args.recipient_sex,
        dob=args.recipient_dob,
        transplant_type=args.transplant_type,
        transplant_date=args.transplant_date,
        donors=donors,
        sample_dates=sample_dates,
    )


def _build_report_params(args: argparse.Namespace) -> dict:
    """Collect the analysis parameters shown in the HTML report footer."""
    return {
        "panel_vcf": args.panel_vcf,
        "admix_vcf": args.admix_vcf,
        "error_table": args.error_table,
        "bias_table": args.bias_table,
        "contamination_table": args.contamination_table,
        "min_dp": args.min_dp,
        "min_gq": args.min_gq,
        "error_rate": args.error_rate,
        "robust": args.robust,
        "robust_k": args.robust_k,
        "marker_type_overdispersion": args.marker_type_overdispersion,
        "no_bias_correction": args.no_bias_correction,
        "estimate_bias": args.estimate_bias,
        "no_error_correction": args.no_error_correction,
        "contamination_correction": args.contamination_correction,
        "host_presence": not args.no_host_presence,
        "artifact_filter": not args.no_artifact_filter,
        "use_sex_chroms": args.use_sex_chroms,
    }


def _report_timestamp(args: argparse.Namespace) -> str:
    """Report generation time: the --report-timestamp override, else local now.

    Pinning the override makes the rendered report byte-reproducible (the only
    other nondeterminism, the trend-chart PNG, is fixed for fixed input data).
    """
    override = getattr(args, "report_timestamp", None)
    if override:
        return override
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def _validate_expected_relatedness(args: argparse.Namespace) -> None:
    """Check --expected-relatedness count matches the number of donors.

    Fails early with a clear message rather than as a later strict-zip error in QC.
    """
    er = args.expected_relatedness
    if er is not None and len(er) != len(args.donor_sample):
        raise SystemExit(
            f"--expected-relatedness given {len(er)} value(s) but there are "
            f"{len(args.donor_sample)} donor(s); provide exactly one per "
            "--donor-sample, in the same order (use NA for no expectation)"
        )


def _validate_sample_names(vcf_path: str, required: list[str]) -> None:
    """Check that all required sample names exist in the VCF header."""
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
    marker_type_overdispersion: bool = True,
    expected_relatedness: list[str] | None = None,
    relatedness_tolerance: int = 1,
    run_unit: RunUnitInfo | None = None,
) -> tuple:
    """Run the chimerism pipeline for one admixture sample.

    Takes pre-parsed host and donor markers to avoid redundant VCF reads, then
    delegates to ``allomix.analysis.analyse_sample`` (shared with the diagnostic
    scripts). Multi-donor estimation kicks in when more than one donor is given.

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
        marker_type_overdispersion=marker_type_overdispersion,
        expected_relatedness=expected_relatedness,
        relatedness_tolerance=relatedness_tolerance,
        run_unit=run_unit,
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


def _add_output_args(parser: argparse.ArgumentParser, *, allow_tsv: bool) -> None:
    """Add the per-artifact output flags shared by monitor and timeline.

    Any combination may be given in one run (e.g. ``--json r.json --html r.html``);
    ``-`` writes to stdout. With no output flag the command falls back to its
    default (TSV for monitor, JSON for timeline) on stdout.
    """
    parser.add_argument(
        "--json",
        metavar="PATH",
        help="Write the structured report JSON (the canonical artifact the HTML "
        "report is rendered from) to PATH ('-' for stdout).",
    )
    parser.add_argument(
        "--html",
        metavar="PATH",
        help="Write the self-contained HTML report to PATH ('-' for stdout).",
    )
    parser.add_argument(
        "--template",
        metavar="DIR",
        help="Directory of HTML report template overrides, searched ahead of the "
        "built-in templates. Drop in a styles.css to restyle, or a report.html / "
        "timeline.html / macros.html to restructure; anything absent falls back to "
        "the built-in. See src/allomix/html/templates/ for the files to override.",
    )
    if allow_tsv:
        parser.add_argument(
            "--tsv",
            metavar="PATH",
            help="Write the TSV summary to PATH ('-' for stdout). The default "
            "when no output flag is given.",
        )


def _write_text(path: str, text: str) -> None:
    """Write text to a file, or to stdout when ``path`` is ``-``."""
    if path == "-":
        sys.stdout.write(text)
    else:
        Path(path).write_text(text, encoding="utf-8")


def _load_calibration(args: argparse.Namespace) -> PanelCalibration:
    """Build the per-marker calibration from the CLI table options.

    Per-marker bias only helps markers that are informative for this run, and
    those are hom in both contributors, so they cannot be measured inline from one
    host/donor pair. Build a reusable table ahead of time with
    ``allomix estimate-bias`` (from a reference cohort, or ``--both-het`` across a
    patient cohort) and pass it with ``--bias-table``.
    """
    estimate_bias = getattr(args, "estimate_bias", False)
    if estimate_bias and args.bias_table:
        raise SystemExit("Use either --bias-table or --estimate-bias, not both")

    if estimate_bias and not args.no_bias_correction:
        samples = list(VCF(args.panel_vcf).samples)
        marker_lists = [parse_vcf(args.panel_vcf, sample=s, min_dp=0, min_gq=0) for s in samples]
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
    contamination_correction = None
    if getattr(args, "contamination_correction", False):
        if not args.contamination_table:
            raise SystemExit("--contamination-correction requires --contamination-table")
        contamination_correction = load_contamination_table(args.contamination_table)
    return PanelCalibration(
        biases=biases,
        errors=errors,
        contamination_correction=contamination_correction,
    )


def cmd_monitor(args: argparse.Namespace) -> int:
    """Run the monitor subcommand."""
    _validate_expected_relatedness(args)
    _validate_sample_names(args.panel_vcf, [args.host_sample] + args.donor_sample)
    _validate_sample_names(args.admix_vcf, args.sample)

    # No output flag defaults to TSV on stdout (the historical default).
    want_json = args.json is not None
    want_html = args.html is not None
    want_tsv = args.tsv is not None
    want_csv = args.marker_csv is not None
    if not (want_json or want_html or want_tsv or want_csv):
        want_tsv = True
        args.tsv = "-"

    # The single-sample report covers one sample; multiple timepoints belong in
    # the timeline report (which adds the trend chart).
    if (want_json or want_html) and len(args.sample) != 1:
        raise SystemExit(
            "monitor --json/--html produce one report for a single --sample; got "
            f"{len(args.sample)}. Use 'allomix timeline' for multiple timepoints."
        )

    calibration = _load_calibration(args)

    # Parse host and donors once (same for every timepoint).
    # gt_ad_consistency=True is the panel-side miscall guard: drops markers where
    # the called GT contradicts the AD VAF (e.g. GATK calling het from 20%-VAF
    # reads in a 2-sample joint call). Without it the wider gnomAD-derived panel
    # recovers markers that bias the estimator toward false host signal (Step 23).
    host = parse_vcf(
        args.panel_vcf, sample=args.host_sample, min_gq=args.min_gq, gt_ad_consistency=True
    )
    donors = [
        parse_vcf(args.panel_vcf, sample=d, min_gq=args.min_gq, gt_ad_consistency=True)
        for d in args.donor_sample
    ]

    # Run-unit metadata stamped on the admix VCF header by the pipeline
    # (index-hopping check, issue #12). Empty when the VCF carries none.
    run_units = read_run_units(args.admix_vcf)

    results: list[tuple[str, object, object]] = []
    marker_rows: list[tuple[str, object]] = []
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
            marker_type_overdispersion=args.marker_type_overdispersion,
            expected_relatedness=args.expected_relatedness,
            relatedness_tolerance=args.relatedness_tolerance,
            run_unit=run_units.get(sample_name),
        )
        results.append((genotypes.sample_name, result, qc))
        marker_rows.append((genotypes.sample_name, result))

    # One envelope feeds both the JSON and the HTML so they always agree.
    if want_json or want_html:
        name, result, qc = results[0]
        data = report_data(
            result,
            qc,
            sample_name=name,
            meta=_build_report_meta(args),
            params=_build_report_params(args),
            timestamp=_report_timestamp(args),
        )
        if want_json:
            _write_text(args.json, json.dumps(data, indent=2) + "\n")
        if want_html:
            _write_text(args.html, render_single(data, template_dir=args.template))

    if want_tsv:
        out = _open_output(args.tsv)
        try:
            for name, result, qc in results:
                to_tsv(result, qc, out, verbose=args.verbose, sample_name=name)
        finally:
            if out is not sys.stdout:
                out.close()

    if want_csv:
        to_marker_csv(marker_rows, args.marker_csv)

    return 0


def cmd_timeline(args: argparse.Namespace) -> int:
    """Run the timeline subcommand."""
    want_json = args.json is not None
    want_html = args.html is not None
    want_csv = args.marker_csv is not None
    if not (want_json or want_html or want_csv):
        want_json = True
        args.json = "-"

    if want_html and importlib.util.find_spec("matplotlib") is None:
        raise SystemExit(
            "timeline --html needs matplotlib for the trend chart: pip install 'allomix[report]'"
        )
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

    run_units = read_run_units(args.admix_vcf)

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
            marker_type_overdispersion=args.marker_type_overdispersion,
            expected_relatedness=args.expected_relatedness,
            relatedness_tolerance=args.relatedness_tolerance,
            run_unit=run_units.get(sample_name),
        )
        results.append((genotypes.sample_name, result, qc))

    if want_json or want_html:
        data = timeline_report_data(
            results,
            meta=_build_report_meta(args),
            params=_build_report_params(args),
            timestamp=_report_timestamp(args),
        )
        if want_json:
            _write_text(args.json, json.dumps(data, indent=2) + "\n")
        if want_html:
            # Deferred import: timeline.py imports matplotlib (the optional
            # `report` extra) at module top, so importing here, after the
            # capability check, keeps the json path working on base deps alone.
            from allomix.html.timeline import render_timeline

            _write_text(
                args.html,
                render_timeline(data, log_scale=args.log_scale, template_dir=args.template),
            )

    if want_csv:
        to_marker_csv([(name, result) for name, result, _ in results], args.marker_csv)

    return 0


def cmd_report(args: argparse.Namespace) -> int:
    """Render an HTML report from a saved report-data JSON file.

    Reads the envelope written by ``monitor --json`` / ``timeline --json`` and
    renders the same HTML, so report generation can be split from analysis (or
    re-run later). The ``kind`` field selects single-sample vs timeline.
    """
    raw = sys.stdin.read() if args.input == "-" else Path(args.input).read_text(encoding="utf-8")
    data = json.loads(raw)

    is_timeline = data.get("kind") == "timeline" or "timepoints" in data
    if is_timeline:
        if importlib.util.find_spec("matplotlib") is None:
            raise SystemExit(
                "rendering a timeline report needs matplotlib for the trend "
                "chart: pip install 'allomix[report]'"
            )
        # Deferred import: keeps the single-sample path matplotlib-free.
        from allomix.html.timeline import render_timeline

        html = render_timeline(data, log_scale=args.log_scale, template_dir=args.template)
    else:
        html = render_single(data, template_dir=args.template)

    _write_text(args.output, html)
    return 0


def cmd_estimate_bias(args: argparse.Namespace) -> int:
    """Run the estimate-bias subcommand."""
    if args.both_het:
        return _cmd_estimate_bias_both_het(args)

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


def _cmd_estimate_bias_both_het(args: argparse.Namespace) -> int:
    """Build a pooled bias table from admix samples at both-het markers.

    A marker heterozygous in both host and every donor has true admix VAF 0.5
    regardless of mixing fraction, so the admix AD there gives the per-marker bias
    directly, from the same caller as the admix (issue #11). Such markers are
    non-informative for that same host/donor pair, so this is a cohort table
    builder: run it across patients and apply the table (``--bias-table``) to
    other patients, whose informative markers it covers.
    """
    if args.vcfs or args.samples:
        raise SystemExit(
            "--both-het uses --vcf with --host-sample/--donor-sample, not --vcfs/--samples"
        )
    if not args.vcf or not args.host_sample or not args.donor_sample or not args.admix_vcfs:
        raise SystemExit(
            "--both-het requires --vcf, --host-sample, --donor-sample, and --admix-vcfs"
        )

    _validate_sample_names(args.vcf, [args.host_sample, *args.donor_sample])
    host = parse_vcf(args.vcf, sample=args.host_sample, min_dp=0, min_gq=0)
    donors = [parse_vcf(args.vcf, sample=d, min_dp=0, min_gq=0) for d in args.donor_sample]

    admix_lists = []
    for vcf_path in args.admix_vcfs:
        for sample in VCF(vcf_path).samples:
            admix_lists.append(parse_vcf(vcf_path, sample=sample, min_dp=0, min_gq=0))

    biases = estimate_biases_both_het(host, donors, admix_lists, min_het=args.min_het)
    save_bias_table(biases, args.output)
    print(
        f"Estimated both-het bias for {len(biases)} markers from "
        f"{len(admix_lists)} admix sample(s) in {len(args.admix_vcfs)} VCF(s) "
        f"-> {args.output}",
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

    if args.homref_vcf and not args.samples:
        raise SystemExit("--homref-vcf requires --samples (the shared sample names)")

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

    # Hom-ref background VCFs (e.g. raw mpileup at force-called amplicon
    # midpoints) carry the same samples; their hom-ref calls add the ref->alt
    # observations the variant-only joint call cannot supply, pooled through the
    # same estimator with the panel hom-ref sites.
    for homref_path in args.homref_vcf:
        _validate_sample_names(homref_path, args.samples)
        for sample in args.samples:
            markers = parse_vcf(homref_path, sample=sample, min_dp=0, min_gq=args.min_gq)
            marker_lists.append(markers)
        n_source += f" + {len(args.samples)} hom-ref samples from {homref_path}"

    errors = estimate_error_rates(
        marker_lists,
        min_reads=args.min_reads,
        max_vaf_homref=args.max_vaf_homref,
        min_vaf_homalt=args.min_vaf_homalt,
    )
    save_error_table(errors, args.output)
    print(
        f"Estimated error rates for {len(errors)} sites from {n_source} -> {args.output}",
        file=sys.stderr,
    )
    return 0


def cmd_build_contamination_table(args: argparse.Namespace) -> int:
    """Run the build-contamination-table subcommand (Step 30, issue #30).

    Builds a per-marker co-pooled contamination correction for one patient on a
    flowcell: carrier counts from the cohort genotypes, a per-flowcell gate from
    the consensus-hom dose response, and a correction slope calibrated on the
    patient's informative donor-hom markers pooled across the supplied timepoints.
    """
    _validate_sample_names(args.vcf, [args.host_sample, *args.donor_sample])
    _validate_sample_names(args.admix_vcf, args.sample)

    # Same panel-side miscall guard (gt_ad_consistency) the analysis uses.
    host = parse_vcf(args.vcf, sample=args.host_sample, min_gq=args.min_gq, gt_ad_consistency=True)
    donors = [
        parse_vcf(args.vcf, sample=d, min_gq=args.min_gq, gt_ad_consistency=True)
        for d in args.donor_sample
    ]

    admix_lists = [parse_vcf(args.admix_vcf, sample=s, min_dp=0) for s in args.sample]

    # Carrier counts come from the co-pooled flowcell individuals; default to
    # every sample in --vcf (which should include host and donors).
    cohort_samples = args.cohort_samples or list(VCF(args.vcf).samples)
    cohort_genotypes = [
        parse_vcf(args.vcf, sample=s, min_gq=args.min_gq, gt_ad_consistency=True)
        for s in cohort_samples
    ]

    correction = estimate_contamination_table(
        host,
        donors,
        admix_lists,
        cohort_genotypes,
        min_dp=args.min_dp,
        min_gq=args.min_gq,
        dose_cap=args.dose_cap,
        alpha=args.alpha,
        min_slope=args.min_slope,
    )
    save_contamination_table(correction, args.output)

    verdict = "gated IN (correcting)" if correction.gated else "gated OUT (no-op)"
    print(
        f"Built contamination table from {len(args.sample)} timepoint(s), "
        f"{len(cohort_samples)} cohort sample(s), {len(correction.carriers)} markers: "
        f"{verdict}; consensus slope {correction.gate_slope * 100:.4f}%/carrier "
        f"(p={correction.gate_p_value:.2e}), correction slope "
        f"{correction.slope * 100:.4f}%/carrier -> {args.output}",
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
    _add_report_meta_args(monitor_parser)
    _add_output_args(monitor_parser, allow_tsv=True)

    timeline_parser = subparsers.add_parser(
        "timeline",
        help="Generate chimerism timeline across timepoints",
    )
    _add_common_args(timeline_parser)
    _add_report_meta_args(timeline_parser)
    _add_output_args(timeline_parser, allow_tsv=False)
    timeline_parser.add_argument(
        "--log-scale",
        action="store_true",
        help="Use a logarithmic y axis on the HTML trend chart.",
    )

    report_parser = subparsers.add_parser(
        "report",
        help="Render an HTML report from a saved monitor/timeline JSON",
    )
    report_parser.add_argument(
        "input",
        help="Report-data JSON from 'monitor --json' or 'timeline --json' ('-' to read stdin)",
    )
    report_parser.add_argument(
        "--output",
        "-o",
        default="-",
        help="HTML output file (default: stdout). A timeline JSON needs the "
        "'report' extra for the trend chart: pip install 'allomix[report]'.",
    )
    report_parser.add_argument(
        "--log-scale",
        action="store_true",
        help="Use a logarithmic y axis on the HTML trend chart (timeline JSON).",
    )
    report_parser.add_argument(
        "--template",
        metavar="DIR",
        help="Directory of HTML report template overrides, searched ahead of the "
        "built-in templates (see src/allomix/html/templates/).",
    )

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
        help="Sample names to extract from --vcf (het-site mode)",
    )
    bias_parser.add_argument(
        "--both-het",
        action="store_true",
        help="Both-het mode: estimate bias from admix samples at markers where "
        "the host and every donor are heterozygous (true VAF 0.5 regardless "
        "of mixing). Use this when you only have admix VCFs (no mpileup'd "
        "reference samples) and need a caller-consistent table. Requires "
        "--vcf (genotypes), --host-sample, --donor-sample, and --admix-vcfs. "
        "A pair's both-het markers are non-informative for that same pair, "
        "so build the table from a cohort and apply it to other patients.",
    )
    bias_parser.add_argument(
        "--host-sample",
        help="Host sample name in --vcf (both-het mode)",
    )
    bias_parser.add_argument(
        "--donor-sample",
        action="append",
        default=[],
        metavar="SAMPLE_NAME",
        help="Donor sample name in --vcf (both-het mode; repeat for multi-donor)",
    )
    bias_parser.add_argument(
        "--admix-vcfs",
        nargs="+",
        metavar="VCF",
        help="Admix VCFs whose samples supply the both-het observations "
        "(both-het mode); all samples in each are pooled",
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
        help="Sample names to extract from --vcf (and from every --homref-vcf)",
    )
    err_parser.add_argument(
        "--homref-vcf",
        nargs="+",
        default=[],
        metavar="VCF",
        help="VCF(s) of force-called hom-ref background positions for the SAME "
        "samples (raw bcftools mpileup, e.g. at amplicon midpoints, NOT "
        "GATK, which strips the minority ALT reads being measured). The "
        "stray ALT reads at these hom-ref calls supply the ref->alt error "
        "background, which a variant-only joint call cannot (it emits no "
        "all-hom-ref sites). Folded into the same estimate. Requires "
        "--samples and uses the same sample names as --vcf/--vcfs.",
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
        help="Minimum total reads per direction to retain a site's estimate (default: 1000)",
    )
    err_parser.add_argument(
        "--max-vaf-homref",
        type=float,
        default=0.10,
        help="Drop hom-ref training observations with vaf > this (default: 0.10)",
    )
    err_parser.add_argument(
        "--min-vaf-homalt",
        type=float,
        default=0.90,
        help="Drop hom-alt training observations with vaf < this (default: 0.90)",
    )
    err_parser.add_argument(
        "--min-gq",
        type=int,
        default=DEFAULT_MIN_GQ,
        help=f"Minimum GQ for training calls (default: {DEFAULT_MIN_GQ})",
    )

    contam_parser = subparsers.add_parser(
        "build-contamination-table",
        help="Build a per-marker co-pooled contamination table (Step 30, issue #30)",
    )
    contam_parser.add_argument(
        "--vcf",
        required=True,
        metavar="VCF",
        help="Joint-called genotype VCF with the host, donor(s), and co-pooled "
        "cohort samples (the flowcell's joint call)",
    )
    contam_parser.add_argument(
        "--host-sample",
        required=True,
        help="Host sample name in --vcf",
    )
    contam_parser.add_argument(
        "--donor-sample",
        action="append",
        default=[],
        required=True,
        metavar="SAMPLE_NAME",
        help="Donor sample name in --vcf (repeat for multi-donor)",
    )
    contam_parser.add_argument(
        "--admix-vcf",
        required=True,
        metavar="VCF",
        help="Admix VCF holding this patient's serial timepoints (forced pileup AD)",
    )
    contam_parser.add_argument(
        "--sample",
        nargs="+",
        required=True,
        metavar="SAMPLE_NAME",
        help="Admix timepoint sample names in --admix-vcf to pool over",
    )
    contam_parser.add_argument(
        "--cohort-samples",
        nargs="+",
        default=None,
        metavar="SAMPLE_NAME",
        help="Co-pooled cohort sample names in --vcf for carrier counts "
        "(default: all samples in --vcf)",
    )
    contam_parser.add_argument(
        "--output",
        "-o",
        default="contamination_table.tsv",
        help="Output contamination table TSV (default: contamination_table.tsv)",
    )
    contam_parser.add_argument(
        "--dose-cap",
        type=int,
        default=DEFAULT_DOSE_CAP,
        help=f"Maximum co-pooled carrier dose (default: {DEFAULT_DOSE_CAP})",
    )
    contam_parser.add_argument(
        "--alpha",
        type=float,
        default=DEFAULT_GATE_ALPHA,
        help=f"Significance level for the consensus-hom slope gate (default: {DEFAULT_GATE_ALPHA})",
    )
    contam_parser.add_argument(
        "--min-slope",
        type=float,
        default=DEFAULT_GATE_MIN_SLOPE,
        help="Minimum consensus-hom slope per carrier worth correcting "
        f"(fraction of depth; default: {DEFAULT_GATE_MIN_SLOPE})",
    )
    contam_parser.add_argument(
        "--min-dp",
        type=int,
        default=DEFAULT_MIN_DP,
        help=f"Minimum admix depth per marker (default: {DEFAULT_MIN_DP})",
    )
    contam_parser.add_argument(
        "--min-gq",
        type=int,
        default=DEFAULT_MIN_GQ,
        help=f"Minimum host/donor GQ (default: {DEFAULT_MIN_GQ})",
    )

    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 1

    if args.command == "monitor":
        return cmd_monitor(args)
    if args.command == "timeline":
        return cmd_timeline(args)
    if args.command == "report":
        return cmd_report(args)
    if args.command == "estimate-bias":
        return cmd_estimate_bias(args)
    if args.command == "estimate-errors":
        return cmd_estimate_errors(args)
    if args.command == "build-contamination-table":
        return cmd_build_contamination_table(args)
    parser.print_help()
    return 1


if __name__ == "__main__":
    sys.exit(main())
