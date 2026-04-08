#!/usr/bin/env python3
"""Verify allomix Pearson overdispersion against statsmodels GLM.

Computes the Pearson dispersion estimate two ways:
  1. allomix's _compute_overdispersion (direct formula)
  2. statsmodels GLM with Binomial family (pearson_chi2 / df_resid)

and confirms they agree numerically.

This is an offline verification script, not a runtime dependency.
statsmodels is NOT required by allomix — only by this script.

Usage:
    pip install statsmodels
    python scripts/verify_overdispersion.py
"""

from __future__ import annotations

import logging
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import numpy as np
import statsmodels.api as sm
from scipy.special import logit

from allomix.chimerism import (
    _compute_overdispersion,
    estimate_single_donor,
    expected_weight,
)
from allomix.genotype import classify_markers, parse_vcf
from allomix.simulate import blend_vcfs, write_vcf

log = logging.getLogger(__name__)


def allomix_overdispersion(markers, f_donor, error_rate):
    """Extract the Pearson chi-squared components for comparison."""
    pearson_residuals_sq = []
    for m in markers:
        dp = m.admix_ad_ref + m.admix_ad_alt
        if dp == 0:
            continue
        w = expected_weight(m.host_gt, m.donor_gts[0], f_donor)
        e = error_rate
        p_alt = (1.0 - w) * (1.0 - e) + w * e / 3.0
        var_vaf = p_alt * (1.0 - p_alt) / dp
        if var_vaf < 1e-12:
            continue
        observed_alt_frac = m.admix_ad_alt / dp
        residual = observed_alt_frac - p_alt
        pearson_residuals_sq.append(residual**2 / var_vaf)
    return pearson_residuals_sq


def statsmodels_dispersion(markers, f_donor, error_rate):
    """Compute dispersion via statsmodels GLM with Binomial family.

    We fit a Binomial GLM with a single offset term (the logit of the
    expected probability from our mixture model) and no free parameters,
    then extract pearson_chi2 / df_resid.
    """
    successes = []
    trials = []
    offsets = []

    for m in markers:
        dp = m.admix_ad_ref + m.admix_ad_alt
        if dp == 0:
            continue
        w = expected_weight(m.host_gt, m.donor_gts[0], f_donor)
        e = error_rate
        p_alt = (1.0 - w) * (1.0 - e) + w * e / 3.0
        p_alt = max(1e-10, min(1.0 - 1e-10, p_alt))

        successes.append(m.admix_ad_alt)
        trials.append(dp)
        offsets.append(logit(p_alt))

    successes = np.array(successes, dtype=float)
    trials = np.array(trials, dtype=float)
    offsets = np.array(offsets, dtype=float)

    # Binomial GLM: endog = successes/trials, offset = logit(expected_p)
    # No covariates (intercept-free) — the offset fully specifies the expected value
    endog = np.column_stack([successes, trials - successes])
    exog = np.zeros((len(successes), 0))  # no covariates

    # Use GLM with offset only — fit with 0 free parameters
    model = sm.GLM(
        endog,
        np.ones((len(successes), 1)),  # dummy intercept column
        family=sm.families.Binomial(),
        offset=offsets,
    )
    # Fix intercept to 0 so offset alone determines expected values
    # Use start_params=[0] and maxiter=0 to avoid re-estimation
    result = model.fit(start_params=[0.0], maxiter=0, disp=False)

    return result.pearson_chi2, result.df_resid, result.pearson_chi2 / result.df_resid


def main() -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    log.info("Verifying allomix overdispersion against statsmodels GLM...")

    host_vcf = "tests/test_data/host.vcf"
    donor_vcf = "tests/test_data/donor.vcf"

    # Test at several donor fractions with realistic noise
    fractions = [0.05, 0.20, 0.50, 0.80, 0.95]
    all_match = True

    for frac in fractions:
        # Generate synthetic chimeric VCF
        result = blend_vcfs(
            host_path=host_vcf,
            donor_path=donor_vcf,
            donor_fraction=frac,
            target_depth=500,
            sample_name=f"test_{int(frac*100)}",
            seed=42,
            depth_cv=0.43,
            locus_dropout_rate=0.016,
        )
        tmp_path = Path(f"/tmp/verify_od_{int(frac*100)}.vcf")
        write_vcf(result, tmp_path)

        host = parse_vcf(host_vcf, min_dp=0, min_gq=0)
        donor = parse_vcf(donor_vcf, min_dp=0, min_gq=0)
        admix = parse_vcf(str(tmp_path), min_dp=0, min_gq=0)
        genotypes = classify_markers(host, [donor], admix, min_dp=0, min_gq=0, pass_only=False)

        est = estimate_single_donor(genotypes.informative, error_rate=0.01)
        f_mle = est.donor_fraction

        # Method 1: allomix
        phi_allomix = _compute_overdispersion(genotypes.informative, f_mle, 0.01)
        resid_sq = allomix_overdispersion(genotypes.informative, f_mle, 0.01)
        chi2_allomix = sum(resid_sq)
        df_allomix = len(resid_sq) - 1

        # Method 2: statsmodels GLM
        chi2_sm, df_sm, phi_sm = statsmodels_dispersion(genotypes.informative, f_mle, 0.01)

        # Compare
        chi2_match = abs(chi2_allomix - chi2_sm) < 0.01
        df_match = df_allomix == df_sm
        phi_match = abs(phi_allomix - phi_sm) < 0.001

        status = "PASS" if (chi2_match and df_match and phi_match) else "FAIL"
        if status == "FAIL":
            all_match = False

        log.info("  f_true=%.2f  f_mle=%.4f  n_markers=%d", frac, f_mle, len(resid_sq))
        log.info(
            "    allomix:     chi2=%.4f  df=%d  phi=%.4f",
            chi2_allomix, df_allomix, phi_allomix,
        )
        log.info("    statsmodels: chi2=%.4f  df=%d  phi=%.4f", chi2_sm, int(df_sm), phi_sm)
        log.info("    [%s]", status)

    if all_match:
        log.info("All fractions match — allomix overdispersion agrees with statsmodels GLM.")
        return 0
    else:
        log.error("MISMATCH detected — investigate differences above.")
        return 1


if __name__ == "__main__":
    sys.exit(main())
