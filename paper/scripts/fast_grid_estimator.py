"""Fast vectorized grid single-donor estimator (paper-only, opt-in).

This is a speed optimisation used only by the paper's LoD sweep
(``run_lod_validation.py``). It lives outside the ``allomix`` package on purpose:
it is not part of the clinical tool and the core package should stay small enough
to scrutinise thoroughly. It reuses the building blocks of the exact estimator in
``allomix.estimate.chimerism`` (some of them private), so a refactor of those internals can
break this helper; that is an accepted cost of keeping it out of core.

The exact estimator (``allomix.estimate.chimerism.estimate_single_donor_bb``) profiles rho
with ``minimize_scalar`` at every grid f and then runs a joint Nelder-Mead refine,
so each call makes hundreds of scalar scipy optimisations. For the LoD sweep we
only need ``donor_fraction``, and that is recovered to well within 1e-4 by
maximising the same beta-binomial log-likelihood on a 2-D (f, rho) grid, then a
single 1-D bounded f-search bracketed around the grid argmax (with rho profiled
out). This path replaces the per-f scipy calls of the grid search with one
vectorized array pass and a short local refine.

Accuracy is pinned against the exact estimator in
``tests/test_fast_grid_estimator.py``: across the fraction, depth, bias,
error-table and two-rho cases the fraction agrees to < 1e-4 (so the resulting
LoD-summary percentages stay well under 0.01 pp). The exact estimator remains
the default; this is opt-in only.
"""

import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import expit, gammaln, logit
from scipy.stats import chi2

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent / "src"))

from allomix.estimate.chimerism import (  # noqa: E402
    _RHO_MAX,
    MIN_CLASS_MARKERS,
    _donor_het_mask,
    _profile_rho_at_f,
    _two_rho_profile_ll,
)
from allomix.constants import (  # noqa: E402
    CI_LEVEL,
    DEFAULT_ERROR_RATE,
    N_OTHER_BASES,
)
from allomix.genotype import InformativeMarker  # noqa: E402
from allomix.estimate.likelihood import (  # noqa: E402
    W_EPS,
    PanelCalibration,
    _ll_from_p_alt,
    _MarkerArrays,
    _p_alt_for_f,
    _precompute_marker_arrays,
)


@dataclass
class GridChimerismResult:
    """Lightweight result of the fast grid single-donor estimator.

    Carries only what the LoD sweep consumes (``donor_fraction``) plus a few
    fields so callers that read ``ChimerismResult`` attributes still work. CI is
    a coarse profile-likelihood bracket on the grid; ``detection_limit`` style
    fields are not computed.
    """

    donor_fraction: float
    donor_fraction_ci: tuple[float, float]
    host_fraction: float
    log_likelihood: float
    n_informative: int
    n_markers_used: int
    error_rate: float
    rho: float = float("inf")
    # Per-marker-type overdispersion (issue #33). None unless the two-rho mode ran
    # with both classes above MIN_CLASS_MARKERS; then rho is the het-class rho.
    rho_hom: float | None = None
    rho_het: float | None = None
    marker_type_overdispersion_fallback: str | None = None


def _p_alt_grid(
    arr: _MarkerArrays,
    f_grid: np.ndarray,
    error_rate: float = DEFAULT_ERROR_RATE,
) -> np.ndarray:
    """Per-marker P(observe ALT) for a whole f-grid at once.

    Vectorized form of ``_p_alt_for_f`` over an ``(n_f,)`` fraction grid; returns
    an ``(n_f, M)`` array whose row i is exactly ``_p_alt_for_f(arr, f_grid[i])``.
    """
    f = f_grid[:, None]  # (n_f, 1)
    host = arr.host_ref_dose[None, :]  # (1, M)
    donor = arr.donor_ref_dose[None, :]
    w = (1.0 - f) * host / 2.0 + f * donor / 2.0  # (n_f, M)

    if arr.has_bias:
        # Inlined apply_bias with the bias-only logit(p) term precomputed.
        bm = arr.bias_mask
        wm = np.clip(w[:, bm], W_EPS, 1.0 - W_EPS)
        w[:, bm] = np.clip(
            expit(logit(wm) - arr.logit_bias_masked[None, :]), W_EPS, 1.0 - W_EPS
        )

    e = error_rate
    e_specific = e / N_OTHER_BASES
    p_alt_raw = (1.0 - w) * (1.0 - e) + w * e_specific
    p_ref_raw = w * (1.0 - e) + (1.0 - w) * e_specific
    p_alt = p_alt_raw / (p_ref_raw + p_alt_raw)
    if arr.has_error:
        p_alt_asym = w * arr.e_refalt[None, :] + (1.0 - w) * (1.0 - arr.e_altref[None, :])
        p_alt = np.where(arr.error_mask[None, :], p_alt_asym, p_alt)
    return np.clip(p_alt, 1e-6, 1.0 - 1e-6)


def _ll_grid_over_rho(
    arr: _MarkerArrays,
    p_alt: np.ndarray,
    rho_grid: np.ndarray,
) -> np.ndarray:
    """Total log-likelihood on the full ``(n_f, n_rho)`` grid.

    ``p_alt`` is the ``(n_f, M)`` array from ``_p_alt_grid``. Loops over the rho-grid
    (cheap, ``n_rho`` is small) so the largest temporary is one ``(n_f, M)`` block per
    rho, bounding memory. Each rho column is the vectorized counterpart of
    ``_ll_from_p_alt`` evaluated for every f at once. Returns an ``(n_f, n_rho)`` array.
    """
    n, k = arr.n[None, :], arr.k[None, :]
    out = np.empty((p_alt.shape[0], rho_grid.shape[0]), dtype=float)
    for j, rho in enumerate(rho_grid):
        a = np.maximum(p_alt * rho, 1e-10)
        b = np.maximum((1.0 - p_alt) * rho, 1e-10)
        ll = (
            gammaln(k + a)
            + gammaln(n - k + b)
            - gammaln(n + rho)
            - gammaln(a)
            - gammaln(b)
            + math.lgamma(rho)
        )
        out[:, j] = ll.sum(axis=1)
    return out


# rho range for the fast grid and its profile. The lower bound matches the exact
# estimator's grid rho-profiling; the upper bound is _RHO_MAX (50000), the same
# ceiling the exact estimator's Nelder-Mead refine allows. Capping the fast
# profile at the lower 10000 (the exact estimator's *grid* bound) leaves a ~0.01
# pp bias at high-concentration samples where the optimum rho sits above 10000;
# matching _RHO_MAX removes it.
_GRID_RHO_LO = 1.0


def estimate_single_donor_bb_grid(
    markers: list[InformativeMarker],
    error_rate: float = DEFAULT_ERROR_RATE,
    calibration: PanelCalibration | None = None,
    n_f: int = 201,
    n_rho: int = 32,
    refine: bool = True,
    marker_type_overdispersion: bool = True,
) -> GridChimerismResult:
    """Fast approximate single-donor MLE via a vectorized (f, rho) grid.

    Opt-in alternative to ``estimate_single_donor_bb`` for parameter sweeps
    where only ``donor_fraction`` is needed. Builds an f-grid on ``[0, 1]`` and a
    log-spaced rho-grid on ``[1, _RHO_MAX]``, evaluates the beta-binomial total
    log-likelihood on the full grid in a few vectorized array passes, and takes
    the grid argmax. When ``refine`` (the default) it then polishes the donor
    fraction with a 1-D bounded search bracketed to +/- 2 f-grid steps around the
    argmax, profiling rho out at each candidate f (the same profile the exact
    estimator uses). Because the total log-likelihood is unimodal in f, this
    local solve lands on the exact estimator's fraction to < 1e-4 across the LoD
    space (pinned in ``tests/test_fast_grid_estimator.py``). The 1-D profiled
    refine is both faster and tighter than a joint Nelder-Mead polish, which can
    stall in the very flat rho direction.

    The exact estimator (``estimate_single_donor_bb``) stays the default and is
    untouched; this path is selected explicitly by the caller.

    Performance and accuracy (measured on the LoD sweep parameter space during
    development; the accuracy bound is pinned in tests/test_fast_grid_estimator.py):

      - About 6.5x faster per call than the exact estimator (~30 ms vs ~191 ms;
        4-7x depending on panel size, with the grid build dominating and the
        refine ~5 ms). The win comes from replacing the exact estimator's
        Python-level f-grid loop (a bounded rho profile per f point) and the
        joint Nelder-Mead refinement with a few vectorized array passes.
      - Donor-fraction agreement with the exact estimator: median 1e-6 pp,
        worst case 0.0115 pp over fractions <= 5% (the whole LoD-sweep regime).
        The only deviations above 0.01 pp are at f=0.5 (outside the LoD range),
        and there the grid finds the strictly higher likelihood, i.e. it is the
        more accurate of the two, not the reverse.
      - End-to-end effect on the reported LoD: across the full 60-cell
        ``lod_summary`` grid (relatedness x depth x panel size), the per-cell
        ``lod_pct`` matches the exact estimator to a median of 0.0000 pp and a
        maximum of 0.0011 pp, comfortably under a 0.01 pp tolerance.

    So this path is appropriate for fast iteration on the LoD sweeps; run the
    exact estimator (the default) for the final publication figures.

    Args:
        error_rate: Fallback when a marker lacks per-direction rates.
        n_f: Number of f-grid points on ``[0, 1]``.
        n_rho: Number of log-spaced rho-grid points on ``[1, _RHO_MAX]``.
        refine: Run the 1-D profiled local polish from the grid argmax
            (default True).
        marker_type_overdispersion: Profile a separate rho per marker class
            (donor-hom vs donor-het) on the grid and sum the per-class
            profiled-over-rho curves before the f argmax (issue #33). On by
            default (matches the exact estimator's default). Set False for the
            legacy single-rho grid. Falls back to the single-rho grid when a class
            has fewer than ``MIN_CLASS_MARKERS`` markers.

    Returns:
        ``GridChimerismResult`` with the donor-fraction estimate and a coarse CI.
    """
    cal = calibration or PanelCalibration()
    n_informative = len(markers)
    if n_informative == 0:
        return GridChimerismResult(
            donor_fraction=0.0,
            donor_fraction_ci=(0.0, 0.0),
            host_fraction=1.0,
            log_likelihood=0.0,
            n_informative=0,
            n_markers_used=0,
            error_rate=error_rate,
        )

    if marker_type_overdispersion:
        het_mask = _donor_het_mask(markers)
        n_het = int(het_mask.sum())
        n_hom = n_informative - n_het
        if n_het >= MIN_CLASS_MARKERS and n_hom >= MIN_CLASS_MARKERS:
            return _estimate_single_donor_bb_grid_two_rho(
                markers, het_mask, error_rate, cal, n_f, n_rho, refine
            )
        fell_back: str | None = (
            f"a marker class is sparse (hom={n_hom}, het={n_het}, "
            f"min={MIN_CLASS_MARKERS}); used shared rho for this sample"
        )
    else:
        fell_back = None

    arr = _precompute_marker_arrays(markers, cal)

    f_grid = np.linspace(0.0, 1.0, n_f)
    rho_grid = np.exp(np.linspace(math.log(_GRID_RHO_LO), math.log(_RHO_MAX), n_rho))

    p_alt = _p_alt_grid(arr, f_grid, error_rate)  # (n_f, M)
    ll = _ll_grid_over_rho(arr, p_alt, rho_grid)  # (n_f, n_rho)

    flat = int(np.argmax(ll))
    fi, ri = np.unravel_index(flat, ll.shape)
    best_f = float(f_grid[fi])
    best_rho = float(rho_grid[ri])
    best_ll = float(ll[fi, ri])

    def profile_ll_f(f_val: float) -> float:
        """Max LL over rho at a fixed f (rho profiled out, as in the exact MLE)."""
        p = _p_alt_for_f(arr, f_val, error_rate)
        opt_rho = minimize_scalar(
            lambda log_r, _p=p: -_ll_from_p_alt(arr, _p, math.exp(log_r)),
            bounds=(math.log(_GRID_RHO_LO), math.log(_RHO_MAX)),
            method="bounded",
        )
        return -float(opt_rho.fun)

    if refine:
        # LL is unimodal in f, so the optimum lies within one grid step of the argmax;
        # bracket +/- 2 steps for safety and solve f with rho profiled out.
        step = 1.0 / (n_f - 1)
        lo = max(0.0, best_f - 2.0 * step)
        hi = min(1.0, best_f + 2.0 * step)
        if hi > lo:
            opt = minimize_scalar(
                lambda f: -profile_ll_f(f),
                bounds=(lo, hi),
                method="bounded",
                options={"xatol": 1e-8},
            )
            f_ref = max(0.0, min(1.0, float(opt.x)))
            ll_ref = -float(opt.fun)
            if ll_ref >= best_ll:
                best_f = f_ref
                best_ll = ll_ref
                # Profile rho at the refined f for the reported concentration.
                p = _p_alt_for_f(arr, best_f, error_rate)
                opt_rho = minimize_scalar(
                    lambda log_r: -_ll_from_p_alt(arr, p, math.exp(log_r)),
                    bounds=(math.log(_GRID_RHO_LO), math.log(_RHO_MAX)),
                    method="bounded",
                )
                best_rho = math.exp(float(opt_rho.x))

    # Coarse profile-likelihood CI: profile rho out (max over the rho axis at each f),
    # then bracket where the profile drops by chi2(0.95, df=1)/2 from the maximum.
    prof = ll.max(axis=1)
    half_threshold = float(chi2.ppf(CI_LEVEL, df=1)) / 2.0
    above = prof >= (best_ll - half_threshold)
    idx = np.nonzero(above)[0]
    if idx.size:
        f_lo = float(f_grid[idx[0]])
        f_hi = float(f_grid[idx[-1]])
    else:
        f_lo = f_hi = best_f

    return GridChimerismResult(
        donor_fraction=best_f,
        donor_fraction_ci=(f_lo, f_hi),
        host_fraction=1.0 - best_f,
        log_likelihood=best_ll,
        n_informative=n_informative,
        n_markers_used=n_informative,
        error_rate=error_rate,
        rho=best_rho,
        marker_type_overdispersion_fallback=fell_back,
    )


def _estimate_single_donor_bb_grid_two_rho(
    markers: list[InformativeMarker],
    het_mask: np.ndarray,
    error_rate: float,
    cal: PanelCalibration,
    n_f: int,
    n_rho: int,
    refine: bool,
) -> GridChimerismResult:
    """Fast grid two-rho estimator (issue #33), the §5.9 counterpart of the core.

    Profiles each class's rho out on the same log-spaced rho-grid the single-rho
    fast path uses (``_ll_grid_over_rho`` -> max over the rho axis), sums the two
    per-class profiled-over-rho curves, and takes the f argmax. The optional
    refine and the reported per-class rhos reuse the exact estimator's 1-D rho
    profile (``_profile_rho_at_f`` / ``_two_rho_profile_ll`` from
    ``allomix.estimate.chimerism``), so the polished fraction matches the exact two-rho MLE.
    """
    hom_markers = [m for m, h in zip(markers, het_mask) if not h]
    het_markers = [m for m, h in zip(markers, het_mask) if h]
    arr_hom = _precompute_marker_arrays(hom_markers, cal)
    arr_het = _precompute_marker_arrays(het_markers, cal)

    f_grid = np.linspace(0.0, 1.0, n_f)
    rho_grid = np.exp(np.linspace(math.log(_GRID_RHO_LO), math.log(_RHO_MAX), n_rho))

    # Profile rho out of each class on the grid, then sum the per-class profiles.
    p_hom = _p_alt_grid(arr_hom, f_grid, error_rate)  # (n_f, M_hom)
    p_het = _p_alt_grid(arr_het, f_grid, error_rate)  # (n_f, M_het)
    ll_hom = _ll_grid_over_rho(arr_hom, p_hom, rho_grid)  # (n_f, n_rho)
    ll_het = _ll_grid_over_rho(arr_het, p_het, rho_grid)  # (n_f, n_rho)
    prof_hom = ll_hom.max(axis=1)  # (n_f,)
    prof_het = ll_het.max(axis=1)
    prof_total = prof_hom + prof_het  # rho profiled, per f

    fi = int(np.argmax(prof_total))
    best_f = float(f_grid[fi])
    best_ll = float(prof_total[fi])
    best_rho_hom = float(rho_grid[int(np.argmax(ll_hom[fi]))])
    best_rho_het = float(rho_grid[int(np.argmax(ll_het[fi]))])

    if refine:
        step = 1.0 / (n_f - 1)
        lo = max(0.0, best_f - 2.0 * step)
        hi = min(1.0, best_f + 2.0 * step)
        if hi > lo:
            opt = minimize_scalar(
                lambda f: -_two_rho_profile_ll(arr_hom, arr_het, f, error_rate),
                bounds=(lo, hi),
                method="bounded",
                options={"xatol": 1e-8},
            )
            f_ref = max(0.0, min(1.0, float(opt.x)))
            ll_ref = -float(opt.fun)
            if ll_ref >= best_ll:
                best_f = f_ref
                best_ll = ll_ref
                # Exact per-class rho profile at the refined f (mirrors the core).
                _, best_rho_hom = _profile_rho_at_f(arr_hom, best_f, error_rate)
                _, best_rho_het = _profile_rho_at_f(arr_het, best_f, error_rate)

    # Coarse CI from the summed profiled-over-rho curve (same bracket as single-rho).
    half_threshold = float(chi2.ppf(CI_LEVEL, df=1)) / 2.0
    above = prof_total >= (best_ll - half_threshold)
    idx = np.nonzero(above)[0]
    if idx.size:
        f_lo = float(f_grid[idx[0]])
        f_hi = float(f_grid[idx[-1]])
    else:
        f_lo = f_hi = best_f

    return GridChimerismResult(
        donor_fraction=best_f,
        donor_fraction_ci=(f_lo, f_hi),
        host_fraction=1.0 - best_f,
        log_likelihood=best_ll,
        n_informative=len(markers),
        n_markers_used=len(markers),
        error_rate=error_rate,
        rho=best_rho_het,  # headline = het class
        rho_hom=best_rho_hom,
        rho_het=best_rho_het,
    )
