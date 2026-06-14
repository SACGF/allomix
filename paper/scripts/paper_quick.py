"""Quick-build switch for the paper pipeline.

Set the environment variable ``ALLOMIX_PAPER_QUICK=1`` (or build with
``snakemake -s paper/Snakefile --config quick=1``) to run the heavy validation
scripts with reduced iteration counts (fewer pairs, replicates, depths, panel
sizes). This trades statistical precision for speed, so the figures are NOT
publication quality.

Two things are provided:

  - ``QUICK`` / ``qval(full, quick)``: scripts pick reduced loop counts.
  - A watermark: when ``QUICK`` is set, importing this module patches
    ``Figure.savefig`` so every figure the build writes is stamped with a
    "QUICK BUILD" banner. That makes a low-iteration figure impossible to
    mistake for a real one.

Every figure-producing paper script imports this module (``import paper_quick``)
so the watermark applies build-wide, regardless of which script drew the figure.
"""

import os

from matplotlib.figure import Figure

_TRUE = {"1", "true", "yes", "on"}

QUICK: bool = os.environ.get("ALLOMIX_PAPER_QUICK", "").strip().lower() in _TRUE


def qval(full, quick):
    """Return ``quick`` in quick-build mode, else ``full``.

    Use at module load or argument-default time to scale loop counts, grids,
    and replicate numbers, e.g. ``N_SEQ_REPS = qval(30, 5)``.
    """
    return quick if QUICK else full


_WATERMARK_TEXT = "QUICK BUILD: LOW-ITERATION, NOT FOR PUBLICATION"


def stamp_quick(fig: Figure) -> None:
    """Draw the quick-build watermark across a figure (no-op unless QUICK)."""
    if not QUICK or getattr(fig, "_quick_stamped", False):
        return
    fig.text(
        0.5, 0.5, _WATERMARK_TEXT,
        transform=fig.transFigure, ha="center", va="center",
        rotation=28, fontsize=30, fontweight="bold", color="red",
        alpha=0.18, zorder=10000,
    )
    fig._quick_stamped = True  # type: ignore[attr-defined]


# Patch savefig once, at import, so any figure written during a quick build is
# stamped without each call site having to remember. plt.savefig() routes
# through Figure.savefig(), so this covers both styles.
if QUICK and not getattr(Figure, "_quick_savefig_patched", False):
    _orig_savefig = Figure.savefig

    def _savefig(self, *args, **kwargs):  # noqa: ANN001, ANN002, ANN003
        stamp_quick(self)
        return _orig_savefig(self, *args, **kwargs)

    Figure.savefig = _savefig  # type: ignore[method-assign]
    Figure._quick_savefig_patched = True  # type: ignore[attr-defined]
