"""Quality control, identity checks, host-presence detection, and run metadata.

The QC assessment (``qc``) plus the independent sub-analyses it aggregates:
host-presence detection (``host_presence``), in-data contamination
(``sample_contamination``), relatedness / sample-swap checks (``relatedness``),
and sequencing run-unit metadata (``runmeta``).

This package's ``__init__`` is intentionally empty: ``qc.qc`` imports
``allomix.results``, which imports ``qc`` submodules in turn, so keeping the
package import side-effect-free avoids a partial-initialisation cycle.
"""
