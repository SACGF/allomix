"""Per-marker calibration tables estimated from cohort data.

Amplification bias, per-site sequencing error rates, and the co-pooled
contamination correction. Each module estimates a table from reference VCFs,
loads/saves it as TSV, and is applied by the estimator through
``allomix.estimate.likelihood.PanelCalibration``.
"""
