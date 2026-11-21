# cytokineR 0.1.0

First release of the rewritten package.

* `cyto` data container (one row per well x analyte) with a processing log
  and stored estimates for every step.
* Import of Bio-Plex Manager exports (multi-analyte and per-analyte layouts,
  Excel/CSV), Luminex xPONENT CSV, and generic long/wide tables.
* Quality control: bead counts, replicate CVs, PCA by plate/batch, plate
  summaries, Levey-Jennings charts with Westgard rules, variance
  decomposition.
* Signal representations: log2 FI, background-subtracted, multiple of blank
  (with empirical-Bayes blank shrinkage), multiple of standard (level or
  reference concentration), multiple of control, vendor concentrations.
* Plate/batch correction: sample centring (within strata), quantile mapping,
  anchor offsets or affine maps (controls, standards), bridging samples
  (offset or Deming), covariate-protected regression, ComBat.
* Calibration: per-curve 5PL/4PL (optionally robust), monotone spline
  interpolation, robust shared-shape 5PL across plates with shrunken
  gain/potency/background shifts; one curve per batch or per study.
* Imputation: fraction of limit, censored-normal (Tobit) with covariates or
  draws, rank-preserving Tobit using the fluorescence signal.
* Pipelines (`pipe_*()`, `cyto_pipeline()`, `expand_pipelines()`), method
  registry for extensions, hold-out benchmarking with scale-free metrics,
  and a simulator with known truth (`simulate_cyto()`,
  `benchmark_simulations()`).
