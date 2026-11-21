# cytokineR (development version)

* `waterfall_lattice()` runs a pipeline with every on/off combination of its
  optional steps and shows how much each step contributes to a result: the
  change along the usual waterfall route, the order-free Shapley share
  (shares add up to the total change), and the range across all combinations
  of the other steps. `plot()` draws the lattice (every pipeline a node, the
  waterfall route in bold) and the per-step contributions; see
  `vignette("waterfall-lattice")`.
* `metric_effect()`: standardised effect of a study variable (Cohen's d, or
  the within-subject change with `adjust = ~ subject`) with its p-value, for
  following the reported result through the pipeline.
* `model_statistic()` turns any model called as `model(formula, data = )`
  (`lm()`, `glm()`, `lme4::lmer()`, ...) into a statistic or metric: a
  coefficient (with SE and p-value) per analyte or from one model on the
  wide table, or out-of-sample prediction scores (AUC, R^2, RMSE,
  correlation, or your own), cross-validated within training samples or
  scored on unseen test samples.
* Unseen test samples: `split_samples()` sets aside samples or whole
  subjects (bridging samples stay in training); the new `test` argument of
  `waterfall_lattice()` and `benchmark_pipelines()` hides them from every
  fitting step and passes them to the statistic.
* `select_pipeline()` chooses the best combination of steps on the training
  samples and only then scores it on the test samples, reporting the
  optimism of the choice; `plot()` shows the training against the test score
  of every pipeline, and `$pipeline` is the chosen pipeline.

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
