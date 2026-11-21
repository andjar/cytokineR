#' cytokineR: transparent, modular pre-processing of multiplex cytokine data
#'
#' cytokineR turns raw bead-based immunoassay exports (Bio-Rad Bio-Plex
#' Manager, Luminex xPONENT, or any long/wide table) into analysis-ready values
#' through a sequence of small, swappable steps:
#'
#' 1. **Import** -- [read_bioplex()], [read_xponent()], [cyto_from_long()],
#'    [cyto_from_wide()] create a [cyto] object (one row per well x analyte).
#' 2. **Quality control** -- [flag_beads()], [qc_replicates()], [qc_pca()],
#'    [qc_plates()], [qc_controls()], [qc_variance()].
#' 3. **Signal representation** -- [cyto_transform()]: log2 FI, background
#'    subtraction, multiple-of-blank, multiple-of-standard, control-anchored
#'    ratios, or the vendor's concentrations.
#' 4. **Plate and batch correction** -- [cyto_normalize()]: sample-centring,
#'    quantile, anchor-based (controls, standards, bridging samples), affine
#'    harmonisation, covariate-protected regression, ComBat.
#' 5. **Calibration** -- [cyto_calibrate()]: per-plate 5PL/4PL, monotone
#'    interpolation, or a shared-shape 5PL that borrows strength across plates.
#' 6. **Censored values** -- [cyto_impute()]: fraction-of-limit, censored-normal
#'    (Tobit) conditional means, or a rank-preserving Tobit that keeps the
#'    ordering information still present in the fluorescence signal.
#'
#' Each step is a registered method (see [cyto_methods()] and
#' [register_cyto_method()]) that records its arguments and estimates in the
#' object ([cyto_steps()], [cyto_fit()]). Pipelines ([cyto_pipeline()]) can be
#' compared with [benchmark_pipelines()] using scale-free metrics on held-out
#' replicate material, and developed against data with known truth from
#' [simulate_cyto()]. [waterfall_lattice()] shows how much each step
#' contributes to a result, in every combination of the other steps.
#'
#' @keywords internal
"_PACKAGE"

#' @import data.table
#' @importFrom stats median mad sd var quantile lm residuals coef optim nlminb
#'   approx splinefun isoreg prcomp complete.cases pnorm dnorm qnorm rnorm
#'   runif rpois rbinom setNames as.formula model.matrix fitted aggregate
#'   na.omit reformulate
#' @importFrom utils head tail read.csv
NULL

# data.table is used internally with non-standard evaluation
.datatable.aware <- TRUE

utils::globalVariables(c(
  ".", ".N", ".SD", ".I", ".GRP", ".BY", ":=", ".row",
  "plate", "batch", "well", "sample_id", "sample_type", "analyte", "fi",
  "level", "expected_conc", "dilution", "bead_count", "vendor_conc",
  "vendor_flag", "value", "censored", "limit_low", "limit_high", "qc",
  "use_for_fit", "signal", "truth", "x", "y", "offset", "slope", "n",
  "material", "consensus", "dev", "group", "ref", "w", "conc", "fit_ok",
  "recovery", "in_range", "lloq", "uloq", "metric", "score", "pipeline",
  "estimate", "PC1", "PC2", "z", "mean_v", "sd_v", "rule", "flagged",
  "median_value", "q25", "q75", "x_true", "x_hat", "curve_id", "shift_x",
  "shift_y", "bg", "ref_fi", "term", "r2", "within_var", "total_var",
  "n_obs", "label", "status", "stat", "lo", "hi", "type", "rank_in",
  "stratum", "pos", "dist", "plate_label", "is_outlier", "flag_reason",
  # simulate
  "gain", "shift", "background", "lot_error", "u", "e", "xw", "fi_true",
  "s1", "mu", "has_effect", "bottom", "top", "b", "g", "mid", "tf", "sid_tmp",
  "nb", "subject", "time",
  # transform / normalize / calibrate / impute
  "lb", "vw", "blank_raw", "blank_used", "gmean", "sw2", "tau2", "weight",
  "lfi", "lconc", "ulfi", "ref_level", "ref_conc", "N", "intercept",
  "n_anchor", "ng", "center", "target", "n_ref", "n_strata", "yref", "a",
  "d", "c", "log_b", "log_g", "borrowed", "shift_bottom", "region", "field",
  "v", ".rows", ".code", ".well", ".fi", ".v", ".analyte", "vendor_conc",
  "Location", "Sample", "description", "type_code", ".wid", "newv",
  # qc / metrics / benchmark / plots
  "cv", "m", "run", "spread", "r4", "var_tech", "ss", "df", "direction",
  "mean_rank", "correction", "k", "component", ".col", "nw", "scenario",
  "fill",
  # waterfall lattice
  "mask", "n_on", "p", "significant", "route", "step", "from", "to", "from_n",
  "from_value", "to_value", "delta", "weight", "waterfall", "shapley",
  "flips", "none", "total", "sig"
))
