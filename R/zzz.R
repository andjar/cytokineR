.register_builtin <- function() {
  r <- function(stage, name, fun, desc) {
    register_cyto_method(stage, name, fun, desc, overwrite = TRUE)
  }
  # transform
  r("transform", "log2", .tf_log2, "log2(FI)")
  r("transform", "bgsub", .tf_bgsub, "log2(FI - plate blank)")
  r("transform", "mob", .tf_mob,
    "Multiple of blank: log2(FI / plate blank), EB-shrunk blank")
  r("transform", "mos", .tf_mos,
    "Multiple of standard: log2(FI / FI of a reference standard or concentration)")
  r("transform", "moc", .tf_moc,
    "Multiple of control: log2(FI / control material level on the plate)")
  r("transform", "vendor", .tf_vendor,
    "log2 of the concentration reported by the instrument software")
  # normalize
  r("normalize", "median", .nm_median,
    "Align plate centres of comparable samples (optionally within strata)")
  r("normalize", "quantile", .nm_quantile,
    "Quantile-map each plate's sample distribution to the pooled one")
  r("normalize", "anchor", .nm_anchor,
    "Offsets (or affine maps) from anchor material on every plate")
  r("normalize", "bridge", .nm_bridge,
    "Batch offsets from samples measured in several batches")
  r("normalize", "regress", .nm_regress,
    "Remove group effects from a linear model protecting covariates")
  r("normalize", "combat", .nm_combat, "Empirical-Bayes ComBat (sva)")
  # calibrate
  r("calibrate", "5pl", .cal_5pl, "Five-parameter logistic per curve")
  r("calibrate", "4pl", .cal_4pl, "Four-parameter logistic per curve")
  r("calibrate", "spline", .cal_spline,
    "Monotone interpolation through standards, no extrapolation")
  r("calibrate", "shared", .cal_shared,
    "Shared-shape 5PL across plates with shrunken plate shifts")
  # impute
  r("impute", "fraction", .im_fraction, "Limit x fraction (e.g. LLOQ/2)")
  r("impute", "tobit", .im_tobit,
    "Censored-normal conditional expectation (or draws)")
  r("impute", "rank_tobit", .im_rank_tobit,
    "Censored-normal quantiles assigned by within-plate signal rank")
  invisible(NULL)
}

.onLoad <- function(libname, pkgname) {
  .register_builtin()
}
