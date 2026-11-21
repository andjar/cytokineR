# Metrics return data.table(metric, analyte, value, n). They are scale-free
# (ratios of variances, R^2, rank correlations, proportions), so pipelines
# that produce different scales (log2 FI, ratios, concentrations) can be
# compared on equal terms.

.metric_direction <- c(reliability = 1, plate_r2 = -1, detected = 1,
                       bio_r2 = 1, truth_rho = 1, truth_rmse = -1)

.direction_of <- function(metric) {
  base <- sub("_(plate|batch|[A-Za-z0-9]+)$", "", metric)
  d <- .metric_direction[metric]
  miss <- is.na(d)
  d[miss] <- .metric_direction[base[miss]]
  d[is.na(d)] <- 1
  unname(d)
}

#' Pipeline quality metrics
#'
#' Scale-free metrics computed on a processed [cyto] object. Used by
#' [benchmark_pipelines()], but also useful on their own.
#'
#' * `metric_reliability()` -- *technical reliability*: one minus the ratio of
#'   the technical variance of repeatedly measured material (pooled variance of
#'   its plate means across `between` groups) to the total variance of the
#'   study samples. 1 = perfect, 0 = technical noise as large as the whole
#'   biological spread. Computed on `materials` (selector; default: every
#'   control/unknown `sample_id` measured in at least two `between` groups).
#'   For an honest score the materials must **not** have been used to fit the
#'   normalisation -- see the `holdout` argument of [benchmark_pipelines()].
#' * `metric_plate_effect()` -- adjusted partial R^2 of plate among the study
#'   samples, after accounting for the strata in `within` (e.g. time point).
#'   Close to 0 is ideal for randomised layouts.
#' * `metric_bio()` -- adjusted partial R^2 of known biology (`formula`, e.g.
#'   `~ time`), after nuisance terms in `adjust` (e.g. `~ subject`). Higher
#'   means more biological signal retained.
#' * `metric_detection()` -- fraction of study-sample values that are
#'   quantified (finite and not censored).
#' * `metric_truth()` -- for simulated data: Spearman correlation between
#'   values and true concentrations (censored values ranked at the extremes,
#'   so censoring is not rewarded) and, on the concentration scale only, the
#'   RMSE of log2 concentrations after removing each analyte's overall bias
#'   (absolute bias reflects the reference lot's assigned standard
#'   concentrations, which no pre-processing can correct).
#'
#' @param x A processed [cyto] object.
#' @param materials Selector of repeated material.
#' @param between Grouping across which repeats are compared (`"plate"` or
#'   `"batch"`).
#' @param exclude_imputed Ignore imputed values?
#' @param types Sample types forming the study population.
#' @return A data.table with columns `metric`, `analyte`, `value`, `n`.
#' @examples
#' x <- cyto_transform(cyto_example, "mob")
#' metric_reliability(x, materials = "ctrl_high")
#' metric_plate_effect(x)
#' metric_bio(x, ~ time, adjust = ~ subject)
#' metric_truth(x)
#' @export
metric_reliability <- function(x, materials = NULL, between = "plate",
                               exclude_imputed = FALSE,
                               types = "unknown") {
  .check_cyto(x)
  d <- x$data
  .check_by(d, between)
  ok <- is.finite(d$value) & d$sample_type %in% c("control", "unknown")
  if (exclude_imputed) ok <- ok & d$censored == "none"
  rep_rows <- if (is.null(materials)) {
    m <- d[ok, .(ng = data.table::uniqueN(get(between))), by = sample_id][
      ng > 1, sample_id]
    ok & d$sample_id %in% m
  } else {
    ok & select_rows(d, materials)
  }
  name <- paste0("reliability_", between)
  if (!any(rep_rows)) {
    return(data.table::data.table(metric = name, analyte = unique(d$analyte),
                                  value = NA_real_, n = 0L))
  }
  ym <- d[rep_rows, .(y = mean(value)), by = c("analyte", "sample_id", between)]
  tech <- ym[, .(ss = sum((y - mean(y))^2), df = .N - 1L),
             by = .(analyte, sample_id)][df > 0,
                                         .(var_tech = sum(ss) / sum(df),
                                           n = sum(df) + .N),
                                         by = analyte]
  tot <- d[ok & d$sample_type %in% types,
           .(var_total = stats::var(value)), by = analyte]
  out <- merge(tot, tech, by = "analyte", all.x = TRUE)
  out[, .(metric = name, analyte, value = 1 - var_tech / var_total,
          n = ifelse(is.na(n), 0L, as.integer(n)))]
}

.partial_r2 <- function(df, test_terms, base_terms) {
  df <- df[stats::complete.cases(df), , drop = FALSE]
  n <- nrow(df)
  if (n < 5) return(c(NA_real_, n))
  mk <- function(tt) {
    if (!length(tt)) return(stats::as.formula("value ~ 1"))
    stats::reformulate(tt, "value")
  }
  red <- tryCatch(stats::lm(mk(base_terms), data = df), error = function(e) NULL)
  full <- tryCatch(stats::lm(mk(c(base_terms, test_terms)), data = df),
                   error = function(e) NULL)
  if (is.null(red) || is.null(full)) return(c(NA_real_, n))
  rss_r <- sum(stats::residuals(red)^2)
  rss_f <- sum(stats::residuals(full)^2)
  k_r <- red$rank
  k_f <- full$rank
  if (n - k_f <= 0 || rss_r <= 0) return(c(NA_real_, n))
  c(1 - (rss_f / (n - k_f)) / (rss_r / (n - k_r)), n)
}

#' @rdname metric_reliability
#' @param by Technical grouping (default `"plate"`).
#' @param within Strata (column names) accounted for before testing `by`.
#' @export
metric_plate_effect <- function(x, by = "plate", within = NULL,
                                types = "unknown") {
  .check_cyto(x)
  d <- x$data[sample_type %in% types & is.finite(value)]
  .check_by(d, c(by, within))
  fac <- function(v) paste0("factor(`", v, "`)")
  d[, {
    r <- .partial_r2(as.data.frame(.SD), fac(by),
                     if (length(within)) fac(within) else character())
    list(metric = paste0("plate_r2", if (by != "plate") paste0("_", by) else ""),
         value = r[1], n = as.integer(r[2]))
  }, by = analyte, .SDcols = c("value", by, within)][
    , .(metric, analyte, value, n)]
}

#' @rdname metric_reliability
#' @param formula Biological terms, e.g. `~ time` or `~ group`.
#' @param adjust Nuisance terms, e.g. `~ subject`.
#' @export
metric_bio <- function(x, formula, adjust = NULL, types = "unknown") {
  .check_cyto(x)
  tt <- attr(stats::terms(formula), "term.labels")
  at <- if (is.null(adjust)) character() else
    attr(stats::terms(adjust), "term.labels")
  vars <- unique(c(all.vars(formula), all.vars(adjust)))
  d <- x$data[sample_type %in% types & is.finite(value)]
  .check_by(d, vars)
  d[, {
    r <- .partial_r2(as.data.frame(.SD), tt, at)
    list(metric = "bio_r2", value = r[1], n = as.integer(r[2]))
  }, by = analyte, .SDcols = c("value", vars)][, .(metric, analyte, value, n)]
}

#' @rdname metric_reliability
#' @export
metric_detection <- function(x, types = "unknown") {
  .check_cyto(x)
  x$data[sample_type %in% types,
         .(metric = "detected",
           value = mean(is.finite(value) & censored == "none"), n = .N),
         by = analyte][, .(metric, analyte, value, n)]
}

#' @rdname metric_reliability
#' @export
metric_truth <- function(x, types = "unknown") {
  .check_cyto(x)
  d <- x$data[sample_type %in% types & is.finite(truth)]
  if (!nrow(d)) {
    return(data.table::data.table(metric = character(), analyte = character(),
                                  value = numeric(), n = integer()))
  }
  is_conc <- identical(x$scale, "log2_conc")
  d[, {
    v <- value
    v[is.na(v) & censored == "left"] <- -Inf
    v[is.na(v) & censored == "right"] <- Inf
    ok <- !is.na(v)
    rv <- rank(v[ok])
    rt <- rank(truth[ok])
    rho <- if (sum(ok) > 3 && stats::sd(rv) > 0 && stats::sd(rt) > 0)
      stats::cor(rv, rt) else NA_real_
    fin <- is.finite(value)
    rmse <- if (is_conc && sum(fin) > 3) {
      err <- value[fin] - truth[fin]
      sqrt(mean((err - stats::median(err))^2))
    } else NA_real_
    list(metric = c("truth_rho", "truth_rmse"), value = c(rho, rmse),
         n = c(sum(ok), sum(fin)))
  }, by = analyte][, .(metric, analyte, value, n)]
}
