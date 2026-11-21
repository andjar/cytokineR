#' Choose a pipeline on training samples, judge it on unseen ones
#'
#' Trying many pre-processing combinations and keeping the one that gives the
#' best result is a way to overfit: the winner is partly the pipeline that
#' happened to fit the noise. `select_pipeline()` keeps the choice honest. It
#' runs every on/off combination of the optional steps (as
#' [waterfall_lattice()]) with the test samples hidden from every fitting
#' step, scores each pipeline on the **training** samples, picks the best, and
#' only then looks at how that choice does on the **test** samples, which no
#' pre-processing step, model fit or choice has seen.
#'
#' The gap between the training and test scores of the chosen pipeline is the
#' optimism of the choice. The test scores of all pipelines are reported too,
#' to show how much the choice mattered, but choosing by them would overfit
#' the test set: they are for reading, not for choosing.
#'
#' @section Criteria:
#' `criteria` is one statistic or a named list of them, typically a mix of
#' * **harmonisation**: residual batch, cohort or plate effects among the study
#'   samples, e.g. `sample_metric(metric_plate_effect, by = "batch", within =
#'   "time")` (lower is better);
#' * **biology kept**: e.g. `sample_metric(metric_bio, ~ time, adjust = ~
#'   subject)` or a [model_statistic()] coefficient or prediction score
#'   (higher is better);
#' * **agreement of repeated material**: e.g.
#'   `default_metrics()$reliability_batch` with `holdout = "ctrl_high"`.
#'
#' Each criterion is evaluated per pipeline on the training samples
#' (`on = "train"`) and on the test samples (`on = "test"`); see
#' [sample_metric()] and [model_statistic()]. A criterion without an `on`
#' argument, such as the reliability of held-out control material, gives the
#' same value for both -- that material is then used for choosing, so its
#' score is not an independent check.
#'
#' Pipelines are ranked per criterion (respecting whether higher or lower is
#' better), and the ranks are averaged with `weights`. The pipeline with the
#' best mean training rank is chosen; ties go to the pipeline with fewer
#' steps.
#'
#' Ranks treat a difference of 0.001 like a difference of 0.1. Many
#' corrections often end up practically equal, and choosing between them by
#' the third decimal chases noise. With `resolution`, differences smaller
#' than the resolution count as ties (values are rounded to it before
#' ranking), so the simpler pipeline wins unless a more complex one is
#' better by a margin you care about -- a guard against overfitting the
#' choice. Give one number for all criteria or a named vector, on each
#' metric's own scale (e.g. `c(batch = 0.01, biology = 0.01)`).
#'
#' @param x A [cyto] object (raw, i.e. before processing).
#' @param steps,fixed,route,holdout,max_steps,verbose As in
#'   [waterfall_lattice()].
#' @param criteria A statistic or a named list of statistics; see Criteria.
#' @param test The test set: a selector or a vector of `sample_id`s, e.g. from
#'   [split_samples()].
#' @param weights Optional named numeric vector of criterion weights (default
#'   1 each).
#' @param maximize Optional named logical vector: is higher better? Default
#'   from the metric name (`plate_r2*`, `pred_rmse` and `truth_rmse` are
#'   minimised, everything else maximised).
#' @param analytes Analytes whose mean score counts (default all).
#' @param resolution Optional smallest difference that counts, per criterion
#'   (named) or for all; see Criteria.
#' @return An object of class `cyto_selection`: a list with
#'   * `scores`: one row per pipeline with each criterion's mean over analytes
#'     on training and test samples (`<criterion>_train`,
#'     `<criterion>_test`), the mean ranks (`rank_train`, `rank_test`) and
#'     `selected`; with a single criterion also `train` and `test`;
#'   * `criteria`: name, metric, direction and weight of each criterion;
#'   * `selected` (name) and `pipeline` (the chosen [cyto_pipeline()]);
#'   * `lattices`: a `cyto_lattice` per criterion and split (`train`/`test`
#'     for a single criterion);
#'   * the test set and settings, used by [finalise_pipeline()] and
#'     [processing_report()].
#'
#'   Has `print()` and `plot()` methods; `plot()` puts the training score
#'   against the test score of every pipeline, per criterion.
#' @seealso [split_samples()], [sample_metric()], [model_statistic()],
#'   [finalise_pipeline()], `vignette("harmonising")`
#' @examples
#' test <- split_samples(cyto_example, prop = 0.4, by = "subject",
#'                       strata = "batch")
#' steps <- list(log2   = pipe_transform("log2"),
#'               median = pipe_normalize("median", within = "time"),
#'               bridge = pipe_correct_batch("bridge"))
#' sel <- select_pipeline(
#'   cyto_example, steps, fixed = "log2", test = test,
#'   criteria = list(
#'     batch   = sample_metric(metric_plate_effect, by = "batch",
#'                             within = "time"),
#'     biology = sample_metric(metric_bio, ~ time, adjust = ~ subject)))
#' sel
#' plot(sel)
#' sel$pipeline
#' @export
select_pipeline <- function(x, steps, criteria, test, fixed = character(),
                            route = NULL, holdout = NULL, weights = NULL,
                            maximize = NULL, analytes = NULL,
                            resolution = NULL, max_steps = 8,
                            verbose = FALSE) {
  .check_cyto(x)
  if (missing(test) || is.null(test)) {
    stop("select_pipeline() needs a test set; see split_samples().",
         call. = FALSE)
  }
  single <- is.function(criteria)
  if (single) criteria <- list(criteria)
  if (!is.list(criteria) || !length(criteria) ||
      !all(vapply(criteria, is.function, TRUE))) {
    stop("'criteria' must be a function or a list of functions, e.g. from ",
         "model_statistic() or sample_metric().", call. = FALSE)
  }
  nms <- names(criteria) %||% rep("", length(criteria))
  nms[!nzchar(nms)] <- paste0("criterion", which(!nzchar(nms)))
  names(criteria) <- make.unique(nms)
  steps <- .lattice_steps(steps)
  bad <- setdiff(fixed, names(steps))
  if (length(bad)) {
    stop("Unknown fixed step(s): ", paste(bad, collapse = ", "), call. = FALSE)
  }
  optional <- setdiff(names(steps), fixed)
  K <- length(optional)
  if (K < 1) stop("At least one step must be optional.", call. = FALSE)
  if (K > max_steps) {
    stop(K, " optional steps would need ", 2^K, " pipelines. Fix some steps ",
         "or raise 'max_steps'.", call. = FALSE)
  }
  route <- route %||% optional
  if (!setequal(route, optional) || length(route) != K) {
    stop("'route' must list each optional step once: ",
         paste(optional, collapse = ", "), call. = FALSE)
  }
  test_ids <- .test_ids(x$data, test)
  xf <- .hide_rows(x, holdout, test_ids)

  # two statistics per criterion: on the training and on the test samples
  split_stat <- function(f, on) {
    has_on <- "on" %in% names(formals(args(f)))
    function(x, holdout = NULL, test = NULL) {
      if (has_on) {
        a <- list(x, on = on)
        fa <- names(formals(args(f)))
        if ("holdout" %in% fa) a["holdout"] <- list(holdout)
        if ("test" %in% fa) a["test"] <- list(test)
        do.call(f, a)
      } else {
        .call_statistic(f, x, holdout, test)
      }
    }
  }
  stats <- list()
  for (cn in names(criteria)) {
    stats[[paste0(cn, "::train")]] <- split_stat(criteria[[cn]], "train")
    stats[[paste0(cn, "::test")]] <- split_stat(criteria[[cn]], "test")
  }
  run <- .lattice_run(xf, steps, optional, stats, holdout, test_ids, verbose)
  lats <- list()
  for (sn in names(stats)) {
    l <- tryCatch(.lattice_build(run, sn, steps, optional, fixed, route, NULL,
                                 0.05),
                  error = function(e) e)
    if (inherits(l, "error")) {
      stop("Criterion '", sub("::.*", "", sn), "' could not be computed on ",
           "the ", sub(".*::", "", sn), "ing samples: ", conditionMessage(l),
           call. = FALSE)
    }
    l$test <- test_ids
    lats[[sn]] <- l
  }

  # per criterion: mean over analytes, direction, ranks
  crit <- data.table::data.table(
    criterion = names(criteria),
    metric = vapply(names(criteria),
                    function(cn) lats[[paste0(cn, "::train")]]$label, ""))
  dir <- vapply(seq_len(nrow(crit)), function(i) {
    cn <- crit$criterion[i]
    if (!is.null(maximize) && cn %in% names(maximize)) {
      return(if (isTRUE(maximize[[cn]])) 1 else -1)
    }
    if (!is.null(maximize) && length(maximize) == 1L && is.null(names(maximize))) {
      return(if (isTRUE(maximize)) 1 else -1)
    }
    .direction_of(crit$metric[i])
  }, 1)
  crit[, direction := dir]
  crit[, weight := if (is.null(weights)) 1 else
    ifelse(criterion %in% names(weights), weights[criterion], 1)]

  agg <- function(l) {
    nd <- l$nodes
    if (!is.null(analytes)) nd <- nd[analyte %in% analytes]
    if (!nrow(nd)) stop("None of 'analytes' was found.", call. = FALSE)
    nd[, .(value = mean(value), n_on = n_on[1]), by = .(mask, pipeline)]
  }
  first <- agg(lats[[1]])
  sc <- first[, .(mask, pipeline, n_on)]
  rk_tr <- rk_te <- matrix(NA_real_, nrow(sc), nrow(crit))
  for (i in seq_len(nrow(crit))) {
    cn <- crit$criterion[i]
    tr <- agg(lats[[paste0(cn, "::train")]])
    te <- agg(lats[[paste0(cn, "::test")]])
    v_tr <- tr$value[match(sc$mask, tr$mask)]
    v_te <- te$value[match(sc$mask, te$mask)]
    data.table::set(sc, j = paste0(cn, "_train"), value = v_tr)
    data.table::set(sc, j = paste0(cn, "_test"), value = v_te)
    res <- if (is.null(resolution)) NA_real_ else
      if (!is.null(names(resolution))) {
        if (cn %in% names(resolution)) resolution[[cn]] else NA_real_
      } else resolution[[1]]
    snap <- function(v) if (is.na(res) || res <= 0) v else round(v / res) * res
    rk_tr[, i] <- rank(-crit$direction[i] * snap(v_tr),
                       ties.method = "average", na.last = "keep")
    rk_te[, i] <- rank(-crit$direction[i] * snap(v_te),
                       ties.method = "average", na.last = "keep")
  }
  w <- crit$weight
  wmean <- function(r) {
    out <- apply(r, 1, function(v) if (anyNA(v)) NA_real_ else
      sum(w * v) / sum(w))
    out
  }
  sc[, mean_rank_train := wmean(rk_tr)]
  sc[, mean_rank_test := wmean(rk_te)]
  if (all(is.na(sc$mean_rank_train))) {
    stop("No pipeline could be scored on the training samples; see ",
         "$errors of waterfall_lattice() with the same arguments.",
         call. = FALSE)
  }
  sc[, rank_train := rank(mean_rank_train, ties.method = "min",
                          na.last = "keep")]
  sc[, rank_test := rank(mean_rank_test, ties.method = "min",
                         na.last = "keep")]
  data.table::setorderv(sc, c("mean_rank_train", "n_on", "mask"),
                        na.last = TRUE)
  best <- sc[1]
  sc[, selected := mask == best$mask]
  if (single) {
    cn <- crit$criterion[1]
    sc[, `:=`(train = get(paste0(cn, "_train")),
              test = get(paste0(cn, "_test")))]
  }
  on_steps <- optional[bitwAnd(best$mask, as.integer(2^(seq_len(K) - 1))) > 0]
  chosen <- steps[names(steps) %in% c(fixed, on_steps)]
  pipe <- do.call(cyto_pipeline, c(unname(chosen),
                                   list(name = best$pipeline)))
  if (single) {
    names(lats) <- sub(".*::", "", names(lats))
  }

  structure(list(
    scores = sc, criteria = crit, selected = best$pipeline, pipeline = pipe,
    lattices = lats,
    train = if (single) lats$train else NULL,
    test = if (single) lats$test else NULL,
    label = if (single) crit$metric[1] else "mean rank",
    maximize = if (single) crit$direction[1] > 0 else NULL,
    single = single, test_ids = test_ids, holdout = holdout,
    steps = steps, fixed = fixed, optional = optional,
    n_train = data.table::uniqueN(
      xf$data[sample_type == "unknown" & !sample_id %in% test_ids, sample_id]),
    errors = data.table::rbindlist(run$errors), analytes = analytes,
    resolution = resolution, weights = weights,
    time = Sys.time()
  ), class = "cyto_selection")
}


#' @export
print.cyto_selection <- function(x, digits = 3, ...) {
  sc <- x$scores
  K <- length(x$optional)
  s <- sc[selected == TRUE]
  r <- function(v) format(round(v, digits), nsmall = digits)
  cat("<cyto_selection> ", nrow(sc), " pipelines (", K, " optional step",
      if (K > 1) "s", ")", if (x$single) paste0(" | statistic: ", x$label),
      "\n", sep = "")
  cat("Training: ", x$n_train, " samples | test: ", length(x$test_ids),
      " unseen samples (not used for fitting or choosing)\n", sep = "")
  if (!is.null(x$analytes)) {
    cat("Scores: mean over ", paste(x$analytes, collapse = ", "), "\n",
        sep = "")
  }
  cat("\nChosen on the training samples: ", s$pipeline, "\n", sep = "")
  if (x$single) {
    cat("  training score: ", r(s$train), "\n", sep = "")
    cat("  test score:     ", r(s$test), "  (rank ", s$rank_test, " of ",
        sum(is.finite(sc$test)), " on test)\n", sep = "")
  } else {
    tab <- x$criteria[, .(criterion, metric,
                          better = ifelse(direction > 0, "higher", "lower"),
                          weight)]
    tab[, train := vapply(criterion, function(cn)
      s[[paste0(cn, "_train")]], 1)]
    tab[, test := vapply(criterion, function(cn) s[[paste0(cn, "_test")]], 1)]
    tab[, `:=`(train = round(train, digits), test = round(test, digits))]
    print(tab, row.names = FALSE)
    cat("  mean rank: ", r(s$mean_rank_train), " in training, ",
        r(s$mean_rank_test), " on test (rank ", s$rank_test, " of ",
        sum(is.finite(sc$mean_rank_test)), ")\n", sep = "")
  }
  bt <- sc[which.min(mean_rank_test)]
  cat("Best on test (for reading, not choosing): ", bt$pipeline,
      if (x$single) paste0(", ", r(bt$test)), "\n", sep = "")
  out <- data.table::copy(sc)
  out <- if (x$single) {
    out[, .(pipeline, train = round(train, digits), test = round(test, digits),
            rank_train, rank_test)]
  } else {
    out[, .(pipeline, mean_rank_train = round(mean_rank_train, 2),
            mean_rank_test = round(mean_rank_test, 2), rank_train,
            rank_test)]
  }
  cat("\n")
  print(utils::head(out, 8), row.names = FALSE)
  if (nrow(out) > 8) cat("... ", nrow(out) - 8, " more in $scores\n", sep = "")
  if (nrow(x$errors)) {
    cat("\n", nrow(x$errors), " error(s); see $errors.\n", sep = "")
  }
  invisible(x)
}

#' @rdname select_pipeline
#' @param ... Unused.
#' @export
plot.cyto_selection <- function(x, ...) {
  sc <- x$scores
  bt <- sc[which.min(mean_rank_test), mask]
  long <- data.table::rbindlist(lapply(seq_len(nrow(x$criteria)), function(i) {
    cn <- x$criteria$criterion[i]
    data.table::data.table(
      criterion = paste0(cn, " (", x$criteria$metric[i], ", ",
                         ifelse(x$criteria$direction[i] > 0, "higher",
                                "lower"), " is better)"),
      mask = sc$mask, pipeline = sc$pipeline, n_on = sc$n_on,
      selected = sc$selected, train = sc[[paste0(cn, "_train")]],
      test = sc[[paste0(cn, "_test")]])
  }))
  long <- long[is.finite(train) & is.finite(test)]
  long[, what := ifelse(selected, "chosen on training",
                        ifelse(mask == bt, "best on test", "other"))]
  long[, what := factor(what, levels = c("chosen on training", "best on test",
                                         "other"))]
  lab <- long[what != "other" | n_on == 0L]
  lab[, vj := ifelse(what == "best on test", 1.8, -1)]
  ggplot2::ggplot(long, ggplot2::aes(x = .data$train, y = .data$test)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey60",
                         linetype = 2) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$what, size = .data$what)) +
    ggplot2::geom_text(data = lab, ggplot2::aes(label = .data$pipeline,
                                                vjust = .data$vj),
                       size = 2.6) +
    ggplot2::scale_colour_manual(values = c(`chosen on training` = "#D55E00",
                                            `best on test` = "#0072B2",
                                            other = "grey45"),
                                 drop = FALSE, name = NULL) +
    ggplot2::scale_size_manual(values = c(`chosen on training` = 3.2,
                                          `best on test` = 2.6, other = 1.8),
                               drop = FALSE, name = NULL) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.15)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.15)) +
    ggplot2::facet_wrap(~criterion, scales = "free") +
    ggplot2::labs(x = "training samples", y = "unseen test samples",
                  subtitle = paste("One point per pipeline; dashed line: test",
                                   "= training")) +
    .theme()
}

#' Evaluate a metric on the training or the test samples only
#'
#' Wraps a `metric_*()` function (or any function of a [cyto] object) so that
#' it can be computed on the training samples, the test samples or all study
#' samples -- the form [select_pipeline()] needs for its criteria. Rows that
#' are not study samples (standards, blanks, controls) are always kept.
#'
#' @param metric A function whose first argument is a [cyto] object, e.g.
#'   [metric_plate_effect()], [metric_bio()], [metric_detection()].
#' @param ... Further arguments to `metric`.
#' @param types Sample types forming the study population.
#' @param pool_test Compute the test value on all study samples instead of
#'   the test samples alone? Use it for effects *between* groups (batch,
#'   cohort) when the test set is a whole group ([split_samples()] with
#'   `values`): a held-out cohort on its own has no cohort effect to measure,
#'   but pooled with the training samples it shows whether it lines up with
#'   them. The test samples were still not used for fitting.
#' @return A function `function(x, holdout = NULL, test = NULL, on = "test")`.
#' @examples
#' batch_r2 <- sample_metric(metric_plate_effect, by = "batch",
#'                           within = "time")
#' x <- cyto_transform(cyto_example, "log2")
#' test <- split_samples(cyto_example, prop = 0.4, by = "subject")
#' batch_r2(x, test = test)               # among the test samples
#' batch_r2(x, test = test, on = "train") # among the training samples
#' @export
sample_metric <- function(metric, ..., types = "unknown", pool_test = FALSE) {
  if (!is.function(metric)) stop("'metric' must be a function.", call. = FALSE)
  args <- list(...)
  function(x, holdout = NULL, test = NULL, on = c("test", "train", "all")) {
    on <- match.arg(on)
    .check_cyto(x)
    if (pool_test && on == "test") on <- "all"
    if (!is.null(test) && on != "all") {
      ids <- .test_ids(x$data, test)
      study <- x$data$sample_type %in% types
      it <- x$data$sample_id %in% ids
      keep <- !study | (if (on == "test") it else !it)
      x <- .cyto_copy(x)
      x$data <- x$data[keep]
    }
    do.call(metric, c(list(x), args))
  }
}
