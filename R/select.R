#' Choose a pipeline on training samples, judge it on unseen ones
#'
#' Trying many pre-processing combinations and keeping the one that gives the
#' best result is a way to overfit: the winner is partly the pipeline that
#' happened to fit the noise. `select_pipeline()` keeps the choice honest. It
#' runs every on/off combination of the optional steps (as
#' [waterfall_lattice()]) with the test samples hidden from every fitting
#' step, scores each pipeline on the **training** samples (cross-validated
#' for prediction models), picks the best, and only then looks at how that
#' choice does on the **test** samples, which no pre-processing step, model
#' fit or choice has seen.
#'
#' The gap between the training and test scores of the chosen pipeline is the
#' optimism of the choice. The test scores of all pipelines are reported too,
#' to show how much the choice mattered, but choosing by them would overfit
#' the test set: they are for reading, not for choosing.
#'
#' @param x A [cyto] object (raw, i.e. before processing).
#' @param steps,fixed,route,holdout,max_steps,verbose As in
#'   [waterfall_lattice()].
#' @param statistic A [model_statistic()] (or any function with arguments
#'   `test` and `on` that follows the same rules).
#' @param test The test set: a selector or a vector of `sample_id`s, e.g. from
#'   [split_samples()].
#' @param maximize Choose the pipeline with the largest training score?
#'   Default: `TRUE`, except for metrics where lower is better (`pred_rmse`).
#' @param analytes Analytes whose mean score decides (default all).
#' @return An object of class `cyto_selection`: a list with `scores` (one row
#'   per pipeline: mean `train` and `test` score over the analytes, their
#'   ranks, and `selected`), `selected` (the chosen pipeline's name),
#'   `pipeline` (the chosen [cyto_pipeline()], ready for [run_pipeline()]),
#'   and `train` and `test` (`cyto_lattice` objects for each). Has `print()`
#'   and `plot()` methods; `plot()` puts the training score against the test
#'   score of every pipeline.
#' @seealso [model_statistic()], [split_samples()], [waterfall_lattice()]
#' @examples
#' test <- split_samples(cyto_example, prop = 0.4, by = "subject")
#' steps <- list(log2   = pipe_transform("log2"),
#'               median = pipe_normalize("median", within = "time"),
#'               bridge = pipe_correct_batch("bridge"))
#' auc <- model_statistic(factor(time) ~ ., model = glm, family = binomial,
#'                        data = "wide", cv_by = "subject")
#' sel <- select_pipeline(cyto_example, steps, auc, test = test,
#'                        fixed = "log2")
#' sel
#' plot(sel)
#' sel$pipeline
#' @export
select_pipeline <- function(x, steps, statistic, test, fixed = character(),
                            route = NULL, holdout = NULL, maximize = NULL,
                            analytes = NULL, max_steps = 8, verbose = FALSE) {
  .check_cyto(x)
  if (missing(test) || is.null(test)) {
    stop("select_pipeline() needs a test set; see split_samples().",
         call. = FALSE)
  }
  if (!is.function(statistic) ||
      !all(c("test", "on") %in% names(formals(args(statistic))))) {
    stop("'statistic' must be a function with arguments 'test' and 'on', ",
         "e.g. from model_statistic().", call. = FALSE)
  }
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
  stats <- list(
    train = function(x, holdout = NULL, test = NULL)
      statistic(x, holdout = holdout, test = test, on = "train"),
    test = function(x, holdout = NULL, test = NULL)
      statistic(x, holdout = holdout, test = test, on = "test"))
  run <- .lattice_run(xf, steps, optional, stats, holdout, test_ids, verbose)
  lat <- lapply(c(train = "train", test = "test"), function(w) {
    l <- .lattice_build(run, w, steps, optional, fixed, route, NULL, 0.05)
    l$test <- test_ids
    l
  })
  label <- lat$train$label
  maximize <- maximize %||% (.direction_of(label) > 0)

  agg <- function(l) {
    nd <- l$nodes
    if (!is.null(analytes)) nd <- nd[analyte %in% analytes]
    if (!nrow(nd)) stop("None of 'analytes' was found.", call. = FALSE)
    nd[, .(value = mean(value), n_on = n_on[1]), by = .(mask, pipeline)]
  }
  tr <- agg(lat$train)
  te <- agg(lat$test)
  sc <- merge(tr[, .(mask, pipeline, n_on, train = value)],
              te[, .(mask, test = value)], by = "mask")
  sgn <- if (maximize) -1 else 1
  sc[, rank_train := rank(sgn * train, ties.method = "min", na.last = "keep")]
  sc[, rank_test := rank(sgn * test, ties.method = "min", na.last = "keep")]
  if (all(is.na(sc$train))) {
    stop("No pipeline could be scored on the training samples; see ",
         "$errors of waterfall_lattice() with the same arguments.",
         call. = FALSE)
  }
  best <- sc[which.min(rank_train)]
  sc[, selected := mask == best$mask]
  data.table::setorderv(sc, "rank_train", na.last = TRUE)
  on_steps <- optional[bitwAnd(best$mask, as.integer(2^(seq_len(K) - 1))) > 0]
  chosen <- steps[names(steps) %in% c(fixed, on_steps)]
  pipe <- do.call(cyto_pipeline, c(unname(chosen),
                                   list(name = best$pipeline)))

  structure(list(
    scores = sc, selected = best$pipeline, pipeline = pipe,
    train = lat$train, test = lat$test, label = label, maximize = maximize,
    test_ids = test_ids,
    n_train = data.table::uniqueN(
      xf$data[sample_type == "unknown" & !sample_id %in% test_ids, sample_id]),
    errors = lat$train$errors, analytes = analytes
  ), class = "cyto_selection")
}

#' @export
print.cyto_selection <- function(x, digits = 3, ...) {
  sc <- x$scores
  K <- length(x$train$optional)
  s <- sc[selected == TRUE]
  r <- function(v) format(round(v, digits), nsmall = digits)
  cat("<cyto_selection> ", nrow(sc), " pipelines (", K, " optional step",
      if (K > 1) "s", ") | statistic: ", x$label, "\n", sep = "")
  cat("Training: ", x$n_train, " samples | test: ", length(x$test_ids),
      " unseen samples (not used for fitting or choosing)\n", sep = "")
  if (!is.null(x$analytes)) {
    cat("Score: mean over ", paste(x$analytes, collapse = ", "), "\n",
        sep = "")
  } else if (data.table::uniqueN(x$train$nodes$analyte) > 1) {
    cat("Score: mean over analytes\n")
  }
  cat("\nChosen on the training samples: ", s$pipeline, "\n", sep = "")
  cat("  training score: ", r(s$train), "\n", sep = "")
  cat("  test score:     ", r(s$test), "  (rank ", s$rank_test, " of ",
      sum(is.finite(sc$test)), " on test)\n", sep = "")
  bt <- sc[which.min(rank_test)]
  cat("Best on test (for reading, not choosing): ", bt$pipeline, ", ",
      r(bt$test), "\n", sep = "")
  out <- data.table::copy(sc)[, .(pipeline, train, test, rank_train,
                                  rank_test)]
  out[, `:=`(train = round(train, digits), test = round(test, digits))]
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
  sc <- data.table::copy(x$scores)[is.finite(train) & is.finite(test)]
  bt <- sc[which.min(rank_test), mask]
  sc[, what := ifelse(selected, "chosen on training",
                      ifelse(mask == bt, "best on test", "other"))]
  sc[, what := factor(what, levels = c("chosen on training", "best on test",
                                       "other"))]
  lab <- sc[what != "other" | n_on == 0L | n_on == max(n_on)]
  rng <- range(c(sc$train, sc$test))
  ggplot2::ggplot(sc, ggplot2::aes(x = .data$train, y = .data$test)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, colour = "grey60",
                         linetype = 2) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$what, size = .data$what)) +
    ggplot2::geom_text(data = lab, ggplot2::aes(label = .data$pipeline),
                       vjust = -1, size = 2.8) +
    ggplot2::scale_colour_manual(values = c(`chosen on training` = "#D55E00",
                                            `best on test` = "#0072B2",
                                            other = "grey45"),
                                 drop = FALSE, name = NULL) +
    ggplot2::scale_size_manual(values = c(`chosen on training` = 3.2,
                                          `best on test` = 2.6, other = 1.8),
                               drop = FALSE, name = NULL) +
    ggplot2::coord_cartesian(xlim = rng, ylim = rng) +
    ggplot2::labs(x = paste0(x$label, " on training samples"),
                  y = paste0(x$label, " on unseen test samples"),
                  subtitle = paste("One point per pipeline; dashed line: test",
                                   "= training (no optimism)")) +
    .theme()
}
