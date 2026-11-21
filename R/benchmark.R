#' Default benchmark metrics
#'
#' Builds the list of metric functions used by [benchmark_pipelines()].
#' Each element is a function `function(x, holdout)` returning the long
#' metric table described in [metric_reliability()].
#'
#' @param bio Optional formula of known biology (e.g. `~ time`); adds
#'   `bio_r2` and makes `plate_r2` account for these variables.
#' @param adjust Optional nuisance formula for `bio_r2` (e.g. `~ subject`).
#' @param within Strata for `plate_r2` (default: the variables in `bio`).
#' @param truth Include truth metrics when the data has a `truth` column
#'   (simulations)?
#' @return A named list of functions.
#' @examples
#' names(default_metrics(bio = ~ time, adjust = ~ subject))
#' @export
default_metrics <- function(bio = NULL, adjust = NULL, within = NULL,
                            truth = TRUE) {
  within <- within %||% (if (!is.null(bio)) all.vars(bio) else NULL)
  m <- list(
    reliability_plate = function(x, holdout)
      metric_reliability(x, materials = holdout, between = "plate"),
    reliability_batch = function(x, holdout) {
      if (data.table::uniqueN(x$data$batch) < 2) return(NULL)
      metric_reliability(x, materials = holdout, between = "batch")
    },
    plate_r2 = function(x, holdout) metric_plate_effect(x, within = within),
    detected = function(x, holdout) metric_detection(x)
  )
  if (!is.null(bio)) {
    m$bio_r2 <- function(x, holdout) metric_bio(x, bio, adjust = adjust)
  }
  if (truth) {
    m$truth <- function(x, holdout) {
      if (all(is.na(x$data$truth))) return(NULL)
      metric_truth(x)
    }
  }
  m
}

#' Benchmark pre-processing pipelines
#'
#' Runs every pipeline on the same data and scores the results with
#' scale-free metrics, so that strategies producing different scales (log2
#' FI, multiple-of-blank, concentrations, ...) can be compared directly.
#'
#' @section Hold-out material:
#' The key to an honest comparison is that material used to *fit* a
#' correction must not be used to *judge* it. Rows selected by `holdout`
#' (e.g. one control material, or the bridging samples) get
#' `use_for_fit = FALSE` before each pipeline runs, so no method can use them
#' for estimation; the reliability metrics are then computed on exactly these
#' rows. Give a named list of selectors to cross-fit (each material held out
#' in turn); results are reported per fold and summarised over folds.
#'
#' @param x A [cyto] object (raw, i.e. before processing).
#' @param pipelines A (named) list of [cyto_pipeline()]s, e.g. from
#'   [expand_pipelines()].
#' @param holdout `NULL`, a selector (see [selectors]), or a named list of
#'   selectors (folds).
#' @param metrics List of metric functions, see [default_metrics()]. A
#'   metric with an argument called `test` receives the test set, e.g. a
#'   [model_statistic()].
#' @param test Optional unseen test set: a selector or a vector of
#'   `sample_id`s (see [split_samples()]). All rows of these samples get
#'   `use_for_fit = FALSE` in every fold, so no pipeline is fitted on them.
#' @param keep Keep the processed objects in the result?
#' @param verbose Print progress?
#' @return An object of class `cyto_benchmark`: a list with `results` (long
#'   table: pipeline, fold, metric, analyte, value, n), `summary` (mean per
#'   pipeline x metric with ranks and an overall score), `errors`,
#'   `warnings` (raised while running each pipeline), and optionally
#'   `objects`. Has `print()`, `summary()` and `plot()` methods.
#' @examples
#' pipes <- list(
#'   log2_raw   = cyto_pipeline(pipe_transform("log2")),
#'   mob        = cyto_pipeline(pipe_transform("mob")),
#'   ctrl_anchor = cyto_pipeline(pipe_transform("log2"),
#'                               pipe_normalize("anchor", anchors = "ctrl_low",
#'                                              scope = "batch"))
#' )
#' bm <- benchmark_pipelines(cyto_example, pipes, holdout = "ctrl_high",
#'                           metrics = default_metrics(bio = ~ time),
#'                           verbose = FALSE)
#' bm
#' plot(bm)
#' @export
benchmark_pipelines <- function(x, pipelines, holdout = NULL,
                                metrics = default_metrics(), keep = FALSE,
                                verbose = TRUE, test = NULL) {
  .check_cyto(x)
  if (!is.null(test)) {
    test <- .test_ids(x$data, test)
    x <- .hide_rows(x, test = test)
  }
  if (inherits(pipelines, c("cyto_pipeline", "cyto_step"))) {
    pipelines <- list(pipelines)
  }
  nms <- names(pipelines)
  if (is.null(nms)) nms <- rep("", length(pipelines))
  for (i in seq_along(pipelines)) {
    if (!nzchar(nms[i])) nms[i] <- attr(pipelines[[i]], "name") %||%
        paste0("pipeline", i)
  }
  names(pipelines) <- make.unique(nms)
  folds <- if (is.null(holdout)) list(all = NULL) else
    if (is.list(holdout) && !inherits(holdout, "formula")) holdout else
      list(holdout = holdout)
  if (is.null(names(folds))) names(folds) <- paste0("fold", seq_along(folds))

  res <- list()
  errs <- list()
  warns <- list()
  objs <- list()
  for (f in names(folds)) {
    sel <- folds[[f]]
    xf <- x
    if (!is.null(sel)) {
      xf <- .cyto_copy(x)
      h <- select_rows(xf$data, sel)
      if (!any(h)) stop("Hold-out '", f, "' selects no rows.", call. = FALSE)
      xf$data[h, use_for_fit := FALSE]
    }
    for (p in names(pipelines)) {
      if (verbose) message("[", f, "] ", p)
      wl <- character()
      out <- tryCatch(withCallingHandlers(
        run_pipeline(xf, pipelines[[p]]),
        warning = function(w) {
          wl <<- c(wl, conditionMessage(w))
          invokeRestart("muffleWarning")
        }),
        error = function(e) e)
      if (length(wl)) {
        warns[[length(warns) + 1]] <- data.table::data.table(
          pipeline = p, fold = f, message = unique(wl))
      }
      if (inherits(out, "error")) {
        errs[[length(errs) + 1]] <- data.table::data.table(
          pipeline = p, fold = f, message = conditionMessage(out))
        next
      }
      if (keep) objs[[paste(f, p, sep = "::")]] <- out
      for (mn in names(metrics)) {
        mt <- tryCatch(.call_metric(metrics[[mn]], out, sel, test),
                       error = function(e) e)
        if (inherits(mt, "error")) {
          errs[[length(errs) + 1]] <- data.table::data.table(
            pipeline = p, fold = f,
            message = paste0("metric ", mn, ": ", conditionMessage(mt)))
          next
        }
        if (is.null(mt) || !nrow(mt)) next
        mt <- data.table::as.data.table(mt)
        mt[, `:=`(pipeline = p, fold = f)]
        res[[length(res) + 1]] <- mt
      }
    }
  }
  results <- data.table::rbindlist(res, use.names = TRUE, fill = TRUE)
  if (nrow(results)) {
    data.table::setcolorder(results, c("pipeline", "fold", "metric", "analyte",
                                       "value", "n"))
  }
  errors <- data.table::rbindlist(errs)
  bm <- structure(list(results = results, errors = errors,
                       warnings = data.table::rbindlist(warns),
                       pipelines = pipelines, holdout = folds, test = test,
                       objects = if (keep) objs else NULL),
                  class = "cyto_benchmark")
  bm$summary <- summary(bm)
  bm
}

# Metrics are called as f(x, holdout); those with a `test` argument also get
# the test set
.call_metric <- function(f, x, holdout, test) {
  if ("test" %in% names(formals(args(f)))) f(x, holdout, test = test) else
    f(x, holdout)
}

#' @export
summary.cyto_benchmark <- function(object, stat = c("mean", "median"),
                                   weights = NULL, ...) {
  stat <- match.arg(stat)
  f <- if (stat == "mean") mean else stats::median
  r <- object$results
  if (!nrow(r)) return(data.table::data.table())
  # average over folds first, then summarise over analytes
  r <- r[is.finite(value), .(value = mean(value)),
         by = .(pipeline, metric, analyte)]
  s <- r[, .(value = f(value), n_analytes = .N), by = .(pipeline, metric)]
  s[, direction := .direction_of(metric)]
  s[, rank := rank(-direction * value, ties.method = "average"), by = metric]
  w <- if (is.null(weights)) NULL else weights
  sc <- s[, .(score = if (is.null(w)) mean(rank) else
    stats::weighted.mean(rank, ifelse(metric %in% names(w), w[metric], 1))),
    by = pipeline]
  wide <- data.table::dcast(s, pipeline ~ metric, value.var = "value")
  wide <- merge(wide, sc, by = "pipeline")
  data.table::setnames(wide, "score", "mean_rank")
  wide[order(mean_rank)]
}

#' @export
print.cyto_benchmark <- function(x, digits = 3, ...) {
  cat("<cyto_benchmark> ", length(x$pipelines), " pipelines, ",
      length(x$holdout), " fold(s)", if (!is.null(x$holdout[[1]]))
        paste0(" (hold-out: ", paste(names(x$holdout), collapse = ", "), ")"),
      "\n", sep = "")
  cat("Mean over analytes (and folds); lower mean_rank is better.\n")
  cat("Directions: reliability/detected/bio_r2/truth_rho higher is better; ",
      "plate_r2/truth_rmse lower is better.\n", sep = "")
  s <- data.table::copy(x$summary)
  num <- names(s)[vapply(s, is.numeric, TRUE)]
  s[, (num) := lapply(.SD, round, digits), .SDcols = num]
  print(s)
  if (nrow(x$errors)) {
    cat("\nErrors:\n")
    print(x$errors)
  }
  invisible(x)
}

#' Benchmark pipelines across simulated scenarios
#'
#' Method development needs data where the answer is known. This runs
#' [benchmark_pipelines()] on data sets generated by [simulate_cyto()] under
#' several scenarios (each a list of arguments to `simulate_cyto()`, e.g.
#' stronger plate effects, a confounded layout, more bad standards, fewer
#' blanks), with `reps` replicate data sets per scenario. The result shows
#' which strategies are robust under which failure modes -- before spending
#' real samples on them.
#'
#' @param pipelines A (named) list of pipelines.
#' @param scenarios A named list; each element is a list of arguments passed
#'   to [simulate_cyto()].
#' @param reps Replicate data sets per scenario.
#' @param holdout,metrics Passed to [benchmark_pipelines()].
#' @param seed Base seed; data set *r* of scenario *s* uses `seed + 1000*s + r`.
#' @param verbose Print progress?
#' @return A `cyto_benchmark` object whose `results` have extra columns
#'   `scenario` and `rep`; `summary` is computed per scenario.
#' @examples
#' \donttest{
#' pipes <- list(
#'   mob = cyto_pipeline(pipe_transform("mob")),
#'   anchor = cyto_pipeline(pipe_transform("log2"),
#'                          pipe_normalize("anchor", anchors = "ctrl_low",
#'                                         scope = "batch"),
#'                          pipe_correct_batch("bridge"))
#' )
#' sc <- list(default = list(),
#'            strong_plate = list(plate_sd = c(gain = 0.5, shift = 0.5,
#'                                             background = 0.5)))
#' bs <- benchmark_simulations(pipes, sc, reps = 1, verbose = FALSE)
#' bs$summary
#' }
#' @export
benchmark_simulations <- function(pipelines,
                                  scenarios = list(default = list()),
                                  reps = 3, holdout = "ctrl_high",
                                  metrics = default_metrics(
                                    bio = ~ time, adjust = ~ subject),
                                  seed = 1, verbose = TRUE) {
  if (is.null(names(scenarios)) || any(!nzchar(names(scenarios)))) {
    stop("'scenarios' must be a named list.", call. = FALSE)
  }
  res <- list()
  errs <- list()
  warns <- list()
  for (si in seq_along(scenarios)) {
    s <- names(scenarios)[si]
    for (r in seq_len(reps)) {
      if (verbose) message("Scenario '", s, "', data set ", r, "/", reps)
      args <- utils::modifyList(scenarios[[si]],
                                list(seed = seed + 1000 * si + r))
      x <- do.call(simulate_cyto, args)
      bm <- benchmark_pipelines(x, pipelines, holdout = holdout,
                                metrics = metrics, verbose = FALSE)
      if (nrow(bm$results)) {
        res[[length(res) + 1]] <- bm$results[, `:=`(scenario = s, rep = r)]
      }
      if (nrow(bm$errors)) {
        errs[[length(errs) + 1]] <- bm$errors[, `:=`(scenario = s, rep = r)]
      }
      if (nrow(bm$warnings)) {
        warns[[length(warns) + 1]] <- bm$warnings[, `:=`(scenario = s, rep = r)]
      }
    }
  }
  results <- data.table::rbindlist(res, fill = TRUE)
  out <- structure(list(results = results,
                        errors = data.table::rbindlist(errs, fill = TRUE),
                        warnings = data.table::rbindlist(warns, fill = TRUE),
                        pipelines = pipelines, holdout = list(holdout),
                        scenarios = scenarios, objects = NULL),
                   class = c("cyto_benchmark_sim", "cyto_benchmark"))
  out$summary <- summary(out)
  out
}

#' @export
summary.cyto_benchmark_sim <- function(object, stat = c("mean", "median"),
                                       ...) {
  stat <- match.arg(stat)
  f <- if (stat == "mean") mean else stats::median
  r <- object$results
  if (!nrow(r)) return(data.table::data.table())
  s <- r[is.finite(value), .(value = f(value)),
         by = .(scenario, pipeline, metric)]
  s[, direction := .direction_of(metric)]
  s[, rank := rank(-direction * value, ties.method = "average"),
    by = .(scenario, metric)]
  sc <- s[, .(mean_rank = mean(rank)), by = .(scenario, pipeline)]
  wide <- data.table::dcast(s, scenario + pipeline ~ metric, value.var = "value")
  wide <- merge(wide, sc, by = c("scenario", "pipeline"))
  wide[order(scenario, mean_rank)]
}

#' @export
print.cyto_benchmark_sim <- function(x, digits = 3, ...) {
  cat("<cyto_benchmark_sim> ", length(x$pipelines), " pipelines x ",
      length(x$scenarios), " scenario(s)\n", sep = "")
  cat("Mean over analytes and replicate data sets; lower mean_rank is better.\n")
  s <- data.table::copy(x$summary)
  num <- names(s)[vapply(s, is.numeric, TRUE)]
  s[, (num) := lapply(.SD, round, digits), .SDcols = num]
  print(s)
  if (nrow(x$errors)) {
    cat("\nErrors (", nrow(x$errors), "):\n", sep = "")
    print(utils::head(x$errors))
  }
  invisible(x)
}

#' @rdname benchmark_simulations
#' @param x A `cyto_benchmark_sim` object.
#' @param metric Metric to show in the heat map (default `"mean_rank"`).
#' @param ... Unused.
#' @export
plot.cyto_benchmark_sim <- function(x, metric = "mean_rank", ...) {
  s <- data.table::as.data.table(x$summary)
  if (!metric %in% names(s)) stop("Unknown metric '", metric, "'.",
                                  call. = FALSE)
  s[, fill := get(metric)]
  ord <- s[, .(m = mean(mean_rank)), by = pipeline][order(-m), pipeline]
  s[, pipeline := factor(pipeline, levels = ord)]
  ggplot2::ggplot(s, ggplot2::aes(x = .data$scenario, y = .data$pipeline,
                                  fill = .data$fill)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = round(.data$fill, 2)), size = 3) +
    ggplot2::scale_fill_gradient(low = if (metric %in% c("mean_rank",
                                                        "plate_r2",
                                                        "truth_rmse"))
      "#2c7fb8" else "#f7fbff",
      high = if (metric %in% c("mean_rank", "plate_r2", "truth_rmse"))
        "#f7fbff" else "#2c7fb8") +
    ggplot2::labs(x = "scenario", y = NULL, fill = metric) +
    .theme()
}

#' Split repeated material into hold-out folds
#'
#' Bridging samples are usually the only anchor for batch correction, so
#' holding all of them out leaves nothing to fit. `holdout_folds()` splits the
#' selected material (by `sample_id`) into `k` groups; each fold holds out one
#' group, so the correction is fitted on the others and judged on unseen
#' material. Pass the result (possibly combined with other selectors) as
#' `holdout` to [benchmark_pipelines()].
#'
#' @param x A [cyto] object.
#' @param materials Selector of the material to split (default `"bridge"`).
#' @param k Number of folds.
#' @param seed Random seed for the split.
#' @param prefix Prefix for fold names.
#' @return A named list of character vectors of `sample_id`s.
#' @examples
#' f <- holdout_folds(cyto_example, "bridge", k = 2)
#' lengths(f)
#' @export
holdout_folds <- function(x, materials = "bridge", k = 2, seed = 1,
                          prefix = "fold") {
  .check_cyto(x)
  ids <- sort(unique(x$data$sample_id[select_rows(x$data, materials)]))
  if (length(ids) < k) stop("Fewer materials (", length(ids), ") than folds.",
                            call. = FALSE)
  old <- if (exists(".Random.seed", envir = globalenv()))
    get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()),
          add = TRUE)
  set.seed(seed)
  grp <- sample(rep_len(seq_len(k), length(ids)))
  out <- split(ids, grp)
  names(out) <- paste0(prefix, seq_len(k))
  out
}
