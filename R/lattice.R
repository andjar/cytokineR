#' Waterfall lattice: what did each pre-processing step contribute?
#'
#' A waterfall chart, or a table of results "after each step", answers *how
#' much did each step change the result?* But the answer depends on the order
#' in which the steps are switched on. When steps interact -- two plate
#' corrections that partly do the same job, an imputation that only matters
#' once values are censored -- the same step can look decisive in one order
#' and irrelevant in another.
#'
#' `waterfall_lattice()` runs the pipeline with every on/off combination of the
#' optional steps (\eqn{2^K} pipelines for \eqn{K} steps, always applied in
#' processing order) and computes a statistic for each. The result is a
#' lattice: one node per pipeline, placed at its statistic and stacked by the
#' number of steps switched on, with one edge per step switched on. The usual
#' waterfall is one route through it (`route`). Each step's contribution is
#' summarised in three ways:
#'
#' * `waterfall` -- its change along the route (what a waterfall chart shows);
#' * `shapley` -- its average change over all orders of switching steps on
#'   (the Shapley value). This share does not depend on any order, and the
#'   shares of all steps add up exactly to the total change from no optional
#'   steps to all of them;
#' * `min`, `max` -- the range of its change across all \eqn{2^{K-1}}
#'   combinations of the other steps. When the range crosses zero
#'   (`flips = TRUE`), whether the step raises or lowers the statistic depends
#'   on what else was done.
#'
#' @section The statistic:
#' `statistic` is a function of the processed [cyto] object. It may return
#' * a single number;
#' * a named numeric vector (one value per analyte);
#' * a data frame with columns `analyte` and `value`, and optionally `p` (a
#'   p-value, used to mark pipelines with p < `alpha`) -- the format of the
#'   `metric_*()` functions, e.g. [metric_effect()], [metric_bio()] or
#'   [metric_reliability()]. A table with several metrics must be filtered to
#'   one.
#'
#' If the function has an argument called `holdout`, the `holdout` selector
#' is passed to it, so the elements of [default_metrics()] can be used
#' directly.
#'
#' @param x A [cyto] object (raw, i.e. before processing).
#' @param steps A named list of [pipe_step()]s (or small [cyto_pipeline()]s)
#'   in processing order. Every combination is applied in this order.
#' @param statistic A function of a processed [cyto] object; see The
#'   statistic.
#' @param fixed Names of steps in `steps` that are always applied (e.g. the
#'   signal transform). The others are switched on and off.
#' @param route Order in which the optional steps are switched on for the
#'   waterfall; defaults to processing order.
#' @param holdout Optional selector (see [selectors]); these rows get
#'   `use_for_fit = FALSE` before every pipeline runs, as in
#'   [benchmark_pipelines()].
#' @param label Name of the statistic used in printing and plots (default:
#'   the `metric` returned by `statistic`, or `"statistic"`).
#' @param alpha Significance level for marking nodes when `statistic` returns
#'   p-values.
#' @param max_steps Maximum number of optional steps (the number of pipelines
#'   doubles with each one).
#' @param verbose Print progress?
#' @return An object of class `cyto_lattice`: a list with
#'   * `nodes`: one row per pipeline x analyte (`mask`, `n_on`, `pipeline`,
#'     `analyte`, `value`, `p`, `significant`, `route`);
#'   * `edges`: one row per step switched on in each context (`analyte`,
#'     `step`, `from`, `to`, `from_n`, `from_value`, `to_value`, `delta`,
#'     `route`);
#'   * `contributions`: one row per analyte x step (`waterfall`, `shapley`,
#'     `min`, `max`, `mean`, `n_contexts`, `flips`);
#'   * `totals`: one row per analyte (value with no optional steps and with
#'     all of them, total change, number of pipelines with p < `alpha`);
#'   * `errors`, `warnings` raised while running, and the settings.
#'
#'   Has `print()`, `summary()` (the `contributions` table) and `plot()`
#'   methods.
#' @seealso [metric_effect()] for the standardised effect of a study variable;
#'   [benchmark_pipelines()] to compare whole pipelines;
#'   `vignette("waterfall-lattice")`.
#' @references
#' Shapley LS (1953). A value for n-person games. In: Kuhn HW, Tucker AW
#' (eds), *Contributions to the Theory of Games II*, 307--317. Princeton
#' University Press.
#'
#' \enc{Grömping}{Groemping} U (2006). Relative importance for linear regression in R: the
#' package relaimpo. *Journal of Statistical Software* 17(1).
#' @examples
#' steps <- list(
#'   log2  = pipe_transform("log2"),
#'   beads = pipe_custom(function(x) flag_beads(x, 35), label = "beads"),
#'   plate = pipe_normalize("median", within = "time"),
#'   batch = pipe_correct_batch("bridge")
#' )
#' lat <- waterfall_lattice(
#'   cyto_filter(cyto_example, analyte %in% c("IL-8", "IP-10")),
#'   steps, fixed = "log2",
#'   statistic = function(x) metric_effect(x, ~ time, adjust = ~ subject)
#' )
#' lat
#' summary(lat)
#' plot(lat)
#' plot(lat, "contributions")
#' @export
waterfall_lattice <- function(x, steps, statistic, fixed = character(),
                              route = NULL, holdout = NULL, label = NULL,
                              alpha = 0.05, max_steps = 8, verbose = FALSE) {
  .check_cyto(x)
  if (!is.function(statistic)) {
    stop("'statistic' must be a function of a cyto object.", call. = FALSE)
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

  xf <- x
  if (!is.null(holdout)) {
    xf <- .cyto_copy(x)
    h <- select_rows(xf$data, holdout)
    if (!any(h)) stop("'holdout' selects no rows.", call. = FALSE)
    xf$data[h, use_for_fit := FALSE]
  }

  run <- .lattice_run(xf, steps, optional, statistic, holdout, verbose)
  .lattice_build(run, steps, optional, fixed, route, label, alpha)
}

# Name the steps and check them
.lattice_steps <- function(steps) {
  if (inherits(steps, c("cyto_step", "cyto_pipeline")) || !is.list(steps)) {
    stop("'steps' must be a named list of pipe_*() steps.", call. = FALSE)
  }
  ok <- vapply(steps, inherits, TRUE, c("cyto_step", "cyto_pipeline"))
  if (!all(ok)) {
    stop("Every element of 'steps' must be a pipe_*() step or a ",
         "cyto_pipeline().", call. = FALSE)
  }
  nms <- names(steps) %||% rep("", length(steps))
  for (i in seq_along(steps)) {
    if (!nzchar(nms[i])) {
      s <- steps[[i]]
      nms[i] <- if (inherits(s, "cyto_step")) s$label else
        attr(s, "name") %||% paste0("step", i)
    }
  }
  names(steps) <- make.unique(nms)
  steps
}

# Normalise what the statistic returned to data.table(analyte, value, p)
.lattice_stat <- function(res) {
  if (is.data.frame(res)) {
    res <- data.table::as.data.table(res)
    if (!"value" %in% names(res)) {
      stop("The statistic returned a table without a 'value' column.",
           call. = FALSE)
    }
    metric <- NULL
    if ("metric" %in% names(res)) {
      metric <- unique(res$metric)
      if (length(metric) > 1) {
        stop("The statistic returned several metrics (",
             paste(metric, collapse = ", "), "); return only one.",
             call. = FALSE)
      }
    }
    an <- if ("analyte" %in% names(res)) as.character(res$analyte) else
      rep("all", nrow(res))
    if (anyDuplicated(an)) {
      stop("The statistic returned more than one value per analyte.",
           call. = FALSE)
    }
    out <- data.table::data.table(
      analyte = an, value = as.numeric(res$value),
      p = if ("p" %in% names(res)) as.numeric(res$p) else NA_real_)
    attr(out, "metric") <- if (length(metric)) as.character(metric) else NULL
    return(out)
  }
  if (is.numeric(res)) {
    if (length(res) == 1L && is.null(names(res))) {
      return(data.table::data.table(analyte = "all", value = as.numeric(res),
                                    p = NA_real_))
    }
    if (is.null(names(res)) || anyDuplicated(names(res))) {
      stop("A numeric statistic with several values needs unique names ",
           "(one per analyte).", call. = FALSE)
    }
    return(data.table::data.table(analyte = names(res),
                                  value = as.numeric(res), p = NA_real_))
  }
  stop("The statistic must return a number, a named numeric vector or a ",
       "table with columns 'analyte' and 'value'.", call. = FALSE)
}

.call_statistic <- function(statistic, x, holdout) {
  if ("holdout" %in% names(formals(args(statistic)))) {
    statistic(x, holdout = holdout)
  } else {
    statistic(x)
  }
}

# Walk the binary tree of include/exclude decisions in processing order, so
# that shared prefixes are computed once (about 2^(K+1) step applications
# instead of K * 2^(K-1)), and evaluate the statistic at every leaf.
.lattice_run <- function(x, steps, optional, statistic, holdout, verbose) {
  K <- length(optional)
  bits <- stats::setNames(as.integer(2^(seq_len(K) - 1)), optional)
  env <- new.env(parent = emptyenv())
  env$stats <- vector("list", 2^K)
  env$metric <- NULL
  env$errors <- list()
  env$warnings <- list()
  env$done <- 0L
  pipe_label <- function(path) {
    if (length(path)) paste(path, collapse = " > ") else "(no steps)"
  }
  apply_step <- function(x, nm, path) {
    wl <- character()
    out <- tryCatch(withCallingHandlers(
      # copy first: a custom step that modifies x$data by reference must not
      # leak into the branch where the step is off
      run_pipeline(.cyto_copy(x), steps[[nm]]),
      warning = function(w) {
        wl <<- c(wl, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
      error = function(e) e)
    lab <- pipe_label(c(path, nm))
    if (length(wl)) {
      env$warnings[[length(env$warnings) + 1L]] <- data.table::data.table(
        pipeline = lab, message = unique(wl))
    }
    if (inherits(out, "error")) {
      env$errors[[length(env$errors) + 1L]] <- data.table::data.table(
        pipeline = lab, message = conditionMessage(out))
    }
    out
  }
  leaf <- function(x, mask, path) {
    env$done <- env$done + 1L
    if (verbose) message("[", env$done, "/", 2^K, "] ", pipe_label(path))
    if (inherits(x, "error")) return(invisible(NULL))
    st <- tryCatch(.lattice_stat(.call_statistic(statistic, x, holdout)),
                   error = function(e) e)
    if (inherits(st, "error")) {
      env$errors[[length(env$errors) + 1L]] <- data.table::data.table(
        pipeline = pipe_label(path),
        message = paste0("statistic: ", conditionMessage(st)))
      return(invisible(NULL))
    }
    env$metric <- env$metric %||% attr(st, "metric")
    env$stats[[mask + 1L]] <- st
    invisible(NULL)
  }
  walk <- function(x, i, mask, path) {
    if (i > length(steps)) return(leaf(x, mask, path))
    nm <- names(steps)[i]
    if (nm %in% optional) {
      walk(x, i + 1L, mask, path)
      y <- if (inherits(x, "error")) x else apply_step(x, nm, path)
      walk(y, i + 1L, bitwOr(mask, bits[[nm]]), c(path, nm))
    } else {
      y <- if (inherits(x, "error")) x else apply_step(x, nm, path)
      walk(y, i + 1L, mask, c(path, nm))
    }
  }
  walk(x, 1L, 0L, character())
  env
}

.lattice_build <- function(run, steps, optional, fixed, route, label, alpha) {
  K <- length(optional)
  bits <- stats::setNames(as.integer(2^(seq_len(K) - 1)), optional)
  masks <- seq_len(2^K) - 1L
  on <- outer(masks, bits, function(m, b) bitwAnd(m, b) > 0L)
  dim(on) <- c(length(masks), K)
  colnames(on) <- optional
  n_on <- as.integer(rowSums(on))
  pipe_names <- apply(on, 1, function(r) {
    if (any(r)) paste(optional[r], collapse = " + ") else "none"
  })
  route_masks <- c(0L, cumsum(bits[route]))

  st <- run$stats
  analytes <- unique(unlist(lapply(st, function(s) s$analyte)))
  if (!length(analytes)) {
    msg <- if (length(run$errors)) {
      data.table::rbindlist(run$errors)$message[1]
    } else "no values"
    stop("No pipeline produced a statistic (first error: ", msg, ").",
         call. = FALSE)
  }
  vals <- data.table::rbindlist(lapply(seq_along(st), function(i) {
    if (is.null(st[[i]])) return(NULL)
    cbind(mask = masks[i], st[[i]])
  }))
  nodes <- data.table::CJ(mask = masks, analyte = analytes, sorted = FALSE)
  nodes <- merge(nodes, vals, by = c("mask", "analyte"), all.x = TRUE,
                 sort = FALSE)
  nodes[, `:=`(n_on = n_on[mask + 1L], pipeline = pipe_names[mask + 1L],
               route = mask %in% route_masks)]
  has_p <- any(is.finite(nodes$p))
  nodes[, significant := if (has_p) is.finite(p) & p < alpha else NA]
  data.table::setcolorder(nodes, c("mask", "n_on", "pipeline", "analyte",
                                   "value", "p", "significant", "route"))
  data.table::setorderv(nodes, c("analyte", "mask"))

  # edges: switch step j on in every context that does not have it
  val_of <- function(a, m) {
    nodes$value[match(paste(a, m), paste(nodes$analyte, nodes$mask))]
  }
  edges <- data.table::rbindlist(lapply(optional, function(s) {
    from <- masks[!on[, s]]
    data.table::data.table(step = s, from = from,
                           to = bitwOr(from, bits[[s]]),
                           from_n = n_on[from + 1L])
  }))
  edges <- edges[, .(analyte = analytes), by = .(step, from, to, from_n)]
  edges[, `:=`(from_value = val_of(analyte, from),
               to_value = val_of(analyte, to))]
  edges[, delta := to_value - from_value]
  route_edges <- paste(route, route_masks[-length(route_masks)])
  edges[, route := paste(step, from) %in% route_edges]
  data.table::setcolorder(edges, c("analyte", "step", "from", "to", "from_n",
                                   "from_value", "to_value", "delta", "route"))

  # Shapley weight of a context with k of the other K - 1 steps on
  edges[, weight := factorial(from_n) * factorial(K - from_n - 1) /
          factorial(K)]
  contr <- edges[, .(
    waterfall = if (any(route)) delta[route] else NA_real_,
    shapley = if (anyNA(delta)) NA_real_ else sum(weight * delta),
    min = if (all(is.na(delta))) NA_real_ else min(delta, na.rm = TRUE),
    max = if (all(is.na(delta))) NA_real_ else max(delta, na.rm = TRUE),
    mean = if (all(is.na(delta))) NA_real_ else mean(delta, na.rm = TRUE),
    n_contexts = sum(is.finite(delta)),
    flips = any(delta > 1e-12, na.rm = TRUE) & any(delta < -1e-12,
                                                   na.rm = TRUE)
  ), by = .(analyte, step)]
  contr[, step := factor(step, levels = route)]
  data.table::setorderv(contr, c("analyte", "step"))
  contr[, step := as.character(step)]
  edges[, weight := NULL]

  totals <- nodes[, .(none = value[mask == 0L],
                      all = value[mask == max(masks)],
                      n_pipelines = .N,
                      n_finite = sum(is.finite(value)),
                      n_significant = if (has_p) sum(significant) else
                        NA_integer_),
                  by = analyte]
  totals[, total := all - none]
  data.table::setcolorder(totals, c("analyte", "none", "all", "total"))

  structure(list(
    nodes = nodes, edges = edges, contributions = contr, totals = totals,
    errors = data.table::rbindlist(run$errors),
    warnings = data.table::rbindlist(run$warnings),
    steps = steps, optional = optional, fixed = fixed, route = route,
    label = label %||% run$metric %||% "statistic", alpha = alpha,
    has_p = has_p
  ), class = "cyto_lattice")
}

#' @export
print.cyto_lattice <- function(x, digits = 3, ...) {
  K <- length(x$optional)
  na <- data.table::uniqueN(x$nodes$analyte)
  cat("<cyto_lattice> ", K, " optional step", if (K > 1) "s", " -> ", 2^K,
      " pipelines | ", na, " analyte", if (na > 1) "s", "\n", sep = "")
  cat("Statistic: ", x$label, "\n", sep = "")
  if (length(x$fixed)) {
    cat("Always applied: ", paste(x$fixed, collapse = ", "), "\n", sep = "")
  }
  cat("Waterfall route: ", paste(x$route, collapse = " > "), "\n", sep = "")
  ct <- data.table::copy(x$contributions)
  if (na == 1) {
    cat("\nContribution of each step (waterfall: along the route; shapley: ",
        "averaged over all orders;\nmin/max: over all combinations of the ",
        "other steps):\n", sep = "")
    out <- ct[, .(step, waterfall, shapley, min, max, flips)]
  } else {
    cat("\nContribution of each step, median over analytes (flips: analytes ",
        "where the sign\ndepends on the other steps):\n", sep = "")
    out <- ct[, .(waterfall = stats::median(waterfall, na.rm = TRUE),
                  shapley = stats::median(shapley, na.rm = TRUE),
                  flips = paste0(sum(flips), "/", .N)), by = step]
  }
  num <- names(out)[vapply(out, is.numeric, TRUE)]
  out[, (num) := lapply(.SD, round, digits), .SDcols = num]
  print(out, row.names = FALSE)
  tt <- x$totals
  if (na == 1) {
    cat("\nTotal change (no optional steps -> all): ",
        round(tt$total, digits), "\n", sep = "")
  }
  if (isTRUE(x$has_p)) {
    ns <- tt$n_significant
    cat("Pipelines with p < ", x$alpha, ": ", sep = "")
    if (na == 1) {
      cat(ns, " of ", 2^K, "\n", sep = "")
    } else {
      cat("median ", stats::median(ns), " of ", 2^K, " per analyte (range ",
          min(ns), "-", max(ns), ")\n", sep = "")
    }
  }
  if (nrow(x$errors)) {
    cat("\n", nrow(x$errors), " error(s); see $errors.\n", sep = "")
  }
  invisible(x)
}

#' @rdname waterfall_lattice
#' @param object A `cyto_lattice` object.
#' @param ... Unused.
#' @export
summary.cyto_lattice <- function(object, ...) {
  data.table::copy(object$contributions)
}

# Okabe-Ito colours, yellow last
.lattice_colours <- function(steps) {
  pal <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9",
           "#000000", "#F0E442")
  stats::setNames(rep_len(pal, length(steps)), steps)
}

#' @rdname waterfall_lattice
#' @param type `"lattice"` (every pipeline as a node, every step as a
#'   coloured edge, the waterfall route in bold), `"contributions"` (each
#'   step's change in every context as dots, its Shapley share as a bar and its
#'   waterfall change as a diamond) or `"both"` side by side (needs the
#'   patchwork package).
#' @param analytes Analytes to show (default: up to 6).
#' @export
plot.cyto_lattice <- function(x, type = c("lattice", "contributions", "both"),
                              analytes = NULL, ...) {
  type <- match.arg(type)
  an <- .pick_analytes(x$nodes, analytes, max_n = 6)
  if (type == "both") {
    .require("patchwork", "to draw both panels")
    return(patchwork::wrap_plots(
      .plot_lattice(x, an), .plot_contributions(x, an), widths = c(3, 2)))
  }
  if (type == "lattice") .plot_lattice(x, an) else .plot_contributions(x, an)
}

.plot_lattice <- function(x, an) {
  K <- length(x$optional)
  cols <- .lattice_colours(x$route)
  nd <- x$nodes[analyte %in% an & is.finite(value)]
  ed <- x$edges[analyte %in% an & is.finite(delta)]
  ed[, step := factor(step, levels = x$route)]
  route_masks <- c(0L, cumsum(as.integer(2^(match(x$route, x$optional) - 1))))
  lab <- nd[route == TRUE]
  idx <- match(lab$mask, route_masks) - 1L
  lab[, label := ifelse(idx == 0L, "none",
                        paste0("+", x$route[pmax(idx, 1L)]))]
  g <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = ed[route == FALSE],
      ggplot2::aes(x = .data$from_value, xend = .data$to_value,
                   y = .data$from_n, yend = .data$from_n + 1,
                   colour = .data$step),
      linewidth = 0.45, alpha = 0.55) +
    ggplot2::geom_segment(
      data = ed[route == TRUE],
      ggplot2::aes(x = .data$from_value, xend = .data$to_value,
                   y = .data$from_n, yend = .data$from_n + 1,
                   colour = .data$step),
      linewidth = 1.6, lineend = "round")
  if (isTRUE(x$has_p)) {
    nd[, sig := factor(ifelse(significant, "yes", "no"),
                       levels = c("yes", "no"))]
    g <- g +
      ggplot2::geom_point(data = nd,
                          ggplot2::aes(x = .data$value, y = .data$n_on,
                                       fill = .data$sig),
                          shape = 21, size = 2.3, stroke = 0.6) +
      ggplot2::scale_fill_manual(
        values = c(yes = "black", no = "white"), drop = FALSE,
        name = paste0("p < ", x$alpha))
  } else {
    g <- g + ggplot2::geom_point(data = nd,
                                 ggplot2::aes(x = .data$value, y = .data$n_on),
                                 shape = 21, size = 2.3, stroke = 0.6,
                                 fill = "grey30")
  }
  g +
    ggplot2::geom_text(data = lab,
                       ggplot2::aes(x = .data$value, y = .data$n_on,
                                    label = .data$label),
                       vjust = -0.9, size = 2.8) +
    ggplot2::scale_colour_manual(values = cols, drop = FALSE,
                                 name = "step switched on") +
    ggplot2::scale_y_continuous(breaks = 0:K,
                                expand = ggplot2::expansion(add = c(0.3,
                                                                    0.6))) +
    ggplot2::facet_wrap(~analyte, scales = "free_x") +
    ggplot2::labs(x = x$label, y = "optional steps switched on",
                  subtitle = "Every pipeline is a node; bold: the waterfall route") +
    .theme()
}

.plot_contributions <- function(x, an) {
  cols <- .lattice_colours(x$route)
  ed <- x$edges[analyte %in% an & is.finite(delta)]
  ct <- x$contributions[analyte %in% an]
  lv <- rev(x$route)
  ed[, step := factor(step, levels = lv)]
  ct[, step := factor(step, levels = lv)]
  ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = 0, colour = "grey55") +
    ggplot2::geom_col(data = ct[is.finite(shapley)],
                      ggplot2::aes(x = .data$shapley, y = .data$step,
                                   fill = .data$step),
                      alpha = 0.35, width = 0.6) +
    ggplot2::geom_point(data = ed,
                        ggplot2::aes(x = .data$delta, y = .data$step,
                                     colour = .data$step),
                        alpha = 0.7, size = 1.3,
                        position = ggplot2::position_jitter(height = 0.12,
                                                            width = 0,
                                                            seed = 1)) +
    ggplot2::geom_point(data = ct[is.finite(waterfall)],
                        ggplot2::aes(x = .data$waterfall, y = .data$step),
                        shape = 23, fill = "white", colour = "black",
                        size = 2.8) +
    ggplot2::scale_colour_manual(values = cols, guide = "none") +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::facet_wrap(~analyte, scales = "free_x") +
    ggplot2::labs(
      x = paste("change in", x$label), y = NULL,
      caption = paste0("Dots: change when the step is switched on, in each ",
                       "combination of the other steps.\nBar: Shapley share ",
                       "(order-free; the shares add up to the total change).",
                       " Diamond: change along the waterfall route.")) +
    .theme()
}

#' Standardised effect of a study variable
#'
#' Fits `value ~ term + adjust` per analyte among the study samples and
#' returns the coefficient of `term` divided by the residual standard
#' deviation: for a two-level group without adjustment this is Cohen's *d*
#' with a pooled SD; with `adjust = ~ subject` it is the within-subject
#' standardised change. A typical `statistic` for [waterfall_lattice()]: how
#' much does each pre-processing step move the effect you are going to report?
#'
#' @param x A processed [cyto] object.
#' @param formula One-sided formula with one term giving a single
#'   coefficient: a numeric variable or a two-level factor, e.g. `~ time` or
#'   `~ group`.
#' @param adjust Optional nuisance terms, e.g. `~ subject`.
#' @param types Sample types forming the study population.
#' @return A data.table with columns `metric` (`"effect"`), `analyte`,
#'   `value` (standardised effect), `n`, `p` (two-sided t-test of the
#'   coefficient) and `estimate` (the unstandardised coefficient, in units of
#'   the current value scale).
#' @examples
#' x <- cyto_transform(cyto_example, "log2")
#' metric_effect(x, ~ time, adjust = ~ subject)
#' @export
metric_effect <- function(x, formula, adjust = NULL, types = "unknown") {
  .check_cyto(x)
  tt <- attr(stats::terms(formula), "term.labels")
  if (length(tt) != 1L) {
    stop("'formula' must have exactly one term, e.g. ~ time.", call. = FALSE)
  }
  at <- if (is.null(adjust)) character() else
    attr(stats::terms(adjust), "term.labels")
  vars <- unique(c(all.vars(formula), all.vars(adjust)))
  d <- x$data[sample_type %in% types & is.finite(value)]
  .check_by(d, vars)
  d[, {
    r <- .effect_fit(as.data.frame(.SD), tt, at)
    list(metric = "effect", value = r[["d"]], n = as.integer(r[["n"]]),
         p = r[["p"]], estimate = r[["estimate"]])
  }, by = analyte, .SDcols = c("value", vars)][
    , .(metric, analyte, value, n, p, estimate)]
}

.effect_fit <- function(df, term, adjust) {
  df <- df[stats::complete.cases(df), , drop = FALSE]
  n <- nrow(df)
  out <- c(d = NA_real_, estimate = NA_real_, p = NA_real_, n = n)
  if (n < 4) return(out)
  fit <- tryCatch(stats::lm(stats::reformulate(c(term, adjust), "value"),
                            data = df),
                  error = function(e) NULL)
  if (is.null(fit)) return(out)
  mm <- stats::model.matrix(fit)
  col <- which(attr(mm, "assign") == 1L)
  if (length(col) != 1L) {
    stop("The term '", term, "' gives ", length(col), " coefficients; use a ",
         "numeric variable or a two-level factor.", call. = FALSE)
  }
  cf <- summary(fit)$coefficients
  nm <- colnames(mm)[col]
  sigma <- summary(fit)$sigma
  if (!nm %in% rownames(cf) || !is.finite(sigma) || sigma <= 0) return(out)
  c(d = unname(cf[nm, 1] / sigma), estimate = unname(cf[nm, 1]),
    p = unname(cf[nm, 4]), n = n)
}
