#' Pipeline steps
#'
#' Steps describe a processing operation without running it, so that
#' pipelines can be stored, printed, combined, and compared with
#' [benchmark_pipelines()]. `pipe_step()` is the general constructor; the
#' `pipe_*()` helpers fix the stage. `pipe_custom()` wraps any function that
#' takes and returns a [cyto] object (e.g. a QC filter).
#'
#' @param stage One of `"transform"`, `"normalize"`, `"calibrate"`,
#'   `"impute"`.
#' @param method Name of a registered method (see [cyto_methods()]).
#' @param ... Arguments passed to the method.
#' @param label Optional short label used when naming pipelines.
#' @return An object of class `cyto_step`.
#' @examples
#' pipe_transform("mob")
#' pipe_normalize("anchor", anchors = "control", scope = "batch")
#' pipe_custom(function(x) flag_beads(x, 35), label = "beads35")
#' @export
pipe_step <- function(stage, method, ..., label = NULL) {
  stage <- match.arg(stage, .stages)
  structure(list(stage = stage, method = method, args = list(...),
                 label = label %||% method), class = "cyto_step")
}

#' @rdname pipe_step
#' @export
pipe_transform <- function(method = "log2", ..., label = NULL) {
  pipe_step("transform", method, ..., label = label)
}

#' @rdname pipe_step
#' @param by Grouping whose effect is removed (normalize).
#' @export
pipe_normalize <- function(method = "median", by = "plate", ..., label = NULL) {
  pipe_step("normalize", method, by = by, ..., label = label)
}

#' @rdname pipe_step
#' @export
pipe_correct_batch <- function(method = "bridge", ..., label = NULL) {
  pipe_step("normalize", method, by = "batch", ..., label = label %||%
         paste0("batch-", method))
}

#' @rdname pipe_step
#' @export
pipe_calibrate <- function(method = "5pl", ..., label = NULL) {
  pipe_step("calibrate", method, ..., label = label)
}

#' @rdname pipe_step
#' @export
pipe_impute <- function(method = "tobit", ..., label = NULL) {
  pipe_step("impute", method, ..., label = label)
}

#' @rdname pipe_step
#' @param fun A function taking and returning a [cyto] object.
#' @export
pipe_custom <- function(fun, label = "custom") {
  stopifnot(is.function(fun))
  structure(list(stage = "custom", method = label, fun = fun, args = list(),
                 label = label), class = "cyto_step")
}

#' @export
print.cyto_step <- function(x, ...) {
  cat("<step> ", x$stage, ": ", x$method,
      if (length(x$args)) paste0(" (", .format_args(x$args), ")") else "",
      "\n", sep = "")
  invisible(x)
}

#' Build and run pipelines
#'
#' A pipeline is an ordered list of [pipe_step()]s. `run_pipeline()` applies it to
#' a [cyto] object; every step is recorded in the object's log. Pipelines can
#' also be applied with `x |> run_pipeline(p)`.
#'
#' @param ... Steps (or pipelines, which are flattened).
#' @param name Optional pipeline name.
#' @return `cyto_pipeline()`: an object of class `cyto_pipeline`.
#' @examples
#' p <- cyto_pipeline(pipe_transform("log2"),
#'                    pipe_normalize("anchor", anchors = "standard",
#'                                   scope = "batch", fit = "affine"),
#'                    pipe_calibrate("5pl", curve_by = "batch"),
#'                    pipe_impute("tobit"),
#'                    name = "harmonise-then-calibrate")
#' p
#' x <- run_pipeline(cyto_example, p)
#' cyto_steps(x)
#' @export
cyto_pipeline <- function(..., name = NULL) {
  steps <- list(...)
  flat <- list()
  for (s in steps) {
    if (inherits(s, "cyto_pipeline")) flat <- c(flat, unclass(s)) else
      if (inherits(s, "cyto_step")) flat <- c(flat, list(s)) else
        stop("Pipelines are built from pipe_*() steps.", call. = FALSE)
  }
  structure(flat, class = "cyto_pipeline", name = name)
}

#' @export
print.cyto_pipeline <- function(x, ...) {
  nm <- attr(x, "name")
  cat("<cyto_pipeline>", if (!is.null(nm)) nm, "\n")
  if (!length(x)) cat("  (no steps)\n")
  for (i in seq_along(x)) {
    s <- x[[i]]
    cat(sprintf("  %d. %-10s %s%s\n", i, s$stage, s$method,
                if (length(s$args)) paste0(" (", .format_args(s$args), ")")
                else ""))
  }
  invisible(x)
}

#' @rdname cyto_pipeline
#' @param x A [cyto] object.
#' @param pipeline A `cyto_pipeline` (or a single step).
#' @export
run_pipeline <- function(x, pipeline) {
  .check_cyto(x)
  if (inherits(pipeline, "cyto_step")) pipeline <- cyto_pipeline(pipeline)
  for (s in pipeline) {
    x <- if (s$stage == "custom") {
      out <- s$fun(x)
      if (!is_cyto(out)) stop("Custom step '", s$label,
                              "' did not return a cyto object.", call. = FALSE)
      out
    } else if (s$stage == "transform") {
      do.call(cyto_transform, c(list(x = x, method = s$method), s$args))
    } else if (s$stage == "calibrate") {
      do.call(cyto_calibrate, c(list(x = x, method = s$method), s$args))
    } else {
      .run_stage(x, s$stage, s$method, s$args)
    }
  }
  x
}

#' Generate a grid of pipelines
#'
#' Builds every combination of the alternatives given for each stage, which is
#' the quickest way to search the space of pre-processing strategies with
#' [benchmark_pipelines()]. Each argument is a (named) list of steps or
#' pipelines; `NULL` entries mean "skip this stage". Combinations are named by
#' joining the element names (or step labels) with `" > "`.
#'
#' @param ... Named or unnamed lists of alternatives, in processing order.
#' @param .filter Optional function `function(pipeline) logical` to drop
#'   invalid combinations.
#' @return A named list of `cyto_pipeline` objects.
#' @examples
#' grid <- expand_pipelines(
#'   signal = list(log2 = pipe_transform("log2"), mob = pipe_transform("mob")),
#'   plate = list(none = NULL, median = pipe_normalize("median"),
#'                ctrl = pipe_normalize("anchor", anchors = "control",
#'                                      scope = "batch")),
#'   batch = list(bridge = pipe_correct_batch("bridge"))
#' )
#' names(grid)
#' @export
expand_pipelines <- function(..., .filter = NULL) {
  alts <- list(...)
  alts <- lapply(alts, function(a) {
    if (inherits(a, c("cyto_step", "cyto_pipeline"))) a <- list(a)
    if (is.null(names(a))) {
      names(a) <- vapply(a, function(s) {
        if (is.null(s)) "none" else if (inherits(s, "cyto_step")) s$label else
          attr(s, "name") %||% "pipeline"
      }, "")
    }
    names(a)[!nzchar(names(a))] <- "none"
    a
  })
  idx <- expand.grid(lapply(alts, seq_along), KEEP.OUT.ATTRS = FALSE)
  out <- list()
  for (i in seq_len(nrow(idx))) {
    parts <- list()
    nms <- character()
    for (j in seq_along(alts)) {
      el <- alts[[j]][[idx[i, j]]]
      nms <- c(nms, names(alts[[j]])[idx[i, j]])
      if (!is.null(el)) parts <- c(parts, list(el))
    }
    nm <- paste(nms[nms != "none"], collapse = " > ")
    if (!nzchar(nm)) nm <- "raw"
    p <- do.call(cyto_pipeline, c(parts, list(name = nm)))
    if (is.null(.filter) || isTRUE(.filter(p))) out[[nm]] <- p
  }
  out
}
