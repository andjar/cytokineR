.registry <- new.env(parent = emptyenv())
.stages <- c("transform", "normalize", "calibrate", "impute")

#' Register a processing method
#'
#' Every processing step in cytokineR is a *method* registered under one of
#' four stages: `"transform"` (choose the signal representation),
#' `"normalize"` (remove plate/batch effects), `"calibrate"` (signal to
#' concentration) and `"impute"` (fill censored values). Registering your own
#' method makes it available to [cyto_transform()], [cyto_normalize()],
#' [cyto_calibrate()], [cyto_impute()], to [pipe_step()] in pipelines, and hence to
#' [benchmark_pipelines()].
#'
#' @section Method contract:
#' A method is a function `function(x, ...)` that receives a [cyto] object and
#' the user's arguments, and returns a list with
#'
#' * `value` (required): numeric vector, one entry per row of `x$data`, the new
#'   working value (always on a log2 scale).
#' * `scale` (optional): new scale label, e.g. `"log2_conc"`.
#' * `censored`, `limit_low`, `limit_high`, `qc`, `use_for_fit` (optional):
#'   replacement vectors for those columns.
#' * `columns` (optional): named list of additional columns to add or replace.
#' * `fit` (optional): any object describing what was estimated (typically a
#'   data frame of per-plate/per-analyte parameters). It is stored and
#'   retrievable with [cyto_fit()], which is what makes each step auditable.
#' * `note` (optional): a short character string for the step log.
#'
#' Methods must not modify `x$data` by reference; work on copies.
#' Methods that *estimate* something should only use rows with
#' `use_for_fit == TRUE`, so that [benchmark_pipelines()] can hold material out.
#'
#' @param stage One of `"transform"`, `"normalize"`, `"calibrate"`, `"impute"`.
#' @param name Method name (character scalar).
#' @param fun The method function (see contract).
#' @param description One-line description shown by [cyto_methods()].
#' @param overwrite Replace an existing method of the same name?
#' @return Invisibly, the registered entry.
#' @examples
#' # A "multiple of the plate's 75th percentile" representation
#' register_cyto_method("transform", "mo_q75", function(x, ...) {
#'   d <- x$data
#'   q <- d[, .(q75 = stats::quantile(fi[sample_type == "unknown"], 0.75,
#'                                     na.rm = TRUE)), by = .(plate, analyte)]
#'   d <- merge(d[, .(.row, plate, analyte, fi)], q, by = c("plate", "analyte"))
#'   data.table::setorder(d, .row)
#'   list(value = log2(d$fi / d$q75), scale = "log2_ratio",
#'        fit = q, note = "FI / plate 75th percentile")
#' }, description = "Ratio to the plate's 75th percentile of unknowns")
#' x <- cyto_transform(cyto_example, "mo_q75")
#' cyto_steps(x)
#' unregister_cyto_method("transform", "mo_q75")
#' @export
register_cyto_method <- function(stage, name, fun, description = "",
                                 overwrite = FALSE) {
  stage <- match.arg(stage, .stages)
  stopifnot(is.character(name), length(name) == 1L, is.function(fun))
  key <- paste(stage, name, sep = "::")
  if (exists(key, envir = .registry, inherits = FALSE) && !overwrite) {
    stop("Method '", name, "' already registered for stage '", stage,
         "'. Use overwrite = TRUE to replace it.", call. = FALSE)
  }
  entry <- list(stage = stage, name = name, fun = fun,
                description = description)
  assign(key, entry, envir = .registry)
  invisible(entry)
}

#' @rdname register_cyto_method
#' @export
unregister_cyto_method <- function(stage, name) {
  stage <- match.arg(stage, .stages)
  key <- paste(stage, name, sep = "::")
  if (exists(key, envir = .registry, inherits = FALSE)) {
    rm(list = key, envir = .registry)
  }
  invisible(NULL)
}

#' List available processing methods
#'
#' @param stage Optional stage to filter on.
#' @return A data frame with columns `stage`, `method`, `description`.
#' @examples
#' cyto_methods()
#' cyto_methods("calibrate")
#' @export
cyto_methods <- function(stage = NULL) {
  keys <- ls(.registry)
  if (!length(keys)) {
    return(data.frame(stage = character(), method = character(),
                      description = character()))
  }
  out <- do.call(rbind, lapply(keys, function(k) {
    e <- get(k, envir = .registry)
    data.frame(stage = e$stage, method = e$name,
               description = e$description, stringsAsFactors = FALSE)
  }))
  out$stage <- factor(out$stage, levels = .stages)
  out <- out[order(out$stage, out$method), , drop = FALSE]
  out$stage <- as.character(out$stage)
  if (!is.null(stage)) out <- out[out$stage %in% stage, , drop = FALSE]
  rownames(out) <- NULL
  out
}

.get_method <- function(stage, name) {
  key <- paste(stage, name, sep = "::")
  if (!exists(key, envir = .registry, inherits = FALSE)) {
    avail <- cyto_methods(stage)$method
    stop("Unknown ", stage, " method '", name, "'. Available: ",
         paste(avail, collapse = ", "), call. = FALSE)
  }
  get(key, envir = .registry)
}
