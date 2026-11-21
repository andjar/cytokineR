.required_cols <- c("plate", "well", "sample_id", "sample_type", "analyte", "fi")
.sample_types <- c("standard", "blank", "control", "unknown")

# Columns that are specific to a well x analyte measurement (everything else
# is treated as well-level sample metadata, e.g. study, subject, time).
.measure_cols <- c(".row", "analyte", "fi", "bead_count", "expected_conc",
                   "vendor_conc", "vendor_flag", "value", "censored",
                   "limit_low", "limit_high", "qc", "use_for_fit", "signal",
                   "truth")

.scales <- c("log2_fi", "log2_ratio", "log2_conc")

#' Create a cyto object
#'
#' A `cyto` object is the container that every cytokineR function takes and
#' returns. It holds one row per *well x analyte* in a long table, plus a log of
#' all processing steps and the estimates each step produced.
#'
#' @section Columns:
#' Required: `plate`, `well`, `sample_id`, `sample_type` (one of `"standard"`,
#' `"blank"`, `"control"`, `"unknown"`), `analyte`, `fi` (raw median
#' fluorescence intensity).
#'
#' Optional (filled with defaults when absent): `batch` (default `"1"`; use it
#' for kit lot / run campaign), `level` (standard level, 1 = most
#' concentrated), `expected_conc` (assigned concentration of standards),
#' `dilution` (sample dilution factor, default 1), `bead_count`, `vendor_conc`
#' and `vendor_flag` (concentrations and flags reported by the instrument
#' software: `""`, `"<"`, `">"`, `"*"`), and `truth` (true log2 concentration,
#' only for simulated data).
#'
#' Working columns maintained by the pipeline: `value` (current value, always
#' log2-scaled; initialised to `log2(fi)`), `censored` (`"none"`, `"left"`,
#' `"right"`), `limit_low`/`limit_high` (censoring limits on the value scale),
#' `qc` (semicolon-separated QC flags) and `use_for_fit` (whether the row may be
#' used to *estimate* normalisation/calibration parameters).
#'
#' Any other column (study, subject, time point, ...) is kept as well-level
#' metadata and can be used for grouping, stratification and plotting.
#'
#' Plate identifiers must be unique across batches; if the same plate id occurs
#' in several batches it is made unique as `"<batch>:<plate>"`.
#'
#' @param data A data frame in long format (see Columns).
#' @param fi_floor Values of `fi` below this are floored before `log2()` when
#'   initialising `value`.
#' @return An object of class `cyto`.
#' @seealso [cyto_from_long()], [cyto_from_wide()], [read_bioplex()]
#' @examples
#' d <- data.frame(
#'   plate = "P1", well = c("A1", "A2", "B1", "B2"),
#'   sample_id = c("S1", "B", "ctrl", "subj1"),
#'   sample_type = c("standard", "blank", "control", "unknown"),
#'   analyte = "IL-6", fi = c(12000, 15, 800, 240),
#'   expected_conc = c(5000, NA, NA, NA), level = c(1, NA, NA, NA)
#' )
#' x <- cyto(d)
#' x
#' @export
cyto <- function(data, fi_floor = 1) {
  dt <- data.table::as.data.table(data)
  dt <- data.table::copy(dt)
  miss <- setdiff(.required_cols, names(dt))
  if (length(miss)) {
    stop("Missing required column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  }
  defaults <- list(batch = "1", level = NA_real_, expected_conc = NA_real_,
                   dilution = 1, bead_count = NA_real_, vendor_conc = NA_real_,
                   vendor_flag = "", truth = NA_real_)
  for (nm in names(defaults)) {
    if (!nm %in% names(dt)) data.table::set(dt, j = nm, value = defaults[[nm]])
  }
  for (nm in c("plate", "batch", "well", "sample_id", "analyte")) {
    data.table::set(dt, j = nm, value = as.character(dt[[nm]]))
  }
  st <- tolower(trimws(as.character(dt$sample_type)))
  bad <- setdiff(unique(st), .sample_types)
  if (length(bad)) {
    stop("Unknown sample_type value(s): ", paste(bad, collapse = ", "),
         ". Use one of: ", paste(.sample_types, collapse = ", "), call. = FALSE)
  }
  data.table::set(dt, j = "sample_type", value = st)
  for (nm in c("fi", "level", "expected_conc", "dilution", "bead_count",
               "vendor_conc", "truth")) {
    data.table::set(dt, j = nm, value = suppressWarnings(as.numeric(dt[[nm]])))
  }
  dt[is.na(dilution), dilution := 1]
  dt[is.na(vendor_flag), vendor_flag := ""]

  # plate ids must be unique across batches
  pb <- unique(dt[, .(plate, batch)])
  dup <- pb[, .N, by = plate][N > 1, plate]
  if (length(dup)) {
    message("Plate id(s) shared across batches made unique as '<batch>:<plate>'.")
    dt[, plate := paste0(batch, ":", plate)]
  }

  if (!"value" %in% names(dt)) {
    dt[, value := log2(pmax(fi, fi_floor))]
  }
  if (!"censored" %in% names(dt)) dt[, censored := "none"]
  if (!"limit_low" %in% names(dt)) dt[, limit_low := NA_real_]
  if (!"limit_high" %in% names(dt)) dt[, limit_high := NA_real_]
  if (!"qc" %in% names(dt)) dt[, qc := ""]
  if (!"use_for_fit" %in% names(dt)) dt[, use_for_fit := !is.na(fi)]
  dt[, .row := seq_len(.N)]
  data.table::setcolorder(dt, intersect(
    c(".row", "batch", "plate", "well", "sample_type", "sample_id", "level",
      "analyte", "fi", "value", "censored"), names(dt)))
  new_cyto(dt)
}

new_cyto <- function(data, scale = "log2_fi", steps = list(), fits = list(),
                     meta = list()) {
  structure(list(data = data, scale = scale, steps = steps, fits = fits,
                 meta = meta), class = "cyto")
}

#' @rdname cyto
#' @param x An object.
#' @export
is_cyto <- function(x) inherits(x, "cyto")

.check_cyto <- function(x) {
  if (!is_cyto(x)) stop("Expected a cyto object.", call. = FALSE)
  invisible(x)
}

# copy-on-modify helper: every verb starts with this
.cyto_copy <- function(x) {
  x$data <- data.table::copy(x$data)
  x
}

#' @export
print.cyto <- function(x, ...) {
  d <- x$data
  cat("<cyto> ", format(nrow(d), big.mark = ","), " measurements | ",
      data.table::uniqueN(d$analyte), " analytes | ",
      data.table::uniqueN(d$plate), " plates | ",
      data.table::uniqueN(d$batch), " batch(es)\n", sep = "")
  tab <- d[, .(wells = data.table::uniqueN(paste(plate, well))),
           by = sample_type]
  cat("Wells: ", paste0(tab$sample_type, "=", tab$wells, collapse = ", "),
      "\n", sep = "")
  cat("Current value scale: ", x$scale, "\n", sep = "")
  nc <- d[censored != "none", .N]
  if (nc > 0) {
    cat("Censored: ", nc, " (", round(100 * nc / nrow(d), 1), "%); ",
        d[censored != "none" & !is.na(value), .N], " imputed\n", sep = "")
  }
  nq <- d[qc != "", .N]
  if (nq > 0) cat("QC-flagged: ", nq, "\n", sep = "")
  if (length(x$steps)) {
    cat("Steps:\n")
    for (i in seq_along(x$steps)) {
      s <- x$steps[[i]]
      cat(sprintf("  %d. %-10s %-14s %s\n", i, s$stage, s$method, s$args_txt))
    }
  } else {
    cat("Steps: none\n")
  }
  invisible(x)
}

#' @export
summary.cyto <- function(object, ...) {
  d <- object$data
  d[, .(n = .N,
        median_fi = stats::median(fi, na.rm = TRUE),
        median_value = stats::median(value, na.rm = TRUE),
        pct_censored = round(100 * mean(censored != "none"), 1),
        pct_qc_flag = round(100 * mean(qc != ""), 1)),
    by = .(analyte, sample_type)][order(analyte, sample_type)]
}

#' Access the data in a cyto object
#'
#' `cyto_data()` returns the long table (one row per well x analyte) as a
#' plain data frame. `cyto_wide()` pivots one value column to a wide table
#' with one column per analyte -- the usual input for downstream statistics.
#'
#' @param x A [cyto] object.
#' @param types Sample types to keep (default all for `cyto_data()`,
#'   `"unknown"` for `cyto_wide()`).
#' @param value Column to pivot (default `"value"`).
#' @param as_data_table Return a data.table instead of a data.frame?
#' @return A data frame.
#' @examples
#' head(cyto_data(cyto_example))
#' w <- cyto_wide(cyto_example)
#' dim(w)
#' @export
cyto_data <- function(x, types = NULL, as_data_table = FALSE) {
  .check_cyto(x)
  d <- data.table::copy(x$data)
  if (!is.null(types)) d <- d[sample_type %in% types]
  if (as_data_table) d else as.data.frame(d)
}

#' @rdname cyto_data
#' @export
cyto_wide <- function(x, value = "value", types = "unknown",
                      as_data_table = FALSE) {
  .check_cyto(x)
  d <- x$data
  if (!is.null(types)) d <- d[sample_type %in% types]
  id_cols <- setdiff(names(d), .measure_cols)
  if (!value %in% names(d)) stop("No column '", value, "'.", call. = FALSE)
  f <- stats::as.formula(paste(paste0("`", id_cols, "`", collapse = " + "),
                               "~ analyte"))
  w <- data.table::dcast(d, f, value.var = value)
  if (as_data_table) w else as.data.frame(w)
}

#' Processing log of a cyto object
#'
#' @param x A [cyto] object.
#' @return `cyto_steps()`: a data frame with one row per step (stage, method,
#'   arguments, scale before/after, number of values changed, note, time).
#' @examples
#' x <- cyto_transform(cyto_example, "mob")
#' cyto_steps(x)
#' head(cyto_fit(x))
#' @export
cyto_steps <- function(x) {
  .check_cyto(x)
  if (!length(x$steps)) {
    return(data.frame(step = integer(), stage = character(),
                      method = character(), args = character()))
  }
  do.call(rbind, lapply(seq_along(x$steps), function(i) {
    s <- x$steps[[i]]
    data.frame(step = i, stage = s$stage, method = s$method,
               args = s$args_txt, scale_before = s$scale_before,
               scale_after = s$scale_after, n_changed = s$n_changed,
               note = s$note %||% "", time = format(s$time),
               stringsAsFactors = FALSE)
  }))
}

#' @rdname cyto_steps
#' @param stage Return the fit of the last step of this stage; if `NULL`, the
#'   last step overall.
#' @param step Step number (overrides `stage`).
#' @return `cyto_fit()`: whatever the method stored as its estimates
#'   (typically a data frame), or `NULL`.
#' @export
cyto_fit <- function(x, stage = NULL, step = NULL) {
  .check_cyto(x)
  if (!length(x$steps)) return(NULL)
  if (is.null(step)) {
    stages <- vapply(x$steps, `[[`, "", "stage")
    idx <- if (is.null(stage)) length(stages) else which(stages == stage)
    if (!length(idx)) return(NULL)
    step <- max(idx)
  }
  x$fits[[step]]
}

`%||%` <- function(a, b) if (is.null(a)) b else a

.format_args <- function(args) {
  if (!length(args)) return("")
  parts <- vapply(names(args), function(nm) {
    v <- args[[nm]]
    txt <- if (is.function(v)) {
      "<function>"
    } else if (inherits(v, "formula")) {
      paste(deparse(v), collapse = "")
    } else if (is.data.frame(v)) {
      sprintf("<data.frame %dx%d>", nrow(v), ncol(v))
    } else if (is.atomic(v) && length(v) <= 6) {
      paste(deparse(v), collapse = "")
    } else {
      sprintf("<%s[%d]>", class(v)[1], length(v))
    }
    if (nchar(txt) > 60) txt <- paste0(substr(txt, 1, 57), "...")
    paste0(nm, "=", txt)
  }, "")
  paste(parts, collapse = ", ")
}

.add_step <- function(x, stage, method, args, fit = NULL, note = NULL,
                      scale_before, n_changed = NA_integer_) {
  i <- length(x$steps) + 1L
  x$steps[[i]] <- list(stage = stage, method = method, args = args,
                       args_txt = .format_args(args),
                       scale_before = scale_before, scale_after = x$scale,
                       n_changed = n_changed, note = note, time = Sys.time())
  x$fits[i] <- list(fit)
  x
}

# Run a registered method and record it
.run_stage <- function(x, stage, method, args) {
  .check_cyto(x)
  m <- .get_method(stage, method)
  res <- do.call(m$fun, c(list(x = x), args))
  if (!is.list(res) || is.null(res$value)) {
    stop("Method '", method, "' did not return a list with a 'value' element.",
         call. = FALSE)
  }
  x <- .cyto_copy(x)
  d <- x$data
  if (length(res$value) != nrow(d)) {
    stop("Method '", method, "' returned ", length(res$value),
         " values for ", nrow(d), " rows.", call. = FALSE)
  }
  old <- d$value
  data.table::set(d, j = "value", value = as.numeric(res$value))
  for (nm in c("censored", "limit_low", "limit_high", "qc", "use_for_fit")) {
    if (!is.null(res[[nm]])) data.table::set(d, j = nm, value = res[[nm]])
  }
  if (!is.null(res$columns)) {
    for (nm in names(res$columns)) {
      data.table::set(d, j = nm, value = res$columns[[nm]])
    }
  }
  x$data <- d
  scale_before <- x$scale
  if (!is.null(res$scale)) x$scale <- res$scale
  changed <- sum(!((is.na(old) & is.na(d$value)) |
                     (!is.na(old) & !is.na(d$value) &
                        abs(old - d$value) < 1e-12)))
  .add_step(x, stage, method, args, fit = res$fit, note = res$note,
            scale_before = scale_before, n_changed = changed)
}

#' Subset a cyto object
#'
#' Keep rows matching a condition evaluated in the long data (e.g.
#' `sample_type != "control"`, `analyte %in% c("IL-6", "IL-8")`,
#' `batch == "2"`). The subsetting is recorded in the step log.
#'
#' @param x A [cyto] object.
#' @param ... Logical conditions (combined with `&`), evaluated in the data.
#' @return A [cyto] object.
#' @examples
#' cyto_filter(cyto_example, analyte %in% c("IL-6", "IL-8"), batch == "1")
#' @export
cyto_filter <- function(x, ...) {
  .check_cyto(x)
  conds <- eval(substitute(alist(...)))
  keep <- rep(TRUE, nrow(x$data))
  for (cnd in conds) {
    k <- eval(cnd, x$data, parent.frame())
    k[is.na(k)] <- FALSE
    keep <- keep & k
  }
  x <- .cyto_copy(x)
  n0 <- nrow(x$data)
  x$data <- x$data[keep]
  x$data[, .row := seq_len(.N)]
  txt <- paste(vapply(conds, function(e) paste(deparse(e), collapse = ""), ""),
               collapse = " & ")
  .add_step(x, "filter", "filter", list(condition = txt),
            note = sprintf("kept %d of %d rows", nrow(x$data), n0),
            scale_before = x$scale, n_changed = n0 - nrow(x$data))
}

#' Combine cyto objects
#'
#' Row-binds several [cyto] objects, typically one per plate or per batch. All
#' objects must be on the same value scale. Step logs are not merged (the
#' result starts a fresh log noting the combination).
#'
#' @param ... [cyto] objects, or a single list of them.
#' @return A [cyto] object.
#' @examples
#' a <- cyto_filter(cyto_example, batch == "1")
#' b <- cyto_filter(cyto_example, batch == "2")
#' cyto_bind(a, b)
#' @export
cyto_bind <- function(...) {
  objs <- list(...)
  if (length(objs) == 1L && is.list(objs[[1]]) && !is_cyto(objs[[1]])) {
    objs <- objs[[1]]
  }
  lapply(objs, .check_cyto)
  scales <- unique(vapply(objs, `[[`, "", "scale"))
  if (length(scales) > 1) {
    stop("Cannot combine objects on different scales: ",
         paste(scales, collapse = ", "), call. = FALSE)
  }
  d <- data.table::rbindlist(lapply(objs, `[[`, "data"), fill = TRUE)
  d[, .row := seq_len(.N)]
  for (nm in c("qc", "vendor_flag")) d[is.na(get(nm)), (nm) := ""]
  d[is.na(censored), censored := "none"]
  d[is.na(use_for_fit), use_for_fit := !is.na(fi)]
  x <- new_cyto(d, scale = scales)
  .add_step(x, "combine", "bind", list(n_objects = length(objs)),
            scale_before = scales, n_changed = 0L)
}

#' Collapse replicate wells
#'
#' Averages replicate wells of the same sample on the same plate (on the log2
#' value scale), returning one row per plate x sample x analyte. Censoring is
#' kept if all replicates are censored; otherwise censored replicates are
#' dropped from the mean.
#'
#' @param x A [cyto] object.
#' @param fun Summary function applied to replicate values (default mean).
#' @return A [cyto] object; `well` holds the comma-joined replicate wells.
#' @examples
#' x <- cyto_collapse_replicates(cyto_example)
#' x
#' @export
cyto_collapse_replicates <- function(x, fun = mean) {
  .check_cyto(x)
  d <- data.table::copy(x$data)
  keys <- setdiff(names(d), c(.measure_cols, "well"))
  keys <- c(keys, "analyte")
  num <- c("fi", "bead_count", "expected_conc", "vendor_conc", "limit_low",
           "limit_high", "signal", "truth")
  num <- intersect(num, names(d))
  out <- d[, {
    obs <- censored == "none"
    use <- if (any(obs)) obs else rep(TRUE, .N)
    res <- list(well = paste(sort(well), collapse = ","),
                value = if (all(is.na(value[use]))) NA_real_ else
                  fun(value[use], na.rm = TRUE),
                censored = if (any(obs)) "none" else censored[1],
                qc = paste(unique(qc[qc != ""]), collapse = ";"),
                use_for_fit = any(use_for_fit),
                vendor_flag = vendor_flag[1])
    for (nm in num) {
      v <- get(nm)
      res[[nm]] <- if (all(is.na(v))) NA_real_ else
        if (nm %in% c("limit_low", "limit_high", "expected_conc", "truth"))
          stats::median(v, na.rm = TRUE) else
            2^mean(log2(pmax(v, 1e-9)), na.rm = TRUE)
    }
    res
  }, by = keys]
  out[, .row := seq_len(.N)]
  x$data <- out
  .add_step(x, "collapse", "replicates", list(fun = fun),
            note = sprintf("%d -> %d rows", nrow(d), nrow(out)),
            scale_before = x$scale, n_changed = nrow(d) - nrow(out))
}
