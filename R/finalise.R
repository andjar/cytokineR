#' Produce the analysis dataset with the chosen pipeline
#'
#' Applies the pipeline chosen by [select_pipeline()] to the raw data and
#' records how it was chosen in the result, so that the processed object --
#' and the files written by [export_dataset()] -- carry their own
#' documentation ([processing_report()]).
#'
#' With `refit = TRUE` (default) every sample, including the test samples, is
#' used to fit the final normalisation and batch corrections: the test set has
#' done its job (judging the choice), and more data give more stable
#' estimates. With `refit = FALSE`, the test samples stay out of the fits
#' exactly as in the evaluation, so their values are what the reported test
#' scores were computed on.
#'
#' @param x The raw [cyto] object the selection was made on (or new raw data
#'   to process the same way, e.g. a later batch).
#' @param selection A `cyto_selection` from [select_pipeline()], or a
#'   [cyto_pipeline()].
#' @param refit Use the test samples when fitting the final pipeline?
#' @return The processed [cyto] object; `x$meta$selection` holds the choice
#'   (pipeline, criteria, scores, test set, `refit`).
#' @seealso [select_pipeline()], [processing_report()], [export_dataset()]
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
#' final <- finalise_pipeline(cyto_example, sel)
#' processing_report(final)
#' head(cyto_wide(final))
#' @export
finalise_pipeline <- function(x, selection, refit = TRUE) {
  .check_cyto(x)
  if (inherits(selection, "cyto_pipeline")) {
    return(run_pipeline(x, selection))
  }
  if (!inherits(selection, "cyto_selection")) {
    stop("'selection' must come from select_pipeline() or cyto_pipeline().",
         call. = FALSE)
  }
  xin <- if (refit) x else .hide_rows(x, test = selection$test_ids)
  out <- run_pipeline(xin, selection$pipeline)
  out$meta$selection <- list(
    pipeline = selection$selected,
    steps = vapply(selection$steps, .step_text, ""),
    fixed = selection$fixed, optional = selection$optional,
    criteria = selection$criteria,
    scores = selection$scores,
    chosen = selection$scores[selected == TRUE],
    n_pipelines = nrow(selection$scores),
    n_train = selection$n_train, test_ids = selection$test_ids,
    holdout = if (is.null(selection$holdout)) NULL else
      paste(format(selection$holdout), collapse = " "),
    analytes = selection$analytes, resolution = selection$resolution,
    refit = refit,
    selected_at = selection$time
  )
  out
}

.step_text <- function(s) {
  if (inherits(s, "cyto_pipeline")) {
    return(paste(vapply(unclass(s), .step_text, ""), collapse = " > "))
  }
  paste0(s$stage, ": ", s$method,
         if (length(s$args)) paste0(" (", .format_args(s$args), ")") else "")
}

#' Processing report: how a dataset was made
#'
#' Writes a plain-language (Markdown) account of a processed [cyto] object:
#' the data, every processing step with its arguments and notes, how the
#' pipeline was chosen (when it came from [finalise_pipeline()]), and how
#' many values were censored, imputed or flagged. Paste it into a methods
#' section, or keep it next to the data ([export_dataset()] does this).
#'
#' @param x A processed [cyto] object.
#' @param file Optional path to write the report to.
#' @param title Report title.
#' @return The report as a character vector (invisibly when written to a
#'   file); printed with line breaks.
#' @examples
#' x <- cyto_example |> cyto_transform("log2") |> cyto_correct_batch("bridge")
#' processing_report(x)
#' @export
processing_report <- function(x, file = NULL,
                              title = "Processing report") {
  .check_cyto(x)
  d <- x$data
  r <- function(v, k = 3) format(round(v, k), nsmall = k)
  pct <- function(a, b) if (b > 0) sprintf("%.1f%%", 100 * a / b) else "0%"
  wells <- d[, .(n = data.table::uniqueN(paste(plate, well))),
             by = sample_type]
  study <- d[sample_type == "unknown"]
  L <- c(paste0("# ", title), "",
         paste0("Generated ", format(Sys.time(), "%Y-%m-%d %H:%M"),
                " with cytokineR ",
                as.character(utils::packageVersion("cytokineR")), "."), "",
         "## Data", "",
         paste0("- ", format(nrow(d), big.mark = ","),
                " measurements of ", data.table::uniqueN(d$analyte),
                " analytes (", paste(sort(unique(d$analyte)),
                                      collapse = ", "), ")"),
         paste0("- ", data.table::uniqueN(d$plate), " plates in ",
                data.table::uniqueN(d$batch), " batch(es)"),
         paste0("- wells: ", paste0(wells$sample_type, " ", wells$n,
                                    collapse = ", ")),
         paste0("- ", data.table::uniqueN(study$sample_id),
                " study samples"),
         paste0("- final value scale: ", x$scale), "",
         "## Processing steps", "")
  if (!length(x$steps)) {
    L <- c(L, "No processing steps were applied.")
  } else {
    for (i in seq_along(x$steps)) {
      s <- x$steps[[i]]
      L <- c(L, paste0(i, ". **", s$stage, "** -- `", s$method, "`",
                       if (nzchar(s$args_txt)) paste0(" (", s$args_txt, ")"),
                       if (!is.null(s$note) && nzchar(s$note))
                         paste0(": ", s$note)))
    }
  }
  sel <- x$meta$selection
  if (!is.null(sel)) {
    ch <- sel$chosen
    L <- c(L, "", "## How the pipeline was chosen", "",
           paste0("All ", sel$n_pipelines, " on/off combinations of the ",
                  "optional steps (", paste(sel$optional, collapse = ", "),
                  ")", if (length(sel$fixed))
                    paste0(", with ", paste(sel$fixed, collapse = ", "),
                           " always applied,"), " were compared. ",
                  length(sel$test_ids), " study samples were set aside as ",
                  "an unseen test set: no normalisation, batch correction, ",
                  "imputation, model fit or choice used them. The pipeline ",
                  "was chosen on the remaining ", sel$n_train,
                  " training samples", if (!is.null(sel$holdout))
                    paste0(" (held-out material: ", sel$holdout, ")"),
                  " and then scored on the test samples. Pipelines were ",
                  "ranked per criterion and the mean rank decided",
                  if (!is.null(sel$resolution)) paste0(
                    "; differences smaller than ",
                    paste(if (is.null(names(sel$resolution)))
                      sel$resolution else paste0(names(sel$resolution), " ",
                                                 sel$resolution),
                      collapse = ", "), " counted as ties"),
                  "; ties went to the pipeline with fewer steps."), "",
           "| criterion | metric | better | weight | training | test |",
           "|---|---|---|---|---|---|")
    for (i in seq_len(nrow(sel$criteria))) {
      cr <- sel$criteria[i]
      L <- c(L, paste0("| ", cr$criterion, " | ", cr$metric, " | ",
                       ifelse(cr$direction > 0, "higher", "lower"), " | ",
                       cr$weight, " | ",
                       r(ch[[paste0(cr$criterion, "_train")]]), " | ",
                       r(ch[[paste0(cr$criterion, "_test")]]), " |"))
    }
    L <- c(L, "",
           paste0("Chosen pipeline: **", sel$pipeline, "** (mean rank ",
                  r(ch$mean_rank_train, 2), " of ", sel$n_pipelines,
                  " in training; ranked ", ch$rank_test, " on the test ",
                  "samples)."),
           paste0("Steps: ", paste(sel$steps[names(sel$steps) %in%
                                               c(sel$fixed,
                                                 strsplit(sel$pipeline,
                                                          " \\+ ")[[1]])],
                                   collapse = "; "), "."),
           paste0("The final dataset was fitted on ",
                  if (isTRUE(sel$refit)) "all samples (training and test)."
                  else "the training samples only; the test samples were processed as new data."))
  }
  n <- nrow(study)
  cens <- study[censored != "none", .N]
  imp <- study[censored != "none" & is.finite(value), .N]
  miss <- study[!is.finite(value), .N]
  qcf <- study[qc != "", .N]
  L <- c(L, "", "## Values in the dataset (study samples)", "",
         paste0("- censored (outside the measurable range): ", cens, " (",
                pct(cens, n), "), of which imputed: ", imp),
         paste0("- missing: ", miss, " (", pct(miss, n), ")"),
         paste0("- with QC flags: ", qcf, " (", pct(qcf, n), ")"))
  per <- study[, .(censored = mean(censored != "none"),
                   imputed = mean(censored != "none" & is.finite(value))),
               by = analyte][censored > 0][order(-censored)]
  if (nrow(per)) {
    L <- c(L, "", "| analyte | censored | imputed |", "|---|---|---|",
           paste0("| ", per$analyte, " | ", sprintf("%.1f%%",
                                                     100 * per$censored),
                  " | ", sprintf("%.1f%%", 100 * per$imputed), " |"))
  }
  out <- structure(L, class = "cyto_report")
  if (!is.null(file)) {
    writeLines(L, file, useBytes = TRUE)
    return(invisible(out))
  }
  out
}

#' @export
print.cyto_report <- function(x, ...) {
  cat(x, sep = "\n")
  invisible(x)
}

#' Write the analysis dataset and its documentation
#'
#' Writes a processed [cyto] object as files that can be analysed anywhere:
#' * `<name>_wide.csv`: one row per study sample well, one column per analyte
#'   (the values to analyse);
#' * `<name>_long.csv`: one row per well x analyte with censoring, limits and
#'   QC flags, for all sample types;
#' * `<name>_report.md`: the [processing_report()];
#' * `<name>.rds`: the cyto object itself, with its full processing log.
#'
#' @param x A processed [cyto] object, e.g. from [finalise_pipeline()].
#' @param dir Output directory (created if needed).
#' @param name File name stem.
#' @param types Sample types in the wide table.
#' @return The paths written, invisibly.
#' @examples
#' x <- cyto_example |> cyto_transform("log2") |> cyto_correct_batch("bridge")
#' paths <- export_dataset(x, tempdir())
#' basename(paths)
#' @export
export_dataset <- function(x, dir, name = "cytokines", types = "unknown") {
  .check_cyto(x)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  p <- file.path(dir, paste0(name, c("_wide.csv", "_long.csv", "_report.md",
                                     ".rds")))
  names(p) <- c("wide", "long", "report", "rds")
  data.table::fwrite(cyto_wide(x, types = types, as_data_table = TRUE),
                     p[["wide"]])
  long <- cyto_data(x, as_data_table = TRUE)
  long[, imputed := censored != "none" & is.finite(value)]
  data.table::fwrite(long, p[["long"]])
  processing_report(x, file = p[["report"]])
  saveRDS(x, p[["rds"]])
  invisible(p)
}
