#' Use a statistical model as a statistic or metric
#'
#' `model_statistic()` turns a model -- [stats::lm()], [stats::glm()],
#' `lme4::lmer()` or anything called as `model(formula, data = , ...)` -- into
#' a function that fits it to a processed [cyto] object. Use it as the
#' `statistic` of [waterfall_lattice()] to see how every combination of
#' pre-processing choices moves the result, with [select_pipeline()] to choose
#' a pipeline without overfitting, or as a metric in [benchmark_pipelines()].
#'
#' @section What is returned:
#' * With `term`: the model's coefficient for `term`, with its standard error
#'   and p-value when the model's `summary()` provides them. Because cytokineR
#'   values are always log2-scaled, a coefficient of `value` is a log2 fold
#'   change whatever the pipeline.
#' * Without `term` (prediction mode): how well the model predicts the
#'   outcome, scored with `score` -- on the test samples after fitting on the
#'   training samples, or by `folds`-fold cross-validation within the training
#'   samples (`on = "train"`).
#' * With `extract`: whatever `extract(fit)` returns, a number or
#'   `c(value = , p = )`.
#'
#' @section Data:
#' With `data = "long"` (default), the model is fitted once per analyte to the
#' study samples, with the analyte's processed values in column `value` and all
#' sample metadata (subject, time, outcome, ...) as further columns, e.g.
#' `value ~ time + subject` or `outcome ~ value + age`. With `data = "wide"`, it
#' is fitted once to one row per well with one column per analyte (named as
#' the analytes; use backticks, e.g. `` outcome ~ `IL-6` + age ``). A `.` in
#' the formula then stands for all analytes plus the variables named in the
#' formula. With `term`, the result has one row per term; otherwise one row
#' (`analyte = "all"`).
#'
#' @section Unseen test samples:
#' Give a test set to [waterfall_lattice()], [select_pipeline()] or
#' [benchmark_pipelines()] (argument `test`; see [split_samples()]). Its
#' samples get `use_for_fit = FALSE` before every pipeline runs, so
#' normalisation, batch correction and imputation are fitted without them,
#' then applied to them as to new data. The test set is passed on to the
#' statistic, and `on` decides which samples the model uses:
#'
#' | `on` | with `term` | prediction mode |
#' |---|---|---|
#' | `"test"` (default) | coefficient estimated on the test samples | fitted on the training samples, scored on the test samples |
#' | `"train"` | estimated on the training samples | cross-validated within the training samples |
#' | `"all"` | estimated on all samples | cross-validated over all samples |
#'
#' Without a test set, every sample counts as training. To call the statistic
#' yourself, use `f(x, test = ..., on = ...)`.
#'
#' Cross-validation within the training samples refits the model, not the
#' pre-processing, in each fold; it is meant for choosing between pipelines.
#' The test samples stay clean: no pre-processing step, model fit or choice
#' sees them. Methods that fit on every sample regardless of `use_for_fit`
#' cannot keep test samples unseen; of the built-in methods this applies to
#' `combat`.
#'
#' @param formula Model formula.
#' @param model Model-fitting function called as `model(formula, data = ,
#'   ...)`, default [stats::lm()].
#' @param ... Further arguments to `model` (e.g. `family = binomial`).
#' @param term Coefficient to report (its full name, e.g. `"groupcase"`, or a
#'   unique prefix such as `"group"`). In wide mode it may be several terms.
#' @param data `"long"` (one model per analyte) or `"wide"` (one model with the
#'   analytes as columns).
#' @param on Which samples the model uses; see Unseen test samples.
#' @param score Prediction mode: `"auto"` (AUC for a two-level outcome,
#'   otherwise out-of-sample R^2), `"r2"`, `"rmse"`, `"cor"`, `"auc"`, or a
#'   function `function(observed, predicted)` returning a number.
#' @param folds Number of cross-validation folds for `on = "train"` or
#'   `"all"` in prediction mode.
#' @param cv_by Column whose values are kept together in cross-validation
#'   (e.g. `"subject"`); default `"sample_id"`.
#' @param seed Seed for the cross-validation folds.
#' @param predict Prediction mode: function `function(fit, newdata)` returning
#'   predictions on the outcome scale; default `predict(fit, newdata, type =
#'   "response")`.
#' @param extract Optional function of the fitted model returning a number or
#'   `c(value = , p = )`; overrides `term` and `score`.
#' @param types Sample types forming the study population.
#' @param label Name of the metric (default `coef_<term>`, `pred_<score>` or
#'   `model`).
#' @return A function `function(x, holdout = NULL, test = NULL, on = NULL)`
#'   returning a data.table with columns `metric`, `analyte`, `value`, `n`,
#'   `p` and `se`.
#' @seealso [split_samples()], [select_pipeline()], [waterfall_lattice()],
#'   [metric_effect()]
#' @examples
#' x <- cyto_transform(cyto_example, "log2")
#' fc <- model_statistic(value ~ time + subject, term = "time")
#' fc(x)
#'
#' # the time effect in unseen subjects, with pre-processing fitted on the rest
#' test <- split_samples(cyto_example, prop = 0.5, by = "subject")
#' steps <- list(log2 = pipe_transform("log2"),
#'               plate = pipe_normalize("median", within = "time"),
#'               batch = pipe_correct_batch("bridge"))
#' lat <- waterfall_lattice(cyto_filter(cyto_example, analyte == "IL-8"),
#'                          steps, fixed = "log2", statistic = fc, test = test)
#' lat
#'
#' # prediction: time point from the cytokine profile
#' auc <- model_statistic(factor(time) ~ ., model = glm, family = binomial,
#'                        data = "wide", cv_by = "subject")
#' auc(x, test = test)               # AUC in the unseen subjects
#' auc(x, test = test, on = "train") # cross-validated in the others
#' @export
model_statistic <- function(formula, model = stats::lm, ..., term = NULL,
                            data = c("long", "wide"),
                            on = c("test", "train", "all"),
                            score = "auto", folds = 5, cv_by = "sample_id",
                            seed = 1, predict = NULL, extract = NULL,
                            types = "unknown", label = NULL) {
  if (!inherits(formula, "formula") || length(formula) != 3L) {
    stop("'formula' must be a two-sided formula.", call. = FALSE)
  }
  if (!is.function(model)) stop("'model' must be a function.", call. = FALSE)
  data <- match.arg(data)
  on_default <- match.arg(on)
  margs <- list(...)
  mode <- if (!is.null(extract)) "extract" else
    if (!is.null(term)) "term" else "predict"
  if (mode == "term" && data == "long" && length(term) != 1L) {
    stop("With data = \"long\", give a single 'term'.", call. = FALSE)
  }
  if (!is.function(score)) {
    score <- match.arg(score, c("auto", "r2", "rmse", "cor", "auc"))
  }
  # in prediction mode with score = "auto" the label is set once the outcome
  # is known
  user_label <- label
  label <- label %||% switch(mode,
                             term = if (length(term) == 1L)
                               paste0("coef_", term) else "coef",
                             predict = paste0("pred_",
                                              if (is.function(score)) "custom"
                                              else score),
                             extract = "model")
  force(predict)

  fit_model <- function(df) {
    do.call(model, c(list(formula, data = df), margs))
  }
  observed <- function(df) eval(formula[[2]], df, environment(formula))
  empty <- function() {
    data.table::data.table(term = NA_character_, value = NA_real_,
                           p = NA_real_, se = NA_real_, n = 0L)
  }
  one_value <- function(v, n, p = NA_real_, name = NA_character_) {
    data.table::data.table(term = NA_character_, value = v, p = p,
                           se = NA_real_, n = as.integer(n), score = name)
  }
  cv_score <- function(df, units) {
    u <- unique(units)
    k <- min(folds, length(u))
    if (k < 2) stop("Too few units for cross-validation.", call. = FALSE)
    fold <- .with_seed(seed, sample(rep_len(seq_len(k), length(u))))
    fold <- fold[match(units, u)]
    pred <- rep(NA_real_, nrow(df))
    for (i in seq_len(k)) {
      fit <- fit_model(df[fold != i, , drop = FALSE])
      pred[fold == i] <- .model_predict(fit, df[fold == i, , drop = FALSE],
                                        predict)
    }
    .model_score(observed(df), pred, score)
  }

  # one model -> data.table(term, value, p, se, n)
  run_one <- function(train, test_df, units_train, units_test, on) {
    if (mode == "predict") {
      if (on == "test") {
        if (is.null(test_df)) {
          stop("Scoring on test samples needs a test set; pass 'test', or ",
               "use on = \"train\" for cross-validation.", call. = FALSE)
        }
        pr <- .model_predict(fit_model(train), test_df, predict)
        sc <- .model_score(observed(test_df), pr, score)
      } else {
        df <- if (on == "all" && !is.null(test_df)) rbind(train, test_df) else
          train
        un <- if (on == "all" && !is.null(test_df))
          c(units_train, units_test) else units_train
        sc <- cv_score(df, un)
      }
      return(one_value(sc$value, sc$n, name = sc$name))
    }
    df <- switch(on, test = test_df %||% train, train = train,
                 all = if (is.null(test_df)) train else rbind(train, test_df))
    fit <- fit_model(df)
    n <- tryCatch(stats::nobs(fit), error = function(e) nrow(df))
    if (mode == "extract") {
      r <- extract(fit)
      if (!is.numeric(r) || !length(r)) {
        stop("'extract' must return a number or c(value = , p = ).",
             call. = FALSE)
      }
      v <- if ("value" %in% names(r)) r[["value"]] else r[[1]]
      p <- if ("p" %in% names(r)) r[["p"]] else NA_real_
      return(one_value(v, n, p))
    }
    ct <- .coef_table(fit)
    rows <- vapply(term, .match_term, "", available = rownames(ct))
    data.table::data.table(term = term, value = unname(ct[rows, "estimate"]),
                           p = unname(ct[rows, "p"]),
                           se = unname(ct[rows, "se"]), n = as.integer(n))
  }

  function(x, holdout = NULL, test = NULL, on = NULL) {
    .check_cyto(x)
    on <- if (is.null(on)) on_default else
      match.arg(on, c("test", "train", "all"))
    d <- x$data[sample_type %in% types]
    test_ids <- if (is.null(test)) character() else .test_ids(x$data, test)
    if (!cv_by %in% names(d)) {
      stop("Unknown 'cv_by' column '", cv_by, "'.", call. = FALSE)
    }
    if (data == "long") {
      res <- list()
      errs <- character()
      for (an in sort(unique(d$analyte))) {
        da <- as.data.frame(d[analyte == an])
        it <- da$sample_id %in% test_ids
        r <- tryCatch(run_one(da[!it, , drop = FALSE],
                              if (any(it)) da[it, , drop = FALSE] else NULL,
                              da[[cv_by]][!it], da[[cv_by]][it], on),
                      error = function(e) e)
        if (inherits(r, "error")) {
          errs <- c(errs, conditionMessage(r))
          r <- empty()
        }
        r[, analyte := an]
        res[[an]] <- r
      }
      if (length(errs) == length(res)) stop(errs[1], call. = FALSE)
      out <- data.table::rbindlist(res, fill = TRUE)
    } else {
      w <- cyto_wide(x, types = types, as_data_table = TRUE)
      an <- sort(unique(d$analyte))
      it <- w$sample_id %in% test_ids
      units <- w[[cv_by]]
      vars <- intersect(setdiff(all.vars(formula), "."), names(w))
      w <- as.data.frame(w[, unique(c(vars, an)), with = FALSE],
                         check.names = FALSE)
      out <- run_one(w[!it, , drop = FALSE],
                     if (any(it)) w[it, , drop = FALSE] else NULL,
                     units[!it], units[it], on)
      out[, analyte := if (mode == "term") term else "all"]
    }
    lab <- label
    if (is.null(user_label) && "score" %in% names(out) &&
        any(!is.na(out$score))) {
      lab <- paste0("pred_", out$score[!is.na(out$score)][1])
    }
    out[, .(metric = lab, analyte, value, n, p, se)]
  }
}

.with_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", envir = globalenv()))
    get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv())
          else if (exists(".Random.seed", envir = globalenv()))
            rm(".Random.seed", envir = globalenv()), add = TRUE)
  set.seed(seed)
  expr
}

# sample_ids of the test set (a sample is in the test set if any of its rows
# is selected)
.test_ids <- function(data, test) {
  if (is.character(test) && !any(c("bridge", .sample_types) %in% test)) {
    return(unique(test))
  }
  unique(data$sample_id[select_rows(data, test)])
}

.coef_table <- function(fit) {
  ct <- tryCatch(stats::coef(summary(fit)), error = function(e) NULL)
  if (is.null(ct) || !is.matrix(ct) || is.null(rownames(ct))) {
    b <- stats::coef(fit)
    if (!is.numeric(b) || is.null(names(b))) {
      stop("Cannot read coefficients from the model; use 'extract'.",
           call. = FALSE)
    }
    ct <- matrix(b, ncol = 1, dimnames = list(names(b), "Estimate"))
  }
  pc <- grep("^Pr\\(", colnames(ct))
  sc <- match("Std. Error", colnames(ct))
  out <- cbind(estimate = ct[, 1],
               se = if (is.na(sc)) NA_real_ else ct[, sc],
               p = if (length(pc)) ct[, pc[1]] else NA_real_)
  rownames(out) <- rownames(ct)
  out
}

.match_term <- function(term, available) {
  clean <- gsub("`", "", available)
  t <- gsub("`", "", term)
  hit <- which(clean == t)
  if (!length(hit)) hit <- which(startsWith(clean, t))
  if (length(hit) != 1L) {
    stop("Term '", term, "' ", if (length(hit)) "is ambiguous" else
      "not found", "; coefficients: ", paste(clean, collapse = ", "),
      call. = FALSE)
  }
  available[hit]
}

.model_predict <- function(fit, newdata, predict) {
  if (!is.null(predict)) return(as.numeric(predict(fit, newdata)))
  pr <- tryCatch(stats::predict(fit, newdata = newdata, type = "response"),
                 error = function(e) NULL)
  if (is.null(pr)) pr <- stats::predict(fit, newdata = newdata)
  as.numeric(pr)
}

.model_score <- function(obs, pred, score) {
  two <- is.logical(obs) || (is.factor(obs) && nlevels(obs) == 2L) ||
    (is.numeric(obs) && all(obs[!is.na(obs)] %in% c(0, 1)))
  if (is.factor(obs)) {
    if (nlevels(obs) == 2L) obs <- as.numeric(obs == levels(obs)[2]) else
      obs <- as.numeric(obs)
  }
  obs <- as.numeric(obs)
  ok <- is.finite(obs) & is.finite(pred)
  obs <- obs[ok]
  pred <- pred[ok]
  n <- length(obs)
  if (is.function(score)) {
    return(list(value = as.numeric(score(obs, pred)), n = n, name = "custom"))
  }
  if (score == "auto") score <- if (two) "auc" else "r2"
  if (n < 3) return(list(value = NA_real_, n = n, name = score))
  v <- switch(score,
              r2 = 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2),
              rmse = sqrt(mean((obs - pred)^2)),
              cor = stats::cor(obs, pred),
              auc = .auc(obs, pred))
  list(value = v, n = n, name = score)
}

# Mann-Whitney AUC; obs coded 0/1
.auc <- function(obs, pred) {
  n1 <- sum(obs == 1)
  n0 <- sum(obs == 0)
  if (!n1 || !n0) return(NA_real_)
  r <- rank(pred)
  (sum(r[obs == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

#' Split study samples into training and test sets
#'
#' Picks a random share of the study's sampling units (samples, or subjects
#' with all their samples) as an unseen test set, for the `test` argument of
#' [waterfall_lattice()] and [benchmark_pipelines()]. Bridging samples are
#' kept out of the test set by default, so that batch correction keeps its
#' anchors.
#'
#' @param x A [cyto] object.
#' @param prop Share of units in the test set.
#' @param by Column defining the unit that is kept together (e.g.
#'   `"subject"`, so that all time points of a person end up on the same
#'   side).
#' @param strata Optional column to stratify the split by (e.g. an outcome),
#'   taken from each unit's first row.
#' @param types Sample types that can be in the test set.
#' @param exclude Selector of rows whose samples must stay in training
#'   (default: bridging samples).
#' @param seed Random seed for the split.
#' @return A character vector of test `sample_id`s, usable as a selector.
#' @examples
#' test <- split_samples(cyto_example, prop = 0.3, by = "subject")
#' length(test)
#' @export
split_samples <- function(x, prop = 0.3, by = "sample_id", strata = NULL,
                          types = "unknown", exclude = "bridge", seed = 1) {
  .check_cyto(x)
  if (!is.numeric(prop) || length(prop) != 1L || prop <= 0 || prop >= 1) {
    stop("'prop' must be between 0 and 1.", call. = FALSE)
  }
  d <- x$data
  .check_by(d, c(by, strata))
  ex <- if (is.null(exclude)) rep(FALSE, nrow(d)) else select_rows(d, exclude)
  bad_units <- unique(d[[by]][ex])
  cand <- d[sample_type %in% types & !is.na(get(by)) &
              !get(by) %in% bad_units]
  units <- cand[, .(stratum = if (is.null(strata)) "all" else
    as.character(get(strata)[1])), by = by]
  if (nrow(units) < 2) stop("Too few units to split.", call. = FALSE)
  pick <- .with_seed(seed, units[, .(unit = get(by)[sample.int(.N)][
    seq_len(round(prop * .N))]), by = stratum]$unit)
  sort(unique(cand$sample_id[cand[[by]] %in% pick]))
}
