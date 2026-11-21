test_that("model_statistic() reports lm coefficients per analyte", {
  s <- small_sim()
  x <- cyto_transform(s, "log2")
  fc <- model_statistic(value ~ time + subject, term = "time")
  r <- fc(x)
  expect_setequal(r$analyte, unique(s$data$analyte))
  expect_equal(unique(r$metric), "coef_time")
  d <- as.data.frame(x$data[sample_type == "unknown" & analyte == "IL-8"])
  m <- summary(stats::lm(value ~ time + subject, data = d))$coefficients
  expect_equal(r[analyte == "IL-8", value], m["time", 1])
  expect_equal(r[analyte == "IL-8", se], m["time", 2])
  expect_equal(r[analyte == "IL-8", p], m["time", 4])
  # term prefix matching and errors
  g <- model_statistic(value ~ factor(time), term = "factor(time)")
  expect_equal(nrow(g(x)), 4)
  bad <- model_statistic(value ~ time, term = "nope")
  expect_error(bad(x), "not found")
  expect_error(model_statistic(value ~ time, term = c("a", "b")), "single")
  expect_error(model_statistic(~ time), "two-sided")
})

test_that("on = 'test' / 'train' use the right samples", {
  s <- small_sim()
  x <- cyto_transform(s, "log2")
  test <- split_samples(s, prop = 0.5, by = "subject")
  fc <- model_statistic(value ~ time + subject, term = "time")
  d <- as.data.frame(x$data[sample_type == "unknown" & analyte == "IL-6"])
  it <- d$sample_id %in% test
  b_test <- stats::coef(stats::lm(value ~ time + subject, d[it, ]))[["time"]]
  b_train <- stats::coef(stats::lm(value ~ time + subject, d[!it, ]))[["time"]]
  expect_equal(fc(x, test = test)[analyte == "IL-6", value], b_test)
  expect_equal(fc(x, test = test, on = "train")[analyte == "IL-6", value],
               b_train)
  expect_equal(fc(x, test = test)[analyte == "IL-6", n], sum(it))
  # without a test set, everything is training data
  expect_equal(fc(x)$value, fc(x, on = "train")$value)
})

test_that("prediction mode scores unseen samples and cross-validates", {
  s <- small_sim()
  x <- cyto_transform(s, "log2")
  test <- split_samples(s, prop = 0.4, by = "subject")
  auc <- model_statistic(factor(time) ~ ., model = stats::glm,
                         family = stats::binomial, data = "wide",
                         cv_by = "subject")
  a <- auc(x, test = test)
  expect_equal(a$metric, "pred_auc")
  expect_equal(a$analyte, "all")
  # manual check
  w <- cyto_wide(x)
  it <- w$sample_id %in% test
  an <- sort(unique(s$data$analyte))
  dd <- w[, c("time", an)]
  fit <- stats::glm(factor(time) ~ ., family = stats::binomial,
                    data = dd[!it, ])
  pr <- stats::predict(fit, dd[it, ], type = "response")
  y <- as.numeric(dd$time[it] == max(dd$time))
  n1 <- sum(y)
  n0 <- sum(1 - y)
  manual <- (sum(rank(pr)[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  expect_equal(a$value, manual)
  expect_equal(a$n, sum(it))
  # the outcome variable is not among the predictors of `.`
  ncoef <- model_statistic(factor(time) ~ ., model = stats::glm,
                           family = stats::binomial, data = "wide",
                           extract = function(f) length(stats::coef(f)))
  expect_equal(ncoef(x)$value, length(an) + 1)
  # cross-validation is reproducible and needs no test set
  cv1 <- auc(x, on = "train")
  cv2 <- auc(x, on = "train")
  expect_equal(cv1$value, cv2$value)
  expect_true(cv1$value > 0.5 && cv1$value <= 1)
  expect_error(auc(x), "needs a test set")
  # R^2 and RMSE for a continuous outcome
  r2 <- model_statistic(value ~ time, data = "long", score = "r2")
  rm <- model_statistic(value ~ time, data = "long", score = "rmse")
  expect_equal(unique(r2(x, test = test)$metric), "pred_r2")
  expect_true(all(rm(x, test = test)$value > 0))
  cs <- model_statistic(value ~ time, score = function(o, p) 42)
  expect_equal(unique(cs(x, test = test)$value), 42)
})

test_that("wide mode reports several terms", {
  s <- small_sim()
  x <- cyto_transform(s, "log2")
  m <- model_statistic(time ~ `IL-6` + `IL-8`, data = "wide",
                       term = c("IL-6", "IL-8"))
  r <- m(x)
  expect_equal(r$analyte, c("IL-6", "IL-8"))
  w <- cyto_wide(x)
  cf <- stats::coef(stats::lm(time ~ `IL-6` + `IL-8`, data = w))
  expect_equal(r$value, unname(cf[c("`IL-6`", "`IL-8`")]))
})

test_that("split_samples() keeps units together and bridges in training", {
  s <- small_sim()
  t1 <- split_samples(s, prop = 0.3, by = "subject")
  expect_identical(t1, split_samples(s, prop = 0.3, by = "subject"))
  d <- s$data[sample_type == "unknown"]
  ts <- unique(d[sample_id %in% t1, subject])
  tr <- unique(d[!sample_id %in% t1, subject])
  expect_length(intersect(ts, tr), 0)
  expect_false(any(d$sample_id[select_rows(d, "bridge")] %in% t1))
  # subjects with any bridging sample stay in training as a whole
  bridged <- unique(d$subject[select_rows(d, "bridge")])
  n_subj <- data.table::uniqueN(d[!subject %in% bridged, subject])
  expect_equal(length(ts), round(0.3 * n_subj))
  st <- split_samples(s, prop = 0.5, by = "subject", strata = "batch")
  expect_true(all(c("1", "2") %in% d[sample_id %in% st, batch]))
  expect_error(split_samples(s, prop = 1), "between 0 and 1")
  expect_error(split_samples(s, by = "nope"), "Unknown")
})

test_that("test samples are hidden from pre-processing in the lattice", {
  s <- small_sim()
  test <- split_samples(s, prop = 0.4, by = "subject")
  seen <- NULL
  stat <- function(x, test) {
    seen <<- test
    expect_false(any(x$data$use_for_fit[x$data$sample_id %in% test]))
    mean(x$data[sample_type == "unknown", value])
  }
  steps <- list(log2 = pipe_transform("log2"),
                plate = pipe_normalize("median", within = "time"))
  lat <- waterfall_lattice(s, steps, fixed = "log2", statistic = stat,
                           test = test)
  expect_identical(sort(seen), sort(test))
  expect_identical(lat$test, test)
  expect_output(print(lat), "Test set")
  # same as normalising with the test samples excluded by hand
  h <- data.table::copy(s$data)
  h[sample_id %in% test, use_for_fit := FALSE]
  sh <- s
  sh$data <- h
  ref <- run_pipeline(sh, cyto_pipeline(steps$log2, steps$plate))
  expect_equal(lat$nodes[mask == 1L, value],
               mean(ref$data[sample_type == "unknown", value]))
  # and different from normalising with them
  all_in <- run_pipeline(s, cyto_pipeline(steps$log2, steps$plate))
  expect_false(isTRUE(all.equal(
    ref$data$value, all_in$data$value)))
  expect_error(waterfall_lattice(s, steps, fixed = "log2", statistic = stat,
                                 test = "no_such_sample"),
               "selects no samples")
})

test_that("select_pipeline() chooses on training and reports test", {
  s <- small_sim()
  test <- split_samples(s, prop = 0.4, by = "subject")
  steps <- list(log2 = pipe_transform("log2"),
                plate = pipe_normalize("median", within = "time"),
                batch = pipe_correct_batch("bridge"))
  auc <- model_statistic(factor(time) ~ ., model = stats::glm,
                         family = stats::binomial, data = "wide",
                         cv_by = "subject")
  sel <- select_pipeline(s, steps, auc, test = test, fixed = "log2")
  expect_s3_class(sel, "cyto_selection")
  sc <- sel$scores
  expect_equal(nrow(sc), 4)
  expect_equal(sc[selected == TRUE, train], max(sc$train))
  expect_equal(sel$selected, sc[selected == TRUE, pipeline])
  expect_true(sel$maximize)
  # the chosen pipeline reproduces its scores
  h <- .hide_rows(s, test = test)
  out <- run_pipeline(h, sel$pipeline)
  expect_equal(auc(out, test = test, on = "train")$value,
               sc[selected == TRUE, train])
  expect_equal(auc(out, test = test)$value, sc[selected == TRUE, test])
  expect_s3_class(sel$train, "cyto_lattice")
  expect_s3_class(sel$test, "cyto_lattice")
  expect_output(print(sel), "Chosen on the training samples")
  expect_s3_class(plot(sel), "ggplot")
  # lower-is-better metrics are minimised
  rm <- model_statistic(value ~ time, score = "rmse", cv_by = "subject")
  sel2 <- select_pipeline(s, steps, rm, test = test, fixed = "log2",
                          analytes = "IL-8")
  expect_false(sel2$maximize)
  expect_equal(sel2$scores[selected == TRUE, train], min(sel2$scores$train))
  expect_error(select_pipeline(s, steps, auc, fixed = "log2"), "test set")
  expect_error(select_pipeline(s, steps, mean, test = test, fixed = "log2"),
               "'test' and 'on'")
})

test_that("benchmark_pipelines() hides the test set and passes it on", {
  s <- small_sim()
  test <- split_samples(s, prop = 0.4, by = "subject")
  got <- NULL
  m <- list(probe = function(x, holdout, test) {
    got <<- test
    data.table::data.table(
      metric = "probe", analyte = "all",
      value = sum(x$data$use_for_fit[x$data$sample_id %in% test]), n = 1L)
  })
  bm <- benchmark_pipelines(s, list(a = cyto_pipeline(pipe_transform("log2"))),
                            metrics = m, test = test, verbose = FALSE)
  expect_identical(sort(got), sort(test))
  expect_equal(bm$results$value, 0)
  expect_identical(bm$test, test)
})
