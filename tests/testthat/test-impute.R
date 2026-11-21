test_that("imputation fills censored values below their limits", {
  x <- cyto_calibrate(small_sim(seed = 8, abundance = c(-3, 3)), "5pl")
  n_left <- x$data[censored == "left" & sample_type == "unknown", .N]
  expect_gt(n_left, 0)
  for (m in c("fraction", "tobit", "rank_tobit")) {
    im <- cyto_impute(x, m)
    l <- im$data[censored == "left" & sample_type == "unknown"]
    expect_true(all(is.finite(l$value)))
    expect_true(all(l$value <= l$limit_low + 1e-8))
    expect_true(all(grepl("imputed", l$qc)))
  }
  fr <- cyto_impute(x, "fraction", frac = 0.5)
  l <- fr$data[censored == "left" & sample_type == "unknown"]
  expect_equal(l$value, l$limit_low - 1)
})

test_that("rank_tobit preserves signal order within plate", {
  x <- cyto_calibrate(small_sim(seed = 9, abundance = c(-3, 3)), "5pl")
  im <- cyto_impute(x, "rank_tobit")
  l <- im$data[censored == "left" & sample_type == "unknown"]
  chk <- l[, .(ok = if (.N > 2) stats::cor(signal, value, method = "spearman") > 0.99
               else TRUE), by = .(plate, analyte)]
  expect_true(all(chk$ok))
  expect_error(cyto_impute(small_sim(), "rank_tobit"), "signal")
})

test_that("tobit with covariates and draws", {
  x <- cyto_calibrate(small_sim(seed = 10, abundance = c(-3, 3)), "5pl")
  a <- cyto_impute(x, "tobit", formula = ~ factor(time), target = "unknown")
  expect_true(all(is.finite(a$data[censored != "none" & sample_type == "unknown",
                                   value])))
  d1 <- cyto_impute(x, "tobit", draw = TRUE, seed = 1)
  d2 <- cyto_impute(x, "tobit", draw = TRUE, seed = 1)
  expect_identical(d1$data$value, d2$data$value)
})

test_that("tobit falls back safely when almost everything is censored", {
  x <- cyto_calibrate(small_sim(seed = 12, abundance = c(-9, -7)), "5pl")
  im <- cyto_impute(x, "rank_tobit")
  f <- cyto_fit(im)
  expect_true(any(f$status == "fallback_too_few_observed"))
  l <- im$data[censored == "left" & sample_type == "unknown"]
  expect_true(all(l$value >= l$limit_low - 2 - 1e-8 & l$value <= l$limit_low + 1e-8))
})
