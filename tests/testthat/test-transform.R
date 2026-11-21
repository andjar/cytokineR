test_that("signal representations behave as defined", {
  x <- small_sim()
  l <- cyto_transform(x, "log2")
  expect_equal(l$data$value, log2(pmax(x$data$fi, 1)))
  m <- cyto_transform(x, "mob", shrink = FALSE)
  expect_equal(m$scale, "log2_ratio")
  # blanks are ~0 on the multiple-of-blank scale
  expect_lt(abs(stats::median(m$data[sample_type == "blank", value])), 0.2)
  ms <- cyto_transform(x, "mob", shrink = TRUE)
  f <- cyto_fit(ms)
  expect_true(all(f$weight >= 0 & f$weight <= 1))
  s <- cyto_transform(x, "mos", level = 4)
  expect_lt(abs(stats::median(s$data[sample_type == "standard" & level == 4,
                                     value])), 0.2)
  s2 <- cyto_transform(x, "mos", by_conc = TRUE)
  expect_true(all(is.finite(s2$data$value)))
  mc <- cyto_transform(x, "moc", control = "ctrl_low")
  expect_lt(abs(stats::median(mc$data[sample_id == "ctrl_low", value])), 0.2)
  b <- cyto_transform(x, "bgsub")
  expect_true(all(is.finite(b$data$value)))
})

test_that("transform after normalize warns", {
  x <- cyto_normalize(small_sim(), "median")
  expect_warning(cyto_transform(x, "log2"), "overridden")
})
