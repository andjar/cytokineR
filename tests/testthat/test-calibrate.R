rmse <- function(x) x$data[sample_type == "unknown" & is.finite(value),
                           sqrt(mean((value - truth - stats::median(value - truth))^2))]

test_that("calibration methods produce accurate concentrations", {
  x <- small_sim(seed = 5, bad_standard_prob = 0)
  for (m in c("5pl", "4pl", "spline", "shared")) {
    cx <- cyto_calibrate(x, m)
    expect_equal(cx$scale, "log2_conc")
    expect_lt(rmse(cx), 0.6)
    expect_true("signal" %in% names(cx$data))
    cen <- cx$data[censored == "left"]
    expect_true(all(is.na(cen$value)))
    expect_true(all(is.finite(cen$limit_low)))
  }
})

test_that("spline never extrapolates and 5pl/asymptote keeps extrapolated", {
  x <- small_sim(seed = 6)
  s <- cyto_calibrate(x, "spline")
  f <- cyto_fit(s)$curves
  ok <- s$data[is.finite(value) & sample_type == "unknown"]
  lim <- merge(ok, f[, .(plate, analyte, lloq, uloq)], by = c("plate", "analyte"))
  expect_true(all(lim$value - log2(lim$dilution) >= lim$lloq - 1e-8))
  a <- cyto_calibrate(x, "5pl", censor = "asymptote")
  expect_true(any(grepl("extrapolated", a$data$qc)))
})

test_that("one curve from a reference lot and shared borrowed curves", {
  x <- small_sim(seed = 7)
  one <- cyto_calibrate(x, "5pl", curve_by = NULL, standards = ~ batch == "1")
  expect_equal(nrow(cyto_fit(one)$curves), 4)
  # remove all standards of one plate: shared calibration borrows the curve
  y <- flag_wells(x, plate == "b1p2" & sample_type == "standard")
  sh <- cyto_calibrate(y, "shared")
  expect_true(any(grepl("borrowed_curve", sh$data[plate == "b1p2", qc])))
  expect_true(any(is.finite(sh$data[plate == "b1p2", value])))
  p5 <- cyto_calibrate(y, "5pl")
  expect_true(all(is.na(p5$data[plate == "b1p2", value])))
})

test_that("calibration refuses concentration input", {
  x <- cyto_calibrate(small_sim(), "5pl")
  expect_error(cyto_calibrate(x, "5pl"), "requires")
})
