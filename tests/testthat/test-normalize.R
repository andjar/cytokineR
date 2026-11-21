test_that("median and anchor normalisation recover known plate shifts", {
  x <- cyto(toy_shifted())
  m <- cyto_normalize(x, "median")
  f <- cyto_fit(m)
  off <- stats::setNames(f$offset, f$plate)
  expect_equal(unname(off["P2"] - off["P1"]), 1, tolerance = 1e-8)
  expect_equal(unname(off["P3"] - off["P1"]), -0.5, tolerance = 1e-8)
  v <- m$data[sample_type == "unknown"]
  expect_lt(max(v[, stats::sd(value), by = sample_id]$V1), 1e-8)
  a <- cyto_normalize(x, "anchor", anchors = "ctrl")
  v <- a$data[sample_type == "unknown"]
  expect_lt(max(v[, stats::sd(value), by = sample_id]$V1), 1e-8)
  q <- cyto_normalize(x, "quantile")
  expect_s3_class(q, "cyto")
})

test_that("regress protects covariates and removes plate effects", {
  x <- cyto(toy_shifted())
  r <- cyto_normalize(x, "regress")
  v <- r$data[sample_type == "unknown"]
  expect_lt(max(v[, stats::sd(value), by = sample_id]$V1), 1e-8)
})

test_that("anchor with no material on a plate warns and leaves it", {
  d <- toy_shifted()
  d <- d[!(d$plate == "P3" & d$sample_type == "control"), ]
  x <- cyto(d)
  expect_warning(a <- cyto_normalize(x, "anchor", anchors = "ctrl"),
                 "no anchor")
  expect_equal(cyto_fit(a)[plate == "P3", status], "no_anchor")
})

test_that("bridge correction aligns batches", {
  x <- small_sim()
  l <- cyto_transform(x, "log2")
  b <- cyto_correct_batch(l)
  f <- cyto_fit(b)
  expect_true(all(f[status == "reference", offset] == 0))
  br <- cytokineR:::.is_bridge(b$data)
  dev <- b$data[br, .(y = mean(value)), by = .(batch, sample_id, analyte)][
    , .(d = diff(range(y))), by = .(sample_id, analyte)]
  dev0 <- l$data[br, .(y = mean(value)), by = .(batch, sample_id, analyte)][
    , .(d = diff(range(y))), by = .(sample_id, analyte)]
  expect_lt(stats::median(dev$d), stats::median(dev0$d))
  bd <- cyto_correct_batch(l, fit = "deming")
  expect_true(all(is.finite(cyto_fit(bd)$slope)))
})

test_that("affine anchor fit with standards", {
  x <- cyto_transform(small_sim(), "log2")
  a <- cyto_normalize(x, "anchor", anchors = "standard", scope = "batch",
                      fit = "affine")
  f <- cyto_fit(a)
  expect_true(all(f$slope >= 0.5 & f$slope <= 2))
  sd_before <- x$data[sample_type == "standard",
                      stats::sd(value), by = .(batch, level, analyte)]$V1
  sd_after <- a$data[sample_type == "standard",
                     stats::sd(value), by = .(batch, level, analyte)]$V1
  expect_lt(stats::median(sd_after), stats::median(sd_before))
})
