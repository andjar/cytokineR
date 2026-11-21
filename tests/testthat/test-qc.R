test_that("flagging excludes rows from fitting", {
  x <- small_sim()
  x <- flag_beads(x, min_count = 70)
  expect_true(all(!x$data[grepl("low_beads", qc), use_for_fit]))
  y <- flag_wells(x, well == "A1", reason = "test", exclude = FALSE)
  expect_true(any(grepl("test", y$data$qc)))
  r <- flag_replicates(small_sim(), max_cv = 0.01, types = "standard")
  expect_true(any(grepl("high_cv", r$data$qc)))
})

test_that("PCA, plate summary, variance decomposition and Levey-Jennings", {
  x <- cyto_transform(small_sim(), "log2")
  p <- qc_pca(x)
  expect_s3_class(p, "cyto_pca")
  expect_equal(nrow(p$by_plate), 4)
  expect_s3_class(plot(p), "ggplot")
  qp <- qc_plates(x)
  expect_true(all(c("z_samples", "z_blank", "pca_dist", "flagged") %in% names(qp)))
  v <- qc_variance(x, covariates = "time")
  tot <- v[, sum(r2), by = analyte]$V1
  expect_equal(tot, rep(1, length(tot)), tolerance = 1e-8)
  # inject a large shift into one control on the last plate
  d <- data.table::copy(x$data)
  last <- tail(unique(d$plate), 1)
  d[sample_id == "ctrl_high" & plate == last, value := value + 5]
  x2 <- x
  x2$data <- d
  lj <- qc_controls(x2)
  expect_true(any(grepl("1_3s", lj[plate == last & material == "ctrl_high", rule])))
  expect_s3_class(plot(lj), "ggplot")
})

test_that("plots return ggplot objects", {
  x <- small_sim()
  expect_s3_class(plot_plates(x), "ggplot")
  expect_s3_class(plot_standards(cyto_calibrate(x, "shared")), "ggplot")
  expect_s3_class(plot_corrections(cyto_normalize(x, "median")), "ggplot")
})
