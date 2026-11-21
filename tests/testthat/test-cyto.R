test_that("cyto() validates and fills defaults", {
  d <- data.frame(plate = 1, well = "A1", sample_id = "x", sample_type = "Unknown",
                  analyte = "IL-6", fi = 100)
  x <- cyto(d)
  expect_s3_class(x, "cyto")
  expect_equal(x$data$value, log2(100))
  expect_equal(x$data$batch, "1")
  expect_equal(x$data$dilution, 1)
  expect_equal(x$data$sample_type, "unknown")
  expect_true(x$data$use_for_fit)
  expect_error(cyto(d[, -6]), "Missing required")
  d$sample_type <- "weird"
  expect_error(cyto(d), "Unknown sample_type")
})

test_that("plate ids are made unique across batches", {
  d <- data.frame(plate = 1, well = "A1", sample_id = "x", sample_type = "unknown",
                  analyte = "A", fi = 10, batch = c("1", "2"))
  expect_message(x <- cyto(d), "unique")
  expect_setequal(x$data$plate, c("1:1", "2:1"))
})

test_that("filter, bind, wide and collapse work", {
  x <- small_sim()
  a <- cyto_filter(x, batch == "1")
  b <- cyto_filter(x, batch == "2")
  expect_equal(nrow(cyto_bind(a, b)$data), nrow(x$data))
  expect_equal(nrow(cyto_steps(a)), 1)
  w <- cyto_wide(x)
  expect_true(all(c("IL-6", "IL-8") %in% names(w)))
  expect_equal(nrow(w), x$data[sample_type == "unknown",
                               data.table::uniqueN(paste(plate, well))])
  cl <- cyto_collapse_replicates(x)
  expect_lt(nrow(cl$data), nrow(x$data))
  expect_equal(cl$data[sample_type == "standard" & analyte == "IL-6", .N],
               x$data[sample_type == "standard" & analyte == "IL-6",
                      data.table::uniqueN(paste(plate, level))])
  expect_output(print(x), "measurements")
})

test_that("steps are logged with fits", {
  x <- cyto_transform(small_sim(), "mob")
  s <- cyto_steps(x)
  expect_equal(s$method, "mob")
  expect_equal(s$scale_after, "log2_ratio")
  expect_true(all(c("plate", "analyte", "blank_used") %in% names(cyto_fit(x))))
})
