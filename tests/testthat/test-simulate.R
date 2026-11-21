test_that("simulation is reproducible and complete", {
  a <- small_sim(seed = 3)
  b <- small_sim(seed = 3)
  expect_identical(a$data$fi, b$data$fi)
  d <- a$data
  expect_setequal(unique(d$sample_type), c("standard", "blank", "control",
                                           "unknown"))
  expect_true(all(d[, .N, by = plate]$N == 96 * 4))
  expect_true(all(is.finite(d[sample_type == "unknown", truth])))
  expect_true(any(cytokineR:::.is_bridge(d)))
  expect_true(all(d$fi >= 1))
})

test_that("confounded design puts time points on separate plates", {
  x <- small_sim(design = "by_time")
  tab <- x$data[sample_type == "unknown" & batch == "1",
                .(n_time = data.table::uniqueN(time)), by = plate]
  expect_true(any(tab$n_time == 1))
})
