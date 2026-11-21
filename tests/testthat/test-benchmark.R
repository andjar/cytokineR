test_that("metrics are scale-free and behave as expected", {
  x <- cyto(toy_shifted())
  r0 <- metric_reliability(x, materials = "ctrl")
  expect_lt(r0$value, 1)  # the control moves with the plate shifts
  m <- cyto_normalize(x, "median")
  expect_equal(metric_reliability(m, materials = "ctrl")$value, 1)
  pe0 <- metric_plate_effect(x)$value
  pe1 <- metric_plate_effect(m)$value
  expect_gt(pe0, 0.1)
  expect_lt(pe1, 0.01)
  s <- small_sim()
  expect_true(all(metric_detection(cyto_transform(s, "log2"))$value == 1))
  tr <- metric_truth(cyto_transform(s, "log2"))
  expect_true(all(tr[metric == "truth_rho", value] > 0.5))
  b <- metric_bio(s, ~ time, adjust = ~ subject)
  expect_true(all(c("metric", "analyte", "value") %in% names(b)))
})

test_that("benchmark holds out material from fitting", {
  s <- small_sim()
  p <- list(
    leak = cyto_pipeline(pipe_transform("log2"),
                         pipe_normalize("anchor", anchors = "ctrl_high",
                                        scope = "batch")),
    other = cyto_pipeline(pipe_transform("log2"),
                          pipe_normalize("anchor", anchors = "ctrl_low",
                                         scope = "batch"))
  )
  bm <- benchmark_pipelines(s, p, holdout = "ctrl_high", verbose = FALSE)
  expect_s3_class(bm, "cyto_benchmark")
  # the pipeline anchored on the held-out material finds no usable anchors
  expect_true(any(grepl("no usable anchor|No usable anchor",
                        bm$errors$message)) ||
                any(grepl("anchor", bm$warnings$message)))
  expect_true("other" %in% bm$summary$pipeline)
  expect_s3_class(plot(bm), "ggplot")
  expect_output(print(bm), "cyto_benchmark")
  cf <- benchmark_pipelines(s, p["other"],
                            holdout = list(hi = "ctrl_high", br = "bridge"),
                            verbose = FALSE)
  expect_setequal(unique(cf$results$fold), c("hi", "br"))
})

test_that("simulation benchmark runs across scenarios", {
  p <- list(mob = cyto_pipeline(pipe_transform("mob")))
  bs <- benchmark_simulations(p, list(a = list(plates_per_batch = 2),
                                      b = list(plates_per_batch = 2,
                                               design = "by_time")),
                              reps = 1, verbose = FALSE)
  expect_setequal(unique(bs$summary$scenario), c("a", "b"))
  expect_s3_class(plot(bs), "ggplot")
})
