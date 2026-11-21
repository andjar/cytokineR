test_that("pipelines run and are logged", {
  p <- cyto_pipeline(pipe_transform("log2"),
                     pipe_normalize("anchor", anchors = "control", scope = "batch"),
                     pipe_correct_batch("bridge"),
                     name = "anchor")
  expect_output(print(p), "anchor")
  x <- run_pipeline(small_sim(), p)
  expect_equal(nrow(cyto_steps(x)), 3)
  q <- cyto_pipeline(p, pipe_custom(function(x) cyto_filter(x, analyte != "IL-6")))
  expect_length(q, 4)
  y <- run_pipeline(small_sim(), q)
  expect_false("IL-6" %in% y$data$analyte)
})

test_that("expand_pipelines builds the full grid", {
  g <- expand_pipelines(a = list(log2 = pipe_transform("log2"),
                                 mob = pipe_transform("mob")),
                        b = list(none = NULL, med = pipe_normalize("median")))
  expect_length(g, 4)
  expect_setequal(names(g), c("log2", "mob", "log2 > med", "mob > med"))
})

test_that("custom methods can be registered and used", {
  register_cyto_method("normalize", "nothing", function(x, ...) {
    list(value = x$data$value, note = "no-op")
  }, description = "test")
  on.exit(unregister_cyto_method("normalize", "nothing"))
  expect_true("nothing" %in% cyto_methods("normalize")$method)
  x <- cyto_normalize(small_sim(), "nothing")
  expect_equal(cyto_steps(x)$note, "no-op")
  expect_error(register_cyto_method("normalize", "nothing", identity),
               "already")
  expect_error(cyto_normalize(small_sim(), "does_not_exist"), "Unknown")
})
