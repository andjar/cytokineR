# Custom steps that change `value` by reference, to check that branches of
# the lattice do not leak into each other.
add_one <- pipe_custom(function(x) {
  x$data[, value := value + 1]
  x
}, label = "add1")
double_it <- pipe_custom(function(x) {
  x$data[, value := value * 2]
  x
}, label = "double")
mean_value <- function(x) mean(x$data$value)

test_that("lattice values, waterfall and Shapley shares are exact", {
  x <- cyto(toy_shifted())
  m <- mean(x$data$value)
  lat <- waterfall_lattice(x, list(a = add_one, b = double_it),
                           statistic = mean_value)
  expect_s3_class(lat, "cyto_lattice")
  v <- lat$nodes[order(mask), value]
  expect_equal(v, c(m, m + 1, 2 * m, 2 * m + 2))
  ct <- lat$contributions
  # route a > b: a adds 1, then b doubles m + 1
  expect_equal(ct[step == "a", waterfall], 1)
  expect_equal(ct[step == "b", waterfall], m + 1)
  # Shapley: average over both orders
  expect_equal(ct[step == "a", shapley], 1.5)
  expect_equal(ct[step == "b", shapley], m + 0.5)
  expect_equal(sum(ct$shapley), lat$totals$total)
  expect_equal(ct[step == "a", c(min, max)], c(1, 2))
  expect_false(any(ct$flips))
  # another route changes the waterfall but not the Shapley shares
  lat2 <- waterfall_lattice(x, list(a = add_one, b = double_it),
                            statistic = mean_value, route = c("b", "a"))
  ct2 <- lat2$contributions
  expect_equal(ct2$step, c("b", "a"))
  expect_equal(ct2[step == "a", waterfall], 2)
  expect_equal(ct2[step == "b", waterfall], m)
  expect_equal(ct2$shapley, ct$shapley[c(2, 1)])
  # the original object is untouched
  expect_equal(mean(x$data$value), m)
})

test_that("fixed steps are always applied and order is respected", {
  x <- cyto(toy_shifted())
  m <- mean(x$data$value)
  lat <- waterfall_lattice(x, list(base = add_one, b = double_it),
                           fixed = "base", statistic = mean_value)
  expect_equal(nrow(lat$nodes), 2)
  expect_equal(lat$nodes[order(mask), value], c(m + 1, 2 * m + 2))
  expect_equal(lat$optional, "b")
})

test_that("lattice works with metric tables, p-values and plots", {
  s <- small_sim()
  steps <- list(log2 = pipe_transform("log2"),
                plate = pipe_normalize("median", within = "time"),
                batch = pipe_correct_batch("bridge"))
  lat <- waterfall_lattice(s, steps, fixed = "log2",
                           statistic = function(x)
                             metric_effect(x, ~ time, adjust = ~ subject))
  na <- data.table::uniqueN(s$data$analyte)
  expect_equal(nrow(lat$nodes), 4 * na)
  expect_equal(nrow(lat$edges), 4 * na)
  expect_equal(lat$label, "effect")
  expect_true(lat$has_p)
  expect_true(is.logical(lat$nodes$significant))
  tot <- lat$contributions[, .(s = sum(shapley)), by = analyte]
  tot <- merge(tot, lat$totals, by = "analyte")
  expect_equal(tot$s, tot$total, tolerance = 1e-10)
  expect_equal(nrow(summary(lat)), 2 * na)
  expect_output(print(lat), "cyto_lattice")
  expect_s3_class(plot(lat), "ggplot")
  expect_s3_class(plot(lat, "contributions"), "ggplot")
  one <- waterfall_lattice(cyto_filter(s, analyte == "IL-6"), steps,
                           fixed = "log2", statistic = function(x)
                             metric_effect(x, ~ time, adjust = ~ subject))
  expect_output(print(one), "Total change")
  skip_if_not_installed("patchwork")
  expect_s3_class(plot(one, "both"), "patchwork")
})

test_that("holdout rows are hidden from fitting and passed to the statistic", {
  s <- small_sim()
  seen <- NULL
  stat <- function(x, holdout) {
    seen <<- holdout
    h <- select_rows(x$data, holdout)
    expect_false(any(x$data$use_for_fit[h]))
    metric_reliability(x, materials = holdout)
  }
  lat <- waterfall_lattice(
    s, list(log2 = pipe_transform("log2"),
            ctrl = pipe_normalize("anchor", anchors = "ctrl_low",
                                  scope = "batch")),
    fixed = "log2", statistic = stat, holdout = "ctrl_high")
  expect_identical(seen, "ctrl_high")
  expect_equal(lat$label, "reliability_plate")
  expect_false(lat$has_p)
  # default_metrics() elements work directly
  lat2 <- waterfall_lattice(
    s, list(log2 = pipe_transform("log2"),
            ctrl = pipe_normalize("anchor", anchors = "ctrl_low",
                                  scope = "batch")),
    fixed = "log2", statistic = default_metrics()$reliability_plate,
    holdout = "ctrl_high")
  expect_equal(lat$nodes$value, lat2$nodes$value)
})

test_that("failing steps give missing nodes and are reported", {
  x <- cyto(toy_shifted())
  boom <- pipe_custom(function(x) stop("boom"), label = "boom")
  lat <- waterfall_lattice(x, list(a = add_one, boom = boom),
                           statistic = mean_value)
  expect_true(any(grepl("boom", lat$errors$message)))
  expect_equal(sum(is.na(lat$nodes$value)), 2)
  expect_true(is.na(lat$contributions[step == "boom", shapley]))
  expect_output(print(lat), "error")
  expect_error(waterfall_lattice(x, list(boom = boom, a = add_one),
                                 fixed = "boom", statistic = mean_value),
               "No pipeline produced")
})

test_that("statistic formats are normalised", {
  x <- cyto(toy_shifted())
  named <- function(x) c(A = mean(x$data$value), B = 1)
  lat <- waterfall_lattice(x, list(a = add_one), statistic = named)
  expect_setequal(unique(lat$nodes$analyte), c("A", "B"))
  expect_equal(lat$contributions[analyte == "B", shapley], 0)
  expect_error(waterfall_lattice(x, list(a = add_one),
                                 statistic = function(x) c(1, 2)),
               "No pipeline produced")
  two <- function(x) data.frame(metric = c("m1", "m2"), analyte = "A",
                                value = 1:2)
  expect_error(waterfall_lattice(x, list(a = add_one), statistic = two),
               "several metrics")
})

test_that("arguments are checked", {
  x <- cyto(toy_shifted())
  st <- list(a = add_one, b = double_it)
  expect_error(waterfall_lattice(x, st, mean_value, fixed = "zz"),
               "Unknown fixed")
  expect_error(waterfall_lattice(x, st, mean_value, fixed = c("a", "b")),
               "optional")
  expect_error(waterfall_lattice(x, st, mean_value, route = "a"), "route")
  expect_error(waterfall_lattice(x, st, mean_value, max_steps = 1),
               "max_steps")
  expect_error(waterfall_lattice(x, add_one, mean_value), "named list")
  expect_error(waterfall_lattice(x, list(1), mean_value), "pipe_")
  expect_error(waterfall_lattice(x, st, "mean"), "function")
  expect_error(waterfall_lattice(x, st, mean_value, holdout = "nothing"),
               "selects no rows")
})

test_that("metric_effect returns Cohen's d and the t-test p-value", {
  set.seed(3)
  g1 <- stats::rnorm(12, 10, 1)
  g2 <- stats::rnorm(12, 11, 1)
  d <- data.frame(plate = "P1", well = paste0("W", 1:24),
                  sample_id = paste0("s", 1:24), sample_type = "unknown",
                  analyte = "A", fi = 2^c(g1, g2),
                  group = rep(c("a", "b"), each = 12))
  x <- cyto(d)
  e <- metric_effect(x, ~ group)
  sp <- sqrt((11 * stats::var(g1) + 11 * stats::var(g2)) / 22)
  expect_equal(e$value, (mean(g2) - mean(g1)) / sp)
  expect_equal(e$estimate, mean(g2) - mean(g1))
  expect_equal(e$p, stats::t.test(g2, g1, var.equal = TRUE)$p.value)
  expect_equal(e$n, 24L)
  expect_equal(e$metric, "effect")
  d$group3 <- rep(c("a", "b", "c"), 8)
  expect_error(metric_effect(cyto(d), ~ group3), "two-level")
  expect_error(metric_effect(x, ~ group + plate), "exactly one term")
  expect_error(metric_effect(x, ~ nope), "Unknown grouping")
})
