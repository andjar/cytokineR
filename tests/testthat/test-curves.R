test_that("5PL inverse is the inverse", {
  p <- c(a = 3, d = 14, c = 8, log_b = log(0.8), log_g = log(1.3))
  xx <- seq(0, 14, by = 0.5)
  expect_equal(as.numeric(logistic5_inverse(logistic5(xx, p), p)), xx,
               tolerance = 1e-6)
  inv <- logistic5_inverse(c(2, 15, 8), p)
  expect_equal(attr(inv, "side"), c("left", "right", "none"))
})

test_that("5PL fit recovers a noiseless curve", {
  p <- c(a = 3, d = 14, c = 6, log_b = log(0.7), log_g = 0)
  xx <- rep(seq(-2, 12, by = 2), 2)
  fit <- fit_logistic5(xx, logistic5(xx, p))
  expect_true(fit$converged)
  expect_lt(max(abs(logistic5(xx, fit$par) - logistic5(xx, p))), 0.01)
  f4 <- fit_logistic5(xx, logistic5(xx, p), fix_g = 1)
  expect_equal(unname(f4$par["log_g"]), 0)
})

test_that("shared fit recovers plate shifts", {
  p <- c(a = 3, d = 14, c = 6, log_b = log(0.7), log_g = 0)
  xx <- seq(-2, 12, by = 2)
  sh <- c(P1 = 0, P2 = 0.6, P3 = -0.4)
  x <- rep(xx, 3)
  cur <- rep(names(sh), each = length(xx))
  y <- logistic5(x, p) + sh[cur]
  f <- cytokineR:::.fit_shared5(x, y, cur, shift = "y", tau = c(y = 10))
  est <- f$shifts$shift_y - f$shifts$shift_y[1]
  expect_equal(est, unname(sh - sh[1]), tolerance = 0.02)
})
