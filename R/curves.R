# Five-parameter logistic in log-log space
#
#   y = a + (d - a) / (1 + exp(-b * (x - c)))^g
#
# x = log2 concentration, y = log2 signal. a = bottom asymptote, d = top
# asymptote, c = inflection location, b > 0 slope, g > 0 asymmetry (g = 1 is
# the 4PL). Parameters are stored as c(a, d, c, log_b, log_g).

#' Five-parameter logistic helpers
#'
#' Low-level functions used by [cyto_calibrate()], exported for users who want
#' to build their own calibration methods. The curve is parameterised in
#' log-log space: `y = a + (d - a) / (1 + exp(-b * (x - c)))^g`, where `x` is
#' log2 concentration and `y` a log2 signal; `par = c(a, d, c, log(b), log(g))`.
#'
#' @param x Log2 concentrations.
#' @param y Log2 signal values.
#' @param par Parameter vector `c(a, d, c, log_b, log_g)`.
#' @return `logistic5()`: predicted signal. `logistic5_inverse()`: log2
#'   concentration (`NA` outside the asymptotes), with attribute `side`
#'   (`"left"`, `"right"` or `"none"`).
#' @examples
#' p <- c(a = 3, d = 14, c = 8, log_b = log(0.8), log_g = 0)
#' y <- logistic5(1:12, p)
#' logistic5_inverse(y, p)
#' fit <- fit_logistic5(rep(1:12, 2), rep(y, 2) + rnorm(24, sd = 0.05))
#' fit$par
#' @export
logistic5 <- function(x, par) {
  a <- par[1]
  d <- par[2]
  cc <- par[3]
  b <- exp(par[4])
  g <- exp(par[5])
  a + (d - a) / (1 + exp(-b * (x - cc)))^g
}

#' @rdname logistic5
#' @export
logistic5_inverse <- function(y, par) {
  a <- par[1]
  d <- par[2]
  cc <- par[3]
  b <- exp(par[4])
  g <- exp(par[5])
  out <- rep(NA_real_, length(y))
  side <- rep("none", length(y))
  ok <- is.finite(y) & y > a & y < d
  side[is.finite(y) & y <= a] <- "left"
  side[is.finite(y) & y >= d] <- "right"
  e <- ((d - a) / (y[ok] - a))^(1 / g) - 1
  out[ok] <- cc - log(e) / b
  bad <- ok
  bad[ok] <- !is.finite(out[ok])
  side[bad & y <= (a + d) / 2] <- "left"
  side[bad & y > (a + d) / 2] <- "right"
  out[bad] <- NA_real_
  attr(out, "side") <- side
  out
}

#' @rdname logistic5
#' @param w Optional weights.
#' @param fix_g Fix the asymmetry `g` (e.g. `1` for a 4PL); `NULL` estimates it.
#' @param lambda_g Ridge penalty on `log(g)`; shrinks towards the symmetric
#'   4PL and stabilises fits with few points. `0` disables.
#' @param start Optional starting values.
#' @param robust Use a Huber loss (scale from the MAD of least-squares
#'   residuals) so that single aberrant standard wells have bounded influence.
#' @param .loss_scale Internal.
#' @return `fit_logistic5()`: a list with `par`, `sse`, `sigma`,
#'   `converged`, `n`.
#' @export
fit_logistic5 <- function(x, y, w = NULL, fix_g = NULL, lambda_g = 0.1,
                          start = NULL, robust = FALSE, .loss_scale = NULL) {
  if (robust) {
    f0 <- fit_logistic5(x, y, w, fix_g, lambda_g, start, robust = FALSE)
    if (anyNA(f0$par)) return(f0)
    ok <- is.finite(x) & is.finite(y)
    sc <- .robust_scale(y[ok] - logistic5(x[ok], f0$par))
    f1 <- fit_logistic5(x, y, w, fix_g, lambda_g, start = f0$par,
                        robust = FALSE, .loss_scale = sc)
    if (anyNA(f1$par)) return(f0)
    return(f1)
  }
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (is.null(w)) w <- rep(1, length(x)) else w <- w[ok]
  n <- length(x)
  fail <- list(par = rep(NA_real_, 5), sse = NA_real_, sigma = NA_real_,
               converged = FALSE, n = n)
  if (n < 4 || data.table::uniqueN(x) < 4) return(fail)
  ry <- diff(range(y))
  if (ry < 1e-6) return(fail)

  # starting values
  if (is.null(start)) {
    a0 <- min(y) - 0.05 * ry
    d0 <- max(y) + 0.05 * ry
    agg <- stats::aggregate(y, list(x = x), mean)
    agg <- agg[order(agg$x), ]
    mid <- (a0 + d0) / 2
    # location where the (monotonised) mean signal crosses the midpoint
    ymono <- cummax(agg[[2]])
    if (any(ymono >= mid) && any(ymono < mid)) {
      i <- which(ymono >= mid)[1]
      c0 <- agg$x[max(1, i - 1)] + (agg$x[i] - agg$x[max(1, i - 1)]) / 2
    } else {
      c0 <- stats::median(agg$x)
    }
    s0 <- stats::coef(stats::lm(agg[[2]] ~ agg$x))[2]
    if (!is.finite(s0) || s0 <= 0) s0 <- 0.5
    b0 <- max(min(4 * s0 / (d0 - a0), 5), 0.02)
    starts <- list(c(a0, d0, c0, log(b0), 0),
                   c(a0, d0, c0, log(b0), log(0.5)),
                   c(a0, d0, c0, log(b0), log(2)))
  } else {
    starts <- list(start)
  }
  if (!is.null(fix_g)) starts <- unique(lapply(starts, function(s) {
    s[5] <- log(fix_g)
    s
  }))

  lower <- c(min(y) - 8, max(y) - 0.25 * ry, min(x) - 12, log(0.01), log(0.1))
  upper <- c(min(y) + 0.25 * ry, max(y) + 4, max(x) + 12, log(10), log(10))

  obj <- function(p) {
    if (!is.null(fix_g)) p[5] <- log(fix_g)
    r <- y - logistic5(x, p)
    s <- if (is.null(.loss_scale)) sum(w * r^2) else
      sum(w * .huber(r, .loss_scale))
    if (is.null(fix_g) && lambda_g > 0) s <- s + lambda_g * p[5]^2
    if (!is.finite(s)) 1e10 else s
  }
  best <- NULL
  for (s in starts) {
    s <- pmin(pmax(s, lower + 1e-8), upper - 1e-8)
    r <- tryCatch(
      stats::nlminb(s, obj, lower = lower, upper = upper,
                    control = list(iter.max = 300, eval.max = 600)),
      error = function(e) NULL)
    if (!is.null(r) && (is.null(best) || r$objective < best$objective)) best <- r
  }
  if (is.null(best)) return(fail)
  p <- best$par
  if (!is.null(fix_g)) p[5] <- log(fix_g)
  names(p) <- c("a", "d", "c", "log_b", "log_g")
  res <- y - logistic5(x, p)
  df <- max(n - (if (is.null(fix_g)) 5 else 4), 1)
  list(par = p, sse = sum(res^2), sigma = sqrt(sum(res^2) / df),
       converged = best$convergence == 0, n = n)
}

# Huber loss (times 2, so it equals r^2 in the quadratic zone)
.huber <- function(r, s, k = 1.345) {
  a <- abs(r) / s
  ifelse(a <= k, r^2, s^2 * (2 * k * a - k^2))
}

.robust_scale <- function(r, floor = 0.03) {
  r <- r[is.finite(r)]
  max(1.4826 * stats::median(abs(r - stats::median(r))), floor)
}

# Shared-shape 5PL across curves with curve-specific shifts.
#
#   y_ij = a + beta_j + delta_j + (d - a - beta_j) /
#          (1 + exp(-b (x_ij - c - eta_j)))^g
#
# i.e. a common shape theta = (a, d, c, b, g) and per-curve shifts:
# eta (horizontal, "x": potency), delta (vertical, "y": gain) and beta
# (bottom asymptote, "bottom": background). Shifts are penalised towards zero
# (Gaussian prior with sd tau; residual sd sigma) and estimated by
# alternating between the shared shape and the per-curve shifts.
.shared_pred <- function(theta, sh, x, ci) {
  e <- sh[ci, 1]
  dd <- sh[ci, 2]
  bt <- sh[ci, 3]
  a <- theta[1] + bt + dd
  d <- theta[2] + dd
  a + (d - a) / (1 + exp(-exp(theta[4]) * (x - theta[3] - e)))^exp(theta[5])
}

.fit_shared5 <- function(x, y, curve, shift = c("x", "y", "bottom"),
                         tau = c(x = 1, y = 1, bottom = 1), fix_g = NULL,
                         lambda_g = 0.1, max_iter = 25, tol = 1e-6,
                         robust = TRUE) {
  if (identical(shift, "both")) shift <- c("x", "y")
  shift <- match.arg(shift, c("x", "y", "bottom"), several.ok = TRUE)
  tau0 <- c(x = 1, y = 1, bottom = 1)
  nt <- intersect(names(tau), names(tau0))
  tau0[nt] <- tau[nt]
  curve <- as.character(curve)
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  curve <- curve[ok]
  ids <- unique(curve)
  ci <- match(curve, ids)
  active <- c("x", "y", "bottom") %in% shift
  sh <- matrix(0, length(ids), 3, dimnames = list(ids, c("x", "y", "bottom")))
  fit <- fit_logistic5(x, y, fix_g = fix_g, lambda_g = lambda_g,
                       robust = robust)
  if (any(is.na(fit$par))) return(NULL)
  theta <- fit$par
  sigma <- if (robust) .robust_scale(y - logistic5(x, theta), 0.05) else
    max(fit$sigma, 0.05)
  loss <- if (robust) function(r) .huber(r, sigma) else function(r) r^2
  ry <- diff(range(y))
  lower <- c(min(y) - 8, max(y) - 0.25 * ry, min(x) - 12, log(0.01), log(0.1))
  upper <- c(min(y) + 0.25 * ry, max(y) + 4, max(x) + 12, log(10), log(10))
  if (!is.null(fix_g)) {
    lower[5] <- upper[5] <- log(fix_g)
    theta[5] <- log(fix_g)
  }
  pen <- 1 / tau0^2
  prev <- Inf
  it <- 0L
  for (it in seq_len(max_iter)) {
    # per-curve shifts given shape
    for (j in seq_along(ids)) {
      k <- ci == j
      xj <- x[k]
      yj <- y[k]
      f2 <- function(s) {
        full <- rep(0, 3)
        full[active] <- s
        r <- yj - .shared_pred(theta, matrix(full, 1), xj, rep(1L, length(xj)))
        sum(loss(r)) / sigma^2 + sum(pen * full^2)
      }
      st <- sh[j, active]
      o <- if (sum(active) == 1) {
        oo <- stats::optimize(f2, c(-6, 6))
        list(par = oo$minimum)
      } else {
        stats::optim(st, f2, method = "BFGS")
      }
      sh[j, active] <- o$par
    }
    # shared shape given shifts
    f1 <- function(p) {
      r <- y - .shared_pred(p, sh, x, ci)
      v <- sum(loss(r))
      if (is.null(fix_g) && lambda_g > 0) v <- v + lambda_g * p[5]^2
      if (is.finite(v)) v else 1e10
    }
    o <- tryCatch(stats::nlminb(pmin(pmax(theta, lower), upper), f1,
                                lower = lower, upper = upper),
                  error = function(e) NULL)
    if (is.null(o)) break
    theta <- o$par
    names(theta) <- c("a", "d", "c", "log_b", "log_g")
    res <- y - .shared_pred(theta, sh, x, ci)
    sigma <- if (robust) .robust_scale(res, 0.05) else max(stats::sd(res), 0.05)
    loss <- if (robust) function(r) .huber(r, sigma) else function(r) r^2
    objv <- sum(loss(res)) / sigma^2 + sum(sweep(sh^2, 2, pen, `*`))
    if (abs(prev - objv) < tol * max(1, abs(objv))) break
    prev <- objv
  }
  list(par = theta,
       shifts = data.table::data.table(curve_id = ids, shift_x = sh[, 1],
                                       shift_y = sh[, 2],
                                       shift_bottom = sh[, 3]),
       sigma = sigma, iterations = it)
}

# Longest contiguous run of acceptable standard levels.
# rec: data.table with columns x_true (log2 conc) and in_range (logical),
# one row per level.
.contiguous_range <- function(x_true, in_range) {
  o <- order(x_true)
  x_true <- x_true[o]
  in_range <- in_range[o]
  in_range[is.na(in_range)] <- FALSE
  if (!any(in_range)) return(c(NA_real_, NA_real_))
  r <- rle(in_range)
  ends <- cumsum(r$lengths)
  starts <- ends - r$lengths + 1
  cand <- which(r$values)
  best <- cand[which.max(r$lengths[cand])]
  c(x_true[starts[best]], x_true[ends[best]])
}
