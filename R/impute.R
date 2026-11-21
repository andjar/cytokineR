#' Fill censored values
#'
#' After calibration (or with vendor concentrations), values outside the
#' working range are `NA` with `censored = "left"`/`"right"` and known limits.
#' `cyto_impute()` replaces them; the rows stay marked as censored (and get QC
#' flag `imputed`) so downstream analyses can treat them specially.
#'
#' @section Methods:
#' * `fraction` -- left-censored values become `limit_low + log2(frac)`
#'   (default `frac = 0.5`, i.e. LLOQ/2); right-censored become
#'   `limit_high + log2(high)` (default `high = 1`). Simple and transparent, but
#'   biased and variance-shrinking.
#' * `tobit` -- fits a censored-normal (Tobit) model per analyte (optionally per
#'   `by` group, optionally with covariates `formula`, e.g. `~ time`) to the
#'   `reference` rows by maximum likelihood, and replaces each censored value
#'   by its conditional expectation given that it lies beyond its limit, or by a
#'   random draw (`draw = TRUE`, for multiple imputation).
#' * `rank_tobit` -- like `tobit`, but uses the information still present in
#'   the fluorescence signal: censored values within a `rank_by` group
#'   (default plate x analyte) are ordered by their pre-calibration `signal`
#'   and assigned evenly spaced quantiles of their truncated distribution.
#'   Values stay below the LLOQ, but their ranking is preserved.
#'
#' The distribution below the limit can only be estimated when enough values
#' are observed: with fewer than `min_obs` (default 5) or `min_frac` (default
#' 10%) quantified reference values in a group, `tobit` falls back to LLOQ/2
#' and `rank_tobit` to values between LLOQ/4 and LLOQ ordered by signal; the
#' fit table reports `status = "fallback_too_few_observed"`. The SD of the
#' censored-normal model is bounded to keep fits with heavy censoring stable.
#'
#' @param x A [cyto] object.
#' @param method A registered impute method (see [cyto_methods()]).
#' @param ... Method arguments: `frac`, `high` (fraction); `reference`
#'   (selector of rows used to fit, default `"unknown"`), `target` (rows to
#'   impute, default `c("unknown", "control")`), `by`, `formula`, `draw`,
#'   `seed`, `min_obs`, `min_frac` (tobit); `rank_by` (rank_tobit).
#' @return A [cyto] object; [cyto_fit()] returns the fitted means/SDs.
#' @examples
#' x <- cyto_calibrate(cyto_example, "5pl")
#' x <- cyto_impute(x, "tobit")
#' head(cyto_fit(x))
#' @export
cyto_impute <- function(x, method = "tobit", ...) {
  .check_cyto(x)
  .run_stage(x, "impute", method, list(...))
}

.target_rows <- function(d, target) {
  select_rows(d, target) & d$censored %in% c("left", "right")
}

.im_fraction <- function(x, frac = 0.5, high = 1,
                         target = c("unknown", "control"), ...) {
  d <- x$data
  t <- .target_rows(d, target)
  v <- d$value
  l <- t & d$censored == "left" & is.finite(d$limit_low)
  r <- t & d$censored == "right" & is.finite(d$limit_high)
  v[l] <- d$limit_low[l] + log2(frac)
  v[r] <- d$limit_high[r] + log2(high)
  list(value = v, qc = .append_qc(d$qc, "imputed", l | r),
       note = sprintf("%d left, %d right imputed", sum(l), sum(r)))
}

# Censored-normal MLE; returns list(beta, sigma, mu (per row of X))
.tobit_mle <- function(y, cens, L, U, X) {
  o <- cens == "none" & is.finite(y)
  l <- cens == "left" & is.finite(L)
  r <- cens == "right" & is.finite(U)
  if (sum(o) < 3) return(NULL)
  m0 <- stats::lm.fit(X[o, , drop = FALSE], y[o])
  b0 <- m0$coefficients
  b0[is.na(b0)] <- 0
  s0 <- max(stats::sd(y[o]), 0.1)
  smax <- max(3 * stats::sd(c(y[o], L[l], U[r])), 3)
  nll <- function(p) {
    b <- p[-length(p)]
    s <- exp(p[length(p)])
    mu <- drop(X %*% b)
    -(sum(stats::dnorm(y[o], mu[o], s, log = TRUE)) +
        sum(stats::pnorm((L[l] - mu[l]) / s, log.p = TRUE)) +
        sum(stats::pnorm((U[r] - mu[r]) / s, lower.tail = FALSE, log.p = TRUE)))
  }
  np <- length(b0) + 1
  lo <- c(rep(-Inf, np - 1), log(0.05))
  hi <- c(rep(Inf, np - 1), log(smax))
  st <- c(b0, min(max(log(s0), lo[np] + 1e-6), hi[np] - 1e-6))
  fit <- tryCatch(stats::optim(st, nll, method = "L-BFGS-B", lower = lo,
                               upper = hi),
                  error = function(e) NULL)
  if (is.null(fit) || !all(is.finite(fit$par))) return(NULL)
  b <- fit$par[-length(fit$par)]
  list(beta = b, sigma = exp(fit$par[length(fit$par)]),
       mu = drop(X %*% b), converged = fit$convergence == 0,
       n_obs = sum(o), n_left = sum(l), n_right = sum(r))
}

.im_tobit_core <- function(x, reference = "unknown",
                           target = c("unknown", "control"), by = NULL,
                           formula = NULL, draw = FALSE, seed = NULL,
                           rank_by = NULL, min_obs = 5, min_frac = 0.1) {
  d <- x$data
  if (!is.null(seed)) set.seed(seed)
  ref <- select_rows(d, reference) & d$use_for_fit
  tgt <- .target_rows(d, target)
  grp <- c("analyte", by)
  .check_by(d, grp)
  gk <- do.call(paste, c(d[, grp, with = FALSE], sep = "|"))
  v <- d$value
  fits <- list()
  for (g in unique(gk[tgt])) {
    kr <- ref & gk == g
    kt <- tgt & gk == g
    dd <- as.data.frame(d[kr | kt])
    if (is.null(formula)) {
      X <- matrix(1, nrow(dd), 1)
    } else {
      mf <- stats::model.frame(formula, dd, na.action = stats::na.pass)
      X <- stats::model.matrix(formula, mf)
      # rows with missing covariates (e.g. controls) get the average design
      bad <- !stats::complete.cases(X)
      if (any(bad)) {
        cm <- colMeans(X[!bad, , drop = FALSE])
        X[bad, ] <- matrix(cm, sum(bad), ncol(X), byrow = TRUE)
      }
    }
    isref <- (kr)[kr | kt]
    n_obs <- sum(dd$censored[isref] == "none" & is.finite(dd$value[isref]))
    enough <- n_obs >= max(min_obs, ceiling(min_frac * sum(isref)))
    f <- if (enough) .tobit_mle(dd$value[isref], dd$censored[isref],
                                dd$limit_low[isref], dd$limit_high[isref],
                                X[isref, , drop = FALSE]) else NULL
    if (is.null(f)) {
      # Too little observed data to estimate the distribution below the
      # limit: LLOQ/2 (or, for rank_tobit, LLOQ/4..LLOQ spread by signal rank)
      l <- kt & d$censored == "left"
      r <- kt & d$censored == "right"
      v[l] <- d$limit_low[l] - 1
      v[r] <- d$limit_high[r]
      if (!is.null(rank_by)) {
        rk <- do.call(paste, c(d[l, rank_by, with = FALSE], sep = "|"))
        il <- which(l)
        for (h in unique(rk)) {
          ii <- il[rk == h]
          if (length(ii) < 2) next
          q <- (rank(d$signal[ii], ties.method = "average") - 0.5) / length(ii)
          v[ii] <- d$limit_low[ii] - 2 + 2 * q
        }
      }
      fits[[g]] <- data.table::data.table(group = g, mean = NA_real_,
                                          sigma = NA_real_, n_obs = n_obs,
                                          status = "fallback_too_few_observed")
      next
    }
    mu <- drop(X %*% f$beta)
    s <- f$sigma
    rows <- which(kr | kt)
    tt <- which(kt[rows])
    for (j in tt) {
      i <- rows[j]
      if (d$censored[i] == "left") {
        a <- (d$limit_low[i] - mu[j]) / s
        pa <- max(stats::pnorm(a), 1e-12)
        if (draw) {
          v[i] <- mu[j] + s * stats::qnorm(stats::runif(1) * pa)
        } else {
          v[i] <- mu[j] - s * stats::dnorm(a) / pa
        }
      } else {
        b <- (d$limit_high[i] - mu[j]) / s
        pb <- max(1 - stats::pnorm(b), 1e-12)
        if (draw) {
          v[i] <- mu[j] + s * stats::qnorm(1 - stats::runif(1) * pb)
        } else {
          v[i] <- mu[j] + s * stats::dnorm(b) / pb
        }
      }
    }
    if (!is.null(rank_by)) {
      # rank-preserving assignment within rank_by groups
      rk <- do.call(paste, c(d[rows[tt], rank_by, with = FALSE], sep = "|"))
      for (side in c("left", "right")) {
        for (h in unique(rk)) {
          jj <- tt[rk == h & d$censored[rows[tt]] == side]
          if (length(jj) < 2) next
          ii <- rows[jj]
          sig <- d$signal[ii]
          q <- (rank(sig, ties.method = "average") - 0.5) / length(ii)
          if (side == "left") {
            pa <- pmax(stats::pnorm((d$limit_low[ii] - mu[jj]) / s), 1e-12)
            v[ii] <- mu[jj] + s * stats::qnorm(q * pa)
          } else {
            pb <- pmax(1 - stats::pnorm((d$limit_high[ii] - mu[jj]) / s), 1e-12)
            v[ii] <- mu[jj] + s * stats::qnorm(1 - (1 - q) * pb)
          }
        }
      }
    }
    fits[[g]] <- data.table::data.table(group = g, mean = mean(mu[isref]),
                                        sigma = s, n_obs = f$n_obs,
                                        n_left = f$n_left, n_right = f$n_right,
                                        status = if (f$converged) "ok" else
                                          "not_converged")
  }
  fit <- data.table::rbindlist(fits, fill = TRUE)
  if (nrow(fit)) {
    parts <- data.table::tstrsplit(fit$group, "|", fixed = TRUE)
    for (i in seq_along(grp)) data.table::set(fit, j = grp[i], value = parts[[i]])
  }
  imputed <- tgt & is.finite(v)
  list(value = v, qc = .append_qc(d$qc, "imputed", imputed), fit = fit,
       note = sprintf("%d values imputed", sum(imputed)))
}

.im_tobit <- function(x, reference = "unknown", target = c("unknown", "control"),
                      by = NULL, formula = NULL, draw = FALSE, seed = NULL,
                      min_obs = 5, min_frac = 0.1, ...) {
  .im_tobit_core(x, reference, target, by, formula, draw, seed, rank_by = NULL,
                 min_obs = min_obs, min_frac = min_frac)
}

.im_rank_tobit <- function(x, reference = "unknown",
                           target = c("unknown", "control"), by = NULL,
                           formula = NULL, rank_by = "plate", min_obs = 5,
                           min_frac = 0.1, ...) {
  if (!"signal" %in% names(x$data)) {
    stop("'rank_tobit' needs the pre-calibration 'signal' column; run ",
         "cyto_calibrate() first.", call. = FALSE)
  }
  .im_tobit_core(x, reference, target, by, formula, draw = FALSE, seed = NULL,
                 rank_by = rank_by, min_obs = min_obs, min_frac = min_frac)
}
