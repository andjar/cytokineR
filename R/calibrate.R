#' Convert signal to concentration
#'
#' `cyto_calibrate()` fits standard curves to the standard wells (their current
#' `value` against `log2(expected_conc)`) and inverts them for every row. It
#' works on any signal representation (`log2_fi` or `log2_ratio`), so it can be
#' run after plate harmonisation ("harmonise first, calibrate once"). The
#' result is the log2 concentration of the original sample (well
#' concentration times `dilution`).
#'
#' @section Methods:
#' * `5pl` / `4pl` -- the conventional approach: one five- (four-) parameter
#'   logistic per curve (default per plate x analyte), fitted in log-log space
#'   by penalised least squares with several starting values. Fragile when a
#'   plate has few usable standards or a bad standard well.
#' * `spline` -- model-free: monotone (isotonic, then monotone Hermite spline)
#'   interpolation through the per-level mean signals. Never extrapolates;
#'   flat ends of the curve (signal steps smaller than `min_step` log2 units
#'   between adjacent levels) are trimmed from the working range.
#' * `shared` -- a shared-shape 5PL: within each `shape_by` group (default
#'   `"batch"`, i.e. kit lot) all plates share the curve *shape*, and each plate
#'   only gets a horizontal (`shift = "x"`, potency) and/or vertical
#'   (`shift = "y"`, gain) shift, shrunk towards zero (prior SD `tau`, log2
#'   units). Borrowing strength across plates makes the curve robust to single
#'   bad standards and lets plates without usable standards borrow the shared
#'   curve (flagged `borrowed_curve`).
#'
#' @section Working range and censoring:
#' Standards are back-calculated; a level is acceptable when its recovery
#' (observed / expected) lies within `recovery` (default 70--130%). The
#' longest contiguous run of acceptable levels defines LLOQ and ULOQ. With
#' `censor = "loq"` (default) values outside this range are set to `NA` and
#' marked `censored = "left"`/`"right"` with the limits stored in
#' `limit_low`/`limit_high` (on the value scale), ready for [cyto_impute()].
#' With `censor = "asymptote"` only values beyond the curve asymptotes are
#' censored and extrapolated values are kept (flagged `extrapolated`). The
#' pre-calibration value is stored in column `signal`.
#'
#' @param x A [cyto] object on a signal scale.
#' @param method A registered calibrate method (see [cyto_methods()]).
#' @section Harmonise first, calibrate once:
#' A plate correction applied *before* per-plate calibration is undone by it,
#' because the standards move with the samples. To keep a correction, fit
#' fewer curves: `curve_by = "batch"` (one curve per batch after aligning its
#' plates) or `curve_by = NULL` (one curve per analyte for the whole study),
#' optionally restricted to the standards of a reference lot with
#' `standards = ~ batch == "1"` -- after aligning batches with bridging
#' samples, this puts every batch on the reference lot's concentration scale.
#'
#' @param ... Method arguments: `curve_by` (default `"plate"`; `NULL` for one
#'   curve per analyte), `standards` (selector restricting the standards used
#'   for fitting, see [selectors]), `recovery`
#'   (default `c(0.7, 1.3)`), `censor` (`"loq"` or `"asymptote"`),
#'   `lambda_g` (5PL asymmetry penalty), `robust` (Huber loss; default
#'   `FALSE` for 5pl/4pl, `TRUE` for shared), `drop_outliers` (max. number of
#'   standard wells to drop per curve when their residual exceeds 3.5 robust
#'   SDs, 5pl/4pl only; default 0), `min_step` (spline), `shape_by`, `shift`,
#'   `tau` (shared).
#' @return A [cyto] object on the `log2_conc` scale. [cyto_fit()] returns a
#'   list with `curves` (parameters, LLOQ/ULOQ per curve) and `standards`
#'   (back-calculated standards with recoveries).
#' @examples
#' x <- cyto_calibrate(cyto_example, "5pl")
#' head(cyto_fit(x)$curves)
#' xs <- cyto_calibrate(cyto_example, "shared", shift = "both")
#' head(cyto_fit(xs)$curves)
#' @export
cyto_calibrate <- function(x, method = "5pl", ...) {
  .check_cyto(x)
  .need_scale(x, c("log2_fi", "log2_ratio"), "cyto_calibrate()")
  .run_stage(x, "calibrate", method, list(...))
}

.curve_key <- function(d, curve_by) {
  .check_by(d, curve_by)
  do.call(paste, c(d[, c(curve_by, "analyte"), with = FALSE], sep = "|"))
}

.std_rows <- function(d, standards = NULL) {
  s <- d$sample_type == "standard" & d$use_for_fit & is.finite(d$value) &
    is.finite(d$expected_conc) & d$expected_conc > 0
  if (!is.null(standards)) s <- s & select_rows(d, standards)
  s
}

# back-calculation & working range for a parametric curve
.param_range <- function(xs, ys, lv, par, recovery) {
  xh <- logistic5_inverse(ys, par)
  rec <- data.table::data.table(level = lv, x_true = xs, x_hat = as.numeric(xh))
  lvt <- rec[, .(x_true = x_true[1],
                 recovery = if (all(is.na(x_hat))) NA_real_ else
                   2^(mean(x_hat, na.rm = TRUE) - x_true[1])),
             by = level]
  rlo <- recovery[1]
  rhi <- recovery[2]
  lvt[, in_range := !is.na(recovery) & recovery >= rlo & recovery <= rhi]
  rng <- .contiguous_range(lvt$x_true, lvt$in_range)
  list(range = rng, levels = lvt, wells = rec)
}

.cal_finish <- function(d, xhat, side, lloq, uloq, censor, failed,
                        extra_qc = NULL) {
  dil <- log2(d$dilution)
  side[is.na(side)] <- "none"
  out_lo <- is.finite(xhat) & is.finite(lloq) & xhat < lloq
  out_hi <- is.finite(xhat) & is.finite(uloq) & xhat > uloq
  if (censor == "loq") {
    left <- side == "left" | out_lo
    right <- side == "right" | out_hi
  } else {
    left <- side == "left"
    right <- side == "right"
  }
  value <- xhat + dil
  value[left | right] <- NA_real_
  cens <- ifelse(left, "left", ifelse(right, "right", "none"))
  cens[failed] <- "none"
  value[failed] <- NA_real_
  qc <- .append_qc(d$qc, "curve_failed", failed)
  if (censor == "asymptote") {
    qc <- .append_qc(qc, "extrapolated", (out_lo | out_hi) & !failed)
  }
  if (!is.null(extra_qc)) {
    for (nm in names(extra_qc)) qc <- .append_qc(qc, nm, extra_qc[[nm]])
  }
  list(value = value, censored = cens, limit_low = lloq + dil,
       limit_high = uloq + dil, qc = qc, scale = "log2_conc",
       columns = list(signal = d$value))
}

# parametric calibration given a table of curve parameters
.invert_param <- function(d, key, pars) {
  xhat <- rep(NA_real_, nrow(d))
  side <- rep("none", nrow(d))
  idx <- split(seq_len(nrow(d)), key)
  for (k in names(idx)) {
    p <- pars[curve_id == k]
    if (!nrow(p) || anyNA(p[, .(a, d, c, log_b, log_g)])) next
    ii <- idx[[k]]
    inv <- logistic5_inverse(d$value[ii],
                             c(p$a, p$d, p$c, p$log_b, p$log_g))
    xhat[ii] <- as.numeric(inv)
    side[ii] <- attr(inv, "side")
  }
  list(xhat = xhat, side = side)
}

.cal_logistic <- function(x, curve_by = "plate", fix_g = NULL, lambda_g = 0.1,
                          recovery = c(0.7, 1.3), censor = c("loq", "asymptote"),
                          drop_outliers = 0, standards = NULL, robust = FALSE,
                          ...) {
  censor <- match.arg(censor)
  d <- x$data
  key <- .curve_key(d, curve_by)
  s <- .std_rows(d, standards)
  if (!any(s)) stop("No usable standards (need expected_conc).", call. = FALSE)
  sidx <- split(which(s), key[s])
  curves <- list()
  stdtab <- list()
  for (k in names(sidx)) {
    ii <- sidx[[k]]
    xs <- log2(d$expected_conc[ii])
    ys <- d$value[ii]
    lv <- d$level[ii]
    keep <- rep(TRUE, length(ii))
    fit <- fit_logistic5(xs, ys, fix_g = fix_g, lambda_g = lambda_g,
                         robust = robust)
    nd <- 0
    while (nd < drop_outliers && !anyNA(fit$par) && sum(keep) > 6) {
      r <- ys - logistic5(xs, fit$par)
      sc <- .robust_scale(r[keep], 0.05)
      r[!keep] <- 0
      j <- which.max(abs(r))
      if (abs(r[j]) < 3.5 * sc) break
      keep[j] <- FALSE
      nd <- nd + 1
      fit <- fit_logistic5(xs[keep], ys[keep], fix_g = fix_g,
                           lambda_g = lambda_g, robust = robust)
    }
    if (anyNA(fit$par)) {
      curves[[k]] <- data.table::data.table(curve_id = k, a = NA_real_,
        d = NA_real_, c = NA_real_, log_b = NA_real_, log_g = NA_real_,
        sigma = NA_real_, n_std = length(ii), n_dropped = nd,
        converged = FALSE, lloq = NA_real_, uloq = NA_real_)
      next
    }
    pr <- .param_range(xs[keep], ys[keep], lv[keep], fit$par, recovery)
    curves[[k]] <- data.table::data.table(curve_id = k, a = fit$par[[1]],
      d = fit$par[[2]], c = fit$par[[3]], log_b = fit$par[[4]],
      log_g = fit$par[[5]], sigma = fit$sigma, n_std = length(ii),
      n_dropped = nd, converged = fit$converged, lloq = pr$range[1],
      uloq = pr$range[2])
    st <- pr$levels
    st[, curve_id := k]
    stdtab[[k]] <- st
  }
  curves <- data.table::rbindlist(curves)
  inv <- .invert_param(d, key, curves)
  lloq <- curves$lloq[match(key, curves$curve_id)]
  uloq <- curves$uloq[match(key, curves$curve_id)]
  failed <- !key %in% curves$curve_id[is.finite(curves$lloq)]
  res <- .cal_finish(d, inv$xhat, inv$side, lloq, uloq, censor, failed)
  res$fit <- list(curves = .split_key(curves, curve_by),
                  standards = .split_key(data.table::rbindlist(stdtab), curve_by))
  res$note <- sprintf("%d curves, %d failed", nrow(curves),
                      sum(!is.finite(curves$lloq)))
  res
}

.cal_5pl <- function(x, ...) .cal_logistic(x, ...)
.cal_4pl <- function(x, ...) .cal_logistic(x, fix_g = 1, ...)

.split_key <- function(tab, curve_by) {
  if (!nrow(tab)) return(tab)
  parts <- data.table::tstrsplit(tab$curve_id, "|", fixed = TRUE)
  nms <- c(curve_by, "analyte")
  for (i in seq_along(nms)) data.table::set(tab, j = nms[i], value = parts[[i]])
  data.table::setcolorder(tab, c("curve_id", nms))
  tab
}

.cal_spline <- function(x, curve_by = "plate", min_step = 0.1,
                        censor = c("loq", "asymptote"), standards = NULL, ...) {
  censor <- match.arg(censor)
  d <- x$data
  key <- .curve_key(d, curve_by)
  s <- .std_rows(d, standards)
  if (!any(s)) stop("No usable standards (need expected_conc).", call. = FALSE)
  sidx <- split(which(s), key[s])
  xhat <- rep(NA_real_, nrow(d))
  side <- rep("none", nrow(d))
  lloq <- uloq <- rep(NA_real_, nrow(d))
  ridx <- split(seq_len(nrow(d)), key)
  curves <- list()
  for (k in names(ridx)) {
    ii <- sidx[[k]]
    if (is.null(ii)) next
    agg <- data.table::data.table(x = log2(d$expected_conc[ii]),
                                  y = d$value[ii])[, .(y = mean(y)), by = x]
    data.table::setorder(agg, x)
    if (nrow(agg) < 2) next
    yf <- stats::isoreg(agg$x, agg$y)$yf
    st <- diff(yf) >= min_step
    if (!any(st)) next
    r <- rle(st)
    ends <- cumsum(r$lengths)
    starts <- ends - r$lengths + 1
    cand <- which(r$values)
    b <- cand[which.max(r$lengths[cand])]
    i1 <- starts[b]
    i2 <- ends[b] + 1
    xs <- agg$x[i1:i2]
    ys <- yf[i1:i2]
    f <- if (length(xs) >= 3) stats::splinefun(ys, xs, method = "monoH.FC") else
      function(v) stats::approx(ys, xs, xout = v)$y
    rr <- ridx[[k]]
    v <- d$value[rr]
    inside <- is.finite(v) & v >= ys[1] & v <= ys[length(ys)]
    xhat[rr[inside]] <- f(v[inside])
    side[rr[is.finite(v) & v < ys[1]]] <- "left"
    side[rr[is.finite(v) & v > ys[length(ys)]]] <- "right"
    lloq[rr] <- xs[1]
    uloq[rr] <- xs[length(xs)]
    curves[[k]] <- data.table::data.table(curve_id = k, n_levels = nrow(agg),
                                          n_used = length(xs), lloq = xs[1],
                                          uloq = xs[length(xs)])
  }
  curves <- data.table::rbindlist(curves)
  failed <- !key %in% curves$curve_id
  res <- .cal_finish(d, xhat, side, lloq, uloq, censor, failed)
  res$fit <- list(curves = .split_key(curves, curve_by), standards = NULL)
  res$note <- sprintf("%d monotone interpolation curves", nrow(curves))
  res
}

.cal_shared <- function(x, shape_by = "batch", curve_by = "plate",
                        shift = c("x", "y", "bottom"),
                        tau = c(x = 1, y = 1, bottom = 1),
                        fix_g = NULL, lambda_g = 0.1, recovery = c(0.7, 1.3),
                        censor = c("loq", "asymptote"), standards = NULL,
                        robust = TRUE, ...) {
  censor <- match.arg(censor)
  if (identical(shift, "both")) shift <- c("x", "y")
  shift <- match.arg(shift, c("x", "y", "bottom"), several.ok = TRUE)
  if (is.null(names(tau))) {
    tau <- stats::setNames(rep(tau, length.out = 3), c("x", "y", "bottom"))
  }
  rlo <- recovery[1]
  rhi <- recovery[2]
  d <- x$data
  .check_by(d, c(shape_by, curve_by))
  key <- .curve_key(d, curve_by)
  gkey <- .curve_key(d, shape_by)
  s <- .std_rows(d, standards)
  if (!any(s)) stop("No usable standards (need expected_conc).", call. = FALSE)
  curves <- list()
  stdtab <- list()
  for (g in unique(gkey)) {
    ig <- which(s & gkey == g)
    plates_g <- unique(key[gkey == g])
    if (length(ig) < 4) next
    sh <- .fit_shared5(log2(d$expected_conc[ig]), d$value[ig], key[ig],
                       shift = shift, tau = tau, fix_g = fix_g,
                       lambda_g = lambda_g, robust = robust)
    if (is.null(sh)) next
    th <- sh$par
    shifts <- merge(data.table::data.table(curve_id = plates_g), sh$shifts,
                    by = "curve_id", all.x = TRUE)
    shifts[, borrowed := is.na(shift_x)]
    shifts[is.na(shift_x), `:=`(shift_x = 0, shift_y = 0, shift_bottom = 0)]
    shifts[, `:=`(group = g, a = th[[1]] + shift_bottom + shift_y,
                  d = th[[2]] + shift_y, c = th[[3]] + shift_x,
                  log_b = th[[4]], log_g = th[[5]], sigma = sh$sigma)]
    # pooled working range from all plates' back-calculated standards
    xh <- .invert_param(d[ig], key[ig], shifts)$xhat
    rec <- data.table::data.table(level = d$level[ig], curve_id = key[ig],
                                  x_true = log2(d$expected_conc[ig]),
                                  x_hat = xh)
    lvt <- rec[, .(x_true = x_true[1],
                   recovery = if (all(is.na(x_hat))) NA_real_ else
                     2^(stats::median(x_hat - x_true, na.rm = TRUE))),
               by = level]
    lvt[, in_range := !is.na(recovery) & recovery >= rlo & recovery <= rhi]
    rng <- .contiguous_range(lvt$x_true, lvt$in_range)
    shifts[, `:=`(lloq = rng[1], uloq = rng[2])]
    curves[[g]] <- shifts
    lvt[, group := g]
    stdtab[[g]] <- lvt
  }
  curves <- data.table::rbindlist(curves)
  if (!nrow(curves)) stop("Shared calibration failed for all groups.",
                          call. = FALSE)
  inv <- .invert_param(d, key, curves)
  lloq <- curves$lloq[match(key, curves$curve_id)]
  uloq <- curves$uloq[match(key, curves$curve_id)]
  failed <- !key %in% curves$curve_id[is.finite(curves$lloq)]
  borrowed <- key %in% curves$curve_id[curves$borrowed]
  res <- .cal_finish(d, inv$xhat, inv$side, lloq, uloq, censor, failed,
                     extra_qc = list(borrowed_curve = borrowed))
  res$fit <- list(curves = .split_key(curves, curve_by),
                  standards = data.table::rbindlist(stdtab))
  res$note <- sprintf("shared shape per %s, shifts: %s",
                      paste(shape_by, collapse = "+"),
                      paste(shift, collapse = "+"))
  res
}
