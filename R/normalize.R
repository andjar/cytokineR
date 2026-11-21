#' Remove plate and batch effects
#'
#' `cyto_normalize()` estimates, for every group (plate or batch) and analyte,
#' a correction `value' = intercept + slope * value` on the log2 scale and
#' applies it to *all* rows of the group (standards, blanks and controls
#' included, so that a later [cyto_calibrate()] sees consistently corrected
#' standards). The estimated corrections are returned by [cyto_fit()].
#' `cyto_correct_batch()` is shorthand for `by = "batch"` (default method
#' `"bridge"`).
#'
#' @section Methods:
#' * `median` -- **assumes comparable samples on every plate** (randomised
#'   layout). Each group's centre (median / mean / trimmed mean) of the
#'   `reference` rows (default: unknowns) is aligned to the overall centre.
#'   With `within = c("study", "time")` the centres are compared *within*
#'   strata and pooled, which protects designs where plates differ in
#'   composition; groups that are alone in a stratum are left unchanged.
#' * `quantile` -- like `median` but aligns the whole distribution (quantile
#'   mapping). Strongest assumption; use with fully randomised layouts.
#' * `anchor` -- uses material measured in several groups (`anchors`: pooled
#'   controls, kit controls, standards, bridging samples). Each material's
#'   consensus level (median across groups, or its level in group `ref`) is the
#'   target; the group offset is the median deviation over materials. With
#'   `fit = "affine"` a robust (Theil-Sen) line through the materials also
#'   corrects compression/expansion of the signal (needs >= 3 materials
#'   spanning the range, e.g. standards). Use `scope = "batch"` to align plates
#'   only within each batch (needed when anchors such as standards differ
#'   between kit lots).
#' * `bridge` -- `anchor` with `by = "batch"`, `anchors = "bridge"` (samples
#'   measured in more than one batch) and the batch with most bridging samples
#'   as reference; `fit = "deming"` estimates a slope with errors-in-variables
#'   regression.
#' * `regress` -- per analyte linear model `value ~ protect + group` on the
#'   reference rows; the group effects are removed while biological covariates
#'   in `protect` (e.g. `~ time + group`) are protected (cf.
#'   `limma::removeBatchEffect`). Handles unbalanced designs if the covariates
#'   are known.
#' * `combat` -- empirical-Bayes location/scale adjustment via `sva::ComBat`
#'   (optional dependency). Only the reference rows are adjusted (standards
#'   are not), so do not calibrate afterwards.
#'
#' @param x A [cyto] object.
#' @param method A registered normalize method (see [cyto_methods()]).
#' @param by Grouping whose effect is removed: `"plate"` (default) or
#'   `"batch"` (any column works).
#' @param ... Method arguments: `reference` (selector, see [selectors]),
#'   `within`, `stat` (`"median"`, `"mean"`, `"trimmed"`), `trim`, `scope`,
#'   `anchors`, `ref`, `fit`, `protect`.
#' @return A [cyto] object.
#' @examples
#' x <- cyto_transform(cyto_example, "log2")
#' x <- cyto_normalize(x, "anchor", anchors = "control", by = "plate",
#'                     scope = "batch")
#' head(cyto_fit(x))
#' x <- cyto_correct_batch(x)   # bridging samples
#' cyto_fit(x)
#' @export
cyto_normalize <- function(x, method = "median", by = "plate", ...) {
  .run_stage(x, "normalize", method, c(list(by = by), list(...)))
}

#' @rdname cyto_normalize
#' @export
cyto_correct_batch <- function(x, method = "bridge", ...) {
  .run_stage(x, "normalize", method, c(list(by = "batch"), list(...)))
}

# Apply per-group affine corrections; returns list(value, limit_low, limit_high)
.apply_affine <- function(d, tab, by) {
  keys <- c(by, "analyte")
  icpt <- .lookup(d, tab, keys, "intercept")
  slp <- .lookup(d, tab, keys, "slope")
  icpt[is.na(icpt)] <- 0
  slp[is.na(slp)] <- 1
  list(value = icpt + slp * d$value,
       limit_low = icpt + slp * d$limit_low,
       limit_high = icpt + slp * d$limit_high)
}

.check_by <- function(d, by) {
  miss <- setdiff(by, names(d))
  if (length(miss)) stop("Unknown grouping column(s): ",
                         paste(miss, collapse = ", "), call. = FALSE)
}

# Material-based offsets (used by moc, anchor, bridge).
# a: logical rows used as anchors. group: column(s) defining groups.
# scope: column(s) within which consensus levels are computed.
.anchor_offsets <- function(d, a, group = "plate", scope = NULL, ref = NULL,
                            fit = c("offset", "affine"), stat = "median") {
  fit <- match.arg(fit)
  ad <- d[a & is.finite(value)]
  if (!nrow(ad)) stop("No usable anchor rows.", call. = FALSE)
  # materials: standards are material per level within scope (lot); others by id
  ad[, material := sample_id]
  sc <- if (is.null(scope)) character() else scope
  ym <- ad[, .(y = mean(value), nw = .N), by = c(group, sc, "material", "analyte")]
  if (is.null(ref)) {
    cons <- ym[, .(consensus = stats::median(y), ng = .N),
               by = c(sc, "material", "analyte")]
  } else {
    g1 <- group[1]
    cons <- ym[get(g1) %in% ref, .(consensus = mean(y), ng = .N),
               by = c(sc, "material", "analyte")]
    if (!nrow(cons)) stop("Reference group '", paste(ref, collapse = ","),
                          "' has no anchor material.", call. = FALSE)
  }
  ym <- merge(ym, cons, by = c(sc, "material", "analyte"))
  ym[, dev := y - consensus]
  if (fit == "offset") {
    off <- ym[, .(offset = .center(dev, stat), n_anchor = .N,
                  ref_level = mean(consensus) + .center(dev, stat)),
              by = c(group, "analyte")]
    off[, `:=`(intercept = -offset, slope = 1)]
  } else {
    off <- ym[, {
      if (.N >= 3 && diff(range(y)) > 1) {
        ts <- .theil_sen(y, consensus)
        b <- min(max(ts[2], 0.5), 2)
        a0 <- stats::median(consensus - b * y)
        list(intercept = a0, slope = b, n_anchor = .N,
             offset = stats::median(y - consensus), ref_level = NA_real_)
      } else {
        o <- .center(dev, stat)
        list(intercept = -o, slope = 1, n_anchor = .N, offset = o,
             ref_level = NA_real_)
      }
    }, by = c(group, "analyte")]
  }
  # groups without anchors
  allg <- unique(d[, c(group, "analyte"), with = FALSE])
  off <- merge(allg, off, by = c(group, "analyte"), all.x = TRUE)
  off[, status := ifelse(is.na(n_anchor), "no_anchor", "ok")]
  if (any(off$status == "no_anchor")) {
    warning(sum(off$status == "no_anchor"), " group x analyte combination(s) ",
            "have no anchor material and were left uncorrected.", call. = FALSE)
    off[status == "no_anchor", `:=`(intercept = 0, slope = 1, offset = 0,
                                    n_anchor = 0L)]
  }
  list(offsets = off, deviations = ym)
}

.nm_anchor <- function(x, by = "plate", anchors, scope = NULL, ref = NULL,
                       fit = c("offset", "affine"), stat = "median", ...) {
  if (missing(anchors)) stop("'anchors' must be given (e.g. \"control\").",
                             call. = FALSE)
  fit <- match.arg(fit)
  d <- x$data
  .check_by(d, c(by, scope))
  a <- select_rows(d, anchors) & d$use_for_fit
  off <- .anchor_offsets(d, a, group = by, scope = scope, ref = ref, fit = fit,
                         stat = stat)
  res <- .apply_affine(d, off$offsets, by)
  c(res, list(fit = off$offsets,
              note = sprintf("%s anchors, %s fit", sum(a), fit)))
}

.nm_bridge <- function(x, by = "batch", anchors = "bridge", ref = NULL,
                       fit = c("offset", "deming"), stat = "median", ...) {
  fit <- match.arg(fit)
  d <- x$data
  a <- select_rows(d, anchors) & d$use_for_fit & is.finite(d$value)
  if (!any(a)) stop("No bridging samples available for fitting (none measured ",
                    "in more than one ", by, ", or all excluded/held out).",
                    call. = FALSE)
  ad <- d[a]
  if (is.null(ref)) {
    cnt <- ad[, .(n = data.table::uniqueN(sample_id)), by = by]
    ref <- cnt[[by]][which.max(cnt$n)]
  }
  ym <- ad[, .(y = mean(value)), by = c(by, "sample_id", "analyte")]
  r <- ym[get(by) == ref, .(sample_id, analyte, yref = y)]
  ym <- merge(ym, r, by = c("sample_id", "analyte"))
  tab <- ym[, {
    if (fit == "deming" && .N >= 3) {
      cf <- .deming(y, yref)
      list(intercept = cf[1], slope = cf[2], n_anchor = .N,
           offset = .center(y - yref, stat))
    } else {
      o <- .center(y - yref, stat)
      list(intercept = -o, slope = 1, n_anchor = .N, offset = o)
    }
  }, by = c(by, "analyte")]
  allg <- unique(d[, c(by, "analyte"), with = FALSE])
  tab <- merge(allg, tab, by = c(by, "analyte"), all.x = TRUE)
  tab[, status := ifelse(is.na(n_anchor), "no_bridge", "ok")]
  tab[get(by) == ref, status := "reference"]
  if (any(tab$status == "no_bridge")) {
    nb <- tab[status == "no_bridge"]
    warning("No usable bridging samples (shared with reference ", by, " '",
            ref, "') for ", nrow(nb), " ", by, " x analyte combination(s), ",
            "e.g. ", paste(utils::head(paste(nb[[by]], nb$analyte, sep = "/"),
                                       4), collapse = ", "),
            "; left uncorrected. Censored bridging values are a common ",
            "cause: correct batches before calibration.", call. = FALSE)
  }
  tab[is.na(intercept), `:=`(intercept = 0, slope = 1, offset = 0,
                             n_anchor = 0L)]
  res <- .apply_affine(d, tab, by)
  c(res, list(fit = tab, note = sprintf("reference %s = %s", by, ref)))
}

.nm_median <- function(x, by = "plate", reference = "unknown", within = NULL,
                       stat = c("median", "mean", "trimmed"), trim = 0.1,
                       scope = NULL, ...) {
  stat <- match.arg(stat)
  d <- x$data
  .check_by(d, c(by, within, scope))
  r <- select_rows(d, reference) & d$use_for_fit & is.finite(d$value)
  rd <- d[r]
  if (!nrow(rd)) stop("No reference rows selected.", call. = FALSE)
  st <- c(within, scope)
  if (length(st)) {
    rd[, stratum := do.call(paste, c(.SD, sep = "|")), .SDcols = st]
  } else {
    rd[, stratum := "all"]
  }
  cg <- rd[, .(center = .center(value, stat, trim), n = .N),
           by = c(by, "stratum", "analyte")]
  tg <- rd[, .(target = .center(value, stat, trim)), by = c("stratum", "analyte")]
  cg <- merge(cg, tg, by = c("stratum", "analyte"))
  # groups alone in their stratum carry no information
  cg[, ng := .N, by = c("stratum", "analyte")]
  cg <- cg[ng > 1]
  tab <- cg[, .(offset = .wmedian(center - target, n), n_ref = sum(n),
                n_strata = .N), by = c(by, "analyte")]
  tab[, `:=`(intercept = -offset, slope = 1)]
  allg <- unique(d[, c(by, "analyte"), with = FALSE])
  tab <- merge(allg, tab, by = c(by, "analyte"), all.x = TRUE)
  tab[is.na(offset), `:=`(offset = 0, intercept = 0, slope = 1, n_ref = 0L,
                          n_strata = 0L)]
  res <- .apply_affine(d, tab, by)
  c(res, list(fit = tab, note = sprintf("%s-centring of %d reference rows",
                                        stat, nrow(rd))))
}

.nm_quantile <- function(x, by = "plate", reference = "unknown", scope = NULL,
                         ...) {
  d <- x$data
  .check_by(d, c(by, scope))
  r <- select_rows(d, reference) & d$use_for_fit & is.finite(d$value)
  sc <- if (is.null(scope)) character() else scope
  newv <- d$value
  grp_cols <- c(by, "analyte")
  keys <- unique(d[, c(sc, "analyte"), with = FALSE])
  fits <- list()
  for (k in seq_len(nrow(keys))) {
    in_k <- rep(TRUE, nrow(d))
    for (cc in names(keys)) in_k <- in_k & d[[cc]] == keys[[cc]][k]
    tv <- d$value[in_k & r]
    if (length(tv) < 5) next
    gids <- unique(d[in_k, by, with = FALSE])
    for (gi in seq_len(nrow(gids))) {
      in_g <- in_k
      for (cc in by) in_g <- in_g & d[[cc]] == gids[[cc]][gi]
      gv <- d$value[in_g & r]
      if (length(gv) < 5) next
      p <- stats::ppoints(min(length(gv), 50))
      qg <- stats::quantile(gv, p, names = FALSE)
      qt <- stats::quantile(tv, p, names = FALSE)
      y <- d$value[in_g]
      out <- stats::approx(qg, qt, xout = y, ties = mean)$y
      lo <- is.finite(y) & y < qg[1]
      hi <- is.finite(y) & y > qg[length(qg)]
      out[lo] <- y[lo] + (qt[1] - qg[1])
      out[hi] <- y[hi] + (qt[length(qt)] - qg[length(qg)])
      newv[in_g] <- out
      fits[[length(fits) + 1]] <- cbind(gids[gi], analyte = keys$analyte[k],
                                        median_shift = stats::median(qt - qg))
    }
  }
  fit <- data.table::rbindlist(fits)
  list(value = newv, fit = fit, note = "quantile mapping")
}

.nm_regress <- function(x, by = "plate", reference = "unknown", protect = NULL,
                        ...) {
  d <- x$data
  .check_by(d, by)
  r <- select_rows(d, reference) & d$use_for_fit & is.finite(d$value)
  gcol <- if (length(by) > 1) {
    do.call(paste, c(d[, by, with = FALSE], sep = "|"))
  } else d[[by]]
  tabs <- list()
  for (an in unique(d$analyte)) {
    k <- r & d$analyte == an
    df <- as.data.frame(d[k])
    df$.g <- factor(gcol[k])
    if (nlevels(df$.g) < 2) next
    rhs <- ".g"
    if (!is.null(protect)) {
      rhs <- paste(c(attr(stats::terms(protect), "term.labels"), ".g"),
                   collapse = " + ")
    }
    f <- stats::as.formula(paste("value ~", rhs))
    m <- tryCatch(stats::lm(f, data = df, contrasts = list(.g = "contr.sum")),
                  error = function(e) NULL)
    if (is.null(m)) next
    cf <- stats::coef(m)
    gl <- levels(df$.g)
    eff <- cf[paste0(".g", seq_len(length(gl) - 1))]
    eff[is.na(eff)] <- 0
    eff <- c(eff, -sum(eff))
    tabs[[an]] <- data.table::data.table(.g = gl, analyte = an,
                                         offset = unname(eff))
  }
  tab <- data.table::rbindlist(tabs)
  key <- data.table::data.table(.row = d$.row, .g = gcol, analyte = d$analyte)
  off <- .lookup(key, tab, c(".g", "analyte"), "offset")
  off[is.na(off)] <- 0
  data.table::setnames(tab, ".g", paste(by, collapse = "|"))
  tab[, `:=`(intercept = -offset, slope = 1)]
  list(value = d$value - off, limit_low = d$limit_low - off,
       limit_high = d$limit_high - off, fit = tab,
       note = if (is.null(protect)) "group effects removed" else
         paste("group effects removed, protecting", deparse(protect)))
}

.nm_combat <- function(x, by = "batch", reference = "unknown", protect = NULL,
                       ...) {
  .require("sva", "for method = 'combat'")
  d <- x$data
  r <- select_rows(d, reference) & is.finite(d$value)
  rd <- d[r]
  rd[, .wid := paste(plate, well, sep = "@")]
  w <- data.table::dcast(rd, .wid ~ analyte, value.var = "value")
  mat <- t(as.matrix(w[, -1]))
  colnames(mat) <- w$.wid
  if (anyNA(mat)) stop("ComBat cannot handle missing values; impute first.",
                       call. = FALSE)
  info <- unique(rd[, c(".wid", by, all.vars(protect)), with = FALSE])
  info <- info[match(w$.wid, info$.wid)]
  mod <- if (is.null(protect)) NULL else
    stats::model.matrix(protect, data = as.data.frame(info))
  adj <- sva::ComBat(mat, batch = info[[by]], mod = mod)
  long <- data.table::as.data.table(as.table(adj))
  data.table::setnames(long, c("analyte", ".wid", "newv"))
  rd2 <- merge(rd[, .(.row, .wid, analyte)], long, by = c(".wid", "analyte"))
  newv <- d$value
  newv[match(rd2$.row, d$.row)] <- rd2$newv
  list(value = newv, note = "ComBat on reference rows only",
       fit = data.table::data.table(n_adjusted = nrow(rd2)))
}
