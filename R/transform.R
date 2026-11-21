#' Choose the signal representation
#'
#' `cyto_transform()` (re)computes the working `value` from the raw
#' fluorescence intensity (`fi`). It defines *what is being measured* and is
#' normally the first processing step. All representations are log2-scaled.
#'
#' | method   | value                                         | removes                              | needs           |
#' |----------|-----------------------------------------------|--------------------------------------|-----------------|
#' | `log2`   | log2(FI)                                      | nothing                              | --              |
#' | `bgsub`  | log2(FI - blank)                              | additive background                  | blanks          |
#' | `mob`    | log2(FI / blank): *multiple of blank*         | plate gain (multiplicative)          | blanks          |
#' | `mos`    | log2(FI / FI of a reference standard)         | plate gain, local potency drift      | standards       |
#' | `moc`    | log2(FI / FI of control material)             | plate gain (Olink-style anchoring)   | controls        |
#' | `vendor` | log2(concentration reported by the software)  | (vendor's own per-plate curves)      | `vendor_conc`   |
#'
#' The ratio representations (`mob`, `mos`, `moc`) need no curve fitting and
#' have no out-of-range values; they trade absolute units for robustness and
#' are often the most powerful choice for relative comparisons (cf. Breen et
#' al. 2016, *Sci Rep* 6:26996). They can still be calibrated later: standards
#' are transformed the same way, so [cyto_calibrate()] works on any of them.
#'
#' @param x A [cyto] object.
#' @param method A registered transform method (see [cyto_methods()]).
#' @param ... Method arguments:
#'
#'   * `log2`: `floor` (default 1) -- FI values below are floored.
#'   * `bgsub`: `blank` selector (default `"blank"`), `floor` (default 1).
#'   * `mob`: `blank` selector (default `"blank"`); `shrink` (default `TRUE`)
#'     -- empirical-Bayes shrinkage of each plate's blank level towards the
#'     mean of its `shrink_by` group (default `"batch"`), which stabilises
#'     plates with only one or two noisy blank wells.
#'   * `mos`: `level` -- standard level used as reference, or `"auto"`
#'     (default: per analyte, the level whose FI is closest to the median FI of
#'     the unknowns); `by_conc` (default `FALSE`) -- if `TRUE`, the reference is
#'     a fixed *concentration* (the auto/selected level's geometric mean
#'     assigned concentration across lots), and each plate's reference FI is
#'     interpolated log-log between its neighbouring standards, which makes the
#'     ratio comparable across kit lots with different standard
#'     concentrations.
#'   * `moc`: `control` selector (default `"control"`); several materials are
#'     combined robustly, missing materials on a plate are handled.
#'   * `vendor`: no arguments.
#' @return A [cyto] object; the estimates (blank levels, reference levels) are
#'   available via [cyto_fit()].
#' @seealso [cyto_normalize()] for plate/batch correction on top of any
#'   representation.
#' @examples
#' x <- cyto_transform(cyto_example, "mob")
#' head(cyto_fit(x))
#' x2 <- cyto_transform(cyto_example, "mos", level = "auto")
#' head(cyto_fit(x2))
#' @export
cyto_transform <- function(x, method = "log2", ...) {
  .check_cyto(x)
  prev <- vapply(x$steps, `[[`, "", "stage")
  if (any(prev %in% c("normalize", "calibrate", "impute"))) {
    warning("cyto_transform() recomputes values from raw FI; earlier ",
            "normalize/calibrate/impute steps are overridden.", call. = FALSE)
  }
  .run_stage(x, "transform", method, list(...))
}

.reset_cols <- function(d) {
  list(censored = rep("none", nrow(d)), limit_low = rep(NA_real_, nrow(d)),
       limit_high = rep(NA_real_, nrow(d)))
}

# per plate x analyte blank level (log2), optionally shrunk
.plate_blank <- function(d, blank = "blank", shrink = TRUE,
                         shrink_by = "batch") {
  isb <- select_rows(d, blank) & d$use_for_fit & is.finite(d$fi) & d$fi > 0
  b <- d[isb, .(lb = mean(log2(fi)), vw = if (.N > 1) stats::var(log2(fi))
                else NA_real_, n = .N), by = c("plate", "analyte")]
  allp <- unique(d[, c("plate", "analyte", shrink_by), with = FALSE])
  b <- merge(allp, b, by = c("plate", "analyte"), all.x = TRUE)
  b[is.na(n), n := 0L]
  b[, blank_raw := lb]
  grp <- c("analyte", shrink_by)
  b[, `:=`(gmean = mean(lb, na.rm = TRUE),
           sw2 = mean(vw, na.rm = TRUE)), by = grp]
  b[!is.finite(sw2), sw2 := 0]
  if (shrink) {
    b[, tau2 := {
      v <- stats::var(lb, na.rm = TRUE)
      if (!is.finite(v)) v <- 0
      max(v - mean(sw2 / pmax(n, 1), na.rm = TRUE), 0)
    }, by = grp]
    b[, weight := ifelse(n > 0, tau2 / (tau2 + sw2 / pmax(n, 1)), 0)]
    b[!is.finite(weight), weight := 1]
  } else {
    b[, weight := ifelse(n > 0, 1, 0)]
  }
  b[, blank_used := weight * ifelse(is.na(lb), 0, lb) + (1 - weight) * gmean]
  if (any(b$n == 0)) {
    warning(sum(b$n == 0), " plate x analyte combination(s) without usable ",
            "blank wells; the group mean blank was used.", call. = FALSE)
  }
  b[, .(plate, analyte, n_blank = n, blank_raw, blank_used, weight)]
}

.tf_log2 <- function(x, floor = 1, ...) {
  d <- x$data
  c(list(value = log2(pmax(d$fi, floor)), scale = "log2_fi",
         note = "log2(FI)"), .reset_cols(d))
}

.tf_bgsub <- function(x, blank = "blank", floor = 1, ...) {
  d <- x$data
  b <- .plate_blank(d, blank, shrink = FALSE)
  bg <- 2^.lookup(d, b, c("plate", "analyte"), "blank_used")
  c(list(value = log2(pmax(d$fi - bg, floor)), scale = "log2_fi", fit = b,
         note = "log2(FI - plate blank)"), .reset_cols(d))
}

.tf_mob <- function(x, blank = "blank", shrink = TRUE, shrink_by = "batch",
                    ...) {
  d <- x$data
  b <- .plate_blank(d, blank, shrink = shrink, shrink_by = shrink_by)
  lb <- .lookup(d, b, c("plate", "analyte"), "blank_used")
  c(list(value = log2(pmax(d$fi, 1e-3)) - lb, scale = "log2_ratio", fit = b,
         note = "log2(FI / plate blank)"), .reset_cols(d))
}

.tf_mos <- function(x, level = "auto", by_conc = FALSE, ...) {
  d <- x$data
  std <- d[sample_type == "standard" & use_for_fit & is.finite(fi) & fi > 0 &
             !is.na(level)]
  if (!nrow(std)) stop("No standard wells available for 'mos'.", call. = FALSE)
  lv <- std[, .(lfi = stats::median(log2(fi)),
                lconc = stats::median(log2(expected_conc))),
            by = .(analyte, level)]
  if (identical(level, "auto")) {
    unk <- d[sample_type == "unknown" & is.finite(fi) & fi > 0,
             .(ulfi = stats::median(log2(fi))), by = analyte]
    lv2 <- merge(lv, unk, by = "analyte")
    sel <- lv2[, .SD[which.min(abs(lfi - ulfi))], by = analyte][
      , .(analyte, ref_level = level)]
  } else {
    sel <- data.table::data.table(analyte = unique(d$analyte),
                                  ref_level = as.numeric(
                                    sub("^S", "", as.character(level))))
  }
  if (!by_conc) {
    r <- merge(std, sel, by = "analyte")[level == ref_level,
                                         .(ref_fi = mean(log2(fi)), n = .N),
                                         by = .(plate, analyte, ref_level)]
  } else {
    if (all(is.na(std$expected_conc))) {
      stop("by_conc = TRUE needs expected_conc for standards.", call. = FALSE)
    }
    rc <- merge(std, sel, by = "analyte")[level == ref_level,
                                          .(ref_conc = mean(log2(expected_conc),
                                                            na.rm = TRUE)),
                                          by = analyte]
    lvp <- std[, .(lfi = mean(log2(fi)), lconc = log2(expected_conc[1])),
               by = .(plate, analyte, level)]
    lvp <- merge(lvp, rc, by = "analyte")
    r <- lvp[, {
      o <- order(lconc)
      v <- stats::approx(lconc[o], lfi[o], xout = ref_conc[1], rule = 2)$y
      list(ref_fi = v, n = .N, ref_conc = ref_conc[1])
    }, by = .(plate, analyte)]
    r <- merge(r, sel, by = "analyte")
  }
  miss <- unique(d[, .(plate, analyte)])[!r, on = c("plate", "analyte")]
  if (nrow(miss)) {
    warning(nrow(miss), " plate x analyte combination(s) lack the reference ",
            "standard; values set to NA.", call. = FALSE)
  }
  ref <- .lookup(d, r, c("plate", "analyte"), "ref_fi")
  c(list(value = log2(pmax(d$fi, 1e-3)) - ref, scale = "log2_ratio", fit = r,
         note = if (by_conc) "log2(FI / FI at reference concentration)" else
           "log2(FI / FI of reference standard)"), .reset_cols(d))
}

.tf_moc <- function(x, control = "control", ...) {
  d <- data.table::copy(x$data)
  d[, value := log2(pmax(fi, 1e-3))]
  a <- select_rows(d, control) & d$use_for_fit
  if (!any(a)) stop("No control rows selected for 'moc'.", call. = FALSE)
  off <- .anchor_offsets(d, a, group = "plate", scope = NULL)
  ref <- .lookup(d, off$offsets, c("plate", "analyte"), "ref_level")
  c(list(value = d$value - ref, scale = "log2_ratio", fit = off$offsets,
         note = "log2(FI / control material level on plate)"), .reset_cols(d))
}

.tf_vendor <- function(x, ...) {
  d <- x$data
  if (all(is.na(d$vendor_conc) & d$vendor_flag == "")) {
    stop("No vendor concentrations in the data.", call. = FALSE)
  }
  v <- ifelse(is.finite(d$vendor_conc) & d$vendor_conc > 0,
              log2(d$vendor_conc), NA_real_)
  cens <- ifelse(d$vendor_flag %in% c("<", "OOR <"), "left",
                 ifelse(d$vendor_flag %in% c(">", "OOR >"), "right", "none"))
  v[cens != "none"] <- NA_real_
  # censoring limits: most extreme reported value per plate x analyte
  tmp <- data.table::data.table(.row = d$.row, plate = d$plate,
                                analyte = d$analyte, v = v, cens = cens)
  lim <- tmp[cens == "none" & is.finite(v), .(lo = min(v), hi = max(v)),
             by = .(plate, analyte)]
  lo <- .lookup(tmp, lim, c("plate", "analyte"), "lo")
  hi <- .lookup(tmp, lim, c("plate", "analyte"), "hi")
  qc <- .append_qc(d$qc, "extrapolated", d$vendor_flag == "*")
  list(value = v, scale = "log2_conc", censored = cens,
       limit_low = ifelse(cens == "left", lo, NA_real_),
       limit_high = ifelse(cens == "right", hi, NA_real_), qc = qc,
       note = "log2(vendor concentration)")
}
