# ---- flagging ----------------------------------------------------------------

#' Flag measurements for quality control
#'
#' Flagging functions add a reason to the `qc` column and (by default) set
#' `use_for_fit = FALSE`, so that flagged rows are not used to *estimate*
#' normalisation or calibration parameters. Values are not removed: flagged
#' unknowns are still transformed and calibrated, and can be dropped later
#' with [cyto_filter()] (e.g. `cyto_filter(x, qc == "")`).
#'
#' * `flag_beads()` -- bead count below `min_count` (Bio-Rad recommends at
#'   least 50 beads per region for reliable medians; 35 is a common floor).
#' * `flag_replicates()` -- replicate groups (same sample on the same plate)
#'   whose FI coefficient of variation exceeds `max_cv`.
#' * `flag_wells()` -- arbitrary condition, e.g. a known pipetting error.
#'
#' @param x A [cyto] object.
#' @param min_count Minimum bead count.
#' @param exclude Set `use_for_fit = FALSE` for flagged rows?
#' @return A [cyto] object.
#' @examples
#' x <- flag_beads(cyto_example, min_count = 35)
#' x <- flag_replicates(x, max_cv = 0.3, types = "standard")
#' x <- flag_wells(x, plate == "b1p1" & well == "A1", reason = "pipetting")
#' x
#' table(x$data$qc)
#' @export
flag_beads <- function(x, min_count = 35, exclude = TRUE) {
  .check_cyto(x)
  x <- .cyto_copy(x)
  d <- x$data
  low <- is.finite(d$bead_count) & d$bead_count < min_count
  d[, qc := .append_qc(qc, "low_beads", low)]
  if (exclude) d[low, use_for_fit := FALSE]
  x$data <- d
  .add_step(x, "qc", "flag_beads", list(min_count = min_count,
                                        exclude = exclude),
            note = sprintf("%d flagged", sum(low)), scale_before = x$scale,
            n_changed = 0L)
}

#' @rdname flag_beads
#' @param max_cv Maximum coefficient of variation of FI within a replicate
#'   group.
#' @param types Sample types to check.
#' @export
flag_replicates <- function(x, max_cv = 0.3,
                            types = c("standard", "blank", "control"),
                            exclude = FALSE) {
  .check_cyto(x)
  tab <- qc_replicates(x, types = types)
  bad <- tab[cv > max_cv]
  x <- .cyto_copy(x)
  d <- x$data
  key <- paste(d$plate, d$sample_id, d$analyte)
  hit <- key %in% paste(bad$plate, bad$sample_id, bad$analyte) &
    d$sample_type %in% types
  d[, qc := .append_qc(qc, "high_cv", hit)]
  if (exclude) d[hit, use_for_fit := FALSE]
  x$data <- d
  .add_step(x, "qc", "flag_replicates",
            list(max_cv = max_cv, types = types, exclude = exclude),
            note = sprintf("%d replicate groups flagged", nrow(bad)),
            scale_before = x$scale, n_changed = 0L)
}

#' @rdname flag_beads
#' @param condition A logical expression evaluated in the long data.
#' @param reason Text added to the `qc` column.
#' @export
flag_wells <- function(x, condition, reason = "manual", exclude = TRUE) {
  .check_cyto(x)
  e <- substitute(condition)
  hit <- eval(e, x$data, parent.frame())
  hit[is.na(hit)] <- FALSE
  x <- .cyto_copy(x)
  d <- x$data
  d[, qc := .append_qc(qc, reason, hit)]
  if (exclude) d[hit, use_for_fit := FALSE]
  x$data <- d
  .add_step(x, "qc", "flag_wells",
            list(condition = paste(deparse(e), collapse = ""),
                 reason = reason, exclude = exclude),
            note = sprintf("%d rows flagged", sum(hit)),
            scale_before = x$scale, n_changed = 0L)
}

# ---- summaries ---------------------------------------------------------------

#' Replicate agreement
#'
#' Coefficient of variation (raw FI) and SD (current log2 value) for every
#' replicate group: wells with the same `sample_id` on the same plate.
#'
#' @param x A [cyto] object.
#' @param types Sample types to include (default all).
#' @return A data.table with one row per plate x sample x analyte with at least
#'   two wells.
#' @examples
#' r <- qc_replicates(cyto_example, types = "standard")
#' head(r[order(-cv)])
#' @export
qc_replicates <- function(x, types = NULL) {
  .check_cyto(x)
  d <- x$data
  if (!is.null(types)) d <- d[sample_type %in% types]
  d[, .(n = .N, mean_fi = mean(fi, na.rm = TRUE),
        cv = stats::sd(fi, na.rm = TRUE) / mean(fi, na.rm = TRUE),
        sd_value = stats::sd(value, na.rm = TRUE)),
    by = .(batch, plate, sample_type, sample_id, analyte)][n > 1]
}

#' PCA-based quality control
#'
#' Principal component analysis of the current values (one row per well,
#' one column per analyte). Missing values are replaced by the analyte median
#' for the PCA only. Besides scores and loadings, the result quantifies how
#' much each component is explained by plate and batch, and how far each
#' plate's centroid lies from the others (robust z-score) -- large values point
#' to plates that need attention.
#'
#' @param x A [cyto] object.
#' @param types Sample types to include (default `"unknown"`).
#' @param analytes Analytes to include (default all).
#' @param scale Scale analytes to unit variance?
#' @param n_pc Number of components to keep in the summaries.
#' @return An object of class `cyto_pca` (list with `scores`, `loadings`,
#'   `var_explained`, `by_plate`, `group_r2`); has a `plot()` method.
#' @examples
#' p <- qc_pca(cyto_example)
#' p$group_r2
#' p$by_plate
#' plot(p, color = "batch")
#' @export
qc_pca <- function(x, types = "unknown", analytes = NULL, scale = TRUE,
                   n_pc = 5) {
  .check_cyto(x)
  w <- cyto_wide(x, types = types, as_data_table = TRUE)
  an <- analytes %||% sort(unique(x$data$analyte))
  an <- intersect(an, names(w))
  m <- as.matrix(w[, an, with = FALSE])
  for (j in seq_len(ncol(m))) {
    v <- m[, j]
    v[!is.finite(v)] <- stats::median(v[is.finite(v)])
    m[, j] <- v
  }
  keep <- apply(m, 2, function(v) all(is.finite(v)) && stats::sd(v) > 0)
  m <- m[, keep, drop = FALSE]
  if (ncol(m) < 2) stop("Need at least two analytes with variation.",
                        call. = FALSE)
  pr <- stats::prcomp(m, center = TRUE, scale. = scale)
  k <- min(n_pc, ncol(pr$x))
  ids <- setdiff(names(w), an)
  scores <- cbind(w[, ids, with = FALSE],
                  data.table::as.data.table(pr$x[, seq_len(k), drop = FALSE]))
  ve <- pr$sdev^2 / sum(pr$sdev^2)
  pcs <- paste0("PC", seq_len(k))
  r2 <- function(v, g) {
    g <- factor(g)
    if (nlevels(g) < 2) return(NA_real_)
    1 - sum(stats::residuals(stats::lm(v ~ g))^2) / sum((v - mean(v))^2)
  }
  group_r2 <- data.table::data.table(
    component = pcs, var_explained = ve[seq_len(k)],
    r2_plate = vapply(pcs, function(p) r2(scores[[p]], scores$plate), 0),
    r2_batch = vapply(pcs, function(p) r2(scores[[p]], scores$batch), 0))
  kk <- min(k, 3)
  cen <- scores[, lapply(.SD, stats::median), by = .(batch, plate),
                .SDcols = pcs[seq_len(kk)]]
  cm <- as.matrix(cen[, pcs[seq_len(kk)], with = FALSE])
  zc <- apply(cm, 2, function(v) {
    s <- stats::mad(v)
    if (!is.finite(s) || s == 0) s <- stats::sd(v)
    (v - stats::median(v)) / s
  })
  zc <- matrix(zc, nrow = nrow(cm))
  cen[, dist := sqrt(rowSums(zc^2))]
  structure(list(scores = scores, loadings = pr$rotation[, seq_len(k),
                                                         drop = FALSE],
                 var_explained = ve, by_plate = cen[order(-dist)],
                 group_r2 = group_r2, scale = x$scale),
            class = "cyto_pca")
}

#' @export
print.cyto_pca <- function(x, ...) {
  cat("<cyto_pca> ", nrow(x$scores), " wells, ", nrow(x$loadings),
      " analytes\n", sep = "")
  print(x$group_r2, digits = 3)
  cat("Plates furthest from the others (robust distance of PC centroids):\n")
  print(utils::head(x$by_plate, 5), digits = 3)
  invisible(x)
}

#' Plate-level quality summary
#'
#' One row per plate combining several robust z-scores (relative to the other
#' plates, using median and MAD):
#' `z_samples` (median of the unknowns), `z_blank` (blank level),
#' `z_standards` (standard curve level), `z_controls` (control materials),
#' and `pca_dist` (PCA centroid distance, see [qc_pca()]). Plates with any
#' |z| (or `pca_dist`) above `k` are flagged.
#'
#' @param x A [cyto] object (typically on the `log2_fi` scale).
#' @param k Threshold for flagging.
#' @return A data.table.
#' @examples
#' qc_plates(cyto_example)
#' @export
qc_plates <- function(x, k = 3) {
  .check_cyto(x)
  d <- x$data[is.finite(value)]
  rz <- function(v) {
    s <- stats::mad(v, na.rm = TRUE)
    if (!is.finite(s) || s == 0) s <- stats::sd(v, na.rm = TRUE)
    if (!is.finite(s) || s == 0) return(rep(0, length(v)))
    (v - stats::median(v, na.rm = TRUE)) / s
  }
  per <- function(type) {
    t <- d[sample_type == type, .(m = stats::median(value)),
           by = .(batch, plate, analyte)]
    if (!nrow(t)) return(NULL)
    t[, z := rz(m), by = analyte]
    t[, .(z = stats::median(z)), by = .(batch, plate)]
  }
  out <- unique(d[, .(batch, plate)])
  for (tp in c("unknown", "blank", "standard", "control")) {
    t <- per(tp)
    nm <- paste0("z_", c(unknown = "samples", blank = "blank",
                         standard = "standards", control = "controls")[[tp]])
    if (is.null(t)) {
      out[, (nm) := NA_real_]
    } else {
      data.table::setnames(t, "z", nm)
      out <- merge(out, t, by = c("batch", "plate"), all.x = TRUE)
    }
  }
  p <- tryCatch(qc_pca(x)$by_plate[, .(plate, pca_dist = dist)],
                error = function(e) NULL)
  if (!is.null(p)) out <- merge(out, p, by = "plate", all.x = TRUE) else
    out[, pca_dist := NA_real_]
  zc <- c(grep("^z_", names(out), value = TRUE), "pca_dist")
  out[, flagged := apply(abs(as.matrix(.SD)) > k, 1, any, na.rm = TRUE),
      .SDcols = zc]
  out[is.na(flagged), flagged := FALSE]
  out[order(-flagged, plate)]
}

#' Levey-Jennings monitoring of control materials
#'
#' For each control material x analyte, plate means are ordered by run order
#' and compared with the material's centre and spread across plates. Westgard
#' rules are evaluated: `1_3s` (one point beyond 3 SD), `2_2s` (two
#' consecutive beyond 2 SD on the same side), `R_4s` (within a plate, one
#' material above +2 SD and another below -2 SD), `4_1s` (four consecutive
#' beyond 1 SD on the same side) and `10_x` (ten consecutive on the same side).
#'
#' @param x A [cyto] object.
#' @param controls Selector of control rows (see [selectors]).
#' @param order_by Column giving run order (default: order of first
#'   appearance of plates).
#' @param robust Use median/MAD (default) instead of mean/SD.
#' @return An object of class `cyto_lj` (a data.table) with a `plot()` method.
#' @examples
#' lj <- qc_controls(cyto_example)
#' lj[rule != ""]
#' plot(lj, analytes = c("IL-6", "IL-8"))
#' @export
qc_controls <- function(x, controls = "control", order_by = NULL,
                        robust = TRUE) {
  .check_cyto(x)
  d <- x$data
  a <- select_rows(d, controls) & is.finite(d$value)
  if (!any(a)) stop("No control rows selected.", call. = FALSE)
  porder <- if (is.null(order_by)) unique(d$plate) else
    unique(d[order(get(order_by)), plate])
  t <- d[a, .(y = mean(value)), by = .(material = sample_id, analyte, batch,
                                       plate)]
  t[, run := match(plate, porder)]
  data.table::setorder(t, material, analyte, run)
  t[, `:=`(center = if (robust) stats::median(y) else mean(y),
           spread = if (robust) stats::mad(y) else stats::sd(y)),
    by = .(material, analyte)]
  t[!is.finite(spread) | spread == 0, spread := NA_real_]
  t[, z := (y - center) / spread]
  runlen <- function(cond) {
    # TRUE where cond holds for the current and preceding points (run length)
    r <- integer(length(cond))
    for (i in seq_along(cond)) r[i] <- if (isTRUE(cond[i]))
      (if (i > 1) r[i - 1] else 0L) + 1L else 0L
    r
  }
  t[, rule := {
    zz <- z
    r13 <- abs(zz) > 3
    r22 <- runlen(zz > 2) >= 2 | runlen(zz < -2) >= 2
    r41 <- runlen(zz > 1) >= 4 | runlen(zz < -1) >= 4
    r10 <- runlen(zz > 0) >= 10 | runlen(zz < 0) >= 10
    out <- character(.N)
    out <- .append_qc(out, "1_3s", r13)
    out <- .append_qc(out, "2_2s", r22)
    out <- .append_qc(out, "4_1s", r41)
    out <- .append_qc(out, "10_x", r10)
    out
  }, by = .(material, analyte)]
  r4 <- t[, .(r4 = any(z > 2, na.rm = TRUE) & any(z < -2, na.rm = TRUE)),
          by = .(analyte, plate)]
  t <- merge(t, r4, by = c("analyte", "plate"))
  t[, rule := .append_qc(rule, "R_4s", r4 & abs(z) > 2)]
  t[, r4 := NULL]
  data.table::setorder(t, material, analyte, run)
  class(t) <- c("cyto_lj", class(t))
  t
}

#' Variance decomposition
#'
#' Sequential (type I) proportion of variance of each analyte explained by
#' known covariates (entered first) and technical factors such as batch and
#' plate (in the order given; with plates nested in batches, the plate term is
#' the within-batch plate effect). Useful before/after normalisation.
#'
#' @param x A [cyto] object.
#' @param factors Technical grouping columns.
#' @param covariates Optional biological covariates (column names) entered
#'   before the factors.
#' @param types Sample types to include.
#' @return An object of class `cyto_variance` (data.table with `analyte`,
#'   `term`, `r2`) with a `plot()` method.
#' @examples
#' v <- qc_variance(cyto_transform(cyto_example, "log2"),
#'                  covariates = "time")
#' v[term == "plate"]
#' plot(v)
#' @export
qc_variance <- function(x, factors = c("batch", "plate"), covariates = NULL,
                        types = "unknown") {
  .check_cyto(x)
  d <- x$data[sample_type %in% types & is.finite(value)]
  terms <- c(covariates, factors)
  .check_by(d, terms)
  out <- d[, {
    df <- as.data.frame(.SD)
    tss <- sum((df$value - mean(df$value))^2)
    rss_prev <- tss
    res <- list()
    for (i in seq_along(terms)) {
      used <- terms[seq_len(i)]
      f <- stats::reformulate(paste0("factor(`", used, "`)"), "value")
      rss <- tryCatch(sum(stats::residuals(stats::lm(f, data = df))^2),
                      error = function(e) NA_real_)
      res[[terms[i]]] <- (rss_prev - rss) / tss
      rss_prev <- rss
    }
    res[["residual"]] <- rss_prev / tss
    list(term = names(res), r2 = unlist(res))
  }, by = analyte, .SDcols = c("value", terms)]
  out[, term := factor(term, levels = c(terms, "residual"))]
  class(out) <- c("cyto_variance", class(out))
  out
}
