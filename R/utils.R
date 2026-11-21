# ---- row selection -------------------------------------------------------

#' Selecting rows (anchors, references, hold-outs)
#'
#' Several functions take a *selector* argument (`anchors`, `reference`,
#' `control`, `holdout`, ...) that identifies a set of rows. A selector can be:
#'
#' * a character vector, matched against `sample_type` (e.g. `"control"`,
#'   `"standard"`) **and** `sample_id` (e.g. `"pool_A"`, `"S3"`); the special
#'   keyword `"bridge"` selects unknown samples whose `sample_id` occurs in
#'   more than one batch;
#' * a one-sided formula evaluated in the long data, e.g.
#'   `~ sample_type == "control" & sample_id != "kit_ctrl"`;
#' * a function of the long data returning a logical vector;
#' * a logical vector with one entry per row.
#'
#' @name selectors
#' @param data The long data (`x$data`).
#' @param sel A selector.
#' @return A logical vector.
#' @examples
#' d <- cyto_example$data
#' table(select_rows(d, "control"))
#' table(select_rows(d, ~ sample_type == "unknown" & batch == "1"))
#' table(select_rows(d, "bridge"))
#' @export
select_rows <- function(data, sel) {
  n <- nrow(data)
  if (is.null(sel)) return(rep(TRUE, n))
  if (is.logical(sel)) {
    if (length(sel) != n) stop("Logical selector has wrong length.", call. = FALSE)
    sel[is.na(sel)] <- FALSE
    return(sel)
  }
  if (inherits(sel, "formula")) {
    rhs <- sel[[length(sel)]]
    out <- eval(rhs, data, environment(sel))
    out[is.na(out)] <- FALSE
    return(as.logical(out))
  }
  if (is.function(sel)) {
    out <- sel(data)
    out[is.na(out)] <- FALSE
    return(as.logical(out))
  }
  if (is.character(sel)) {
    out <- data$sample_type %in% sel | data$sample_id %in% sel
    if ("bridge" %in% sel) out <- out | .is_bridge(data)
    return(out)
  }
  stop("Unsupported selector of class ", class(sel)[1], call. = FALSE)
}

.is_bridge <- function(data) {
  u <- data[sample_type == "unknown", .(nb = data.table::uniqueN(batch)),
            by = sample_id][nb > 1, sample_id]
  data$sample_type == "unknown" & data$sample_id %in% u
}

# ---- robust summaries ------------------------------------------------------

.center <- function(v, stat = c("median", "mean", "trimmed"), trim = 0.1) {
  stat <- match.arg(stat)
  v <- v[is.finite(v)]
  if (!length(v)) return(NA_real_)
  switch(stat,
         median = stats::median(v),
         mean = mean(v),
         trimmed = mean(v, trim = trim))
}

# weighted median
.wmedian <- function(v, w = rep(1, length(v))) {
  ok <- is.finite(v) & is.finite(w) & w > 0
  v <- v[ok]
  w <- w[ok]
  if (!length(v)) return(NA_real_)
  o <- order(v)
  v <- v[o]
  w <- w[o]
  cw <- cumsum(w) / sum(w)
  v[which(cw >= 0.5)[1]]
}

# Deming / orthogonal-type regression of y on x with variance ratio lambda
# (error variance of y over error variance of x). Returns c(intercept, slope).
.deming <- function(x, y, lambda = 1) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 3) return(c(NA_real_, NA_real_))
  mx <- mean(x)
  my <- mean(y)
  sxx <- mean((x - mx)^2)
  syy <- mean((y - my)^2)
  sxy <- mean((x - mx) * (y - my))
  if (abs(sxy) < 1e-12) return(c(my - mx, 1))
  b <- (syy - lambda * sxx + sqrt((syy - lambda * sxx)^2 +
                                    4 * lambda * sxy^2)) / (2 * sxy)
  c(my - b * mx, b)
}

# Theil-Sen robust line
.theil_sen <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  if (n < 2) return(c(NA_real_, NA_real_))
  ij <- utils::combn(n, 2)
  dx <- x[ij[2, ]] - x[ij[1, ]]
  keep <- abs(dx) > 1e-12
  if (!any(keep)) return(c(stats::median(y - x), 1))
  b <- stats::median((y[ij[2, keep]] - y[ij[1, keep]]) / dx[keep])
  c(stats::median(y - b * x), b)
}

.append_qc <- function(old, new, where) {
  out <- old
  add <- where & !is.na(where)
  out[add] <- ifelse(out[add] == "", new, paste(out[add], new, sep = ";"))
  out
}

.require <- function(pkg, why) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is needed ", why, ". Install it with ",
         "install.packages('", pkg, "') or BiocManager::install('", pkg, "').",
         call. = FALSE)
  }
}

# Merge a per-group table back onto rows (keeps original row order)
.lookup <- function(d, tab, by, col) {
  key <- d[, c(".row", by), with = FALSE]
  m <- merge(key, tab[, c(by, col), with = FALSE], by = by, all.x = TRUE,
             sort = FALSE)
  data.table::setorderv(m, ".row")
  m[[col]]
}

.need_scale <- function(x, allowed, what) {
  if (!x$scale %in% allowed) {
    stop(what, " requires the value scale to be one of ",
         paste(allowed, collapse = ", "), " (current: ", x$scale, ").",
         call. = FALSE)
  }
}
