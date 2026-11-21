.theme <- function() {
  ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(legend.position = "bottom",
                   strip.background = ggplot2::element_rect(fill = "grey95"),
                   panel.grid.minor = ggplot2::element_blank())
}

.pick_analytes <- function(d, analytes, max_n = 12) {
  an <- analytes %||% sort(unique(d$analyte))
  if (is.null(analytes) && length(an) > max_n) an <- an[seq_len(max_n)]
  an
}

#' Plot values by plate
#'
#' Box plots of the current values per plate, one panel per analyte; the
#' quickest visual check of plate effects before and after normalisation.
#'
#' @param x A [cyto] object.
#' @param types Sample types to show.
#' @param analytes Analytes to show (default: up to 12).
#' @param color Column used to colour the boxes (default `"batch"`).
#' @return A ggplot object.
#' @examples
#' plot_plates(cyto_transform(cyto_example, "log2"))
#' @export
plot_plates <- function(x, types = "unknown", analytes = NULL,
                        color = "batch") {
  .check_cyto(x)
  d <- x$data[sample_type %in% types & is.finite(value)]
  an <- .pick_analytes(d, analytes)
  d <- d[analyte %in% an]
  ggplot2::ggplot(d, ggplot2::aes(x = .data[["plate"]], y = .data[["value"]],
                                  fill = factor(.data[[color]]))) +
    ggplot2::geom_boxplot(outlier.size = 0.6, linewidth = 0.3) +
    ggplot2::facet_wrap(~analyte, scales = "free_y") +
    ggplot2::labs(x = "Plate", y = x$scale, fill = color) +
    .theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5,
                                                       hjust = 1))
}

#' Plot standard curves
#'
#' Standard wells (current value against log2 expected concentration), one
#' line per plate. If the object has been calibrated with a parametric method
#' (`5pl`, `4pl`, `shared`), the fitted curves are drawn and the working range
#' is shaded. Call it on the calibrated object to see the fits.
#'
#' @param x A [cyto] object.
#' @param analytes Analytes to show (default: up to 12).
#' @param highlight Plates to highlight.
#' @return A ggplot object.
#' @examples
#' plot_standards(cyto_example)
#' plot_standards(cyto_calibrate(cyto_example, "shared"), analytes = "IL-6")
#' @export
plot_standards <- function(x, analytes = NULL, highlight = NULL) {
  .check_cyto(x)
  sig <- if (x$scale == "log2_conc" && "signal" %in% names(x$data)) "signal" else
    "value"
  d <- x$data[sample_type == "standard" & is.finite(expected_conc)]
  d[, y := get(sig)]
  an <- .pick_analytes(d, analytes)
  d <- d[analyte %in% an & is.finite(y)]
  d[, x := log2(expected_conc)]
  d[, plate_label := if (is.null(highlight)) plate else
    ifelse(plate %in% highlight, plate, "other")]
  g <- ggplot2::ggplot(d, ggplot2::aes(x = .data$x, y = .data$y,
                                       group = .data$plate,
                                       colour = .data$plate_label)) +
    ggplot2::geom_point(size = 1, alpha = 0.7)
  fit <- cyto_fit(x, "calibrate")
  cur <- if (is.list(fit) && !is.null(fit$curves)) fit$curves else NULL
  if (!is.null(cur) && all(c("a", "d", "c", "log_b", "log_g") %in% names(cur)) &&
      "plate" %in% names(cur)) {
    cur <- cur[analyte %in% an & is.finite(a)]
    rng <- d[, .(lo = min(x) - 1, hi = max(x) + 1), by = analyte]
    cur <- merge(cur, rng, by = "analyte")
    lines <- cur[, {
      xx <- seq(lo, hi, length.out = 80)
      list(x = xx, y = logistic5(xx, c(a, d, c, log_b, log_g)))
    }, by = .(analyte, plate)]
    lines[, plate_label := if (is.null(highlight)) plate else
      ifelse(plate %in% highlight, plate, "other")]
    g <- g + ggplot2::geom_line(data = lines, linewidth = 0.4)
    wr <- cur[, .(lloq = stats::median(lloq), uloq = stats::median(uloq)),
              by = analyte]
    g <- g + ggplot2::geom_rect(data = wr, inherit.aes = FALSE,
                                ggplot2::aes(xmin = .data$lloq,
                                             xmax = .data$uloq,
                                             ymin = -Inf, ymax = Inf),
                                alpha = 0.06, fill = "steelblue")
  } else {
    g <- g + ggplot2::geom_line(linewidth = 0.3, alpha = 0.6)
  }
  g + ggplot2::facet_wrap(~analyte, scales = "free") +
    ggplot2::labs(x = "log2 expected concentration",
                  y = if (sig == "signal") "signal (pre-calibration)" else x$scale,
                  colour = "plate") +
    .theme()
}

#' Plot estimated corrections
#'
#' Shows the per-group corrections (offsets, and slopes for affine fits)
#' estimated by a normalisation step, per analyte -- a direct view of the
#' plate/batch effects that were removed.
#'
#' @param x A [cyto] object.
#' @param step Step number (default: last normalize step).
#' @return A ggplot object.
#' @examples
#' x <- cyto_normalize(cyto_example, "anchor", anchors = "control",
#'                     scope = "batch")
#' plot_corrections(x)
#' @export
plot_corrections <- function(x, step = NULL) {
  .check_cyto(x)
  f <- if (is.null(step)) cyto_fit(x, "normalize") else cyto_fit(x, step = step)
  if (is.null(f) || !"intercept" %in% names(f)) {
    stop("No normalisation estimates found.", call. = FALSE)
  }
  f <- data.table::as.data.table(f)
  grp <- setdiff(names(f), c("analyte", "intercept", "slope", "offset",
                             "n_anchor", "ref_level", "status", "n_ref",
                             "n_strata"))[1]
  f[, correction := -offset]
  ggplot2::ggplot(f, ggplot2::aes(x = .data[[grp]], y = .data$correction,
                                  colour = .data$analyte,
                                  group = .data$analyte)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
    ggplot2::geom_point() + ggplot2::geom_line(alpha = 0.5) +
    ggplot2::labs(x = grp, y = "correction (log2 units)", colour = "analyte") +
    .theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5,
                                                       hjust = 1))
}

#' @rdname qc_pca
#' @param components Components to plot.
#' @param color Column used for colour (default `"plate"`).
#' @param highlight Plates to highlight (others grey).
#' @param ... Unused.
#' @export
plot.cyto_pca <- function(x, components = c(1, 2), color = "plate",
                          highlight = NULL, ...) {
  s <- data.table::copy(x$scores)
  pcs <- paste0("PC", components)
  if (!is.null(highlight)) {
    s[, .col := ifelse(plate %in% highlight, plate, "other")]
    color <- ".col"
  }
  lab <- sprintf("%s (%.0f%%)", pcs, 100 * x$var_explained[components])
  ggplot2::ggplot(s, ggplot2::aes(x = .data[[pcs[1]]], y = .data[[pcs[2]]],
                                  colour = factor(.data[[color]]))) +
    ggplot2::geom_point(alpha = 0.7, size = 1.3) +
    ggplot2::labs(x = lab[1], y = lab[2],
                  colour = if (color == ".col") "plate" else color) +
    .theme()
}

#' @rdname qc_controls
#' @param analytes Analytes to show.
#' @param material Control materials to show.
#' @param ... Unused.
#' @export
plot.cyto_lj <- function(x, analytes = NULL, material = NULL, ...) {
  d <- data.table::as.data.table(x)
  mats <- material
  if (!is.null(mats)) d <- d[material %in% mats]
  an <- .pick_analytes(d, analytes, 8)
  d <- d[analyte %in% an]
  d[, flagged := rule != ""]
  lines <- unique(d[, .(material, analyte, center, spread)])
  lines <- lines[, .(k = -3:3), by = .(material, analyte, center, spread)]
  lines[, y := center + k * spread]
  ggplot2::ggplot(d, ggplot2::aes(x = .data$run, y = .data$y)) +
    ggplot2::geom_hline(data = lines,
                        ggplot2::aes(yintercept = .data$y,
                                     linetype = abs(.data$k) == 3),
                        colour = "grey60", linewidth = 0.3,
                        show.legend = FALSE) +
    ggplot2::geom_line(colour = "grey40") +
    ggplot2::geom_point(ggplot2::aes(colour = .data$flagged)) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "black", `TRUE` = "red")) +
    ggplot2::facet_grid(material ~ analyte, scales = "free_y") +
    ggplot2::labs(x = "Run (plate order)", y = "value",
                  colour = "Westgard violation") +
    .theme()
}

#' @rdname qc_variance
#' @param ... Unused.
#' @export
plot.cyto_variance <- function(x, ...) {
  d <- data.table::as.data.table(x)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$analyte, y = .data$r2,
                                  fill = .data$term)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Proportion of variance", fill = NULL) +
    .theme()
}

#' @rdname benchmark_pipelines
#' @param metrics Metrics to show (default all).
#' @param ... Unused.
#' @export
plot.cyto_benchmark <- function(x, metrics = NULL, ...) {
  r <- x$results[is.finite(value)]
  if (!is.null(metrics)) r <- r[metric %in% metrics]
  ord <- x$summary$pipeline
  r[, pipeline := factor(pipeline, levels = rev(ord))]
  med <- r[, .(value = mean(value)), by = .(pipeline, metric)]
  ggplot2::ggplot(r, ggplot2::aes(x = .data$value, y = .data$pipeline)) +
    ggplot2::geom_point(alpha = 0.35, size = 1,
                        position = ggplot2::position_jitter(height = 0.15,
                                                            width = 0)) +
    ggplot2::geom_point(data = med, colour = "firebrick", size = 2.5,
                        shape = 18) +
    ggplot2::facet_wrap(~metric, scales = "free_x") +
    ggplot2::labs(x = "metric value (points = analytes, diamond = mean)",
                  y = NULL) +
    .theme()
}
