#' Simulate a multiplex cytokine study with known truth
#'
#' Generates plates with standards, blanks, controls and study samples under a
#' physically motivated model of bead-based immunoassays, including the
#' imperfections that make pre-processing hard: plate-specific gain,
#' potency shifts and background, batch (kit lot) effects, errors in the
#' assigned standard concentrations of each lot, occasional bad standard
#' wells, low bead counts, and analytes whose biological levels sit below the
#' lowest standard. Because the true concentration of every well is stored in
#' column `truth`, any pipeline can be scored against reality with
#' [metric_truth()] -- the fastest way to develop and compare methods.
#'
#' @section Model:
#' For analyte *a*, plate *p* and well concentration *C* (pg/mL in the well),
#' \deqn{FI = 2^{\gamma_p} [ B_a 2^{\beta_p} + (T_a - B_a) /
#'   (1 + \exp(-b_a (\log_2 C + \eta_p - c_a)))^{g_a} ]}
#' with multiplicative noise (CV `noise_cv`), additive noise (SD `noise_add`
#' FI) and extra noise for low bead counts. \eqn{\gamma_p} (gain),
#' \eqn{\eta_p} (potency shift) and \eqn{\beta_p} (background) are the sum of a
#' batch and a plate component. Study samples have log2 concentration
#' `mu_a + subject + time effect + noise` and are diluted `dilution`-fold in
#' the well; `truth` is the log2 concentration of the undiluted sample (for
#' standards: of the well).
#'
#' @param n_batches Number of batches (kit lots).
#' @param plates_per_batch Plates per batch.
#' @param analytes Analyte names.
#' @param timepoints Number of time points per subject.
#' @param time_effect Log2 fold change at the last time point (linearly
#'   interpolated for intermediate ones) for every second analyte.
#' @param design `"randomized"` (samples shuffled across plates within batch)
#'   or `"by_time"` (each plate holds mostly one time point -- a confounded
#'   layout that breaks sample-based normalisation).
#' @param unknown_replicates Wells per study sample (1 = singlets).
#' @param std_replicates Wells per standard level.
#' @param n_levels Number of standard levels (serial dilution).
#' @param std_dilution Dilution factor between standard levels.
#' @param n_blanks Blank wells per plate.
#' @param controls Named numeric vector: control materials (run in duplicate on
#'   every plate) and their log2 offset from the typical sample level.
#' @param n_bridge Number of samples from batch 1 re-measured in every other
#'   batch.
#' @param plate_sd,batch_sd Named vectors with SDs (log2) of `gain`, `shift`
#'   and (plate only) `background` effects.
#' @param lot_error_sd SD (log2) of the error in the assigned standard
#'   concentrations of each lot, per analyte.
#' @param noise_cv Multiplicative measurement noise (CV).
#' @param noise_add Additive measurement noise (SD, FI units).
#' @param bad_standard_prob Probability that a standard well is corrupted.
#' @param low_bead_prob Probability of a low bead count in a well x analyte.
#' @param dilution Dilution factor of study samples and controls.
#' @param subject_sd,residual_sd Between-subject and within-subject SDs (log2).
#' @param abundance Range (log2 units, relative to the lowest standard in the
#'   well) from which each analyte's typical sample level is drawn; lower
#'   values give more samples below the standard range.
#' @param seed Optional random seed.
#' @return A [cyto] object with extra columns `subject`, `time` and `truth`;
#'   the generating parameters are in `x$meta$simulation`.
#' @examples
#' x <- simulate_cyto(n_batches = 2, plates_per_batch = 2, seed = 1)
#' x
#' head(x$meta$simulation$analytes)
#' @export
simulate_cyto <- function(n_batches = 2, plates_per_batch = 3,
                          analytes = c("IL-1b", "IL-6", "IL-8", "IL-10",
                                       "IP-10", "MCP-1", "TNF-a", "VEGF"),
                          timepoints = 2, time_effect = 0.7,
                          design = c("randomized", "by_time"),
                          unknown_replicates = 1, std_replicates = 2,
                          n_levels = 8, std_dilution = 4, n_blanks = 2,
                          controls = c(ctrl_low = -1, ctrl_high = 1.5),
                          n_bridge = 8,
                          plate_sd = c(gain = 0.25, shift = 0.25,
                                       background = 0.3),
                          batch_sd = c(gain = 0.3, shift = 0.3),
                          lot_error_sd = 0.2, noise_cv = 0.08, noise_add = 2,
                          bad_standard_prob = 0.03, low_bead_prob = 0.01,
                          dilution = 4, subject_sd = 1, residual_sd = 0.4,
                          abundance = c(-1.5, 6), seed = NULL) {
  design <- match.arg(design)
  dil <- dilution
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv()))
      get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()),
            add = TRUE)
    set.seed(seed)
  }
  na <- length(analytes)

  # ---- analyte curves -------------------------------------------------------
  top_std <- 2^stats::runif(na, log2(2000), log2(30000))     # S1, pg/mL
  low_std <- top_std / std_dilution^(n_levels - 1)
  par_a <- data.table::data.table(
    analyte = analytes,
    bottom = exp(stats::runif(na, log(8), log(60))),
    top = stats::runif(na, 12000, 30000),
    mid = log2(top_std) - stats::runif(na, 4, 7),
    b = stats::runif(na, 0.55, 0.85),
    g = exp(stats::rnorm(na, 0, 0.25)),
    s1 = top_std,
    # biological level relative to the lowest standard (some analytes low)
    mu = log2(low_std) + log2(dil) + stats::runif(na, abundance[1],
                                                  abundance[2]),
    has_effect = rep(c(TRUE, FALSE), length.out = na)
  )

  # ---- layout ---------------------------------------------------------------
  wells96 <- paste0(rep(LETTERS[1:8], 12), rep(1:12, each = 8))
  n_ctrl_w <- 2L * length(controls)
  fixed <- n_levels * std_replicates + n_blanks + n_ctrl_w
  cap <- floor((96 - fixed) / unknown_replicates)
  if (cap < 1) stop("No room for study samples on the plate.", call. = FALSE)
  batches <- as.character(seq_len(n_batches))
  bridge_per_batch <- if (n_batches > 1) n_bridge else 0
  # capacity in batch 1 is all new samples; other batches lose bridge slots
  cap_b <- c(cap * plates_per_batch,
             rep(cap * plates_per_batch - bridge_per_batch, n_batches - 1))
  subj_b <- floor(cap_b / timepoints)

  samples <- data.table::rbindlist(lapply(seq_len(n_batches), function(b) {
    ns <- subj_b[b]
    data.table::data.table(
      batch = batches[b],
      subject = paste0("b", b, "_s", rep(seq_len(ns), each = timepoints)),
      time = rep(seq_len(timepoints), times = ns))
  }))
  samples[, sample_id := paste0(subject, "_t", time)]
  if (bridge_per_batch > 0) {
    br <- samples[batch == "1" & time == 1][seq_len(min(n_bridge, .N))]
    for (b in batches[-1]) {
      samples <- rbind(samples, data.table::copy(br)[, batch := b])
    }
  }

  # assign samples to plates within batch
  lay <- data.table::rbindlist(lapply(batches, function(b) {
    s <- samples[batch == b]
    s <- if (design == "randomized") {
      s[sample(.N)]
    } else {
      s[order(time, stats::runif(.N))]
    }
    plates <- paste0("b", b, "p", seq_len(plates_per_batch))
    s[, plate := rep(plates, each = cap, length.out = .N)]
    if (unknown_replicates > 1) {
      s <- s[rep(seq_len(.N), each = unknown_replicates)]
    }
    s[, sample_type := "unknown"]
    s
  }))

  fixed_rows <- data.table::rbindlist(lapply(batches, function(b) {
    data.table::rbindlist(lapply(paste0("b", b, "p", seq_len(plates_per_batch)),
                                 function(p) {
      rbind(
        data.table::data.table(batch = b, plate = p, sample_type = "standard",
                               sample_id = paste0("S", rep(seq_len(n_levels),
                                                           each = std_replicates)),
                               level = rep(seq_len(n_levels), each = std_replicates)),
        data.table::data.table(batch = b, plate = p, sample_type = "blank",
                               sample_id = "blank", level = NA_real_)[
                                 rep(1, n_blanks)],
        data.table::data.table(batch = b, plate = p, sample_type = "control",
                               sample_id = rep(names(controls), each = 2),
                               level = NA_real_)
      )
    }))
  }))
  layout <- rbind(fixed_rows, lay, fill = TRUE)
  layout[, well := wells96[seq_len(.N)], by = plate]
  layout <- layout[!is.na(well)]

  # ---- effects ----------------------------------------------------------------
  plates_tab <- unique(layout[, .(batch, plate)])
  beff <- data.table::data.table(batch = batches,
                                 bg = stats::rnorm(n_batches, 0, batch_sd[["gain"]]),
                                 bs = stats::rnorm(n_batches, 0, batch_sd[["shift"]]))
  plates_tab <- merge(plates_tab, beff, by = "batch")
  plates_tab[, gain := bg + stats::rnorm(.N, 0, plate_sd[["gain"]])]
  plates_tab[, shift := bs + stats::rnorm(.N, 0, plate_sd[["shift"]])]
  plates_tab[, background := stats::rnorm(.N, 0, plate_sd[["background"]])]
  plates_tab[, c("bg", "bs") := NULL]
  lot_err <- data.table::CJ(batch = batches, analyte = analytes)
  lot_err[, lot_error := stats::rnorm(.N, 0, lot_error_sd)]

  # subject-level biology
  subj <- unique(samples[, .(subject)])
  bio <- data.table::CJ(subject = subj$subject, analyte = analytes)
  bio[, u := stats::rnorm(.N, 0, subject_sd)]
  bio_t <- data.table::CJ(sample_id = unique(samples$sample_id),
                          analyte = analytes)
  bio_t[, e := stats::rnorm(.N, 0, residual_sd)]

  # ---- assemble long table -----------------------------------------------------
  d <- layout[, .(analyte = analytes), by = names(layout)]
  d <- merge(d, par_a, by = "analyte")
  d <- merge(d, plates_tab, by = c("batch", "plate"))
  d <- merge(d, lot_err, by = c("batch", "analyte"))
  d[, sid_tmp := sample_id]
  d <- merge(d, bio, by = c("subject", "analyte"), all.x = TRUE)
  d <- merge(d, bio_t, by.x = c("sid_tmp", "analyte"),
             by.y = c("sample_id", "analyte"), all.x = TRUE)
  d[, sid_tmp := NULL]

  d[, tf := if (timepoints > 1) (time - 1) / (timepoints - 1) else 0]
  d[, dilution := ifelse(sample_type %in% c("unknown", "control"), dil, 1)]
  # nominal (assigned) and true standard concentrations
  d[, expected_conc := ifelse(sample_type == "standard",
                              s1 / std_dilution^(level - 1), NA_real_)]
  d[, truth := NA_real_]
  d[sample_type == "standard", truth := log2(expected_conc) + lot_error]
  d[sample_type == "unknown",
    truth := mu + u + e + has_effect * time_effect * tf]
  d[sample_type == "control", truth := mu + controls[sample_id]]

  # well concentration (log2) and signal
  d[, xw := truth - log2(dilution)]
  d[sample_type == "blank", xw := -Inf]
  d[, fi_true := 2^gain * (bottom * 2^background +
                             (top - bottom) /
                             (1 + exp(-b * (xw + shift - mid)))^g)]
  nr <- nrow(d)
  d[, bead_count := stats::rpois(nr, 80)]
  low <- stats::runif(nr) < low_bead_prob
  d[low, bead_count := stats::rpois(sum(low), 15)]
  bead_noise <- 0.6 / sqrt(pmax(d$bead_count, 1))
  d[, fi := fi_true * exp(stats::rnorm(nr, 0, sqrt(noise_cv^2 + bead_noise^2))) +
      stats::rnorm(nr, 0, noise_add)]
  bad <- d$sample_type == "standard" & stats::runif(nr) < bad_standard_prob
  d[bad, fi := fi * 2^stats::rnorm(sum(bad), 0, 1.5)]
  d[, fi := round(pmax(fi, 1), 1)]

  d[sample_type != "unknown", subject := NA_character_]
  d[sample_type != "unknown", time := NA_integer_]
  out <- d[, .(batch, plate, well, sample_type, sample_id, level, subject,
               time, analyte, fi, bead_count, expected_conc, dilution, truth)]
  data.table::setorder(out, batch, plate, well, analyte)
  x <- cyto(out)
  x$meta$simulation <- list(
    analytes = par_a, plates = plates_tab, lot_error = lot_err,
    call = match.call(), bad_standard_wells = sum(bad))
  x
}
