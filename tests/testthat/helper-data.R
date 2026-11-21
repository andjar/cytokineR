small_sim <- function(seed = 11, ...) {
  simulate_cyto(n_batches = 2, plates_per_batch = 2,
                analytes = c("IL-6", "IL-8", "IP-10", "TNF-a"), seed = seed,
                ...)
}

# Plate-shifted toy data with a known offset per plate
toy_shifted <- function() {
  set.seed(1)
  base <- data.frame(
    plate = rep(c("P1", "P2", "P3"), each = 12),
    well = rep(paste0("W", 1:12), 3),
    sample_type = rep(c(rep("control", 2), rep("unknown", 10)), 3),
    sample_id = rep(c("ctrl", "ctrl", paste0("s", 1:10)), 3),
    analyte = "A"
  )
  bio <- stats::setNames(stats::rnorm(10, 10, 1), paste0("s", 1:10))
  lvl <- ifelse(base$sample_type == "control", 9, bio[base$sample_id])
  shift <- c(P1 = 0, P2 = 1, P3 = -0.5)[base$plate]
  base$fi <- 2^(lvl + shift)
  base
}
