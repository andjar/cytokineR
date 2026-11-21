# Simulation study: which pre-processing strategies are robust, and when?
#
# Runs a set of representative pipelines on data simulated under several
# failure modes and scores them with the hold-out benchmark. Reproduce with
#   Rscript system.file("scripts", "simulation_study.R", package = "cytokineR")
# Results depend on the simulator's model of plate effects (gain, potency
# shift, background, lot errors); confirm conclusions on real data with
# benchmark_pipelines() and hold-out control material.

library(cytokineR)
library(data.table)

pipelines <- list(
  "log2 FI (raw)" =
    cyto_pipeline(pipe_transform("log2")),
  "log2 FI > sample median > bridge" =
    cyto_pipeline(pipe_transform("log2"),
                  pipe_normalize("median", within = "time"),
                  pipe_correct_batch("bridge")),
  "MoB > bridge" =
    cyto_pipeline(pipe_transform("mob"), pipe_correct_batch("bridge")),
  "MoS (ref. conc) > bridge" =
    cyto_pipeline(pipe_transform("mos", by_conc = TRUE),
                  pipe_correct_batch("bridge")),
  "MoC (control) > bridge" =
    cyto_pipeline(pipe_transform("moc", control = "ctrl_low"),
                  pipe_correct_batch("bridge")),
  "log2 FI > std anchor (affine) > bridge" =
    cyto_pipeline(pipe_transform("log2"),
                  pipe_normalize("anchor", anchors = "standard",
                                 scope = "batch", fit = "affine"),
                  pipe_correct_batch("bridge")),
  "5PL/plate > Tobit" =
    cyto_pipeline(pipe_calibrate("5pl"), pipe_impute("tobit")),
  "5PL/plate > Tobit > bridge" =
    cyto_pipeline(pipe_calibrate("5pl"), pipe_impute("tobit"),
                  pipe_correct_batch("bridge")),
  "spline/plate > Tobit > bridge" =
    cyto_pipeline(pipe_calibrate("spline"), pipe_impute("tobit"),
                  pipe_correct_batch("bridge")),
  "shared 5PL > rank-Tobit > bridge" =
    cyto_pipeline(pipe_calibrate("shared"), pipe_impute("rank_tobit"),
                  pipe_correct_batch("bridge")),
  "harmonise > one curve (ref. lot) > rank-Tobit" =
    cyto_pipeline(pipe_transform("log2"),
                  pipe_normalize("anchor", anchors = "standard",
                                 scope = "batch", fit = "affine"),
                  pipe_correct_batch("bridge"),
                  pipe_calibrate("5pl", curve_by = NULL,
                                 standards = ~ batch == "1"),
                  pipe_impute("rank_tobit"))
)

scenarios <- list(
  default      = list(),
  strong_plate = list(plate_sd = c(gain = 0.5, shift = 0.5, background = 0.5)),
  confounded   = list(design = "by_time"),
  bad_standards = list(bad_standard_prob = 0.15),
  low_abundance = list(abundance = c(-4, 3)),
  noisy_blanks = list(n_blanks = 1, noise_add = 6)
)

reps <- as.integer(Sys.getenv("CYTOKINER_REPS", "5"))
res <- benchmark_simulations(pipelines, scenarios, reps = reps,
                             holdout = "ctrl_high",
                             metrics = default_metrics(bio = ~ time,
                                                       adjust = ~ subject),
                             seed = 2026)
print(res)
out <- Sys.getenv("CYTOKINER_OUT", "")
if (nzchar(out)) saveRDS(res, out)
