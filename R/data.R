#' Example multiplex cytokine study (simulated)
#'
#' A small simulated study created with
#' `simulate_cyto(n_batches = 2, plates_per_batch = 3, seed = 2026)`: two
#' batches (kit lots) of three 96-well plates, eight analytes, duplicate
#' 8-level standard curves, two blanks and two control materials
#' (`ctrl_low`, `ctrl_high`, in duplicate) per plate, study samples from
#' subjects measured at two time points (singlets), and eight bridging samples
#' from batch 1 re-measured in batch 2. Plate gain/potency/background effects,
#' lot errors in the standard concentrations, bad standard wells and analytes
#' below the lowest standard are included. The true log2 concentration of each
#' sample is in column `truth`; generating parameters are in
#' `cyto_example$meta$simulation`.
#'
#' @format A [cyto] object.
#' @source `data-raw/make_examples.R`
#' @examples
#' cyto_example
#' head(cyto_data(cyto_example))
"cyto_example"
