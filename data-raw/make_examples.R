# Builds the bundled example data:
#   data/cyto_example.rda                          simulated study (cyto object)
#   inst/extdata/bioplex_multi_example.csv         Bio-Plex Manager multi-analyte layout
#   inst/extdata/xponent_example.csv               Luminex xPONENT batch CSV layout
#   inst/extdata/bioplex_per_analyte/*.csv         per-analyte layout (one file per sheet)
# The per-analyte Excel workbook is assembled from the CSVs by
# data-raw/make_xlsx.py (needs Python + openpyxl).
#
# Run from the package root:  Rscript data-raw/make_examples.R

pkgload::load_all(".")
library(data.table)

cyto_example <- simulate_cyto(n_batches = 2, plates_per_batch = 3, seed = 2026)
cyto_example$meta$simulation$call <- quote(simulate_cyto(n_batches = 2,
                                                          plates_per_batch = 3,
                                                          seed = 2026))
usethis_ok <- dir.exists("data") || dir.create("data")
save(cyto_example, file = "data/cyto_example.rda", compress = "xz")

# ---- one plate as an instrument would export it ---------------------------------
one <- cyto_filter(cyto_example, plate == "b1p1")
vend <- cyto_calibrate(one, "5pl", censor = "asymptote")$data
d <- copy(one$data)
d[, obs := vend$value[match(.row, vend$.row)]]
d[, cens := vend$censored[match(.row, vend$.row)]]
d[, extrap := grepl("extrapolated", vend$qc[match(.row, vend$.row)])]
d[, obs_txt := ifelse(cens == "left", "OOR <", ifelse(cens == "right", "OOR >",
                 ifelse(is.finite(obs), paste0(ifelse(extrap, "*", ""),
                                               formatC(2^obs, format = "f",
                                                       digits = 2)), "***")))]
ctrls <- sort(unique(d[sample_type == "control", sample_id]))
d[, type_code := fcase(sample_type == "standard", paste0("S", level),
                       sample_type == "blank", "B",
                       sample_type == "control", paste0("C", match(sample_id, ctrls)),
                       default = "X")]
unk <- unique(d[sample_type == "unknown", sample_id])
d[sample_type == "unknown", type_code := paste0("X", match(sample_id, unk))]
d[, description := ifelse(sample_type %in% c("unknown", "control"), sample_id, "")]
regions <- setNames(c(39, 19, 54, 56, 48, 51, 36, 45),
                    sort(unique(d$analyte)))
d[, an_label := paste0(analyte, " (", regions[analyte], ")")]
d[, exp_txt := ifelse(sample_type == "standard",
                      formatC(expected_conc, format = "f", digits = 2), "")]
wells <- unique(d[, .(well, type_code, description, sample_type, sample_id)])
an <- sort(unique(d$an_label))

wide_block <- function(col) {
  w <- dcast(d, well + type_code + description ~ an_label, value.var = col)
  w[match(wells$well, w$well)]
}

# Bio-Plex multi-analyte export
lines <- c("Bio-Plex Manager(TM) Report", "File Name: example_plate.rbx",
           "Acquisition Date: 21-Nov-2026, 10:00", "Reader Serial Number: LX00000000",
           "")
add_block <- function(lines, col, dtype) {
  w <- wide_block(col)
  c(lines,
    paste(c("", "", "", an), collapse = ","),
    paste(c("Type", "Well", "Description", rep(dtype, length(an))), collapse = ","),
    apply(w, 1, function(r) paste(c(r[["type_code"]], r[["well"]],
                                    r[["description"]], r[an]), collapse = ",")),
    "")
}
d[, fi_txt := formatC(fi, format = "f", digits = 1)]
d[, bc_txt := as.character(bead_count)]
lines <- add_block(lines, "fi_txt", "FI")
lines <- add_block(lines, "obs_txt", "Obs Conc")
lines <- add_block(lines, "exp_txt", "Exp Conc")
lines <- add_block(lines, "bc_txt", "Bead Count")
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
writeLines(lines, "inst/extdata/bioplex_multi_example.csv")

# xPONENT batch CSV
q <- function(v) paste0("\"", v, "\"")
xn <- sort(unique(d$analyte))
loc <- paste0(seq_len(nrow(wells)), "(1,", wells$well, ")")
xname <- ifelse(wells$sample_type == "standard", sub("^S", "Standard", wells$type_code),
         ifelse(wells$sample_type == "blank", "Background0",
         ifelse(wells$sample_type == "control", sub("^C", "Control", wells$type_code),
                wells$sample_id)))
xsec <- function(col, dtype, extra = NULL) {
  w <- dcast(d, well ~ analyte, value.var = col)
  w <- w[match(wells$well, w$well)]
  rows <- vapply(seq_len(nrow(w)), function(i)
    paste(q(c(loc[i], xname[i], unlist(w[i, xn, with = FALSE]), extra)),
          collapse = ","), "")
  c(paste(q(c("DataType:", dtype)), collapse = ","),
    paste(q(c("Location", "Sample", xn, if (!is.null(extra)) "Total Events")),
          collapse = ","), rows, "")
}
d[, res_txt := ifelse(cens == "none" & is.finite(obs),
                      formatC(2^obs, format = "f", digits = 2), "NaN")]
xl <- c(paste(q(c("Program", "xPONENT")), collapse = ","),
        paste(q(c("Build", "4.3")), collapse = ","),
        paste(q(c("Date", "11/21/2026", "10:00 AM")), collapse = ","), "",
        paste(q(c("Samples", nrow(wells))), collapse = ","), "",
        q("Results"), "",
        xsec("fi_txt", "Median", extra = "1500"),
        xsec("bc_txt", "Count"),
        xsec("res_txt", "Result"))
writeLines(xl, "inst/extdata/xponent_example.csv")

# per-analyte layout: one CSV per analyte (assembled into xlsx by make_xlsx.py)
pa <- "data-raw/per_analyte"
dir.create(pa, showWarnings = FALSE)
blank_fi <- d[sample_type == "blank", .(bg = mean(fi)), by = analyte]
d <- merge(d, blank_fi, by = "analyte")
for (a in sort(unique(d$analyte))) {
  s <- d[analyte == a][match(wells$well, well)]
  out <- c("Bio-Plex Manager(TM) Report", paste0("Analyte: ", s$an_label[1]), "",
           paste(c("Type", "Well", "Outlier", "Description", "FI", "FI - Bkgd",
                   "Obs Conc", "Exp Conc", "Dilution"), collapse = ","),
           paste(s$type_code, s$well, "", s$description,
                 formatC(s$fi, format = "f", digits = 1),
                 formatC(s$fi - s$bg, format = "f", digits = 1),
                 s$obs_txt, s$exp_txt, s$dilution, sep = ","))
  writeLines(out, file.path(pa, paste0(s$an_label[1], ".csv")))
}
message("Done. Now run: python3 data-raw/make_xlsx.py")
