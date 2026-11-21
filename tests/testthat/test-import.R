test_that("value parser handles Bio-Plex codes", {
  p <- parse_bioplex_value(c("12.5", "*3.1", "OOR <", "OOR >", "***", "", "1,234.5"))
  expect_equal(p$value, c(12.5, 3.1, NA, NA, NA, NA, 1234.5))
  expect_equal(p$flag, c("", "*", "<", ">", "", "", ""))
  expect_equal(parse_bioplex_value("1,5", dec = ",")$value, 1.5)
})

test_that("sample types are inferred from codes", {
  expect_equal(cytokineR:::.infer_type(c("S1", "B", "C2", "X14", "Standard3",
                                         "Background0", "patient 7")),
               c("standard", "blank", "control", "unknown", "standard",
                 "blank", "unknown"))
})

test_that("all example export formats give identical data", {
  f1 <- system.file("extdata", "bioplex_multi_example.csv", package = "cytokineR")
  f2 <- system.file("extdata", "xponent_example.csv", package = "cytokineR")
  f3 <- system.file("extdata", "bioplex_per_analyte_example.xlsx",
                    package = "cytokineR")
  x1 <- read_bioplex(f1, plate = "P1")
  x2 <- read_xponent(f2, plate = "P1")
  key <- function(x) x$data[order(well, analyte), .(well, analyte, sample_type,
                                                    fi)]
  expect_equal(key(x1), key(x2))
  # unknowns keep their names in both formats; controls are named in Bio-Plex
  ids <- function(x) x$data[sample_type == "unknown"][order(well, analyte),
                                                      sample_id]
  expect_equal(ids(x1), ids(x2))
  expect_setequal(unique(x1$data[sample_type == "control", sample_id]),
                  c("ctrl_low", "ctrl_high"))
  expect_setequal(unique(x2$data[sample_type == "control", sample_id]),
                  c("Control1", "Control2"))
  expect_equal(data.table::uniqueN(x1$data$analyte), 8)
  expect_true(any(x1$data$vendor_flag == "<"))
  expect_true(all(is.finite(x1$data[sample_type == "standard", expected_conc])))
  expect_true(all(is.finite(x2$data$bead_count)))
  skip_if_not_installed("readxl")
  x3 <- read_bioplex(f3, plate = "P1")
  expect_equal(key(x1), key(x3))
  expect_equal(ids(x1), ids(x3))
  expect_true(all(x3$data[sample_type == "unknown", dilution] == 4))
})

test_that("several files are combined", {
  f1 <- system.file("extdata", "bioplex_multi_example.csv", package = "cytokineR")
  x <- read_bioplex(c(f1, f1), plate = c("A", "B"))
  expect_setequal(unique(x$data$plate), c("A", "B"))
})

test_that("generic wide/long import and annotation", {
  d <- data.frame(Plate = 1, Well = c("A1", "A2", "A3"),
                  Type = c("S1", "B", "X1"), Name = c("", "", "subj7"),
                  IL6 = c(15000, 12, 300), IL8 = c(12000, 9, 950))
  x <- cyto_from_wide(d, analytes = c("IL6", "IL8"), plate = "Plate",
                      well = "Well", sample_id = "Name", sample_type = NULL,
                      type_from = "Type")
  expect_equal(nrow(x$data), 6)
  expect_setequal(unique(x$data$sample_id), c("blank", "S1", "subj7"))
  expect_equal(unique(x$data[sample_type == "standard", level]), 1)
  std <- standards_series(c(IL6 = 4000, IL8 = 8000), dilution = 4, levels = 8)
  expect_equal(nrow(std), 16)
  x <- set_standards(x, std)
  expect_equal(x$data[sample_type == "standard" & analyte == "IL6",
                      expected_conc], 4000)
  x <- annotate_samples(x, data.frame(sample_id = "subj7", group = "case"))
  expect_equal(x$data[sample_id == "subj7", unique(group)], "case")
  expect_error(annotate_samples(x, data.frame(sample_id = "subj7", fi = 1)),
               "core")
  l <- cyto_data(x)
  y <- cyto_from_long(l)
  expect_equal(nrow(y$data), 6)
})

test_that("vendor concentrations are imported with censoring", {
  f1 <- system.file("extdata", "bioplex_multi_example.csv", package = "cytokineR")
  v <- cyto_transform(read_bioplex(f1, plate = "P1"), "vendor")
  expect_equal(v$scale, "log2_conc")
  expect_true(all(is.na(v$data[censored == "left", value])))
  expect_true(all(is.finite(v$data[censored == "left", limit_low])))
})

test_that("sample types can be reassigned", {
  f2 <- system.file("extdata", "xponent_example.csv", package = "cytokineR")
  x <- read_xponent(f2, plate = "P1")
  u <- x$data[sample_type == "unknown", sample_id][1]
  y <- set_sample_types(x, control = u)
  expect_true(all(y$data[sample_id == u, sample_type] == "control"))
  expect_error(set_sample_types(x, foo = "a"), "Unknown")
})
