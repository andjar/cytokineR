# ---- low-level helpers ---------------------------------------------------------

#' Parse values from instrument exports
#'
#' Converts the text values found in Bio-Plex Manager / xPONENT exports to
#' numbers and flags: `"OOR <"` and `"OOR >"` (out of range below/above the
#' curve) become `NA` with flag `"<"`/`">"`; a leading asterisk (`"*12.3"`,
#' extrapolated beyond the standards) gives the number with flag `"*"`;
#' `"***"`, `"---"`, `"NaN"`, `""` and similar become `NA`. Decimal commas are
#' handled with `dec = ","`.
#'
#' @param v Character vector.
#' @param dec Decimal separator.
#' @return A data frame with columns `value` (numeric) and `flag` (character).
#' @examples
#' parse_bioplex_value(c("12.5", "*3.1", "OOR <", "OOR >", "***", "", "1,5"))
#' parse_bioplex_value("1,5", dec = ",")
#' @export
parse_bioplex_value <- function(v, dec = ".") {
  s <- trimws(as.character(v))
  s[is.na(s)] <- ""
  flag <- rep("", length(s))
  up <- toupper(gsub("\\s+", "", s))
  lo <- grepl("^OOR<|^<OOR|^<", up)
  hi <- grepl("^OOR>|^>OOR|^>", up)
  flag[lo] <- "<"
  flag[hi] <- ">"
  star <- grepl("^\\*[^*]", s)
  flag[star] <- "*"
  num <- sub("^\\*", "", s)
  if (dec == ",") {
    num <- gsub(".", "", num, fixed = TRUE)
    num <- sub(",", ".", num, fixed = TRUE)
  } else {
    thousands <- grepl("^-?\\d{1,3}(,\\d{3})+(\\.\\d+)?$", num)
    num[thousands] <- gsub(",", "", num[thousands], fixed = TRUE)
  }
  val <- suppressWarnings(as.numeric(num))
  val[lo | hi] <- NA_real_
  val[!is.finite(val)] <- NA_real_
  data.frame(value = val, flag = flag, stringsAsFactors = FALSE)
}

#' Rules for recognising sample types
#'
#' Regular expressions (case-insensitive) applied to the instrument's sample
#' type code or name (`S1`, `B`, `C1`, `X12`, `Standard3`, `Background0`, ...).
#' Codes matching none of the rules are treated as unknowns.
#'
#' @return A named list of regular expressions.
#' @examples
#' default_type_rules()
#' @export
default_type_rules <- function() {
  list(standard = "^(S|STD|Standard)\\s*[0-9]+$",
       blank = "^(B|BLK|Blank|Background|Bkgd)\\s*[0-9]*$",
       control = "^(C|CTRL|Control)\\s*[0-9]*$",
       unknown = "^(X|U|Unknown|Sample|Unk)\\s*[0-9]*$")
}

.infer_type <- function(code, rules = default_type_rules()) {
  code <- trimws(as.character(code))
  out <- rep(NA_character_, length(code))
  for (nm in names(rules)) {
    hit <- is.na(out) & grepl(rules[[nm]], code, ignore.case = TRUE)
    out[hit] <- nm
  }
  out[is.na(out)] <- "unknown"
  out
}

.file_ext <- function(path) tolower(sub(".*\\.([^.]+)$", "\\1", basename(path)))
.file_stem <- function(path) sub("\\.[^.]+$", "", basename(path))

# Read a file into a list of character matrices (one per sheet)
.read_cells <- function(path, sheet = NULL) {
  ext <- .file_ext(path)
  if (ext %in% c("xlsx", "xls", "xlsm")) {
    .require("readxl", "to read Excel files")
    sheets <- readxl::excel_sheets(path)
    if (!is.null(sheet)) {
      sheets <- if (is.numeric(sheet)) sheets[sheet] else intersect(sheets, sheet)
    }
    out <- lapply(sheets, function(s) {
      m <- suppressMessages(readxl::read_excel(path, sheet = s, col_names = FALSE,
                                               col_types = "text",
                                               .name_repair = "minimal"))
      m <- as.matrix(m)
      m[is.na(m)] <- ""
      m
    })
    names(out) <- sheets
    return(out)
  }
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- sub("^\ufeff", "", lines)
  nonempty <- lines[nzchar(trimws(lines))]
  cnt <- vapply(c("\t", ";", ","), function(s)
    sum(lengths(regmatches(nonempty, gregexpr(s, nonempty, fixed = TRUE)))), 0)
  sep <- names(which.max(cnt))
  nf <- suppressWarnings(utils::count.fields(textConnection(lines), sep = sep,
                                             quote = "\"",
                                             blank.lines.skip = FALSE))
  nmax <- max(nf, na.rm = TRUE)
  df <- utils::read.table(text = lines, sep = sep, quote = "\"", header = FALSE,
                          colClasses = "character", fill = TRUE,
                          col.names = paste0("V", seq_len(nmax)),
                          blank.lines.skip = FALSE, comment.char = "",
                          na.strings = character(), strip.white = TRUE)
  m <- as.matrix(df)
  m[is.na(m)] <- ""
  stats::setNames(list(m), .file_stem(path))
}

.clean_analyte <- function(a) {
  a <- trimws(a)
  region <- suppressWarnings(as.integer(sub(".*\\((\\d+)\\)\\s*$", "\\1", a)))
  region[!grepl("\\(\\d+\\)\\s*$", a)] <- NA_integer_
  name <- trimws(sub("\\s*\\(\\d+\\)\\s*$", "", a))
  list(name = name, region = region)
}

.is_blank_row <- function(m, i) all(!nzchar(trimws(m[i, ])))

.header_info <- function(m) {
  first <- which(tolower(trimws(m[, 1])) == "type")[1]
  if (is.na(first) || first <= 1) return(list())
  h <- m[seq_len(first - 1), 1]
  h <- h[grepl(":", h)]
  keys <- trimws(sub(":.*$", "", h))
  vals <- trimws(sub("^[^:]*:", "", h))
  stats::setNames(as.list(vals), keys)
}

.known_fields <- c("type", "well", "outlier", "description", "fi",
                   "fi - bkgd", "std dev", "%cv", "obs conc", "exp conc",
                   "(obs/exp) * 100", "(obs/exp)*100", "group", "ratio",
                   "dilution", "sampling errors", "bead count", "count",
                   "conc in range")

.field_map <- function(h) {
  h <- tolower(trimws(h))
  h <- gsub("\\s+", " ", h)
  out <- rep(NA_character_, length(h))
  out[h == "type"] <- "type"
  out[h == "well"] <- "well"
  out[h == "description"] <- "description"
  out[h %in% c("fi", "mfi", "median", "median fi")] <- "fi"
  out[h %in% c("fi - bkgd", "fi-bkgd", "net mfi")] <- "fi_bkgd"
  out[h %in% c("obs conc", "observed conc", "result", "obs. conc")] <- "obs_conc"
  out[h %in% c("exp conc", "expected conc", "exp. conc")] <- "exp_conc"
  out[h %in% c("dilution", "dilution factor")] <- "dilution"
  out[h %in% c("bead count", "count", "beads")] <- "bead_count"
  out
}

# ---- Bio-Plex Manager -----------------------------------------------------------

.parse_bioplex_multi <- function(m) {
  tr <- which(tolower(trimws(m[, 1])) == "type" &
                tolower(trimws(m[, 2])) == "well")
  if (!length(tr)) return(NULL)
  blocks <- list()
  for (bi in seq_along(tr)) {
    r <- tr[bi]
    stop_at <- if (bi < length(tr)) tr[bi + 1] - 2 else nrow(m)
    rows <- integer()
    i <- r + 1
    while (i <= stop_at && !.is_blank_row(m, i) && nzchar(trimws(m[i, 1]))) {
      rows <- c(rows, i)
      i <- i + 1
    }
    if (!length(rows) || r < 2) next
    cols <- which(nzchar(trimws(m[r - 1, ])))
    cols <- cols[cols >= 4]
    if (!length(cols)) next
    dtype <- trimws(m[r, 4])
    if (!nzchar(dtype)) dtype <- trimws(m[r, cols[1]])
    field <- .field_map(dtype)
    if (is.na(field)) next
    an <- .clean_analyte(m[r - 1, cols])
    blk <- data.table::data.table(
      type = rep(trimws(m[rows, 1]), times = length(cols)),
      well = rep(trimws(m[rows, 2]), times = length(cols)),
      description = rep(trimws(m[rows, 3]), times = length(cols)),
      analyte = rep(an$name, each = length(rows)),
      region = rep(an$region, each = length(rows)),
      v = as.vector(m[rows, cols, drop = FALSE]))
    blk[, field := field]
    blocks[[length(blocks) + 1]] <- blk
  }
  if (!length(blocks)) return(NULL)
  long <- data.table::rbindlist(blocks)
  long <- unique(long, by = c("type", "well", "description", "analyte", "field"))
  data.table::dcast(long, type + well + description + analyte + region ~ field,
                    value.var = "v")
}

.parse_bioplex_per_analyte <- function(m, sheet_name) {
  tr <- which(tolower(trimws(m[, 1])) == "type" &
                tolower(trimws(m[, 2])) == "well")
  if (!length(tr)) return(NULL)
  blocks <- list()
  for (bi in seq_along(tr)) {
    r <- tr[bi]
    hdr <- .field_map(m[r, ])
    # analyte name: an "Analyte: ..." line just above, else the sheet name
    an_txt <- sheet_name
    for (k in rev(seq_len(min(5, r - 1)))) {
      cell <- m[r - k, ]
      hit <- grep("analyte", cell, ignore.case = TRUE, value = TRUE)
      if (length(hit)) {
        an_txt <- trimws(sub("^.*analyte\\s*:?\\s*", "", hit[1],
                             ignore.case = TRUE))
        if (!nzchar(an_txt)) {
          nxt <- cell[which(grepl("analyte", cell, ignore.case = TRUE))[1] + 1]
          if (!is.na(nxt)) an_txt <- trimws(nxt)
        }
      }
    }
    stop_at <- if (bi < length(tr)) tr[bi + 1] - 1 else nrow(m)
    rows <- integer()
    i <- r + 1
    while (i <= stop_at && !.is_blank_row(m, i) && nzchar(trimws(m[i, 1])) &&
           !grepl("analyte", m[i, 1], ignore.case = TRUE)) {
      rows <- c(rows, i)
      i <- i + 1
    }
    if (!length(rows)) next
    an <- .clean_analyte(an_txt)
    blk <- data.table::data.table(analyte = an$name, region = an$region,
                                  .rows = rows)
    for (f in unique(stats::na.omit(hdr))) {
      j <- which(hdr == f)[1]
      data.table::set(blk, j = f, value = trimws(m[rows, j]))
    }
    blk[, .rows := NULL]
    blocks[[length(blocks) + 1]] <- blk
  }
  if (!length(blocks)) return(NULL)
  data.table::rbindlist(blocks, fill = TRUE)
}

.detect_bioplex_format <- function(m) {
  tr <- which(tolower(trimws(m[, 1])) == "type" &
                tolower(trimws(m[, 2])) == "well")
  if (!length(tr)) return(NA_character_)
  r <- tr[1]
  h <- tolower(trimws(m[r, ]))
  h4 <- h[-(1:3)]
  h4 <- h4[nzchar(h4)]
  above <- if (r > 1) trimws(m[r - 1, -(1:3)]) else character()
  # multi-analyte: analyte names in the row above, one data type repeated
  if (length(h4) && data.table::uniqueN(h4) == 1L && sum(nzchar(above)) >= 1) {
    return("multi")
  }
  if ("fi" %in% h && sum(h %in% .known_fields) >= 3) "per_analyte" else "multi"
}

# long parsed table -> cyto-ready data.table
.finish_import <- function(tab, plate, batch, dec, type_rules) {
  for (f in c("fi", "fi_bkgd", "obs_conc", "exp_conc", "dilution",
              "bead_count", "description")) {
    if (!f %in% names(tab)) data.table::set(tab, j = f, value = NA_character_)
  }
  tab <- tab[nzchar(type) & !is.na(type)]
  # drop aggregated replicate rows ("A1,A2") when single-well rows exist
  agg <- grepl(",", tab$well)
  if (any(agg) && any(!agg)) tab <- tab[!agg]
  fi <- parse_bioplex_value(tab$fi, dec)$value
  oc <- parse_bioplex_value(tab$obs_conc, dec)
  ec <- parse_bioplex_value(tab$exp_conc, dec)$value
  dil <- parse_bioplex_value(tab$dilution, dec)$value
  bc <- parse_bioplex_value(tab$bead_count, dec)$value
  st <- .infer_type(tab$type, type_rules)
  lev <- ifelse(st == "standard",
                suppressWarnings(as.numeric(gsub("\\D", "", tab$type))), NA)
  desc <- ifelse(is.na(tab$description), "", tab$description)
  sid <- ifelse(st == "standard", paste0("S", lev),
         ifelse(st == "blank", "blank",
         ifelse(nzchar(desc), desc,
                ifelse(st == "unknown", paste0(plate, "_", tab$type),
                       tab$type))))
  data.table::data.table(
    batch = batch, plate = plate, well = tab$well, sample_type = st,
    sample_id = sid, level = lev, type_code = tab$type, description = desc,
    analyte = tab$analyte, fi = fi,
    bead_count = bc, expected_conc = ifelse(st == "standard", ec, NA_real_),
    dilution = ifelse(is.finite(dil) & dil > 0, dil, 1),
    vendor_conc = oc$value, vendor_flag = oc$flag)
}

#' Read Bio-Plex Manager exports
#'
#' Reads the Excel (or CSV/TXT) files exported from Bio-Rad Bio-Plex Manager
#' into a [cyto] object. Two layouts are recognised automatically:
#'
#' * **multi-analyte** (`"multi"`): one or more blocks, each starting with a row
#'   of analyte names (e.g. `IL-6 (19)`) above a `Type | Well | Description |
#'   <data type>` row, where the data type is e.g. `FI`, `FI - Bkgd`,
#'   `Obs Conc`, `Exp Conc` or `Bead Count`;
#' * **per-analyte** (`"per_analyte"`): one block per analyte (one sheet per
#'   analyte, or several blocks introduced by `Analyte: ...`), with columns such
#'   as `Type`, `Well`, `Description`, `FI`, `FI - Bkgd`, `Obs Conc`,
#'   `Exp Conc`, `Dilution`.
#'
#' Sample types are inferred from the type codes (`S1`..`S8` standards, `B`
#' blank, `C1`.. controls, `X1`.. unknowns; see [default_type_rules()]). The
#' sample identifier is the description when given (so name your wells in
#' Bio-Plex Manager!), otherwise plate + type code. Concentrations reported by
#' the software are kept as `vendor_conc`/`vendor_flag`; the assigned
#' standard concentrations (`Exp Conc`) become `expected_conc`.
#'
#' Several files (one per plate) can be given at once; plate ids default to the
#' file names.
#'
#' @param path File path(s).
#' @param plate Plate id(s); default: file name without extension.
#' @param batch Batch id(s) (recycled over files).
#' @param format `"auto"`, `"multi"` or `"per_analyte"`.
#' @param sheet Optional sheet names/indices (Excel).
#' @param dec Decimal separator in the file.
#' @param type_rules Rules for sample types, see [default_type_rules()].
#' @return A [cyto] object; import details are in `x$meta$import`.
#' @seealso [read_xponent()], [cyto_from_long()], [cyto_from_wide()],
#'   [set_standards()], [annotate_samples()]
#' @examples
#' f <- system.file("extdata", "bioplex_multi_example.csv", package = "cytokineR")
#' x <- read_bioplex(f, plate = "P1")
#' x
#' @export
read_bioplex <- function(path, plate = NULL, batch = "1",
                         format = c("auto", "multi", "per_analyte"),
                         sheet = NULL, dec = ".",
                         type_rules = default_type_rules()) {
  format <- match.arg(format)
  if (length(path) > 1) {
    plate <- plate %||% .file_stem(path)
    batch <- rep_len(batch, length(path))
    objs <- lapply(seq_along(path), function(i)
      read_bioplex(path[i], plate = plate[i], batch = batch[i],
                   format = format, sheet = sheet, dec = dec,
                   type_rules = type_rules))
    return(cyto_bind(objs))
  }
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  plate <- plate %||% .file_stem(path)
  sheets <- .read_cells(path, sheet)
  parsed <- list()
  info <- list()
  regions <- list()
  for (s in names(sheets)) {
    m <- sheets[[s]]
    f <- if (format == "auto") .detect_bioplex_format(m) else format
    if (is.na(f)) next
    tab <- if (f == "multi") .parse_bioplex_multi(m) else
      .parse_bioplex_per_analyte(m, s)
    if (is.null(tab) || !nrow(tab)) next
    regions[[s]] <- unique(tab[, .(analyte, region)])
    parsed[[s]] <- tab
    info[[s]] <- c(list(format = f), .header_info(m))
  }
  if (!length(parsed)) {
    stop("No Bio-Plex data blocks (a row starting 'Type', 'Well') found in ",
         path, ".", call. = FALSE)
  }
  tab <- data.table::rbindlist(parsed, fill = TRUE)
  tab <- unique(tab, by = intersect(c("type", "well", "analyte"), names(tab)))
  long <- .finish_import(tab, plate, batch, dec, type_rules)
  long <- long[is.finite(fi) | is.finite(vendor_conc)]
  x <- cyto(long)
  x$meta$import <- list(file = path, sheets = info,
                        analytes = unique(data.table::rbindlist(regions)))
  .add_step(x, "import", "read_bioplex",
            list(file = basename(path), plate = plate, batch = batch),
            note = paste(unique(vapply(info, `[[`, "", "format")),
                         collapse = ","), scale_before = x$scale,
            n_changed = nrow(x$data))
}

# ---- Luminex xPONENT ---------------------------------------------------------------

#' Read Luminex xPONENT CSV exports
#'
#' Reads the batch CSV written by xPONENT (Luminex 200, FLEXMAP 3D, MAGPIX),
#' which consists of sections introduced by `"DataType:","<name>"`. The
#' `Median` (or `Net MFI`) section gives the FI, `Count` the bead counts,
#' `Result` the software's concentrations and `Dilution Factor` the dilutions;
#' expected concentrations of standards are taken from a section whose name
#' contains `Expected Concentration` if present. Sample types are inferred
#' from the sample names (`Standard1`, `Background0`, `Control1`, `Unknown1`,
#' `S1`, `B`, ...; see [default_type_rules()]).
#'
#' @inheritParams read_bioplex
#' @param fi Which section to use as FI: `"Median"` (default) or `"Net MFI"`.
#' @return A [cyto] object.
#' @examples
#' f <- system.file("extdata", "xponent_example.csv", package = "cytokineR")
#' x <- read_xponent(f)
#' x
#' @export
read_xponent <- function(path, plate = NULL, batch = "1", fi = "Median",
                         type_rules = default_type_rules()) {
  if (length(path) > 1) {
    plate <- plate %||% .file_stem(path)
    batch <- rep_len(batch, length(path))
    return(cyto_bind(lapply(seq_along(path), function(i)
      read_xponent(path[i], plate[i], batch[i], fi, type_rules))))
  }
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
  plate <- plate %||% .file_stem(path)
  m <- .read_cells(path)[[1]]
  dt_rows <- which(trimws(m[, 1]) == "DataType:")
  if (!length(dt_rows)) stop("No 'DataType:' sections found; not an xPONENT ",
                             "CSV?", call. = FALSE)
  sections <- list()
  for (r in dt_rows) {
    name <- trimws(m[r, 2])
    h <- trimws(m[r + 1, ])
    rows <- integer()
    i <- r + 2
    while (i <= nrow(m) && !.is_blank_row(m, i) &&
           trimws(m[i, 1]) != "DataType:") {
      rows <- c(rows, i)
      i <- i + 1
    }
    if (!length(rows)) next
    keep <- which(nzchar(h))
    sec <- data.table::as.data.table(m[rows, keep, drop = FALSE])
    data.table::setnames(sec, make.unique(h[keep]))
    sections[[name]] <- sec
  }
  med <- sections[[fi]]
  if (is.null(med)) stop("Section '", fi, "' not found. Available: ",
                         paste(names(sections), collapse = ", "), call. = FALSE)
  id_cols <- intersect(c("Location", "Sample"), names(med))
  skip <- c(id_cols, "Total Events", "Notes", "Number of Beads")
  analytes <- setdiff(names(med), skip)
  melt_sec <- function(sec, nm) {
    if (is.null(sec)) return(NULL)
    if (!all(c("Location", "Sample") %in% names(sec))) return(NULL)
    a <- intersect(analytes, names(sec))
    l <- data.table::melt(sec, id.vars = c("Location", "Sample"),
                          measure.vars = a, variable.name = "analyte",
                          value.name = nm, variable.factor = FALSE)
    l
  }
  long <- melt_sec(med, "fi")
  for (p in list(c("Count", "bead_count"), c("Result", "obs_conc"),
                 c("Dilution Factor", "dilution"))) {
    s <- melt_sec(sections[[p[1]]], p[2])
    if (!is.null(s)) long <- merge(long, s, by = c("Location", "Sample",
                                                   "analyte"), all.x = TRUE)
  }
  long[, well := sub("^.*,\\s*([A-Pa-p][0-9]+)\\)?.*$", "\\1", Location)]
  long[!grepl("^[A-Pa-p][0-9]+$", well), well := Location]
  ex_name <- grep("Expected Conc", names(sections), value = TRUE)[1]
  if (!is.na(ex_name)) {
    ex <- sections[[ex_name]]
    if (all(c("Location", "Sample") %in% names(ex))) {
      ex <- melt_sec(ex, "exp_conc")
      long <- merge(long, ex, by = c("Location", "Sample", "analyte"),
                    all.x = TRUE)
    } else if ("Sample" %in% names(ex)) {
      a <- intersect(analytes, names(ex))
      ex <- data.table::melt(ex, id.vars = "Sample", measure.vars = a,
                             variable.name = "analyte", value.name = "exp_conc",
                             variable.factor = FALSE)
      long <- merge(long, ex, by = c("Sample", "analyte"), all.x = TRUE)
    }
  }
  code <- sub("\\s*\\(.*$", "", long$Sample)
  long[, type := gsub("\\s+", "", code)]
  long[, description := Sample]
  st <- .infer_type(long$type, type_rules)
  # named samples (not matching any code) are unknowns with their own id
  long[, description := ifelse(st == "unknown" &
                                 !grepl(type_rules$unknown, type,
                                        ignore.case = TRUE), Sample, "")]
  tab <- long[, intersect(c("type", "well", "description", "analyte", "fi",
                            "obs_conc", "exp_conc", "dilution", "bead_count"),
                          names(long)), with = FALSE]
  out <- .finish_import(tab, plate, batch, ".", type_rules)
  x <- cyto(out)
  x$meta$import <- list(file = path, sections = names(sections))
  .add_step(x, "import", "read_xponent",
            list(file = basename(path), plate = plate, batch = batch),
            scale_before = x$scale, n_changed = nrow(x$data))
}

# ---- generic tables ----------------------------------------------------------------

#' Create a cyto object from a generic long or wide table
#'
#' `cyto_from_long()` takes one row per well x analyte; `cyto_from_wide()` one
#' row per well with one column per analyte (FI or concentrations). Column
#' names are mapped with arguments (`plate = "Plate"`, ...). If there is no
#' sample-type column, types are inferred from a type-code column
#' (`type_from`) with [default_type_rules()].
#'
#' @param data A data frame.
#' @param plate,well,sample_id,sample_type,analyte,fi Column names in `data`.
#'   `well` may be `NULL` (row numbers within plate are used);
#'   `sample_type` may be `NULL` (then `type_from` is used).
#' @param type_from Column with type codes (`S1`, `B`, `C1`, `X3`, ...) used to
#'   infer `sample_type` (and `level` for standards) when `sample_type` is
#'   `NULL`.
#' @param rename Named character vector of further renames,
#'   `c(new_name = "old_name")`, e.g. `c(batch = "Batch", expected_conc =
#'   "conc")`.
#' @param type_rules See [default_type_rules()].
#' @return A [cyto] object.
#' @examples
#' d <- data.frame(Plate = 1, Well = c("A1", "A2", "A3"),
#'                 Type = c("S1", "B", "X1"), Name = c("", "", "subj7"),
#'                 IL6 = c(15000, 12, 300), IL8 = c(12000, 9, 950))
#' x <- cyto_from_wide(d, analytes = c("IL6", "IL8"), plate = "Plate",
#'                     well = "Well", sample_id = "Name", sample_type = NULL,
#'                     type_from = "Type")
#' x$data[, .(well, sample_type, sample_id, level, analyte, fi)]
#' @export
cyto_from_long <- function(data, plate = "plate", well = "well",
                           sample_id = "sample_id", sample_type = "sample_type",
                           analyte = "analyte", fi = "fi", type_from = NULL,
                           rename = NULL, type_rules = default_type_rules()) {
  d <- data.table::as.data.table(data)
  d <- data.table::copy(d)
  .map_cols(d, plate, well, sample_id, sample_type, analyte, fi, type_from,
            rename, type_rules)
}

#' @rdname cyto_from_long
#' @param analytes Names of the analyte columns (wide format).
#' @param value What the analyte columns contain: `"fi"` or `"vendor_conc"`.
#' @param oor_codes For `value = "vendor_conc"`: numeric codes used for
#'   out-of-range values, e.g. `c("<" = -1, ">" = -2)`.
#' @export
cyto_from_wide <- function(data, analytes, plate = "plate", well = "well",
                           sample_id = "sample_id", sample_type = "sample_type",
                           value = c("fi", "vendor_conc"), type_from = NULL,
                           rename = NULL, oor_codes = NULL,
                           type_rules = default_type_rules()) {
  value <- match.arg(value)
  d <- data.table::as.data.table(data)
  d <- data.table::copy(d)
  miss <- setdiff(analytes, names(d))
  if (length(miss)) stop("Analyte column(s) not found: ",
                         paste(miss, collapse = ", "), call. = FALSE)
  if (is.null(well)) {
    d[, .well := paste0("r", seq_len(.N)), by = c(plate)]
    well <- ".well"
  }
  for (a in analytes) data.table::set(d, j = a, value = as.character(d[[a]]))
  long <- data.table::melt(d, measure.vars = analytes, variable.name = ".analyte",
                           value.name = ".v", variable.factor = FALSE)
  if (value == "fi") {
    long[, .fi := parse_bioplex_value(.v)$value]
  } else {
    pv <- parse_bioplex_value(long$.v)
    long[, `:=`(vendor_conc = pv$value, vendor_flag = pv$flag, .fi = NA_real_)]
    if (!is.null(oor_codes)) {
      for (fl in names(oor_codes)) {
        hit <- is.finite(long$vendor_conc) & long$vendor_conc == oor_codes[[fl]]
        long[hit, `:=`(vendor_conc = NA_real_, vendor_flag = fl)]
      }
    }
  }
  long[, .v := NULL]
  .map_cols(long, plate, well, sample_id, sample_type, ".analyte", ".fi",
            type_from, rename, type_rules)
}

.map_cols <- function(d, plate, well, sample_id, sample_type, analyte, fi,
                      type_from, rename, type_rules) {
  if (is.null(well)) {
    # rows assumed to be in the same well order for every analyte
    d[, .well := paste0("r", seq_len(.N)), by = c(plate, analyte)]
    well <- ".well"
  }
  map <- c(plate = plate, well = well, analyte = analyte, fi = fi)
  if (!is.null(rename)) map <- c(map, rename)
  if (!is.null(sample_id)) map <- c(map, sample_id = sample_id)
  if (!is.null(sample_type)) map <- c(map, sample_type = sample_type)
  miss <- setdiff(unname(map), names(d))
  if (length(miss)) stop("Column(s) not found: ", paste(miss, collapse = ", "),
                         call. = FALSE)
  out <- data.table::copy(d)
  for (nm in names(map)) {
    data.table::set(out, j = nm, value = d[[map[[nm]]]])
  }
  drop <- setdiff(unname(map), names(map))
  if (length(drop)) out[, (drop) := NULL]
  if (is.null(sample_type)) {
    if (is.null(type_from)) stop("Give either 'sample_type' or 'type_from'.",
                                 call. = FALSE)
    out[, .code := as.character(d[[type_from]])]
    out[, sample_type := .infer_type(.code, type_rules)]
    if (!"level" %in% names(out)) {
      out[, level := ifelse(sample_type == "standard",
                            suppressWarnings(as.numeric(gsub("\\D", "", .code))),
                            NA_real_)]
    }
    if (is.null(sample_id) || !"sample_id" %in% names(out)) {
      out[, sample_id := .code]
    }
    out[, sample_id := as.character(sample_id)]
    out[sample_type == "standard" & !is.na(level), sample_id := paste0("S", level)]
    out[sample_type == "blank", sample_id := "blank"]
    out[sample_type != "standard" & (is.na(sample_id) | !nzchar(sample_id)),
        sample_id := paste0(plate, "_", .code)]
    out[, .code := NULL]
  }
  cyto(out)
}

# ---- annotation ----------------------------------------------------------------------

#' Expected concentrations of standards
#'
#' `standards_series()` builds the table of assigned standard concentrations
#' for a serial dilution (Bio-Plex Pro standards are typically 4-fold
#' dilutions of S1, whose concentration per analyte is lot-specific).
#' `set_standards()` writes such a table into the `expected_conc` column of
#' the standard wells, matched by analyte and level (and batch, if the table
#' has a `batch` column -- use that for different kit lots).
#'
#' @param top Named numeric vector of S1 concentrations per analyte, or a data
#'   frame with columns `analyte`, `top` and optionally `batch`.
#' @param dilution Dilution factor between levels.
#' @param levels Number of levels.
#' @return `standards_series()`: a data frame with `analyte`, `level`,
#'   `expected_conc` (and `batch`).
#' @examples
#' std <- standards_series(c("IL-6" = 5000, "IL-8" = 9000), dilution = 4)
#' head(std)
#' @export
standards_series <- function(top, dilution = 4, levels = 8) {
  if (is.numeric(top)) {
    top <- data.frame(analyte = names(top), top = unname(top))
  }
  top <- data.table::as.data.table(top)
  out <- top[, .(level = seq_len(levels),
                 expected_conc = top / dilution^(seq_len(levels) - 1)),
             by = setdiff(names(top), "top")]
  as.data.frame(out)
}

#' @rdname standards_series
#' @param x A [cyto] object.
#' @param expected Data frame with `analyte`, `level`, `expected_conc` and
#'   optionally `batch`.
#' @export
set_standards <- function(x, expected) {
  .check_cyto(x)
  e <- data.table::as.data.table(expected)
  need <- c("analyte", "level", "expected_conc")
  if (!all(need %in% names(e))) stop("'expected' needs columns ",
                                     paste(need, collapse = ", "), call. = FALSE)
  by <- intersect(c("batch", "analyte", "level"), names(e))
  e[, analyte := as.character(analyte)]
  e[, level := as.numeric(level)]
  if ("batch" %in% by) e[, batch := as.character(batch)]
  x <- .cyto_copy(x)
  d <- x$data
  key <- d[, c(".row", by), with = FALSE]
  m <- merge(key, e[, c(by, "expected_conc"), with = FALSE], by = by,
             all.x = TRUE)
  data.table::setorder(m, .row)
  hit <- d$sample_type == "standard" & is.finite(m$expected_conc)
  d[hit, expected_conc := m$expected_conc[hit]]
  x$data <- d
  .add_step(x, "annotate", "set_standards", list(expected = expected),
            note = sprintf("%d standard rows updated", sum(hit)),
            scale_before = x$scale, n_changed = 0L)
}

#' Add sample metadata
#'
#' Joins well- or sample-level metadata (study, subject, time point, group,
#' ...) onto the measurements. Existing core columns cannot be overwritten.
#'
#' @param x A [cyto] object.
#' @param metadata A data frame.
#' @param by Join column(s), e.g. `"sample_id"` or `c("plate", "well")`.
#' @return A [cyto] object.
#' @examples
#' meta <- data.frame(sample_id = unique(cyto_example$data[sample_type ==
#'                    "unknown", sample_id]))
#' meta$group <- sample(c("case", "control"), nrow(meta), replace = TRUE)
#' x <- annotate_samples(cyto_example, meta, by = "sample_id")
#' table(x$data$group, useNA = "ifany")
#' @export
annotate_samples <- function(x, metadata, by = "sample_id") {
  .check_cyto(x)
  m <- data.table::as.data.table(metadata)
  m <- data.table::copy(m)
  if (!all(by %in% names(m))) stop("'metadata' lacks join column(s).",
                                   call. = FALSE)
  for (b in by) data.table::set(m, j = b, value = as.character(m[[b]]))
  protected <- c(.required_cols, .measure_cols, "batch", "level", "dilution")
  clash <- setdiff(intersect(names(m), protected), by)
  if (length(clash)) stop("Metadata would overwrite core column(s): ",
                          paste(clash, collapse = ", "), call. = FALSE)
  if (anyDuplicated(m, by = by)) stop("Metadata has duplicated keys.",
                                      call. = FALSE)
  x <- .cyto_copy(x)
  d <- x$data
  over <- setdiff(intersect(names(m), names(d)), by)
  if (length(over)) d[, (over) := NULL]
  d <- merge(d, m, by = by, all.x = TRUE, sort = FALSE)
  data.table::setorder(d, .row)
  x$data <- d
  .add_step(x, "annotate", "annotate_samples",
            list(by = by, columns = setdiff(names(m), by)),
            scale_before = x$scale, n_changed = 0L)
}

#' Reassign sample types
#'
#' Sets the sample type of rows identified by `sample_id` (or another column),
#' e.g. when controls were given custom names in the instrument software and
#' were imported as unknowns. Standards get their `level` from the digits in
#' the identifier if it is missing.
#'
#' @param x A [cyto] object.
#' @param ... Named character vectors: `control = c("ctrl_low", "ctrl_high")`,
#'   `blank = "PBS"`, `standard = paste0("Std", 1:8)`, `unknown = ...`.
#' @param by Column matched against the identifiers.
#' @return A [cyto] object.
#' @examples
#' x <- set_sample_types(cyto_example, unknown = "ctrl_high")
#' table(x$data$sample_type)
#' @export
set_sample_types <- function(x, ..., by = "sample_id") {
  .check_cyto(x)
  sets <- list(...)
  bad <- setdiff(names(sets), .sample_types)
  if (length(bad) || is.null(names(sets))) {
    stop("Unknown sample type(s): ", paste(bad, collapse = ", "),
         ". Use: ", paste(.sample_types, collapse = ", "), call. = FALSE)
  }
  .check_by(x$data, by)
  x <- .cyto_copy(x)
  d <- x$data
  n <- 0L
  for (tp in names(sets)) {
    hit <- d[[by]] %in% sets[[tp]]
    n <- n + sum(hit)
    d[hit, sample_type := tp]
    if (tp == "standard") {
      d[hit & is.na(level),
        level := suppressWarnings(as.numeric(gsub("\\D", "", get(by))))]
    }
  }
  x$data <- d
  .add_step(x, "annotate", "set_sample_types", sets,
            note = sprintf("%d rows reassigned", n), scale_before = x$scale,
            n_changed = 0L)
}
