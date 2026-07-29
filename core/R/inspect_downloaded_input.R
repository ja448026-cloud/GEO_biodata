#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(1L, 2L)) {
  stop("Usage: inspect_downloaded_input.R <run-dir> [input-file]", call. = FALSE)
}

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Missing required package: yaml", call. = FALSE)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
truthy <- function(x) tolower(trimws(as.character(x))) %in% c("true", "t", "yes", "y", "1")

run_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
raw_dir <- file.path(run_dir, "raw")
tables_dir <- file.path(run_dir, "tables")
resources_dir <- file.path(run_dir, "resources")
for (path in c(tables_dir, resources_dir)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

input_files <- if (length(args) == 2L) {
  normalizePath(args[[2L]], winslash = "/", mustWork = TRUE)
} else if (dir.exists(raw_dir)) {
  list.files(raw_dir, full.names = TRUE, recursive = FALSE)
} else {
  character()
}
input_files <- input_files[file.info(input_files)$isdir == FALSE]

detect_format <- function(path) {
  lname <- tolower(basename(path))
  if (grepl("\\.(tsv|txt)(\\.gz)?$", lname)) return("delimited_text")
  if (grepl("\\.csv(\\.gz)?$", lname)) return("csv")
  if (grepl("\\.rds$|\\.rdata$|\\.rda$", lname)) return("r_object")
  if (grepl("\\.h5ad$", lname)) return("h5ad")
  if (grepl("\\.h5$", lname)) return("h5")
  if (grepl("\\.mtx(\\.gz)?$", lname)) return("mtx")
  "unknown"
}

read_head_table <- function(path, format) {
  if (!format %in% c("delimited_text", "csv")) return(NULL)
  sep <- if (identical(format, "csv")) "," else "\t"
  tryCatch(utils::read.table(path, sep = sep, header = TRUE, quote = "\"",
    comment.char = "", check.names = FALSE, stringsAsFactors = FALSE,
    nrows = 2000), error = function(e) NULL)
}

classify_matrix <- function(path) {
  format <- detect_format(path)
  tab <- read_head_table(path, format)
  if (is.null(tab) || ncol(tab) < 2L) {
    return(list(format = format, n_rows_checked = NA_integer_, n_cols = NA_integer_,
      feature_id_column = "", duplicate_feature_ids = NA_integer_,
      numeric_fraction = NA_real_, finite_fraction = NA_real_,
      nonnegative_fraction = NA_real_, integer_like_fraction = NA_real_,
      zero_fraction = NA_real_, min_value = NA_real_, max_value = NA_real_,
      route_candidate = if (format %in% c("r_object", "h5ad", "h5", "mtx")) "scrna_author_object" else "unknown",
      scale_candidate = "unknown", sample_columns = character()))
  }
  feature_ids <- as.character(tab[[1L]])
  values <- suppressWarnings(as.matrix(tab[, -1L, drop = FALSE]))
  storage.mode(values) <- "numeric"
  vals <- as.vector(values)
  finite <- is.finite(vals)
  finite_vals <- vals[finite]
  numeric_fraction <- mean(!is.na(vals))
  finite_fraction <- mean(finite)
  nonnegative_fraction <- if (length(finite_vals)) mean(finite_vals >= 0) else NA_real_
  integer_like_fraction <- if (length(finite_vals)) mean(abs(finite_vals - round(finite_vals)) < 1e-8) else NA_real_
  zero_fraction <- if (length(finite_vals)) mean(abs(finite_vals) < 1e-12) else NA_real_
  max_value <- if (length(finite_vals)) max(finite_vals) else NA_real_
  min_value <- if (length(finite_vals)) min(finite_vals) else NA_real_

  route_candidate <- "bulk_normalized"
  scale_candidate <- "normalized_expression"
  if (isTRUE(integer_like_fraction > 0.95) && isTRUE(nonnegative_fraction > 0.99)) {
    route_candidate <- "bulk_raw_counts"
    scale_candidate <- "raw_integer_counts"
  } else if (isTRUE(nonnegative_fraction > 0.99) && isTRUE(max_value > 100)) {
    scale_candidate <- "tpm_fpkm_cpm_or_unlogged_normalized"
  } else if (isTRUE(min_value >= 0) && isTRUE(max_value < 30)) {
    scale_candidate <- "log_or_transformed_expression"
  }

  list(format = format, n_rows_checked = nrow(tab), n_cols = ncol(tab),
    feature_id_column = names(tab)[[1L]],
    duplicate_feature_ids = sum(duplicated(feature_ids) | !nzchar(feature_ids)),
    numeric_fraction = numeric_fraction, finite_fraction = finite_fraction,
    nonnegative_fraction = nonnegative_fraction,
    integer_like_fraction = integer_like_fraction, zero_fraction = zero_fraction,
    min_value = min_value, max_value = max_value,
    route_candidate = route_candidate, scale_candidate = scale_candidate,
    sample_columns = names(tab)[-1L])
}

file_inventory <- data.frame(
  file = input_files,
  file_name = basename(input_files),
  size_bytes = if (length(input_files)) file.info(input_files)$size else numeric(),
  format_candidate = vapply(input_files, detect_format, character(1)),
  stringsAsFactors = FALSE
)
utils::write.table(file_inventory, file.path(tables_dir, "file_inventory.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

audits <- lapply(input_files, classify_matrix)
matrix_audit <- do.call(rbind, Map(function(path, audit) {
  data.frame(
    file = path,
    file_name = basename(path),
    format_candidate = audit$format,
    n_rows_checked = audit$n_rows_checked,
    n_cols = audit$n_cols,
    feature_id_column = audit$feature_id_column,
    duplicate_or_empty_feature_ids = audit$duplicate_feature_ids,
    numeric_fraction = audit$numeric_fraction,
    finite_fraction = audit$finite_fraction,
    nonnegative_fraction = audit$nonnegative_fraction,
    integer_like_fraction = audit$integer_like_fraction,
    zero_fraction = audit$zero_fraction,
    min_value = audit$min_value,
    max_value = audit$max_value,
    route_candidate = audit$route_candidate,
    scale_candidate = audit$scale_candidate,
    stringsAsFactors = FALSE
  )
}, input_files, audits))
if (is.null(matrix_audit)) {
  matrix_audit <- data.frame(file = character(), file_name = character())
}
utils::write.table(matrix_audit, file.path(tables_dir, "matrix_intake_audit.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

sample_ids <- character()
sample_index_path <- file.path(resources_dir, "sample_index.tsv")
if (file.exists(sample_index_path)) {
  sample_index <- utils::read.delim(sample_index_path, stringsAsFactors = FALSE, check.names = FALSE)
  sample_col <- intersect(c("sample_id", "geo_accession", "gsm", "gsm_accession"), names(sample_index))[1L]
  if (!is.na(sample_col)) sample_ids <- as.character(sample_index[[sample_col]])
}

overlap_rows <- do.call(rbind, Map(function(path, audit) {
  cols <- audit$sample_columns
  data.frame(
    file = path,
    n_matrix_sample_columns = length(cols),
    n_geo_samples = length(sample_ids),
    n_overlap = length(intersect(cols, sample_ids)),
    overlap_fraction = if (length(cols) > 0L) length(intersect(cols, sample_ids)) / length(cols) else NA_real_,
    stringsAsFactors = FALSE
  )
}, input_files, audits))
if (is.null(overlap_rows)) {
  overlap_rows <- data.frame(file = character(), n_matrix_sample_columns = integer(),
    n_geo_samples = integer(), n_overlap = integer(), overlap_fraction = numeric())
}
utils::write.table(overlap_rows, file.path(tables_dir, "sample_overlap.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

best_idx <- if (nrow(matrix_audit) > 0L) {
  valid <- which(matrix_audit$route_candidate != "unknown")
  if (length(valid)) valid[[1L]] else 1L
} else NA_integer_
selected <- if (!is.na(best_idx)) matrix_audit[best_idx, , drop = FALSE] else NULL
state <- if (is.null(selected)) "BLOCKED_INPUT" else "REVIEW_REQUIRED"
next_module <- if (is.null(selected)) "none" else if (selected$route_candidate %in% c("bulk_raw_counts", "bulk_normalized")) {
  "geo-biodata-bulk"
} else if (selected$route_candidate == "scrna_author_object") {
  "geo-biodata-scrna"
} else {
  "review_required"
}

handoff <- list(
  schema_version = "1.0",
  accession = "",
  run_dir = run_dir,
  intake_state = state,
  recommended_next_module = next_module,
  selected_input = list(
    path = if (is.null(selected)) "" else selected$file[[1L]],
    format = if (is.null(selected)) "" else selected$format_candidate[[1L]],
    route = if (is.null(selected)) "" else selected$route_candidate[[1L]],
    input_type = if (is.null(selected)) "" else selected$scale_candidate[[1L]],
    species = ""
  ),
  scale = list(status = "unconfirmed", transformed = NA, transform = "",
    evidence_source = "file_intake_audit", evidence_note = "Audit suggests candidates; human review is required."),
  sample_mapping = list(status = "unreviewed", file = "", biological_unit = "",
    sample_id_column = "", group_column = ""),
  design_readiness = list(status = "not_assessed", suggested_formula = "",
    suggested_contrast = list(factor = "", numerator = "", denominator = "")),
  unresolved_questions = c("Confirm input type, scale, sample mapping, design, and contrast before analysis."),
  warnings = character(),
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
yaml::write_yaml(handoff, file.path(run_dir, "intake_handoff.yaml"))

cat(state, "\n", sep = "")
cat(file.path(run_dir, "intake_handoff.yaml"), "\n", sep = "")
