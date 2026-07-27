#!/usr/bin/env Rscript
#
# run_bulk_normalized.R
# Manifest-driven limma DE driver for bulk_normalized route.
# Accepts: log_normalized, normalized_expression, tpm, fpkm, cpm
# Rejects: raw_integer_counts (use run_bulk_counts.R instead)
# Does NOT use DESeq2 or voom — data is already normalized.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: run_bulk_normalized.R /path/to/run_manifest.yaml", call. = FALSE)
}

manifest_path <- args[[1L]]
if (!file.exists(manifest_path)) stop("Manifest file does not exist: ", manifest_path, call. = FALSE)

required <- c("yaml", "limma", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages({ library(limma) })

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(file.path("skills", "geo-biodata-workflow", "scripts", "drivers", "run_bulk_normalized.R"), mustWork = TRUE)
}
driver_dir <- dirname(script_path)
script_dir <- normalizePath(file.path(driver_dir, ".."), winslash = "/", mustWork = TRUE)
manifest_dir <- dirname(normalizePath(manifest_path, winslash = "/", mustWork = TRUE))

# Source shared limma functions
limma_common <- file.path(script_dir, "bulk_limma_common.R")
if (!file.exists(limma_common)) stop("Could not find bulk_limma_common.R.", call. = FALSE)
source(limma_common)

# Validate manifest before proceeding
validator <- file.path(script_dir, "validate_manifest.R")
if (!file.exists(validator)) stop("Could not find validate_manifest.R.", call. = FALSE)
validation_output <- system2("Rscript", c(shQuote(validator), shQuote(manifest_path)),
  stdout = TRUE, stderr = TRUE)
validation_status <- attr(validation_output, "status") %||% 0L
if (!identical(as.integer(validation_status), 0L)) {
  cat(paste(validation_output, collapse = "\n"), "\n")
  stop("Manifest validation failed; bulk normalized analysis did not run.", call. = FALSE)
}

manifest <- yaml::read_yaml(manifest_path)
if (!identical(manifest$route, "bulk_normalized")) {
  stop("run_bulk_normalized.R only supports route: bulk_normalized", call. = FALSE)
}

allowed_inputs <- c("log_normalized", "normalized_expression", "tpm", "fpkm", "cpm")
input_type <- manifest$input$input_type %||% ""
if (!input_type %in% allowed_inputs) {
  stop("bulk_normalized route requires input_type in: ", paste(allowed_inputs, collapse = ", "),
    ", got: ", input_type, call. = FALSE)
}

resolve_manifest_path <- function(path) {
  if (!nzchar(path %||% "")) return("")
  if (grepl("^[A-Za-z]:[\\\\/]|^/", path)) return(path)
  file.path(manifest_dir, path)
}

input_path <- resolve_manifest_path(manifest$input$file)
sample_path <- resolve_manifest_path(manifest$sample_mapping$file)
sample_id_col <- manifest$sample_mapping$sample_id_column
contrast <- manifest$design$contrast
contrast_factor <- contrast$factor
design_formula <- stats::as.formula(manifest$design$formula)

# Setup output directories
tables_dir <- file.path(manifest_dir, "tables")
figures_dir <- file.path(manifest_dir, "figures")
logs_dir <- file.path(manifest_dir, "logs")
for (path in c(tables_dir, figures_dir, logs_dir)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# Load data
mat <- read_matrix_input(input_path)
sample_map <- read_sample_mapping(sample_path, sample_id_col, contrast_factor)
aligned <- align_samples(mat, sample_map, sample_id_col)
mat <- aligned$matrix
sample_map <- aligned$metadata

# Reject raw-count data: if values look integer-like and non-negative, warn
if (all(mat >= 0) && all(abs(mat - round(mat)) < 1e-8, na.rm = TRUE)) {
  int_ratio <- sum(abs(mat - round(mat)) < 1e-8, na.rm = TRUE) / length(mat)
  if (int_ratio > 0.95) {
    message("NOTE: input matrix appears integer-like (", round(int_ratio * 100),
      "%). If these are raw counts, use route bulk_raw_counts with run_bulk_counts.R instead.")
  }
}

# Write sample mapping
utils::write.table(sample_map, file.path(tables_dir, "sample_mapping_used.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# Low-expression filter
mat_filtered <- filter_low_expression(mat)
message(sprintf("Genes after low-expression filter: %d / %d", nrow(mat_filtered), nrow(mat)))

# Run limma DE
fit_result <- run_limma_de(mat_filtered, sample_map, design_formula, contrast)
result_df <- fit_result$result

# Write outputs
outputs <- write_limma_outputs(result_df, fit_result$ebayes_fit, contrast_factor,
  contrast$numerator, contrast$denominator,
  sample_map, mat_filtered, tables_dir, figures_dir, logs_dir)

# QC checks
qc <- run_limma_qc(fit_result, mat_filtered, sample_map, contrast_factor)
utils::write.table(qc$qc_table, file.path(tables_dir, "qc_checks.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# Determine status
critical_outputs <- c(
  file.path(tables_dir, "sample_mapping_used.tsv"),
  file.path(tables_dir, "design_matrix.tsv"),
  file.path(tables_dir, "library_sizes.tsv"),
  file.path(tables_dir, "pca_coordinates.tsv"),
  file.path(tables_dir, "qc_checks.tsv"),
  outputs["de_path"],
  file.path(figures_dir, "bulk_library_sizes.pdf"),
  file.path(figures_dir, "bulk_pca.pdf"),
  file.path(figures_dir, paste0("bulk_pvalue_histogram_", outputs["contrast_name"], ".pdf")),
  file.path(figures_dir, paste0("bulk_meanvar_", outputs["contrast_name"], ".pdf"))
)

n_total <- nrow(result_df)
n_de <- sum(result_df$adj.P.Val < 0.05, na.rm = TRUE)
status_out <- determine_limma_status(critical_outputs, qc$qc_flags, n_total, n_de)

status <- data.frame(
  state = status_out$state,
  updated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  note = status_out$note,
  stringsAsFactors = FALSE
)
utils::write.table(status, file.path(manifest_dir, "workflow_status.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

session_lines <- utils::capture.output(utils::sessionInfo())
writeLines(session_lines, file.path(logs_dir, "sessionInfo_bulk_normalized.txt"), useBytes = TRUE)

cat(status_out$state, "\n", sep = "")
cat(normalizePath(manifest_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")
