#!/usr/bin/env Rscript
#
# run_bulk_normalized.R
# Manifest-driven driver for bulk_normalized route.
# Two independent code paths: EDA-only (no DE) and Differential Expression (limma).
# Enforces input scale contract.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: run_bulk_normalized.R /path/to/run_manifest.yaml", call. = FALSE)
}

manifest_path <- args[[1L]]
if (!file.exists(manifest_path)) stop("Manifest file does not exist: ", manifest_path, call. = FALSE)

required <- c("yaml", "ggplot2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

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

# Source shared functions
limma_common <- file.path(script_dir, "bulk_limma_common.R")
if (!file.exists(limma_common)) stop("Could not find bulk_limma_common.R.", call. = FALSE)
source(limma_common)

# Validate manifest
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
group_col <- manifest$sample_mapping$group_column

# Scale contract
scale_config <- manifest$input$scale %||% list(transformed = FALSE, transform = "", evidence_source = "")
analysis_intent <- manifest$analysis$intent %||% "differential_expression"

tables_dir <- file.path(manifest_dir, "tables")
figures_dir <- file.path(manifest_dir, "figures")
logs_dir <- file.path(manifest_dir, "logs")
for (path in c(tables_dir, figures_dir, logs_dir)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

fallback_events <- init_fallback_events()

# Load data
mat <- read_matrix_input(input_path)
sample_map <- read_sample_mapping(sample_path, sample_id_col, group_col)
aligned <- align_samples(mat, sample_map, sample_id_col, contrast_factor = group_col)
mat <- aligned$matrix
sample_map <- aligned$metadata

utils::write.table(sample_map, file.path(tables_dir, "sample_mapping_used.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# ── Scale contract validation ────────────────────────────────────────────────

scale_check <- validate_scale_contract(input_type, scale_config, mat)

# ── EDA-only path ────────────────────────────────────────────────────────────

if (identical(analysis_intent, "eda_only")) {

  message("EDA-only mode: skipping differential expression.")
  mat_filtered <- filter_low_expression(mat)

  # Expression summary
  expr_summary <- data.frame(
    sample_id = colnames(mat_filtered),
    mean = colMeans(mat_filtered, na.rm = TRUE),
    median = apply(mat_filtered, 2L, stats::median, na.rm = TRUE),
    sd = apply(mat_filtered, 2L, stats::sd, na.rm = TRUE),
    min = apply(mat_filtered, 2L, min, na.rm = TRUE),
    max = apply(mat_filtered, 2L, max, na.rm = TRUE),
    group = as.character(sample_map[[group_col]]),
    stringsAsFactors = FALSE
  )
  utils::write.table(expr_summary, file.path(tables_dir, "expression_summary.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  # Expression distribution plot
  grDevices::pdf(file.path(figures_dir, "expression_distributions.pdf"), width = 7, height = 5)
  boxplot(mat_filtered, las = 2, cex.axis = 0.7, main = "Expression distributions", ylab = "Expression")
  grDevices::dev.off()

  # PCA
  pca <- stats::prcomp(t(mat_filtered), scale. = TRUE)
  pca_df <- data.frame(
    sample_id = rownames(pca$x),
    PC1 = pca$x[, 1L],
    PC2 = if (ncol(pca$x) >= 2L) pca$x[, 2L] else 0,
    group = as.character(sample_map[rownames(pca$x), group_col]),
    stringsAsFactors = FALSE
  )
  utils::write.table(pca_df, file.path(tables_dir, "pca_coordinates.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  pca_plot <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = group, label = sample_id)) +
    ggplot2::geom_point(size = 3) + ggplot2::geom_text(vjust = -0.8, show.legend = FALSE) +
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(figures_dir, "bulk_pca.pdf"), pca_plot, width = 5.5, height = 4.5)

  # Sample correlation
  cor_mat <- stats::cor(mat_filtered, use = "pairwise.complete.obs")
  utils::write.table(
    data.frame(sample = rownames(cor_mat), cor_mat, check.names = FALSE),
    file.path(tables_dir, "sample_correlation.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  # Correlation heatmap
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    ann <- data.frame(group = sample_map[[group_col]], row.names = rownames(sample_map))
    grDevices::pdf(file.path(figures_dir, "sample_correlation_heatmap.pdf"), width = 7, height = 6)
    pheatmap::pheatmap(cor_mat, annotation_col = ann, main = "Sample correlation")
    grDevices::dev.off()
  }

  critical_outputs <- c(
    file.path(tables_dir, "sample_mapping_used.tsv"),
    file.path(tables_dir, "expression_summary.tsv"),
    file.path(tables_dir, "pca_coordinates.tsv"),
    file.path(tables_dir, "sample_correlation.tsv"),
    file.path(figures_dir, "expression_distributions.pdf"),
    file.path(figures_dir, "bulk_pca.pdf")
  )
  outputs_exist <- all(file.exists(critical_outputs) & file.info(critical_outputs)$size > 0)

  status <- data.frame(
    execution_state = if (outputs_exist) "EXECUTION_COMPLETE" else "EXECUTION_FAILED",
    contract_state = if (!scale_check$allow_de) "NOT_ASSESSED" else "VALID",
    technical_qc = "NOT_ASSESSED",
    result_signal = "NOT_ASSESSED",
    analysis_intent = "EDA_ONLY",
    updated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    note = sprintf("EDA-only completed. %d samples, %d features.", ncol(mat_filtered), nrow(mat_filtered)),
    stringsAsFactors = FALSE
  )
  utils::write.table(status, file.path(manifest_dir, "workflow_status.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  session_lines <- utils::capture.output(utils::sessionInfo())
  writeLines(session_lines, file.path(logs_dir, "sessionInfo_bulk_normalized_eda.txt"), useBytes = TRUE)

  cat("EXECUTION_COMPLETE\n", normalizePath(manifest_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")
  quit(status = 0L)
}

# ── Differential Expression path ─────────────────────────────────────────────

# Enforce scale contract for DE
if (identical(analysis_intent, "differential_expression")) {
  if (!scale_check$allow_de) {
    status <- data.frame(
      execution_state = "EXECUTION_COMPLETE",
      contract_state = "BLOCKED",
      technical_qc = "NOT_ASSESSED",
      result_signal = "NOT_ASSESSED",
      analysis_intent = "DIFFERENTIAL_EXPRESSION",
      updated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      note = paste("DE blocked by scale contract:", paste(scale_check$conditions, collapse = "; ")),
      stringsAsFactors = FALSE
    )
    utils::write.table(status, file.path(manifest_dir, "workflow_status.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    writeLines(c("Scale contract blocks DE:", paste("-", scale_check$conditions)),
      file.path(logs_dir, "scale_blocked.txt"), useBytes = TRUE)
    cat("BLOCKED\n", normalizePath(manifest_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")
    quit(status = 0L)
  }
}

# -- DE requires limma
if (!requireNamespace("limma", quietly = TRUE)) {
  stop("limma package required for DE. Install with: BiocManager::install('limma')", call. = FALSE)
}
suppressPackageStartupMessages({ library(limma) })

contrast <- manifest$design$contrast
contrast_factor <- contrast$factor
design_formula <- stats::as.formula(manifest$design$formula)

mat_filtered <- filter_low_expression(mat)
message(sprintf("Genes after low-expression filter: %d / %d", nrow(mat_filtered), nrow(mat)))

fit_result <- run_limma_de(mat_filtered, sample_map, design_formula, contrast)
result_df <- fit_result$result

outputs <- write_limma_outputs(result_df, fit_result$ebayes_fit,
  fit_result$design, fit_result$contrast_matrix,
  contrast_factor, contrast$numerator, contrast$denominator,
  sample_map, mat_filtered,
  fit_result$factor_levels_before, fit_result$factor_levels_after, fit_result$factor_reference,
  tables_dir, figures_dir, logs_dir)

qc <- run_limma_qc(fit_result, mat_filtered, sample_map, contrast_factor)
utils::write.table(qc$qc_table, file.path(tables_dir, "qc_checks.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

if (nrow(fallback_events) > 0L) {
  utils::write.table(fallback_events, file.path(tables_dir, "fallback_events.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

critical_outputs <- c(
  file.path(tables_dir, "sample_mapping_used.tsv"),
  file.path(tables_dir, "design_matrix_used.tsv"),
  file.path(tables_dir, "contrast_matrix_used.tsv"),
  file.path(tables_dir, "factor_levels_before_relevel.tsv"),
  file.path(tables_dir, "factor_levels_used.tsv"),
  file.path(tables_dir, "library_sizes.tsv"),
  file.path(tables_dir, "pca_coordinates.tsv"),
  file.path(tables_dir, "qc_checks.tsv"),
  outputs["de_path"],
  file.path(figures_dir, "bulk_library_sizes.pdf"),
  file.path(figures_dir, "bulk_pca.pdf"),
  file.path(figures_dir, paste0("bulk_pvalue_histogram_", outputs["contrast_name"], ".pdf")),
  file.path(figures_dir, paste0("bulk_meanvar_", outputs["contrast_name"], ".pdf"))
)

status_out <- determine_limma_status(critical_outputs, qc, fallback_events)

status <- data.frame(
  execution_state = status_out$execution_state,
  contract_state = "VALID",
  technical_qc = status_out$technical_qc,
  result_signal = status_out$result_signal,
  analysis_intent = "DIFFERENTIAL_EXPRESSION",
  updated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  note = status_out$note,
  stringsAsFactors = FALSE
)
utils::write.table(status, file.path(manifest_dir, "workflow_status.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

session_lines <- utils::capture.output(utils::sessionInfo())
writeLines(session_lines, file.path(logs_dir, "sessionInfo_bulk_normalized_de.txt"), useBytes = TRUE)

cat(status_out$execution_state, "\n", sep = "")
cat(normalizePath(manifest_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")
