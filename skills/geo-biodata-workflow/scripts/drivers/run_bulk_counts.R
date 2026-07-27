#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: run_bulk_counts.R /path/to/run_manifest.yaml", call. = FALSE)
}

manifest_path <- args[[1L]]
if (!file.exists(manifest_path)) stop("Manifest file does not exist: ", manifest_path, call. = FALSE)
required <- c("yaml", "DESeq2", "ggplot2", "SummarizedExperiment")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  stop("Missing required packages for bulk raw-count driver: ", paste(missing, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
find_up <- function(starts, relative_path) {
  for (start in starts) {
    current <- normalizePath(start, winslash = "/", mustWork = FALSE)
    for (i in seq_len(8L)) {
      candidate <- file.path(current, relative_path)
      if (file.exists(candidate)) return(candidate)
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  NA_character_
}
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE) else normalizePath(file.path("skills", "geo-biodata-workflow", "scripts", "drivers", "run_bulk_counts.R"), mustWork = TRUE)
driver_dir <- dirname(script_path)
script_dir <- normalizePath(file.path(driver_dir, ".."), winslash = "/", mustWork = TRUE)
manifest_dir <- dirname(normalizePath(manifest_path, winslash = "/", mustWork = TRUE))
validator <- file.path(script_dir, "validate_manifest.R")
if (!file.exists(validator)) stop("Could not find validate_manifest.R beside the driver.", call. = FALSE)

validation_output <- system2("Rscript", c(shQuote(validator), shQuote(manifest_path)), stdout = TRUE, stderr = TRUE)
validation_status <- attr(validation_output, "status") %||% 0L
if (!identical(as.integer(validation_status), 0L)) {
  cat(paste(validation_output, collapse = "\n"), "\n")
  stop("Manifest validation failed; bulk analysis did not run.", call. = FALSE)
}

manifest <- yaml::read_yaml(manifest_path)
if (!identical(manifest$route, "bulk_raw_counts")) {
  stop("run_bulk_counts.R only supports route: bulk_raw_counts", call. = FALSE)
}
if (!identical(manifest$input$input_type, "raw_integer_counts")) {
  stop("run_bulk_counts.R requires input.input_type: raw_integer_counts", call. = FALSE)
}

resolve_manifest_path <- function(path) {
  if (!nzchar(path %||% "")) return("")
  if (grepl("^[A-Za-z]:[\\\\/]|^/", path)) return(path)
  file.path(manifest_dir, path)
}
input_path <- resolve_manifest_path(manifest$input$file)
sample_path <- resolve_manifest_path(manifest$sample_mapping$file)

tables_dir <- file.path(manifest_dir, "tables")
figures_dir <- file.path(manifest_dir, "figures")
logs_dir <- file.path(manifest_dir, "logs")
for (path in c(tables_dir, figures_dir, logs_dir)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

counts_df <- utils::read.delim(input_path, stringsAsFactors = FALSE, check.names = FALSE)
feature_ids <- as.character(counts_df[[1L]])
if (any(!nzchar(feature_ids)) || any(duplicated(feature_ids))) {
  stop("Feature IDs must be non-empty and unique.", call. = FALSE)
}
count_mat <- as.matrix(counts_df[, -1L, drop = FALSE])
suppressWarnings(storage.mode(count_mat) <- "integer")
rownames(count_mat) <- feature_ids

sample_map <- utils::read.delim(sample_path, stringsAsFactors = FALSE, check.names = FALSE)
sample_id_col <- manifest$sample_mapping$sample_id_column
sample_ids <- as.character(sample_map[[sample_id_col]])
sample_map <- sample_map[match(colnames(count_mat), sample_ids), , drop = FALSE]
rownames(sample_map) <- as.character(sample_map[[sample_id_col]])
utils::write.table(sample_map, file.path(tables_dir, "sample_mapping_used.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

design_formula <- stats::as.formula(manifest$design$formula)
design_matrix <- stats::model.matrix(design_formula, data = sample_map)
utils::write.table(
  data.frame(sample_id = rownames(design_matrix), design_matrix, check.names = FALSE),
  file.path(tables_dir, "design_matrix.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

contrast <- manifest$design$contrast
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(count_mat),
  colData = sample_map,
  design = design_formula
)
dds <- dds[rowSums(DESeq2::counts(dds)) > 0L, ]
dds <- DESeq2::estimateSizeFactors(dds)

# ── Dispersion estimation with fallback ──
dispersion_fallback <- FALSE
dispersion_fallback_reason <- ""
dds <- tryCatch(
  DESeq2::estimateDispersions(dds, quiet = TRUE),
  error = function(e) {
    dispersion_fallback <<- TRUE
    dispersion_fallback_reason <<- conditionMessage(e)
    message("Dispersion trend fit failed; using gene-wise dispersion estimates. Reason: ", dispersion_fallback_reason)
    fallback <- DESeq2::estimateDispersionsGeneEst(dds, quiet = TRUE)
    dispersions(fallback) <- S4Vectors::mcols(fallback)$dispGeneEst
    fallback
  }
)
dds <- DESeq2::nbinomWaldTest(dds, quiet = TRUE)
result <- DESeq2::results(
  dds,
  contrast = c(contrast$factor, contrast$numerator, contrast$denominator)
)
result_df <- as.data.frame(result)
result_df$feature_id <- rownames(result_df)
result_df <- result_df[, c("feature_id", setdiff(names(result_df), "feature_id")), drop = FALSE]
contrast_name <- paste(contrast$factor, contrast$numerator, "vs", contrast$denominator, sep = "_")
de_path <- file.path(tables_dir, paste0("de_results_", contrast_name, ".tsv"))
utils::write.table(result_df, de_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# ── Library sizes ──
lib_sizes <- data.frame(
  sample_id = colnames(count_mat),
  library_size = colSums(count_mat),
  group = as.character(sample_map[[contrast$factor]]),
  stringsAsFactors = FALSE
)
utils::write.table(lib_sizes, file.path(tables_dir, "library_sizes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
lib_plot <- ggplot2::ggplot(lib_sizes, ggplot2::aes(x = sample_id, y = library_size, fill = group)) +
  ggplot2::geom_col() +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
ggplot2::ggsave(file.path(figures_dir, "bulk_library_sizes.pdf"), lib_plot, width = 7, height = 4)

# ── PCA ──
transformed <- if (dispersion_fallback) {
  DESeq2::normTransform(dds)
} else {
  DESeq2::varianceStabilizingTransformation(dds, blind = FALSE)
}
pca <- stats::prcomp(t(SummarizedExperiment::assay(transformed)), scale. = FALSE)
pca_df <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1L],
  PC2 = if (ncol(pca$x) >= 2L) pca$x[, 2L] else 0,
  group = as.character(sample_map[rownames(pca$x), contrast$factor]),
  stringsAsFactors = FALSE
)
utils::write.table(pca_df, file.path(tables_dir, "pca_coordinates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")
pca_plot <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = group, label = sample_id)) +
  ggplot2::geom_point(size = 3) +
  ggplot2::geom_text(vjust = -0.8, show.legend = FALSE) +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(figures_dir, "bulk_pca.pdf"), pca_plot, width = 5.5, height = 4.5)

# ── P-value histogram ──
pval_df <- data.frame(pvalue = result_df$pvalue)
pval_plot <- ggplot2::ggplot(pval_df, ggplot2::aes(x = pvalue)) +
  ggplot2::geom_histogram(bins = 30, boundary = 0, color = "white") +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(figures_dir, paste0("bulk_pvalue_histogram_", contrast_name, ".pdf")), pval_plot, width = 5.5, height = 4)

# ── QC metrics ──

qc_flags <- character()

# 1. Dispersion fallback check
if (dispersion_fallback) {
  qc_flags <- c(qc_flags, paste0("dispersion_fallback: ", dispersion_fallback_reason))
}

# 2. DE gene coverage
n_total_genes <- nrow(result_df)
n_non_na_pval <- sum(!is.na(result_df$pvalue))
n_na_pval <- sum(is.na(result_df$pvalue))
prop_na <- n_na_pval / n_total_genes
if (prop_na > 0.5) {
  qc_flags <- c(qc_flags, sprintf("high_NA_proportion: %.1f%% (%d/%d genes have NA p-value)", prop_na * 100, n_na_pval, n_total_genes))
}
if (n_non_na_pval < 100L) {
  qc_flags <- c(qc_flags, sprintf("low_effective_genes: only %d genes with non-NA p-value", n_non_na_pval))
}

# 3. P-value distribution
non_na_pvals <- result_df$pvalue[!is.na(result_df$pvalue)]
if (length(non_na_pvals) >= 100L) {
  pval_hist <- hist(non_na_pvals, breaks = 20, plot = FALSE)
  low_pval_frac <- sum(pval_hist$counts[1:2]) / sum(pval_hist$counts)
  if (low_pval_frac < 0.02) {
    qc_flags <- c(qc_flags, "pvalue_distribution_uniform: very few low p-values suggest no differential signal")
  }
  if (sum(result_df$padj < 0.05, na.rm = TRUE) == 0L) {
    qc_flags <- c(qc_flags, "zero_DE_genes_at_padj0.05: no gene passes FDR=5%")
  }
} else {
  qc_flags <- c(qc_flags, "insufficient_pvalues_for_diagnostics")
}

# 4. Library-size outlier detection
lib_median <- stats::median(lib_sizes$library_size)
lib_mad <- stats::mad(lib_sizes$library_size, constant = 1L)
lib_outliers <- character()
if (lib_mad > 0) {
  lib_outliers <- lib_sizes$sample_id[abs(lib_sizes$library_size - lib_median) / lib_mad > 4]
  if (length(lib_outliers) > 0L) {
    qc_flags <- c(qc_flags, paste0("library_size_outlier: ", paste(lib_outliers, collapse = ", ")))
  }
}

# 5. PCA outlier detection
pca_scores <- pca$x[, 1:min(2, ncol(pca$x)), drop = FALSE]
pca_center <- colMeans(pca_scores)
pca_dist <- sqrt(rowSums((pca_scores - matrix(pca_center, nrow = nrow(pca_scores), ncol = ncol(pca_scores), byrow = TRUE))^2))
pca_mad <- stats::mad(pca_dist, constant = 1L)
pca_center_dist <- stats::median(pca_dist) + 4 * pca_mad
pca_outliers <- character()
if (pca_mad > 0 && any(pca_dist > pca_center_dist)) {
  pca_outliers <- rownames(pca_scores)[pca_dist > pca_center_dist]
  qc_flags <- c(qc_flags, paste0("pca_outlier: ", paste(pca_outliers, collapse = ", ")))
}

# 6. Cook's distance warning
cooks_warning <- FALSE
if (any(grepl("cooks", tolower(names(result_df))))) {
  cooks_col <- grep("cooks", tolower(names(result_df)), value = TRUE)[[1L]]
  if (any(!is.na(result_df[[cooks_col]]) & result_df[[cooks_col]] > stats::qf(0.99, length(all.vars(design_formula)), ncol(count_mat) - length(all.vars(design_formula))))) {
    cooks_warning <- TRUE
    qc_flags <- c(qc_flags, "cooks_distance_outliers_detected")
  }
}

qc_table <- data.frame(
  check = c(
    "dispersion_fallback",
    "na_proportion",
    "n_genes_tested",
    "n_genes_non_na_pval",
    "n_de_genes_padj05",
    "pvalue_low_fraction",
    "library_size_outliers",
    "pca_outliers",
    "cooks_warning"
  ),
  value = c(
    as.character(dispersion_fallback),
    sprintf("%.3f", prop_na),
    as.character(n_total_genes),
    as.character(n_non_na_pval),
    as.character(sum(result_df$padj < 0.05, na.rm = TRUE)),
    sprintf("%.3f", if (length(non_na_pvals) >= 100L) low_pval_frac else NA_real_),
    as.character(length(lib_outliers)),
    as.character(length(pca_outliers)),
    as.character(cooks_warning)
  ),
  flag = c(
    if (dispersion_fallback) "REVIEW" else "PASS",
    if (prop_na > 0.5) "REVIEW" else "PASS",
    if (n_total_genes < 100L) "REVIEW" else "PASS",
    if (n_non_na_pval < 100L) "REVIEW" else "PASS",
    if (sum(result_df$padj < 0.05, na.rm = TRUE) == 0L) "REVIEW" else "PASS",
    if (length(non_na_pvals) >= 100L && low_pval_frac < 0.02) "REVIEW" else "PASS",
    if (length(lib_outliers) > 0L) "REVIEW" else "PASS",
    if (length(pca_outliers) > 0L) "REVIEW" else "PASS",
    if (cooks_warning) "REVIEW" else "PASS"
  ),
  stringsAsFactors = FALSE
)
qc_table <- qc_table[qc_table$flag == "REVIEW" | qc_table$check %in% c("n_genes_tested", "n_genes_non_na_pval", "n_de_genes_padj05", "na_proportion"), ]
utils::write.table(qc_table, file.path(tables_dir, "qc_checks.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# ── Status determination ──

critical_outputs <- c(
  file.path(tables_dir, "sample_mapping_used.tsv"),
  file.path(tables_dir, "design_matrix.tsv"),
  file.path(tables_dir, "library_sizes.tsv"),
  file.path(tables_dir, "pca_coordinates.tsv"),
  file.path(tables_dir, "qc_checks.tsv"),
  de_path,
  file.path(figures_dir, "bulk_library_sizes.pdf"),
  file.path(figures_dir, "bulk_pca.pdf"),
  file.path(figures_dir, paste0("bulk_pvalue_histogram_", contrast_name, ".pdf"))
)
outputs_complete <- all(file.exists(critical_outputs) & file.info(critical_outputs)$size > 0)

qc_pass <- length(qc_flags) == 0L

if (!outputs_complete) {
  state <- "EXECUTION_COMPLETE"
  note <- "Bulk raw-count DESeq2 driver ran but one or more required outputs are missing or empty."
} else if (qc_pass) {
  state <- "BASIC_ANALYSIS_COMPLETE"
  note <- sprintf(
    "Bulk raw-count DESeq2 driver completed. %d genes tested, %d DE genes at padj<0.05. All QC checks passed.",
    n_total_genes, sum(result_df$padj < 0.05, na.rm = TRUE)
  )
} else {
  state <- "QC_REVIEW_REQUIRED"
  note <- paste("QC checks flagged for review:", paste(qc_flags, collapse = "; "))
}

status <- data.frame(
  state = state,
  updated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  note = note,
  stringsAsFactors = FALSE
)
utils::write.table(status, file.path(manifest_dir, "workflow_status.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

session_lines <- utils::capture.output(utils::sessionInfo())
writeLines(session_lines, file.path(logs_dir, "sessionInfo_bulk_counts.txt"), useBytes = TRUE)
cat(state, "\n", sep = "")
cat(normalizePath(manifest_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")
