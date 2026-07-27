# bulk_limma_common.R
# Shared limma-based DE functions for bulk_normalized and microarray_series_matrix routes.
# Source this file; do not run standalone.
#
# Requires: limma, ggplot2, stats, utils

stop_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(missing) > 0L) {
    stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

# ── Matrix and metadata loading ──────────────────────────────────────────────

read_matrix_input <- function(input_path) {
  mat_df <- utils::read.delim(input_path, stringsAsFactors = FALSE, check.names = FALSE)
  feature_ids <- as.character(mat_df[[1L]])
  if (any(!nzchar(feature_ids)) || any(duplicated(feature_ids))) {
    stop("Feature IDs must be non-empty and unique.", call. = FALSE)
  }
  mat <- as.matrix(mat_df[, -1L, drop = FALSE])
  suppressWarnings(storage.mode(mat) <- "numeric")
  if (any(!is.finite(mat))) {
    stop("Input matrix contains NA, NaN, or Inf values.", call. = FALSE)
  }
  rownames(mat) <- feature_ids
  mat
}

read_sample_mapping <- function(sample_path, sample_id_col, contrast_factor) {
  sample_map <- utils::read.delim(sample_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!sample_id_col %in% names(sample_map)) {
    stop("Sample ID column '", sample_id_col, "' not found in sample mapping.", call. = FALSE)
  }
  if (!contrast_factor %in% names(sample_map)) {
    stop("Contrast factor '", contrast_factor, "' not found in sample mapping.", call. = FALSE)
  }
  sample_map
}

align_samples <- function(mat, sample_map, sample_id_col) {
  sample_ids <- as.character(sample_map[[sample_id_col]])
  missing_from_matrix <- setdiff(sample_ids, colnames(mat))
  extra_in_matrix <- setdiff(colnames(mat), sample_ids)
  if (length(missing_from_matrix) > 0L) {
    stop("Samples in mapping not found in matrix columns: ", paste(missing_from_matrix, collapse = ", "), call. = FALSE)
  }
  if (length(extra_in_matrix) > 0L) {
    stop("Matrix columns not found in sample mapping: ", paste(extra_in_matrix, collapse = ", "), call. = FALSE)
  }
  mat <- mat[, sample_ids, drop = FALSE]
  sample_map <- sample_map[match(sample_ids, sample_map[[sample_id_col]]), , drop = FALSE]
  rownames(sample_map) <- sample_ids
  list(matrix = mat, metadata = sample_map)
}

# ── Low-expression filtering for normalized matrices ────────────────────────

filter_low_expression <- function(mat, min_expr = 0.1, min_samples_prop = 0.5) {
  min_samples <- max(1L, floor(ncol(mat) * min_samples_prop))
  keep <- rowSums(mat > min_expr, na.rm = TRUE) >= min_samples
  if (sum(keep) < 10L) {
    warning("Very few genes pass expression filter; using all genes.", call. = FALSE)
    keep <- rep(TRUE, nrow(mat))
  }
  mat[keep, , drop = FALSE]
}

# ── Core limma DE pipeline ───────────────────────────────────────────────────

run_limma_de <- function(mat, sample_map, design_formula, contrast) {
  stop_if_missing(c("limma"))
  design <- stats::model.matrix(design_formula, data = sample_map)

  # Identify the coefficient for the contrast of interest
  contrast_coef <- paste0(contrast$factor, contrast$numerator)
  if (!contrast_coef %in% colnames(design)) {
    all_coefs <- colnames(design)
    coef_candidates <- grep(paste0("^", contrast$factor), all_coefs, value = TRUE)
    if (length(coef_candidates) > 0L) {
      contrast_coef <- coef_candidates[length(coef_candidates)]
    } else {
      stop("Could not find coefficient for ", contrast$factor, contrast$numerator,
        " in design matrix columns: ", paste(all_coefs, collapse = ", "), call. = FALSE)
    }
  }

  fit <- limma::lmFit(mat, design)
  fit2 <- limma::contrasts.fit(fit, coefficients = contrast_coef)
  fit2 <- limma::eBayes(fit2, trend = TRUE)
  result <- limma::topTable(fit2, coef = contrast_coef, number = Inf, sort.by = "none")
  list(fit = fit, ebayes_fit = fit2, result = result)
}

# ── Output writers ───────────────────────────────────────────────────────────

write_limma_outputs <- function(result, ebayes_fit, contrast_factor, numerator, denominator,
                                sample_map, mat, tables_dir, figures_dir, logs_dir) {
  contrast_name <- paste(contrast_factor, numerator, "vs", denominator, sep = "_")

  # Design matrix
  design_formula <- stats::as.formula(paste("~", contrast_factor))
  design_matrix <- stats::model.matrix(design_formula, data = sample_map)
  utils::write.table(
    data.frame(sample_id = rownames(design_matrix), design_matrix, check.names = FALSE),
    file.path(tables_dir, "design_matrix.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )

  # DE table
  de_out <- result
  de_out$feature_id <- rownames(de_out)
  de_out <- de_out[, c("feature_id", setdiff(names(de_out), "feature_id")), drop = FALSE]
  de_path <- file.path(tables_dir, paste0("de_results_", contrast_name, ".tsv"))
  utils::write.table(de_out, de_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  # Library sizes (sum of expression per sample as proxy for signal scale)
  lib_sizes <- data.frame(
    sample_id = colnames(mat),
    total_signal = colSums(mat, na.rm = TRUE),
    group = as.character(sample_map[[contrast_factor]]),
    stringsAsFactors = FALSE
  )
  utils::write.table(lib_sizes, file.path(tables_dir, "library_sizes.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  lib_plot <- ggplot2::ggplot(lib_sizes, ggplot2::aes(x = sample_id, y = total_signal, fill = group)) +
    ggplot2::geom_col() +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  ggplot2::ggsave(file.path(figures_dir, "bulk_library_sizes.pdf"), lib_plot, width = 7, height = 4)

  # PCA
  pca <- stats::prcomp(t(mat), scale. = TRUE)
  pca_df <- data.frame(
    sample_id = rownames(pca$x),
    PC1 = pca$x[, 1L],
    PC2 = if (ncol(pca$x) >= 2L) pca$x[, 2L] else 0,
    group = as.character(sample_map[rownames(pca$x), contrast_factor]),
    stringsAsFactors = FALSE
  )
  utils::write.table(pca_df, file.path(tables_dir, "pca_coordinates.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  pca_plot <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = group, label = sample_id)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(vjust = -0.8, show.legend = FALSE) +
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(figures_dir, "bulk_pca.pdf"), pca_plot, width = 5.5, height = 4.5)

  # P-value histogram
  pval_df <- data.frame(pvalue = de_out$P.Value)
  pval_plot <- ggplot2::ggplot(pval_df, ggplot2::aes(x = pvalue)) +
    ggplot2::geom_histogram(bins = 30, boundary = 0, color = "white") +
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(figures_dir, paste0("bulk_pvalue_histogram_", contrast_name, ".pdf")),
    pval_plot, width = 5.5, height = 4)

  # Mean-variance trend (limma diagnostic)
  grDevices::pdf(file.path(figures_dir, paste0("bulk_meanvar_", contrast_name, ".pdf")), width = 6, height = 5)
  limma::plotSA(ebayes_fit, main = "Mean-variance trend")
  grDevices::dev.off()

  c(tables_dir = tables_dir, figures_dir = figures_dir, logs_dir = logs_dir,
    de_path = de_path, contrast_name = contrast_name)
}

# ── QC checks ────────────────────────────────────────────────────────────────

run_limma_qc <- function(fit_result, mat, sample_map, contrast_factor) {
  qc_flags <- character()

  de_df <- fit_result$result
  n_total <- nrow(de_df)
  n_na_pval <- sum(is.na(de_df$P.Value))
  n_non_na <- n_total - n_na_pval
  n_de <- sum(de_df$adj.P.Val < 0.05, na.rm = TRUE)

  if (n_na_pval / n_total > 0.5) {
    qc_flags <- c(qc_flags, sprintf("high_NA_proportion: %.1f%% (%d/%d genes have NA P.Value)",
      n_na_pval / n_total * 100, n_na_pval, n_total))
  }
  if (n_non_na < 100L) {
    qc_flags <- c(qc_flags, sprintf("low_effective_genes: only %d genes tested", n_non_na))
  }
  if (n_de == 0L) {
    qc_flags <- c(qc_flags, "zero_DE_genes_at_padj0.05: no gene passes FDR=5%")
  }

  non_na_pvals <- de_df$P.Value[!is.na(de_df$P.Value)]
  if (length(non_na_pvals) >= 100L) {
    pval_hist <- hist(non_na_pvals, breaks = 20, plot = FALSE)
    low_pval_frac <- sum(pval_hist$counts[1:2]) / sum(pval_hist$counts)
    if (low_pval_frac < 0.02) {
      qc_flags <- c(qc_flags, "pvalue_distribution_uniform: very few low p-values suggest no differential signal")
    }
  } else {
    qc_flags <- c(qc_flags, "insufficient_pvalues_for_diagnostics")
  }

  # Library-size outlier detection
  lib_totals <- colSums(mat, na.rm = TRUE)
  lib_median <- stats::median(lib_totals)
  lib_mad <- stats::mad(lib_totals, constant = 1L)
  lib_outliers <- character()
  if (lib_mad > 0) {
    lib_outliers <- colnames(mat)[abs(lib_totals - lib_median) / lib_mad > 4]
    if (length(lib_outliers) > 0L) {
      qc_flags <- c(qc_flags, paste0("library_size_outlier: ", paste(lib_outliers, collapse = ", ")))
    }
  }

  # PCA outlier detection
  pca <- stats::prcomp(t(mat), scale. = TRUE)
  pca_scores <- pca$x[, 1:min(2L, ncol(pca$x)), drop = FALSE]
  pca_center <- colMeans(pca_scores)
  pca_dist <- sqrt(rowSums((t(t(pca_scores) - pca_center))^2))
  pca_mad <- stats::mad(pca_dist, constant = 1L)
  pca_outliers <- character()
  if (pca_mad > 0) {
    pca_threshold <- stats::median(pca_dist) + 4 * pca_mad
    pca_outliers <- rownames(pca_scores)[pca_dist > pca_threshold]
    if (length(pca_outliers) > 0L) {
      qc_flags <- c(qc_flags, paste0("pca_outlier: ", paste(pca_outliers, collapse = ", ")))
    }
  }

  qc_table <- data.frame(
    check = c("n_genes_tested", "n_genes_non_na", "n_de_genes_padj05",
      "na_proportion", "pvalue_low_fraction", "library_size_outliers", "pca_outliers"),
    value = c(
      as.character(n_total), as.character(n_non_na), as.character(n_de),
      sprintf("%.3f", n_na_pval / n_total),
      sprintf("%.3f", if (length(non_na_pvals) >= 100L) low_pval_frac else NA_real_),
      as.character(length(lib_outliers)), as.character(length(pca_outliers))
    ),
    flag = c(
      if (n_total < 100L) "REVIEW" else "PASS",
      if (n_non_na < 100L) "REVIEW" else "PASS",
      if (n_de == 0L) "REVIEW" else "PASS",
      if (n_na_pval / n_total > 0.5) "REVIEW" else "PASS",
      if (length(non_na_pvals) >= 100L && low_pval_frac < 0.02) "REVIEW" else "PASS",
      if (length(lib_outliers) > 0L) "REVIEW" else "PASS",
      if (length(pca_outliers) > 0L) "REVIEW" else "PASS"
    ),
    stringsAsFactors = FALSE
  )
  qc_table <- qc_table[qc_table$flag == "REVIEW" |
    qc_table$check %in% c("n_genes_tested", "n_genes_non_na", "n_de_genes_padj05", "na_proportion"), ]

  list(qc_table = qc_table, qc_flags = qc_flags)
}

# ── Status determination ─────────────────────────────────────────────────────

determine_limma_status <- function(output_paths, qc_flags, n_genes, n_de) {
  outputs_complete <- all(file.exists(output_paths) & file.info(output_paths)$size > 0)
  qc_pass <- length(qc_flags) == 0L

  if (!outputs_complete) {
    list(state = "EXECUTION_COMPLETE",
      note = "Limma driver ran but one or more required outputs are missing or empty.")
  } else if (qc_pass) {
    list(state = "BASIC_ANALYSIS_COMPLETE",
      note = sprintf("Limma DE completed. %d genes tested, %d DE genes at padj<0.05. All QC checks passed.", n_genes, n_de))
  } else {
    list(state = "QC_REVIEW_REQUIRED",
      note = paste("QC checks flagged for review:", paste(qc_flags, collapse = "; ")))
  }
}
