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

# ── Interaction detection ─────────────────────────────────────────────────────

detect_interaction <- function(design_formula) {
  terms_str <- as.character(design_formula)
  grepl(":", terms_str[length(terms_str)], fixed = TRUE) ||
    grepl("\\*", terms_str[length(terms_str)])
}

# ── Contrast factor preparation (relevel denominator as reference) ──────────

prepare_contrast_factor <- function(sample_map, factor_name, denominator) {
  x <- factor(sample_map[[factor_name]])
  levels_before <- levels(x)

  if (!denominator %in% levels_before) {
    stop("Denominator '", denominator, "' not found in factor '", factor_name,
      "' levels: ", paste(levels_before, collapse = ", "), call. = FALSE)
  }

  sample_map[[factor_name]] <- stats::relevel(x, ref = denominator)
  levels_after <- levels(sample_map[[factor_name]])

  list(
    sample_map = sample_map,
    levels_before = levels_before,
    levels_after = levels_after,
    reference = denominator
  )
}

# ── Explicit contrast matrix builder ──────────────────────────────────────────

build_contrast_matrix <- function(design, factor, numerator, denominator, sample_map) {
  # After prepare_contrast_factor() has releveled the factor,
  # denominator is the reference level (intercept).
  # The numerator coefficient is simply paste0(factor, numerator).

  if (!factor %in% names(sample_map)) {
    stop("Contrast factor '", factor, "' not found in sample mapping columns: ",
      paste(names(sample_map), collapse = ", "), call. = FALSE)
  }

  factor_values <- as.character(sample_map[[factor]])
  factor_levels <- levels(factor(sample_map[[factor]]))

  if (!numerator %in% factor_levels) {
    stop("Contrast numerator '", numerator, "' not found in factor '",
      factor, "' levels: ", paste(factor_levels, collapse = ", "), call. = FALSE)
  }
  if (!denominator %in% factor_levels) {
    stop("Contrast denominator '", denominator, "' not found in factor '",
      factor, "' levels: ", paste(factor_levels, collapse = ", "), call. = FALSE)
  }

  # After relevel, denominator IS the reference (first level), so its coefficient
  # is part of the intercept. The contrast is just the numerator coefficient.
  numerator_coef <- paste0(factor, numerator)
  design_cols <- colnames(design)

  if (!numerator_coef %in% design_cols) {
    stop("Numerator coefficient '", numerator_coef,
      "' not found in design matrix columns: ", paste(design_cols, collapse = ", "),
      ". This should not happen after prepare_contrast_factor(). ",
      "Check that the design formula includes the contrast factor.",
      call. = FALSE)
  }

  contrast_vec <- rep(0, ncol(design))
  names(contrast_vec) <- design_cols
  contrast_vec[numerator_coef] <- 1

  contrast_name <- paste(numerator, "vs", denominator, sep = "_")
  cm <- matrix(contrast_vec, ncol = 1L, dimnames = list(design_cols, contrast_name))
  cm
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

align_samples <- function(mat, sample_map, sample_id_col, contrast_factor = NULL) {
  sample_ids <- as.character(sample_map[[sample_id_col]])
  missing_from_matrix <- setdiff(sample_ids, colnames(mat))
  extra_in_matrix <- setdiff(colnames(mat), sample_ids)
  if (length(missing_from_matrix) > 0L) {
    stop("Samples in mapping not found in matrix columns: ",
      paste(missing_from_matrix, collapse = ", "), call. = FALSE)
  }
  if (length(extra_in_matrix) > 0L) {
    stop("Matrix columns not found in sample mapping: ",
      paste(extra_in_matrix, collapse = ", "), call. = FALSE)
  }
  mat <- mat[, sample_ids, drop = FALSE]
  sample_map <- sample_map[match(sample_ids, sample_map[[sample_id_col]]), , drop = FALSE]
  rownames(sample_map) <- sample_ids

  # Ensure contrast factor is a factor with explicit levels
  if (!is.factor(sample_map[[contrast_factor]])) {
    sample_map[[contrast_factor]] <- factor(sample_map[[contrast_factor]])
  }

  list(matrix = mat, metadata = sample_map)
}

# ── Scale contract validation ────────────────────────────────────────────────

validate_scale_contract <- function(input_type, scale_config, mat) {
  # Per the normalized input scale contract:
  # - normalized_expression + unknown scale -> EDA only
  # - raw-like integer matrix declared log-normalized -> fail closed
  # - TPM/FPKM/CPM without confirming log -> no limma DE

  transformed <- isTRUE(scale_config$transformed)
  transform_type <- scale_config$transform %||% ""
  evidence_source <- scale_config$evidence_source %||% ""

  # Detect raw-like integer matrix
  mat_vals <- as.vector(mat[!is.na(mat)])
  if (length(mat_vals) > 0L) {
    integer_ratio <- sum(abs(mat_vals - round(mat_vals)) < 1e-8, na.rm = TRUE) / length(mat_vals)
    all_nonneg <- all(mat_vals >= -1e-8, na.rm = TRUE)
  } else {
    integer_ratio <- 0
    all_nonneg <- TRUE
  }

  raw_like <- integer_ratio > 0.95 && all_nonneg

  if (raw_like && transformed && grepl("log", transform_type, ignore.case = TRUE)) {
    stop("INPUT_SCALE_CONFLICT: Matrix appears to be raw integer counts (",
      round(integer_ratio * 100), "% integer-like, all non-negative) ",
      "but manifest declares log-transformed normalized data. ",
      "Raw counts should use route bulk_raw_counts with run_bulk_counts.R. ",
      "If this IS normalized data, set scale.transformed=false or provide evidence.",
      call. = FALSE)
  }

  allow_de <- TRUE
  conditions <- character()

  if (input_type == "normalized_expression" && !nzchar(evidence_source)) {
    allow_de <- FALSE
    conditions <- c(conditions, "normalized_expression with unknown scale: EDA only")
  }

  if (input_type %in% c("tpm", "fpkm", "cpm")) {
    if (!transformed || !nzchar(transform_type)) {
      allow_de <- FALSE
      conditions <- c(conditions,
        sprintf("%s without confirmed log-transform: no limma DE allowed", input_type))
    }
  }

  list(
    allow_de = allow_de,
    conditions = conditions,
    raw_like = raw_like,
    integer_ratio = integer_ratio
  )
}

# ── Low-expression filtering ─────────────────────────────────────────────────

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

  if (detect_interaction(design_formula)) {
    stop("Interaction formulas are not supported by the current limma driver. ",
      "Formula: ", deparse(design_formula),
      ". Use a main-effects-only design or implement interaction support.",
      call. = FALSE)
  }

  # Step 1: Relevel so denominator is reference — makes contrast independent of original levels
  prep <- prepare_contrast_factor(
    sample_map,
    factor_name = contrast$factor,
    denominator = contrast$denominator
  )
  sample_map <- prep$sample_map

  # Step 2: Build design matrix
  design <- stats::model.matrix(design_formula, data = sample_map)

  # Step 3: Build explicit contrast matrix (numerator coefficient = 1)
  contrast_matrix <- build_contrast_matrix(
    design = design,
    factor = contrast$factor,
    numerator = contrast$numerator,
    denominator = contrast$denominator,
    sample_map = sample_map
  )

  # Step 4: Fit
  fit <- limma::lmFit(mat, design)
  fit2 <- limma::contrasts.fit(fit, contrast_matrix)
  fit2 <- limma::eBayes(fit2, trend = TRUE)

  contrast_name <- colnames(contrast_matrix)[[1L]]
  result <- limma::topTable(fit2, coef = contrast_name, number = Inf, sort.by = "none")

  list(
    fit = fit,
    ebayes_fit = fit2,
    result = result,
    design = design,
    contrast_matrix = contrast_matrix,
    contrast_name = contrast_name,
    factor_levels_before = prep$levels_before,
    factor_levels_after = prep$levels_after,
    factor_reference = prep$reference
  )
}

# ── Fallback events logger ───────────────────────────────────────────────────

init_fallback_events <- function() {
  data.frame(
    stage = character(),
    fallback_type = character(),
    trigger = character(),
    original_method = character(),
    fallback_method = character(),
    effect_on_interpretation = character(),
    requires_review = logical(),
    timestamp = character(),
    stringsAsFactors = FALSE
  )
}

log_fallback <- function(events, stage, fallback_type, trigger, original_method,
                          fallback_method, effect_on_interpretation, requires_review = TRUE) {
  rbind(events, data.frame(
    stage = stage,
    fallback_type = fallback_type,
    trigger = trigger,
    original_method = original_method,
    fallback_method = fallback_method,
    effect_on_interpretation = effect_on_interpretation,
    requires_review = requires_review,
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  ))
}

# ── Output writers ───────────────────────────────────────────────────────────

write_limma_outputs <- function(result, ebayes_fit, design, contrast_matrix,
                                 contrast_factor, numerator, denominator,
                                 sample_map, mat,
                                 factor_levels_before, factor_levels_after, factor_reference,
                                 tables_dir, figures_dir, logs_dir) {
  contrast_name <- paste(numerator, "vs", denominator, sep = "_")

  # design_matrix_used.tsv — the actual model.matrix used for fitting
  utils::write.table(
    data.frame(sample_id = rownames(design), design, check.names = FALSE),
    file.path(tables_dir, "design_matrix_used.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )

  # contrast_matrix_used.tsv — the explicit contrast
  utils::write.table(
    data.frame(coefficient = rownames(contrast_matrix),
               contrast_matrix, check.names = FALSE),
    file.path(tables_dir, "contrast_matrix_used.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )

  # factor_levels_before_relevel.tsv
  utils::write.table(
    data.frame(level = factor_levels_before, stringsAsFactors = FALSE),
    file.path(tables_dir, "factor_levels_before_relevel.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )

  # factor_levels_used.tsv (after relevel)
  lvls <- levels(sample_map[[contrast_factor]])
  n_samp <- as.integer(table(sample_map[[contrast_factor]]))
  factor_levels <- data.frame(
    factor = contrast_factor,
    level = lvls,
    n_samples = n_samp,
    reference = lvls == factor_reference,
    stringsAsFactors = FALSE
  )
  utils::write.table(factor_levels,
    file.path(tables_dir, "factor_levels_used.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )

  # contrast_resolution.txt
  resolution_lines <- c(
    paste("Contrast factor:", contrast_factor),
    paste("Numerator:", numerator),
    paste("Denominator (releveled as reference):", denominator),
    "",
    "Original factor levels:",
    paste(" ", factor_levels_before),
    "",
    "After stats::relevel(ref = denominator):",
    paste(" ", factor_levels_after),
    "",
    paste("Contrast:", numerator, "-", denominator),
    paste("Design coefficient:", paste0(contrast_factor, numerator)),
    ""
  )
  writeLines(resolution_lines, file.path(logs_dir, "contrast_resolution.txt"), useBytes = TRUE)

  # DE table
  de_out <- result
  de_out$feature_id <- rownames(de_out)
  de_out <- de_out[, c("feature_id", setdiff(names(de_out), "feature_id")), drop = FALSE]
  de_path <- file.path(tables_dir, paste0("de_results_", contrast_name, ".tsv"))
  utils::write.table(de_out, de_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  # Library sizes
  lib_sizes <- data.frame(
    sample_id = colnames(mat),
    total_signal = colSums(mat, na.rm = TRUE),
    group = as.character(sample_map[[contrast_factor]]),
    stringsAsFactors = FALSE
  )
  utils::write.table(lib_sizes, file.path(tables_dir, "library_sizes.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  lib_plot <- ggplot2::ggplot(lib_sizes, ggplot2::aes(x = sample_id, y = total_signal, fill = group)) +
    ggplot2::geom_col() + ggplot2::theme_bw() +
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
    ggplot2::geom_point(size = 3) + ggplot2::geom_text(vjust = -0.8, show.legend = FALSE) +
    ggplot2::theme_bw()
  ggplot2::ggsave(file.path(figures_dir, "bulk_pca.pdf"), pca_plot, width = 5.5, height = 4.5)

  # P-value histogram
  pval_df <- data.frame(pvalue = de_out$P.Value)
  pval_plot <- ggplot2::ggplot(pval_df, ggplot2::aes(x = pvalue)) +
    ggplot2::geom_histogram(bins = 30, boundary = 0, color = "white") + ggplot2::theme_bw()
  ggplot2::ggsave(file.path(figures_dir, paste0("bulk_pvalue_histogram_", contrast_name, ".pdf")),
    pval_plot, width = 5.5, height = 4)

  # Mean-variance trend
  grDevices::pdf(file.path(figures_dir, paste0("bulk_meanvar_", contrast_name, ".pdf")), width = 6, height = 5)
  limma::plotSA(ebayes_fit, main = "Mean-variance trend")
  grDevices::dev.off()

  # model_specification.txt
  model_spec <- c(
    paste("Contrast factor:", contrast_factor),
    paste("Numerator:", numerator),
    paste("Denominator:", denominator),
    paste("Contrast:", paste(numerator, "-", denominator)),
    "",
    "Factor levels:",
    paste("  ", factor_levels$level, " (n=", factor_levels$n_samples, ")", sep = ""),
    "",
    "Design matrix columns:",
    paste("  ", colnames(design)),
    "",
    "Contrast matrix:",
    utils::capture.output(print(contrast_matrix))
  )
  writeLines(model_spec, file.path(logs_dir, "model_specification.txt"), useBytes = TRUE)

  c(tables_dir = tables_dir, figures_dir = figures_dir, logs_dir = logs_dir,
    de_path = de_path, contrast_name = contrast_name)
}

# ── QC checks (split: technical QC vs result signal) ─────────────────────────

run_limma_qc <- function(fit_result, mat, sample_map, contrast_factor) {
  de_df <- fit_result$result
  n_total <- nrow(de_df)
  n_na_pval <- sum(is.na(de_df$P.Value))
  n_non_na <- n_total - n_na_pval
  n_de <- sum(de_df$adj.P.Val < 0.05, na.rm = TRUE)

  # ── Technical QC (data integrity, model convergence) ──
  technical_flags <- character()

  if (n_na_pval / n_total > 0.5) {
    technical_flags <- c(technical_flags,
      sprintf("high_NA_proportion: %.1f%%", n_na_pval / n_total * 100))
  }
  if (n_non_na < 10L) {
    technical_flags <- c(technical_flags, sprintf("very_few_tested_genes: %d", n_non_na))
  }

  # Library-size outlier detection
  lib_totals <- colSums(mat, na.rm = TRUE)
  lib_median <- stats::median(lib_totals)
  lib_mad <- stats::mad(lib_totals, constant = 1L)
  lib_outliers <- character()
  if (lib_mad > 0) {
    lib_outliers <- colnames(mat)[abs(lib_totals - lib_median) / lib_mad > 4]
    if (length(lib_outliers) > 0L) {
      technical_flags <- c(technical_flags,
        paste0("library_size_outlier: ", paste(lib_outliers, collapse = ", ")))
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
      technical_flags <- c(technical_flags,
        paste0("pca_outlier: ", paste(pca_outliers, collapse = ", ")))
    }
  }

  technical_qc <- if (length(technical_flags) > 0L) "REVIEW_REQUIRED" else "PASS"

  # ── Result signal (biological interpretation, NOT technical failure) ──
  signal_level <- "NOT_ASSESSED"
  non_na_pvals <- de_df$P.Value[!is.na(de_df$P.Value)]

  if (length(non_na_pvals) >= 100L) {
    pval_hist <- hist(non_na_pvals, breaks = 20, plot = FALSE)
    low_pval_frac <- sum(pval_hist$counts[1:2]) / sum(pval_hist$counts)

    signal_level <- if (n_de >= 50L && low_pval_frac > 0.05) {
      "STRONG_SIGNAL"
    } else if (n_de >= 1L && low_pval_frac > 0.02) {
      "WEAK_SIGNAL"
    } else {
      "NO_STRONG_SIGNAL"
    }
  } else {
    signal_level <- "INCONCLUSIVE"
  }

  # QC table — keep for transparency
  qc_table <- data.frame(
    check = c("execution_state", "technical_qc", "result_signal",
      "n_genes_tested", "n_genes_non_na", "n_de_genes_padj05",
      "na_proportion", "pvalue_low_fraction",
      "library_size_outlier_count", "pca_outlier_count"),
    value = c(
      "EXECUTION_COMPLETE", technical_qc, signal_level,
      as.character(n_total), as.character(n_non_na), as.character(n_de),
      sprintf("%.3f", n_na_pval / n_total),
      sprintf("%.3f", if (length(non_na_pvals) >= 100L) low_pval_frac else NA_real_),
      as.character(length(lib_outliers)), as.character(length(pca_outliers))
    ),
    stringsAsFactors = FALSE
  )

  list(
    qc_table = qc_table,
    technical_qc = technical_qc,
    result_signal = signal_level,
    technical_flags = technical_flags,
    n_de = n_de,
    n_total = n_total
  )
}

# ── Status determination (with split dimensions) ─────────────────────────────

determine_limma_status <- function(output_paths, qc_result, fallback_events) {
  outputs_complete <- all(file.exists(output_paths) & file.info(output_paths)$size > 0)

  if (!outputs_complete) {
    return(list(
      execution_state = "EXECUTION_COMPLETE",
      technical_qc = "NOT_ASSESSED",
      result_signal = "NOT_ASSESSED",
      note = "Limma driver ran but one or more required outputs are missing or empty."
    ))
  }

  execution_state <- "EXECUTION_COMPLETE"
  technical_qc <- qc_result$technical_qc
  result_signal <- qc_result$result_signal

  has_fallback_review <- any(fallback_events$requires_review)
  if (has_fallback_review && technical_qc == "PASS") {
    technical_qc <- "REVIEW_REQUIRED"
  }

  note <- sprintf(
    "Limma DE completed. %d genes tested, %d DE genes. technical_qc=%s, result_signal=%s.",
    qc_result$n_total, qc_result$n_de, technical_qc, result_signal
  )

  list(
    execution_state = execution_state,
    technical_qc = technical_qc,
    result_signal = result_signal,
    note = note
  )
}
