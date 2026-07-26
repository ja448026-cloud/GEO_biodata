#!/usr/bin/env Rscript

# ── Bulk RNA-seq / Microarray Analysis Template ──────────────────────────
# SKELETAL REFERENCE — customize thresholds, groups, contrast, and figures per dataset.
# Source this file from a dataset-specific driver script; do not edit the template.

# Expected inputs (set before sourcing):
#   counts     — gene-by-sample matrix (raw counts, normalized matrix, or ExpressionSet)
#   sample_map — data.frame with columns: sample_id, group, [batch, pair, ...]
#   out_dir    — output directory path for derived results
#   contrast   — character vector of length 2: c("numerator_group", "denominator_group")

# Load marker utilities — set template_dir before sourcing, or detect from working directory
if (!exists("template_dir")) {
  template_dir <- getwd()
}
utils_path <- file.path(template_dir, "marker_utilities.R")
if (file.exists(utils_path)) {
  source(utils_path, local = FALSE)
} else {
  cat("Note: marker_utilities.R not found at", utils_path, "- continuing without it.\n")
}

# ── Setup ────────────────────────────────────────────────────────────────

dir.create(file.path(out_dir, "derived"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "logs"), showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(out_dir, "logs", paste0("bulk_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("=== Bulk Analysis Template ===\n")
cat("Started:", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n\n")

# ── Step 0: Package check ────────────────────────────────────────────────

pkg_status <- data.frame(
  package = c("DESeq2", "edgeR", "limma", "ggplot2", "pheatmap"),
  installed = vapply(c("DESeq2", "edgeR", "limma", "ggplot2", "pheatmap"),
                     requireNamespace, logical(1), quietly = TRUE),
  stringsAsFactors = FALSE
)
cat("Package status:\n")
print(pkg_status)

if (!any(pkg_status$installed[1:3])) {
  stop("At least one of DESeq2, edgeR, or limma must be installed for bulk analysis.")
}
cat("Note: if using edgeR in a dataset-specific driver, prefer filterByExpr() and robust quasi-likelihood testing.\n")

# ── Step 1: Intake validation ────────────────────────────────────────────

cat("\n--- Intake validation ---\n")

# Merge counts and sample_map, ensuring alignment
sample_map <- sample_map[match(colnames(counts), sample_map$sample_id), , drop = FALSE]
stopifnot(identical(as.character(sample_map$sample_id), colnames(counts)))

if (!all(contrast %in% sample_map$group)) {
  stop("Contrast groups not found in sample_map$group: ", paste(contrast, collapse = ", "))
}

# Determine value type
is_integer_like <- function(x) {
  all(x == round(x), na.rm = TRUE) && all(x >= 0, na.rm = TRUE)
}
vals <- as.numeric(unlist(counts[1:min(100, nrow(counts)), 1:min(10, ncol(counts))]))
is_count <- isTRUE(all.equal(vals, round(vals), tolerance = 1e-8)) && all(vals >= 0, na.rm = TRUE)
is_log <- any(vals < 0, na.rm = TRUE) || max(vals, na.rm = TRUE) < 50

cat(sprintf(
  "Matrix: %d genes x %d samples\nCount-like: %s\nLog-transformed: %s\n",
  nrow(counts), ncol(counts),
  if (is_count) "likely yes" else "likely no (normalized or microarray)",
  if (is_log) "likely yes" else "likely no"
))

table_dir <- file.path(out_dir, "tables")
utils::write.table(
  sample_map,
  file.path(table_dir, "sample_mapping_used.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ── Step 2: Low-expression filtering ─────────────────────────────────────

cat("\n--- Low-expression filtering ---\n")

if (is_count && !is_log) {
  # DESeq2-style independent filtering
  keep <- rowSums(counts >= 10) >= min(table(sample_map$group))
  cat(sprintf("Filter: >= 10 counts in >= %d samples. Kept %d/%d genes.\n",
              min(table(sample_map$group)), sum(keep), length(keep)))
  filtered <- counts[keep, ]
} else {
  # Variance-based filtering for normalized/microarray data
  row_var <- apply(counts, 1, var, na.rm = TRUE)
  keep <- row_var > quantile(row_var, 0.1, na.rm = TRUE)
  cat(sprintf("Filter: variance > 10th percentile. Kept %d/%d genes.\n",
              sum(keep), length(keep)))
  filtered <- counts[keep, ]
}

# ── Step 3: Sample-level QC figures ──────────────────────────────────────

cat("\n--- Sample QC ---\n")

fig_dir <- file.path(out_dir, "figures")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  # Library size / total expression
  lib_df <- data.frame(
    sample = colnames(counts),
    total = colSums(counts, na.rm = TRUE),
    group = sample_map$group,
    stringsAsFactors = FALSE
  )

  p_lib <- ggplot2::ggplot(lib_df, ggplot2::aes(x = sample, y = total / 1e6, fill = group)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, size = 6)) +
    ggplot2::labs(x = "Sample", y = "Total expression (millions)", title = "Library sizes")

  ggplot2::ggsave(
    file.path(fig_dir, "bulk_library_sizes.pdf"),
    p_lib, width = max(8, ncol(counts) * 0.15), height = 5
  )
}

# ── Step 4: Variance-stabilizing transform for visualization ─────────────

cat("\n--- Transform for visualization ---\n")

if (is_count && !is_log && requireNamespace("DESeq2", quietly = TRUE)) {
  vsd <- DESeq2::vst(as.matrix(round(filtered)))
  vis_mat <- SummarizedExperiment::assay(vsd)
  cat("Applied DESeq2 vst transform.\n")
} else if (!is_log) {
  vis_mat <- log2(filtered + 1)
  cat("Applied log2(x + 1) transform.\n")
} else {
  vis_mat <- filtered
  cat("Using values as-is (appear pre-transformed).\n")
}

# ── Step 5: PCA / sample correlation ─────────────────────────────────────

cat("\n--- PCA and correlation ---\n")

pca <- prcomp(t(vis_mat), center = TRUE, scale. = TRUE)
pca_df <- data.frame(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2],
  sample = colnames(vis_mat),
  group = sample_map$group,
  stringsAsFactors = FALSE
)
utils::write.table(
  pca_df, file.path(table_dir, "pca_coordinates.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p_pca <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = group, label = sample)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_text(size = 2, vjust = -1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "PCA", subtitle = sprintf("PC1: %.1f%%, PC2: %.1f%%",
      summary(pca)$importance[2, 1] * 100, summary(pca)$importance[2, 2] * 100))

  ggplot2::ggsave(file.path(fig_dir, "bulk_pca.pdf"), p_pca, width = 8, height = 6)
  cat("Saved PCA plot.\n")
}

# Sample correlation
if (requireNamespace("pheatmap", quietly = TRUE)) {
  cor_mat <- cor(vis_mat, method = "spearman")
  pheatmap::pheatmap(
    cor_mat,
    annotation_col = data.frame(group = sample_map$group, row.names = sample_map$sample_id),
    filename = file.path(fig_dir, "bulk_sample_correlation.pdf"),
    width = max(8, ncol(vis_mat) * 0.3), height = max(6, ncol(vis_mat) * 0.3)
  )
  cat("Saved correlation heatmap.\n")
}

# ── Step 6: Differential analysis ────────────────────────────────────────

cat("\n--- Differential analysis ---\n")

contrast_str <- paste(contrast[1], "vs", contrast[2])
sample_map$group <- factor(sample_map$group)
sample_map$group <- relevel(sample_map$group, ref = contrast[2])

if (is_count && !is_log && requireNamespace("DESeq2", quietly = TRUE)) {
  cat("Running DESeq2...\n")
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(filtered),
    colData = sample_map,
    design = ~ group
  )
  dds <- DESeq2::DESeq(dds)
  res <- DESeq2::results(dds, contrast = c("group", contrast[1], contrast[2]))
  res <- res[order(res$pvalue), ]

  de_df <- data.frame(
    gene = rownames(res),
    baseMean = res$baseMean,
    log2FC = res$log2FoldChange,
    pvalue = res$pvalue,
    padj = res$padj,
    stringsAsFactors = FALSE
  )

} else if (requireNamespace("limma", quietly = TRUE)) {
  cat("Running limma...\n")
  design <- model.matrix(~ 0 + sample_map$group)
  colnames(design) <- levels(sample_map$group)
  fit <- limma::lmFit(vis_mat, design)
  cm <- limma::makeContrasts(
    contrasts = paste(contrast[1], "-", contrast[2]),
    levels = design
  )
  fit2 <- limma::contrasts.fit(fit, cm)
  fit2 <- limma::eBayes(fit2)
  res <- limma::topTable(fit2, number = Inf, sort.by = "p")

  de_df <- data.frame(
    gene = rownames(res),
    log2FC = res$logFC,
    AveExpr = res$AveExpr,
    pvalue = res$P.Value,
    padj = res$adj.P.Val,
    stringsAsFactors = FALSE
  )
} else {
  stop(
    "No implemented DE route is available. This template implements DESeq2 for raw counts ",
    "and limma for normalized/microarray data; implement edgeR in a dataset-specific driver ",
    "when edgeR is the selected route."
  )
}

cat("Effect-size note: for final DESeq2 count analyses, consider lfcShrink() when stable effect-size ranking is required.\n")

de_df$significant <- !is.na(de_df$padj) & de_df$padj < 0.05
cat(sprintf(
  "DE complete: %d significant genes (padj < 0.05), %d up, %d down.\n",
  sum(de_df$significant, na.rm = TRUE),
  sum(de_df$significant & de_df$log2FC > 0, na.rm = TRUE),
  sum(de_df$significant & de_df$log2FC < 0, na.rm = TRUE)
))

de_path <- file.path(table_dir, paste0("de_results_", contrast_str, ".tsv"))
utils::write.table(de_df, de_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("Saved DE results to", de_path, "\n")

# ── Step 7: Diagnostic figures ───────────────────────────────────────────

cat("\n--- Diagnostic figures ---\n")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  palette_okabe <- c("TRUE" = "#D55E00", "FALSE" = "#999999")

  if ("pvalue" %in% names(de_df)) {
    p_hist <- ggplot2::ggplot(de_df[!is.na(de_df$pvalue), ], ggplot2::aes(x = pvalue)) +
      ggplot2::geom_histogram(bins = 50, fill = "#56B4E9", color = "white") +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = paste("P-value histogram:", contrast_str),
        x = "Nominal p-value", y = "Genes"
      )
    ggplot2::ggsave(
      file.path(fig_dir, paste0("bulk_pvalue_histogram_", contrast_str, ".pdf")),
      p_hist, width = 7, height = 5
    )
  }

  # Volcano plot
  de_plot <- de_df[!is.na(de_df$padj), ]
  p_volc <- ggplot2::ggplot(de_plot, ggplot2::aes(x = log2FC, y = -log10(pvalue), color = significant)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.6) +
    ggplot2::scale_color_manual(values = palette_okabe) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = paste("Volcano:", contrast_str),
      x = "log2 fold change", y = "-log10 p-value"
    ) +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40")

  ggplot2::ggsave(
    file.path(fig_dir, paste0("bulk_volcano_", contrast_str, ".pdf")),
    p_volc, width = 7, height = 6
  )

  # MA plot — use whichever mean-expression column is available
  mean_col <- if ("baseMean" %in% names(de_plot)) "baseMean" else if ("AveExpr" %in% names(de_plot)) "AveExpr" else NULL
  if (!is.null(mean_col)) {
    p_ma <- ggplot2::ggplot(de_plot, ggplot2::aes(
      x = log10(.data[[mean_col]] + 1 - min(.data[[mean_col]], na.rm = TRUE) + 0.1),
      y = log2FC,
      color = significant
    )) +
    ggplot2::geom_point(size = 0.5, alpha = 0.6) +
    ggplot2::scale_color_manual(values = palette_okabe) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste("MA plot:", contrast_str))

    ggplot2::ggsave(
      file.path(fig_dir, paste0("bulk_ma_", contrast_str, ".pdf")),
      p_ma, width = 7, height = 6
    )
  }
  cat("Saved volcano and MA plots.\n")

  # Top-feature heatmap
  top_genes <- head(de_df$gene[de_df$significant], 30)
  if (length(top_genes) >= 10 && requireNamespace("pheatmap", quietly = TRUE)) {
    hm_mat <- vis_mat[top_genes, , drop = FALSE]
    hm_mat <- t(scale(t(hm_mat)))
    clip <- as.numeric(stats::quantile(abs(hm_mat), 0.99, na.rm = TRUE))
    if (is.finite(clip) && clip > 0) {
      hm_mat[hm_mat > clip] <- clip
      hm_mat[hm_mat < -clip] <- -clip
    }

    pheatmap::pheatmap(
      hm_mat,
      annotation_col = data.frame(group = sample_map$group, row.names = sample_map$sample_id),
      clustering_method = "ward.D2",
      filename = file.path(fig_dir, paste0("bulk_top_de_heatmap_", contrast_str, ".pdf")),
      width = max(8, ncol(vis_mat) * 0.3),
      height = max(6, length(top_genes) * 0.2),
      fontsize_row = 6
    )
    cat("Saved top-DE heatmap.\n")
  }
}

# ── Finish ───────────────────────────────────────────────────────────────

status <- data.frame(
  state = "BASIC_ANALYSIS_COMPLETE",
  updated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  note = paste("Bulk template finished for contrast", contrast_str),
  stringsAsFactors = FALSE
)
utils::write.table(
  status,
  file.path(out_dir, "workflow_status.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

writeLines(
  capture.output(utils::sessionInfo()),
  file.path(out_dir, "logs", "session_info.txt")
)

cat("\n=== Bulk analysis complete ===\n")
cat("Finished:", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n")
