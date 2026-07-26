#!/usr/bin/env Rscript

# Lightweight preranked GSEA template.
# Expected inputs when sourced:
#   ranks     - named numeric vector sorted or unsorted; higher means numerator/up
#   gene_sets - named list of character vectors, or a path to a GMT file
#   out_dir   - run directory
# Optional inputs:
#   collection_name - short file-name label; default "custom"
#   min_size, max_size, seed_value
#
# Rank guidance:
#   Prefer DESeq2 Wald stat, limma moderated t, or another signed test statistic.
#   Do not rank by raw p-value alone. Positive ranks should mean numerator/up.

if (!exists("collection_name")) collection_name <- "custom"
if (!exists("min_size")) min_size <- 10L
if (!exists("max_size")) max_size <- 500L
if (!exists("seed_value")) seed_value <- 123L

if (!exists("ranks") || !is.numeric(ranks) || is.null(names(ranks))) {
  stop("Provide ranks as a named numeric vector.", call. = FALSE)
}
if (!exists("gene_sets")) {
  stop("Provide gene_sets as a named list or GMT file path.", call. = FALSE)
}

dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "logs", paste0("gsea_", collection_name, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("=== Preranked GSEA Template ===\n")
cat("Started:", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n")
cat("Collection:", collection_name, "\n")

if (is.character(gene_sets) && length(gene_sets) == 1L) {
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("fgsea is required to read GMT files and run this template.", call. = FALSE)
  }
  gene_sets <- fgsea::gmtPathways(gene_sets)
}
if (!is.list(gene_sets) || is.null(names(gene_sets))) {
  stop("gene_sets must be a named list or a GMT file path.", call. = FALSE)
}
if (!requireNamespace("fgsea", quietly = TRUE)) {
  stop("fgsea is required for this lightweight GSEA template.", call. = FALSE)
}

ranks <- ranks[is.finite(ranks) & !is.na(names(ranks)) & nzchar(names(ranks))]
ranks <- ranks[!duplicated(names(ranks))]
ranks <- sort(ranks, decreasing = TRUE)

set_sizes <- vapply(gene_sets, length, integer(1))
cat(sprintf("Ranked genes: %d\nGene sets: %d\n", length(ranks), length(gene_sets)))
cat(sprintf("Gene-set size range before filtering: %d-%d\n", min(set_sizes), max(set_sizes)))

set.seed(seed_value)
res <- fgsea::fgsea(
  pathways = gene_sets,
  stats = ranks,
  minSize = min_size,
  maxSize = max_size,
  eps = 0
)
res <- res[order(res$padj, -abs(res$NES)), ]
res$leadingEdgeCount <- vapply(res$leadingEdge, length, integer(1))
res$leadingEdge <- vapply(res$leadingEdge, paste, character(1), collapse = "/")

out_path <- file.path(out_dir, "tables", paste0("gsea_results_", collection_name, ".tsv"))
utils::write.table(res, out_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
cat("Saved:", out_path, "\n")

if (requireNamespace("ggplot2", quietly = TRUE) && nrow(res) > 0L) {
  top <- head(res, 20L)
  top$pathway <- factor(top$pathway, levels = rev(top$pathway))
  p <- ggplot2::ggplot(top, ggplot2::aes(x = NES, y = pathway, color = padj, size = size)) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Normalized enrichment score", y = NULL, color = "BH padj", size = "Set size")
  ggplot2::ggsave(
    file.path(out_dir, "figures", paste0("gsea_top_", collection_name, ".pdf")),
    p, width = 8, height = max(4, 0.25 * nrow(top))
  )
  cat("Saved top GSEA diagnostic plot.\n")
}

cat("GSEA_COMPLETE\n")
