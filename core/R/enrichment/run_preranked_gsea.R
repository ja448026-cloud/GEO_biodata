#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
usage <- paste(
  "Usage:",
  "  run_preranked_gsea.R --rank-table de.tsv --gene-column feature_id --rank-column t --gmt sets.gmt --out-dir run_dir [--collection name] [--min-size 10] [--max-size 500]",
  sep = "\n"
)
if (length(args) < 10L) stop(usage, call. = FALSE)

get_opt <- function(flag, default = NA_character_) {
  hit <- which(args == flag)
  if (length(hit) == 0L || hit[[1L]] == length(args)) return(default)
  args[[hit[[1L]] + 1L]]
}

rank_table_path <- get_opt("--rank-table")
gene_column <- get_opt("--gene-column")
rank_column <- get_opt("--rank-column")
gmt_path <- get_opt("--gmt")
out_dir <- get_opt("--out-dir")
collection <- get_opt("--collection", "custom")
min_size <- as.integer(get_opt("--min-size", "10"))
max_size <- as.integer(get_opt("--max-size", "500"))

for (required_path in c(rank_table_path, gmt_path)) {
  if (is.na(required_path) || !file.exists(required_path)) {
    stop("Missing input file: ", required_path, call. = FALSE)
  }
}
if (is.na(out_dir) || !nzchar(out_dir)) stop("--out-dir is required.", call. = FALSE)
if (!is.finite(min_size) || min_size < 1L) min_size <- 10L
if (!is.finite(max_size) || max_size < min_size) max_size <- 500L

required <- c("fgsea")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

tables_dir <- file.path(out_dir, "tables")
figures_dir <- file.path(out_dir, "figures")
logs_dir <- file.path(out_dir, "logs")
for (path in c(tables_dir, figures_dir, logs_dir)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

rank_df <- utils::read.delim(rank_table_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!gene_column %in% names(rank_df)) {
  stop("Gene column not found in rank table: ", gene_column, call. = FALSE)
}
if (!rank_column %in% names(rank_df)) {
  stop("Rank column not found in rank table: ", rank_column, call. = FALSE)
}

genes <- as.character(rank_df[[gene_column]])
ranks <- suppressWarnings(as.numeric(rank_df[[rank_column]]))
keep <- is.finite(ranks) & !is.na(genes) & nzchar(genes)
genes <- genes[keep]
ranks <- ranks[keep]
if (length(ranks) < min_size) {
  stop("Too few finite ranked genes after filtering: ", length(ranks), call. = FALSE)
}

rank_order <- order(abs(ranks), decreasing = TRUE)
genes <- genes[rank_order]
ranks <- ranks[rank_order]
dedup <- !duplicated(genes)
ranks <- ranks[dedup]
names(ranks) <- genes[dedup]
ranks <- sort(ranks, decreasing = TRUE)

gene_sets <- fgsea::gmtPathways(gmt_path)
if (!length(gene_sets)) stop("No gene sets were read from GMT: ", gmt_path, call. = FALSE)

set.seed(123L)
res <- fgsea::fgsea(pathways = gene_sets, stats = ranks, minSize = min_size, maxSize = max_size, eps = 0)
res <- res[order(res$padj, -abs(res$NES)), ]
res$leadingEdgeCount <- vapply(res$leadingEdge, length, integer(1))
res$leadingEdge <- vapply(res$leadingEdge, paste, character(1), collapse = "/")

result_path <- file.path(tables_dir, paste0("gsea_results_", collection, ".tsv"))
utils::write.table(res, result_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

manifest <- data.frame(
  collection = collection,
  rank_table = normalizePath(rank_table_path, winslash = "/", mustWork = TRUE),
  gene_column = gene_column,
  rank_column = rank_column,
  gmt = normalizePath(gmt_path, winslash = "/", mustWork = TRUE),
  gmt_sha256 = if (requireNamespace("digest", quietly = TRUE)) digest::digest(file = gmt_path, algo = "sha256") else "",
  ranked_genes = length(ranks),
  gene_sets_total = length(gene_sets),
  min_size = min_size,
  max_size = max_size,
  positive_rank_interpretation = "originating contrast numerator/up",
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE
)
utils::write.table(manifest, file.path(tables_dir, paste0("gsea_manifest_", collection, ".tsv")),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# ── Pathway redundancy reduction (leading-edge Jaccard overlap) ───────────
sig_mask <- !is.na(res$padj) & res$padj < 0.05
if (sum(sig_mask) > 1L) {
  sig_res <- res[sig_mask, ]
  n_sig <- nrow(sig_res)
  redundancy_table <- data.frame(
    pathway = sig_res$pathway,
    padj = sig_res$padj,
    NES = sig_res$NES,
    size = sig_res$size,
    leadingEdgeCount = sig_res$leadingEdgeCount,
    leadingEdge = sig_res$leadingEdge,
    redundancy_cluster = NA_integer_,
    is_representative = FALSE,
    redundant_pathways = "",
    stringsAsFactors = FALSE
  )

  # Leading-edge Jaccard similarity matrix
  le_sets <- strsplit(sig_res$leadingEdge, "/")
  names(le_sets) <- sig_res$pathway

  cluster_id <- 0L
  assigned <- rep(FALSE, n_sig)
  for (i in seq_len(n_sig)) {
    if (assigned[i]) next
    cluster_id <- cluster_id + 1L
    cluster_members <- i
    for (j in seq_len(n_sig)) {
      if (i == j || assigned[j]) next
      overlap <- length(intersect(le_sets[[i]], le_sets[[j]]))
      union_size <- length(union(le_sets[[i]], le_sets[[j]]))
      jaccard <- if (union_size > 0L) overlap / union_size else 0
      if (jaccard > 0.5) {
        cluster_members <- c(cluster_members, j)
        assigned[j] <- TRUE
      }
    }
    assigned[i] <- TRUE
    redundancy_table$redundancy_cluster[cluster_members] <- cluster_id
    # Representative: smallest padj in cluster
    rep_idx <- cluster_members[which.min(redundancy_table$padj[cluster_members])]
    redundancy_table$is_representative[rep_idx] <- TRUE
    non_rep <- setdiff(cluster_members, rep_idx)
    if (length(non_rep) > 0L) {
      redundancy_table$redundant_pathways[rep_idx] <- paste(
        redundancy_table$pathway[non_rep], collapse = " | "
      )
    }
  }

  utils::write.table(redundancy_table,
    file.path(tables_dir, paste0("gsea_pathway_redundancy_", collection, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  n_clusters <- length(unique(redundancy_table$redundancy_cluster))
  n_dedup <- sum(redundancy_table$is_representative)
  cat(sprintf("Pathway redundancy: %d significant pathways -> %d clusters (%d representative).\n",
    n_sig, n_clusters, n_dedup))
}

if (requireNamespace("ggplot2", quietly = TRUE) && nrow(res) > 0L) {
  top <- head(res, 20L)
  top$pathway <- factor(top$pathway, levels = rev(top$pathway))
  plot <- ggplot2::ggplot(top, ggplot2::aes(x = NES, y = pathway, color = padj, size = size)) +
    ggplot2::geom_point() +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Normalized enrichment score", y = NULL, color = "BH padj", size = "Set size")
  ggplot2::ggsave(file.path(figures_dir, paste0("gsea_top_", collection, ".pdf")),
    plot, width = 8, height = max(4, 0.25 * nrow(top)))
}

status <- data.frame(
  execution_state = "GSEA_COMPLETE",
  collection = collection,
  ranked_genes = length(ranks),
  gene_sets_tested = nrow(res),
  result_file = result_path,
  updated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE
)
utils::write.table(status, file.path(out_dir, "gsea_status.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

writeLines(utils::capture.output(utils::sessionInfo()),
  file.path(logs_dir, paste0("sessionInfo_gsea_", collection, ".txt")), useBytes = TRUE)

cat("GSEA_COMPLETE\n", normalizePath(out_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")
