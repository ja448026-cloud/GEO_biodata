#!/usr/bin/env Rscript
#
# run_ora_enrichment.R
# Over-representation analysis (ORA) driver using MSigDB GMT files.
#
# Usage:
#   Rscript run_ora_enrichment.R \
#     --gene-list de_genes.tsv \
#     --gene-column gene_symbol \
#     --universe universe.txt \
#     --gmt-dir path/to/msigdb \
#     --out-dir run_dir \
#     [--fdr-threshold 0.05] \
#     [--min-set-size 10] \
#     [--max-set-size 500] \
#     [--collections h.all,c2.cp.kegg_medicus,c5.go]

args <- commandArgs(trailingOnly = TRUE)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
usage <- paste(
  "Usage:",
  "  run_ora_enrichment.R --gene-list de.tsv --gene-column symbol --universe bg.txt",
  "    --gmt-dir msigdb/ --out-dir run_dir [--collections h.all] [--fdr 0.05]",
  sep = "\n"
)
if (length(args) < 10L) stop(usage, call. = FALSE)

get_opt <- function(flag, default = NA_character_) {
  hit <- which(args == flag)
  if (length(hit) == 0L || hit[[1L]] == length(args)) return(default)
  args[[hit[[1L]] + 1L]]
}

gene_list_path <- get_opt("--gene-list")
gene_column <- get_opt("--gene-column", "gene_symbol")
universe_path <- get_opt("--universe")
gmt_dir <- get_opt("--gmt-dir")
out_dir <- get_opt("--out-dir")
collections_str <- get_opt("--collections", "h.all,c2.cp.kegg_medicus,c5.go")
fdr_threshold <- as.numeric(get_opt("--fdr-threshold", "0.05"))
min_set_size <- as.integer(get_opt("--min-set-size", "10"))
max_set_size <- as.integer(get_opt("--max-set-size", "500"))

for (req_path in c(gene_list_path, universe_path)) {
  if (is.na(req_path) || !file.exists(req_path)) stop("Missing input: ", req_path)
}
if (is.na(gmt_dir) || !dir.exists(gmt_dir)) stop("GMT directory not found: ", gmt_dir)
if (is.na(out_dir) || !nzchar(out_dir)) stop("--out-dir is required.")

collections <- trimws(strsplit(collections_str, ",")[[1L]])
collections <- collections[nzchar(collections)]

tables_dir <- file.path(out_dir, "tables")
figures_dir <- file.path(out_dir, "figures")
logs_dir <- file.path(out_dir, "logs")
for (d in c(tables_dir, figures_dir, logs_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

write_blocked_status <- function(note, n_query_input = NA_integer_, n_query_in_universe = NA_integer_,
                                 n_query_outside_universe = NA_integer_, n_universe = NA_integer_) {
  status <- data.frame(
    execution_state = "ORA_BLOCKED",
    contract_state = "BLOCKED",
    universe_contract_status = "FAIL",
    n_query_genes_input = n_query_input,
    n_query_genes = n_query_in_universe,
    n_query_genes_outside_universe = n_query_outside_universe,
    n_universe = n_universe,
    gene_sets_tested = 0L,
    n_significant = 0L,
    result_file = "",
    updated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    note = note,
    stringsAsFactors = FALSE
  )
  utils::write.table(status, file.path(out_dir, "ora_status.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  stop(note, call. = FALSE)
}

# ── Load inputs ────────────────────────────────────────────────────────────────

gene_df <- utils::read.delim(gene_list_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!gene_column %in% names(gene_df)) {
  stop("Gene column '", gene_column, "' not found. Available: ", paste(names(gene_df), collapse = ", "))
}
query_genes_input <- unique(as.character(gene_df[[gene_column]]))
query_genes_input <- query_genes_input[nzchar(query_genes_input) & !is.na(query_genes_input)]

universe_lines <- readLines(universe_path, warn = FALSE)
universe <- unique(universe_lines[nzchar(universe_lines)])
# If universe is a one-column table, read first column
if (length(universe) > 0L && (length(universe) == 1L || grepl("\t", universe[1L]))) {
  universe <- unique(utils::read.delim(universe_path, stringsAsFactors = FALSE, check.names = FALSE)[[1L]])
  universe <- as.character(universe[nzchar(universe) & !is.na(universe)])
}
universe <- unique(trimws(as.character(universe)))
universe <- universe[nzchar(universe) & !is.na(universe)]

if (length(universe) < 100L) {
  write_blocked_status(
    sprintf("Universe has %d genes; ORA requires at least 100 background genes.", length(universe)),
    n_query_input = length(query_genes_input),
    n_universe = length(universe)
  )
}

query_outside_universe <- setdiff(query_genes_input, universe)
outside_fraction <- length(query_outside_universe) / max(1L, length(query_genes_input))
query_genes <- intersect(query_genes_input, universe)
if (outside_fraction > 0.10) {
  write_blocked_status(
    sprintf("Query/universe mismatch: %.1f%% of query genes are outside the universe.", outside_fraction * 100),
    n_query_input = length(query_genes_input),
    n_query_in_universe = length(query_genes),
    n_query_outside_universe = length(query_outside_universe),
    n_universe = length(universe)
  )
}

if (length(query_genes) < 5L) {
  write_blocked_status(
    "Fewer than 5 query genes overlap with universe.",
    n_query_input = length(query_genes_input),
    n_query_in_universe = length(query_genes),
    n_query_outside_universe = length(query_outside_universe),
    n_universe = length(universe)
  )
}

# ── GMT loader ─────────────────────────────────────────────────────────────────

load_gmt_collections <- function(gmt_dir, collections) {
  all_sets <- list()
  collection_map <- list()
  for (col in collections) {
    pattern <- paste0("^", col, ".*\\.gmt$")
    files <- list.files(gmt_dir, pattern = pattern, full.names = TRUE)
    if (length(files) == 0L) {
      cat(sprintf("WARNING: no GMT files matched collection '%s'.\n", col))
      next
    }
    for (f in files) {
      gs <- fgsea::gmtPathways(f)
      all_sets <- c(all_sets, gs)
      col_name <- gsub("\\.v\\d+\\.\\d+.*$", "", basename(f))
      col_name <- gsub("\\.Hs\\.symbols$", "", col_name)
      for (nm in names(gs)) {
        collection_map[[nm]] <- col_name
      }
    }
  }
  list(sets = all_sets, collection_map = collection_map)
}

if (!requireNamespace("fgsea", quietly = TRUE)) {
  stop("fgsea package required. Install: BiocManager::install('fgsea')")
}

gmt_data <- load_gmt_collections(gmt_dir, collections)
gene_sets <- gmt_data$sets
collection_map <- gmt_data$collection_map
if (!length(gene_sets)) stop("No gene sets loaded from GMT directory.")

cat(sprintf("Loaded %d gene sets from %d collection(s).\n", length(gene_sets), length(collections)))
cat(sprintf("Query: %d genes. Universe: %d genes.\n", length(query_genes), length(universe)))

# ── Fisher's exact test ORA ────────────────────────────────────────────────────

run_ora_fisher <- function(query_genes, universe, gene_sets, min_size, max_size) {
  results <- data.frame(
    pathway = character(),
    collection = character(),
    overlap_n = integer(),
    set_size = integer(),
    query_size = integer(),
    universe_size = integer(),
    overlap_genes = character(),
    p_value = numeric(),
    odds_ratio = numeric(),
    stringsAsFactors = FALSE
  )

  for (pathway_name in names(gene_sets)) {
    gs <- gene_sets[[pathway_name]]
    gs_in_universe <- intersect(gs, universe)
    gs_size <- length(gs_in_universe)
    if (gs_size < min_size || gs_size > max_size) next

    overlap <- intersect(query_genes, gs_in_universe)
    k <- length(overlap)
    if (k == 0L) next

    # Contingency table:
    #                 In gene set   Not in gene set
    # DE genes             k           n_query - k
    # Not DE          gs_size - k    N - n_query - gs_size + k
    n_query <- length(query_genes)
    N <- length(universe)
    a <- k
    b <- n_query - k
    c <- gs_size - k
    d <- N - n_query - gs_size + k

    mat <- matrix(c(a, b, c, d), nrow = 2L)
    ft <- stats::fisher.test(mat, alternative = "greater")
    or_val <- if (c > 0 && b > 0) (a * d) / (b * c) else Inf

    results <- rbind(results, data.frame(
      pathway = pathway_name,
      collection = collection_map[[pathway_name]] %||% "unknown",
      overlap_n = k,
      set_size = gs_size,
      query_size = n_query,
      universe_size = N,
      overlap_genes = paste(overlap, collapse = "/"),
      p_value = ft$p.value,
      odds_ratio = or_val,
      stringsAsFactors = FALSE
    ))
  }

  if (nrow(results) > 0L) {
    results$padj <- stats::p.adjust(results$p_value, method = "BH")
    results <- results[order(results$padj, -results$odds_ratio), ]
  } else {
    results$padj <- numeric()
  }
  results
}

ora_results <- run_ora_fisher(query_genes, universe, gene_sets, min_set_size, max_set_size)
cat(sprintf("ORA tested %d gene sets; %d with padj < %s.\n",
  nrow(ora_results), sum(ora_results$padj < fdr_threshold, na.rm = TRUE), fdr_threshold))

# ── Save ORA results ──────────────────────────────────────────────────────────

ora_path <- file.path(tables_dir, "ora_results.tsv")
utils::write.table(ora_results, ora_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# ── Per-collection results ────────────────────────────────────────────────────

for (col in unique(ora_results$collection)) {
  col_res <- ora_results[ora_results$collection == col, ]
  if (nrow(col_res) == 0L) next
  col_file <- file.path(tables_dir, paste0("ora_results_", col, ".tsv"))
  utils::write.table(col_res, col_file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

# ── Pathway redundancy reduction ──────────────────────────────────────────────

sig_mask <- !is.na(ora_results$padj) & ora_results$padj < fdr_threshold
if (sum(sig_mask) > 1L) {
  sig_res <- ora_results[sig_mask, ]
  n_sig <- nrow(sig_res)

  redundancy_table <- data.frame(
    pathway = sig_res$pathway,
    collection = sig_res$collection,
    padj = sig_res$padj,
    odds_ratio = sig_res$odds_ratio,
    overlap_n = sig_res$overlap_n,
    set_size = sig_res$set_size,
    overlap_genes = sig_res$overlap_genes,
    redundancy_cluster = NA_integer_,
    is_representative = FALSE,
    redundant_pathways = "",
    stringsAsFactors = FALSE
  )

  le_sets <- strsplit(sig_res$overlap_genes, "/")
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
    file.path(tables_dir, "ora_pathway_redundancy.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  n_clusters <- length(unique(redundancy_table$redundancy_cluster))
  n_dedup <- sum(redundancy_table$is_representative)
  cat(sprintf("Pathway redundancy: %d significant -> %d clusters (%d representative).\n",
    n_sig, n_clusters, n_dedup))
}

# ── Diagnostic plots ──────────────────────────────────────────────────────────

if (requireNamespace("ggplot2", quietly = TRUE) && nrow(ora_results) > 0L) {
  # Top 30 dot plot across all collections
  top <- head(ora_results[order(ora_results$padj), ], 30L)
  top$pathway <- factor(top$pathway, levels = rev(top$pathway))
  top$neg_log10_padj <- -log10(top$padj + .Machine$double.xmin)

  p <- ggplot2::ggplot(top, ggplot2::aes(x = odds_ratio, y = pathway, size = overlap_n, color = neg_log10_padj)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_viridis_c(option = "plasma", direction = -1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Odds ratio", y = NULL, size = "Overlap", color = "-log10(padj)")
  ggplot2::ggsave(file.path(figures_dir, "ora_top_dot.pdf"), p, width = 10, height = max(6, 0.25 * nrow(top)))

  # Per-collection top-20 bar plots
  for (col_name in unique(ora_results$collection)) {
    col_res <- ora_results[ora_results$collection == col_name, ]
    if (nrow(col_res) == 0L) next
    col_top <- head(col_res[order(col_res$padj), ], 20L)
    col_top$pathway <- factor(col_top$pathway, levels = rev(col_top$pathway))
    col_top$neg_log10_padj <- -log10(col_top$padj + .Machine$double.xmin)

    p2 <- ggplot2::ggplot(col_top, ggplot2::aes(x = neg_log10_padj, y = pathway, fill = odds_ratio)) +
      ggplot2::geom_col() +
      ggplot2::scale_fill_viridis_c(option = "plasma") +
      ggplot2::theme_minimal() +
      ggplot2::labs(x = "-log10(padj)", y = NULL, fill = "OR",
                    title = paste("ORA:", col_name))
    safe_col <- gsub("[^A-Za-z0-9_]", "_", col_name)
    ggplot2::ggsave(file.path(figures_dir, paste0("ora_bar_", safe_col, ".pdf")),
      p2, width = 9, height = max(5, 0.25 * nrow(col_top)))
  }
}

# ── Manifest ───────────────────────────────────────────────────────────────────

manifest <- data.frame(
  gene_list = normalizePath(gene_list_path, winslash = "/", mustWork = TRUE),
  gene_column = gene_column,
  universe = normalizePath(universe_path, winslash = "/", mustWork = TRUE),
  n_query_genes_input = length(query_genes_input),
  n_query_genes = length(query_genes),
  n_query_genes_outside_universe = length(query_outside_universe),
  query_outside_universe_fraction = round(outside_fraction, 4),
  n_universe = length(universe),
  universe_contract_status = "PASS",
  gmt_dir = normalizePath(gmt_dir, winslash = "/", mustWork = TRUE),
  collections = collections_str,
  fdr_threshold = fdr_threshold,
  min_set_size = min_set_size,
  max_set_size = max_set_size,
  gene_sets_loaded = length(gene_sets),
  gene_sets_tested = nrow(ora_results),
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE
)
utils::write.table(manifest, file.path(tables_dir, "ora_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# ── Status ────────────────────────────────────────────────────────────────────

status <- data.frame(
  execution_state = "ORA_COMPLETE",
  contract_state = "VALID",
  universe_contract_status = "PASS",
  n_query_genes_input = length(query_genes_input),
  n_query_genes = length(query_genes),
  n_query_genes_outside_universe = length(query_outside_universe),
  n_universe = length(universe),
  gene_sets_tested = nrow(ora_results),
  n_significant = sum(ora_results$padj < fdr_threshold, na.rm = TRUE),
  result_file = ora_path,
  updated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE
)
utils::write.table(status, file.path(out_dir, "ora_status.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

writeLines(utils::capture.output(utils::sessionInfo()),
  file.path(logs_dir, "sessionInfo_ora.txt"), useBytes = TRUE)

cat("ORA_COMPLETE\n", normalizePath(out_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")
