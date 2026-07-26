#!/usr/bin/env Rscript

# ── Single-cell RNA-seq Analysis Template ────────────────────────────────
# SKELETAL REFERENCE — customize thresholds, species, and parameters per dataset.
# Source this file from a dataset-specific driver script; do not edit the template.

# Expected inputs (set before sourcing):
#   scrna_obj — Seurat object, SingleCellExperiment, or h5ad file path
#   out_dir   — output directory path for derived results
#   species   — "human" or "mouse"
#   project_name — short label for file naming
# Optional inputs:
#   apply_qc_filter — TRUE only after reviewing QC distributions; default FALSE

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

for (d in c("derived", "tables", "figures", "logs")) {
  dir.create(file.path(out_dir, d), showWarnings = FALSE, recursive = TRUE)
}

log_file <- file.path(out_dir, "logs", paste0("scrna_analysis_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("=== scRNA Analysis Template ===\n")
cat("Started:", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n")
cat("Project:", project_name, "\n")
cat("Species:", species, "\n\n")

# ── Step 0: Package check ────────────────────────────────────────────────

pkg_status <- data.frame(
  package = c("Seurat", "ggplot2", "patchwork", "Matrix"),
  installed = vapply(c("Seurat", "ggplot2", "patchwork", "Matrix"),
                     requireNamespace, logical(1), quietly = TRUE),
  stringsAsFactors = FALSE
)
cat("Package status:\n")
print(pkg_status)

if (!pkg_status$installed[1]) stop("Seurat is required for scRNA analysis.")

get_assay_data_compat <- function(object, assay = "RNA", layer = "counts") {
  tryCatch(
    Seurat::GetAssayData(object, assay = assay, layer = layer),
    error = function(e) Seurat::GetAssayData(object, assay = assay, slot = layer)
  )
}

# ── Step 1: Intake ───────────────────────────────────────────────────────

cat("\n--- Intake ---\n")

# Determine input type and load if needed
if (is.character(scrna_obj) && endsWith(scrna_obj, ".h5ad")) {
  cat("Loading h5ad...\n")
  if (requireNamespace("Seurat", quietly = TRUE) &&
      exists("ReadH5AD", where = asNamespace("Seurat"))) {
    sobj <- Seurat::ReadH5AD(scrna_obj)
  } else if (requireNamespace("SeuratDisk", quietly = TRUE)) {
    h5seurat_path <- sub("\\.h5ad$", ".h5seurat", scrna_obj)
    if (!file.exists(h5seurat_path)) {
      tryCatch({
        SeuratDisk::Convert(scrna_obj, dest = "h5seurat", overwrite = TRUE)
      }, error = function(e) {
        stop("SeuratDisk::Convert failed (may be incompatible with Seurat v5): ",
             conditionMessage(e))
      })
    }
    sobj <- tryCatch({
      SeuratDisk::LoadH5Seurat(h5seurat_path)
    }, error = function(e) {
      stop("SeuratDisk::LoadH5Seurat failed: ", conditionMessage(e))
    })
  } else if (requireNamespace("reticulate", quietly = TRUE)) {
    cat("Trying reticulate bridge to anndata...\n")
    ad <- reticulate::import("anndata")
    adata <- ad$read_h5ad(scrna_obj)
    counts_raw <- tryCatch(adata$raw$X, error = function(e) adata$X)
    counts_mat <- t(as.matrix(counts_raw))
    rownames(counts_mat) <- adata$var_names$values
    colnames(counts_mat) <- adata$obs_names$values
    sobj <- Seurat::CreateSeuratObject(counts = counts_mat)
  } else {
    stop("Need one of: Seurat::ReadH5AD, SeuratDisk, or reticulate+anndata to read h5ad.")
  }
} else if (inherits(scrna_obj, "Seurat")) {
  sobj <- scrna_obj
} else {
  stop("Input must be a Seurat object or .h5ad file path.")
}
cat(sprintf("Loaded: %d features x %d cells\n", nrow(sobj), ncol(sobj)))

# Determine if raw counts are available
has_raw <- tryCatch({
  !is.null(get_assay_data_compat(sobj, assay = "RNA", layer = "counts"))
}, error = function(e) FALSE)
if (has_raw) {
  cat("Raw counts detected. Proceeding with full QC.\n")
} else {
  cat("No raw counts layer — will reproduce metadata/embeddings only; no QC filtering.\n")
}

# ── Step 2: QC metrics ───────────────────────────────────────────────────

cat("\n--- QC metrics ---\n")

qc_sets <- common_qc_sets(species)

if (has_raw) {
  counts <- get_assay_data_compat(sobj, assay = "RNA", layer = "counts")

  sobj[["percent_mito"]] <- Seurat::PercentageFeatureSet(
    sobj, features = grep(
      if (species == "human") "^MT-" else "^mt-",
      rownames(sobj), value = TRUE
    )
  )
  sobj[["percent_ribo"]] <- Seurat::PercentageFeatureSet(
    sobj, features = grep(
      if (species == "human") "^RP[SL]" else "^Rp[sl]",
      rownames(sobj), value = TRUE
    )
  )
  sobj$nCount_RNA <- sobj$nCount_RNA %||% colSums(counts)
  sobj$nFeature_RNA <- sobj$nFeature_RNA %||% colSums(counts > 0)
} else {
  cat("Skipping QC metrics (no raw counts).\n")
}

# ── Step 3: QC threshold diagnostics ─────────────────────────────────────

cat("\n--- QC diagnostics ---\n")

fig_dir <- file.path(out_dir, "figures")
table_dir <- file.path(out_dir, "tables")

if (has_raw && requireNamespace("ggplot2", quietly = TRUE)) {
  metrics <- c("nCount_RNA", "nFeature_RNA", "percent_mito")
  for (metric in metrics) {
    df <- data.frame(
      value = sobj[[metric]][, 1],
      sample = sobj$orig.ident %||% rep("sample", ncol(sobj)),
      stringsAsFactors = FALSE
    )
    p <- ggplot2::ggplot(df, ggplot2::aes(x = sample, y = value)) +
      ggplot2::geom_violin(ggplot2::aes(fill = sample), alpha = 0.5) +
      ggplot2::stat_summary(fun = median, geom = "crossbar", width = 0.3, color = "black") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(title = paste("Pre-filter:", metric))

    ggplot2::ggsave(
      file.path(fig_dir, paste0("scrna_qc_prefilter_", metric, ".pdf")),
      p, width = max(6, length(unique(df$sample)) * 1.2), height = 5
    )
    cat("Saved pre-filter QC plot for", metric, "\n")
  }

  # Write QC table
  qc_tbl <- sobj[[c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo")]]
  qc_tbl$cell <- colnames(sobj)
  qc_tbl$sample <- sobj$orig.ident %||% "sample"
  utils::write.table(
    qc_tbl, file.path(table_dir, "cell_qc_metrics.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  cat("Saved cell QC metrics table.\n")

  # Proposed thresholds — MUST be reviewed and adjusted per dataset
  cat(sprintf(
    "\n=== Proposed QC thresholds (REVIEW BEFORE APPLYING) ===\n
    nFeature_RNA > %.0f & nFeature_RNA < %.0f\n
    nCount_RNA > %.0f & nCount_RNA < %.0f\n
    percent_mito < %.1f\n",
    max(200, quantile(sobj$nFeature_RNA, 0.01, na.rm = TRUE)),
    quantile(sobj$nFeature_RNA, 0.99, na.rm = TRUE),
    max(500, quantile(sobj$nCount_RNA, 0.01, na.rm = TRUE)),
    quantile(sobj$nCount_RNA, 0.99, na.rm = TRUE),
    min(25, quantile(sobj$percent_mito, 0.95, na.rm = TRUE))
  ))

  if (!exists("apply_qc_filter")) apply_qc_filter <- FALSE
  if (isTRUE(apply_qc_filter)) {
    sobj <- subset(
      sobj,
      subset = nFeature_RNA > 200 &
        nFeature_RNA < quantile(sobj$nFeature_RNA, 0.99) &
        percent_mito < 25
    )
    cat(sprintf("After default filtering: %d cells retained.\n", ncol(sobj)))
  } else {
    cat("QC filtering was not applied. Set apply_qc_filter <- TRUE only after reviewing distributions.\n")
  }

  cat(
    "Doublet note: run scDblFinder, Scrublet, or another sample-aware doublet method ",
    "when raw droplet counts and sample identities are available; save score tables before filtering.\n",
    sep = ""
  )
}

# ── Step 4: Normalization and HVG ────────────────────────────────────────

cat("\n--- Normalization ---\n")

if (has_raw) {
  sobj <- Seurat::NormalizeData(sobj, normalization.method = "LogNormalize")
  sobj <- Seurat::FindVariableFeatures(sobj, selection.method = "vst", nfeatures = 2000)
  cat(sprintf("Normalized. Top 10 variable features: %s\n",
              paste(head(Seurat::VariableFeatures(sobj), 10), collapse = ", ")))
}

# ── Step 5: PCA, clustering, UMAP ────────────────────────────────────────

cat("\n--- Dimensionality reduction ---\n")

if (has_raw) {
  sobj <- Seurat::ScaleData(sobj, features = rownames(sobj))
  sobj <- Seurat::RunPCA(sobj, features = Seurat::VariableFeatures(sobj))
} else {
  sobj <- Seurat::ScaleData(sobj)
  sobj <- Seurat::RunPCA(sobj)
}

# Elbow plot
if (requireNamespace("ggplot2", quietly = TRUE)) {
  elbow <- Seurat::ElbowPlot(sobj, ndims = 50) +
    ggplot2::labs(title = paste(project_name, "PCA elbow"))
  ggplot2::ggsave(file.path(fig_dir, "scrna_pca_elbow.pdf"), elbow, width = 6, height = 5)
}

# Use a conservative default for PCs; agent should adjust based on elbow
n_pcs <- min(30, ncol(sobj@reductions$pca))
cat(sprintf("Using %d PCs.\n", n_pcs))

sobj <- Seurat::FindNeighbors(sobj, dims = 1:n_pcs)
sobj <- tryCatch({
  Seurat::FindClusters(sobj, resolution = 0.8, algorithm = 4)
}, error = function(e) {
  cat("Leiden clustering failed; falling back to Seurat default clustering: ", conditionMessage(e), "\n", sep = "")
  Seurat::FindClusters(sobj, resolution = 0.8)
})
sobj <- Seurat::RunUMAP(sobj, dims = 1:n_pcs)

cat(sprintf("Clusters: %d\n", length(unique(sobj$seurat_clusters))))
cat("Resolution note: review at least one lower and one higher resolution before final labeling.\n")

cluster_resolution_summary <- data.frame(
  resolution = "0.8",
  method = "Seurat FindClusters; Leiden algorithm=4 if available, otherwise default fallback",
  n_clusters = length(unique(sobj$seurat_clusters)),
  stringsAsFactors = FALSE
)
utils::write.table(
  cluster_resolution_summary,
  file.path(table_dir, "cluster_resolution_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

# ── Step 6: UMAP plots ───────────────────────────────────────────────────

cat("\n--- UMAP plots ---\n")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  p_umap_cluster <- Seurat::DimPlot(sobj, group.by = "seurat_clusters", label = TRUE) +
    ggplot2::labs(title = paste(project_name, "UMAP by cluster"))
  ggplot2::ggsave(file.path(fig_dir, "scrna_umap_cluster.pdf"), p_umap_cluster, width = 8, height = 7)

  if ("orig.ident" %in% names(sobj@meta.data)) {
    p_umap_sample <- Seurat::DimPlot(sobj, group.by = "orig.ident") +
      ggplot2::labs(title = paste(project_name, "UMAP by sample"))
    ggplot2::ggsave(file.path(fig_dir, "scrna_umap_sample.pdf"), p_umap_sample, width = 8, height = 7)
  }
  cat("Saved UMAP plots.\n")
}

# ── Step 6b: Generic marker-panel diagnostics ───────────────────────────

cat("\n--- Generic marker panels ---\n")

if (exists("generic_scrna_marker_panels") && exists("marker_presence_table")) {
  marker_panels <- generic_scrna_marker_panels(species)
  expr_for_markers <- tryCatch(
    get_assay_data_compat(sobj, assay = "RNA", layer = "data"),
    error = function(e) get_assay_data_compat(sobj, assay = "RNA", layer = "counts")
  )
  marker_presence <- marker_presence_table(expr_for_markers, marker_panels)
  utils::write.table(
    marker_presence,
    file.path(table_dir, "generic_marker_presence.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  cat("Saved generic marker presence table.\n")

  marker_features <- present_marker_genes(expr_for_markers, marker_panels, max_per_panel = 3L)
  if (length(marker_features) > 0L && requireNamespace("ggplot2", quietly = TRUE)) {
    p_generic_dot <- Seurat::DotPlot(sobj, features = marker_features) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(title = paste(project_name, "generic marker panel"))
    ggplot2::ggsave(
      file.path(fig_dir, "scrna_generic_marker_dotplot.pdf"),
      p_generic_dot,
      width = max(10, length(marker_features) * 0.25),
      height = 6
    )
    cat("Saved generic marker dot plot.\n")
  }
}

# ── Step 7: Cluster markers ──────────────────────────────────────────────

cat("\n--- Cluster markers ---\n")

all_markers <- Seurat::FindAllMarkers(sobj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write_marker_table(all_markers, file.path(table_dir, "cluster_markers.tsv"), label = "seurat_findall")

top_markers <- all_markers[all_markers$p_val_adj < 0.05, ]
top_per_cluster <- do.call(rbind, by(top_markers, top_markers$cluster, function(x) head(x, 5)))
cat(sprintf("Top markers per cluster written to %s\n", file.path(table_dir, "cluster_markers.tsv")))

# Marker dot plot
if (requireNamespace("ggplot2", quietly = TRUE) && nrow(top_per_cluster) > 0) {
  n_top <- min(50, nrow(top_per_cluster))
  p_dot <- Seurat::DotPlot(sobj, features = unique(top_per_cluster$gene[1:n_top])) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = paste(project_name, "cluster markers (top per cluster)"))

  ggplot2::ggsave(
    file.path(fig_dir, "scrna_marker_dotplot.pdf"),
    p_dot,
    width = max(10, n_top * 0.3),
    height = 5
  )
  cat("Saved marker dot plot.\n")
}

# ── Step 8: If author labels exist, compare ──────────────────────────────

author_label_cols <- grep(
  "author|cell_type|celltype|annotation|label",
  names(sobj@meta.data), value = TRUE, ignore.case = TRUE
)
if (length(author_label_cols) > 0L) {
  cat("\n--- Author label comparison ---\n")
  for (alab in author_label_cols) {
    cat(sprintf("\nAuthor label column '%s':\n", alab))
    agreement <- table(
      Cluster = sobj$seurat_clusters,
      AuthorLabel = sobj[[alab]][, 1]
    )
    print(agreement)
  }
}

# ── Step 9: Save outputs ─────────────────────────────────────────────────

cat("\n--- Saving ---\n")

saveRDS(sobj, file.path(out_dir, "derived", paste0(project_name, "_processed.rds")))
cat("Saved processed Seurat object.\n")

# Cluster composition
cluster_comp <- as.data.frame(table(
  cluster = sobj$seurat_clusters,
  sample = sobj$orig.ident %||% "sample"
))
utils::write.table(
  cluster_comp,
  file.path(table_dir, "cluster_composition.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# Session info
writeLines(
  capture.output(utils::sessionInfo()),
  file.path(out_dir, "logs", "session_info.txt")
)

status <- data.frame(
  state = "BASIC_ANALYSIS_COMPLETE",
  updated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  note = paste("scRNA template finished for project", project_name),
  stringsAsFactors = FALSE
)
utils::write.table(
  status,
  file.path(out_dir, "workflow_status.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

cat("\n=== scRNA analysis complete ===\n")
cat("Finished:", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n")
