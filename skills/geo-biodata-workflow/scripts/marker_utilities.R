# Marker gene utility functions - sourced by analysis templates.
# Provides reusable helpers for gene signature scoring, generic marker panels,
# marker presence checks, marker tables, and QC reporting.

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

#' Score samples/cells using a gene signature (mean expression of present genes)
#'
#' @param mat Expression matrix (genes x samples/cells)
#' @param genes Character vector of gene names to score
#' @param name Name for the output score column
#' @return Named numeric vector of scores
score_signature <- function(mat, genes, name = "score") {
  present <- intersect(genes, rownames(mat))
  if (length(present) == 0L) {
    warning("No genes from signature '", name, "' found in matrix.")
    return(setNames(rep(NA_real_, ncol(mat)), colnames(mat)))
  }
  if (length(present) < length(genes)) {
    message(sprintf(
      "Signature '%s': %d/%d genes present in data.",
      name, length(present), length(genes)
    ))
  }
  sub <- as.matrix(mat[present, , drop = FALSE])
  setNames(colMeans(sub, na.rm = TRUE), colnames(mat))
}

#' Score per cluster using average expression of a gene set
#'
#' @param mat Expression matrix (genes x cells)
#' @param clusters Named factor of cluster assignments (names = cell barcodes)
#' @param genes Character vector of gene names
#' @return Matrix (clusters x genes), average expression per cluster
cluster_mean_signature <- function(mat, clusters, genes) {
  cells <- intersect(names(clusters), colnames(mat))
  present_genes <- intersect(genes, rownames(mat))
  if (length(cells) == 0L || length(present_genes) == 0L) {
    return(matrix(NA_real_, nrow = 0L, ncol = 0L))
  }
  mat <- as.matrix(mat[present_genes, cells, drop = FALSE])
  levs <- unique(clusters[cells])
  result <- do.call(rbind, lapply(levs, function(lev) {
    idx <- cells[clusters[cells] == lev]
    if (length(idx) > 0L) colMeans(mat[, idx, drop = FALSE]) else rep(NA_real_, length(present_genes))
  }))
  rownames(result) <- levs
  colnames(result) <- present_genes
  result
}

#' Write a marker table with standard columns
#'
#' @param markers Data frame with gene, cluster, and statistics columns
#' @param path Output file path
#' @param label Optional label to record in the table
write_marker_table <- function(markers, path, label = NULL) {
  if (!is.null(label)) markers$signature_label <- label
  markers$generated_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  utils::write.table(
    markers,
    file = path,
    sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )
  invisible(markers)
}

#' Store commonly used QC gene sets for human and mouse
#'
#' Returns a list of named character vectors usable with score_signature().
#' Add project-specific gene sets to the returned list.
common_qc_sets <- function(species = c("human", "mouse")) {
  species <- match.arg(species)
  if (species == "human") {
    list(
      mito = c(
        "MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2", "MT-ATP8", "MT-ATP6",
        "MT-CO3", "MT-ND3", "MT-ND4L", "MT-ND4", "MT-ND5", "MT-ND6", "MT-CYB"
      ),
      ribo = c("RPL5", "RPL11", "RPL23A", "RPS3", "RPS14", "RPS27A"),
      stress = c("HSP90AA1", "HSPA1A", "HSPA8", "HSPB1", "DNAJB1"),
      dissociation = c("FOS", "JUN", "JUNB", "FOSB", "EGR1", "ZFP36", "NR4A1", "ATF3")
    )
  } else {
    list(
      mito = c(
        "mt-Nd1", "mt-Nd2", "mt-Co1", "mt-Co2", "mt-Atp8", "mt-Atp6",
        "mt-Co3", "mt-Nd3", "mt-Nd4l", "mt-Nd4", "mt-Nd5", "mt-Nd6", "mt-Cytb"
      ),
      ribo = c("Rpl5", "Rpl11", "Rpl23a", "Rps3", "Rps14", "Rps27a"),
      stress = c("Hsp90aa1", "Hspa1a", "Hspa8", "Hspb1", "Dnajb1"),
      dissociation = c("Fos", "Jun", "Junb", "Fosb", "Egr1", "Zfp36", "Nr4a1", "Atf3")
    )
  }
}

#' Generic marker panels for first-pass scRNA annotation review
#'
#' These panels are broad sanity checks for GEO scRNA datasets. They are not a
#' final annotation system and should be reviewed against author labels and
#' dataset-specific marker literature.
#'
#' @param species "human" or "mouse"
#' @return Named list of marker gene vectors
generic_scrna_marker_panels <- function(species = c("human", "mouse")) {
  species <- match.arg(species)
  if (species == "human") {
    list(
      epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "KRT7"),
      t_cell = c("CD3D", "CD3E", "TRAC", "IL7R"),
      cd8_t_cell = c("CD8A", "CD8B", "GZMK"),
      treg = c("FOXP3", "IL2RA", "CTLA4"),
      nk_cytotoxic = c("NKG7", "GNLY", "PRF1", "GZMB", "NCAM1"),
      b_cell = c("MS4A1", "CD79A", "CD79B", "CD19"),
      plasma_cell = c("MZB1", "JCHAIN", "XBP1", "IGKC"),
      myeloid_monocyte = c("LST1", "LYZ", "S100A8", "S100A9", "FCN1", "FCGR3A", "MS4A7"),
      macrophage = c("C1QA", "C1QB", "APOE", "CD68"),
      dendritic_cell = c("FCER1A", "CLEC10A", "CST3", "LILRA4"),
      fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM"),
      endothelial = c("PECAM1", "VWF", "KDR", "PLVAP"),
      platelet = c("PPBP", "PF4"),
      mast_cell = c("TPSAB1", "TPSB2", "KIT", "CPA3"),
      cycling = c("MKI67", "TOP2A", "STMN1", "HMGB2"),
      stress_dissociation = c("FOS", "JUN", "JUNB", "DUSP1", "ATF3", "HSPA1A")
    )
  } else {
    list(
      epithelial = c("Epcam", "Krt8", "Krt18", "Krt19", "Krt7"),
      t_cell = c("Cd3d", "Cd3e", "Trac", "Il7r"),
      cd8_t_cell = c("Cd8a", "Cd8b1", "Gzmk"),
      treg = c("Foxp3", "Il2ra", "Ctla4"),
      nk_cytotoxic = c("Nkg7", "Gzma", "Prf1", "Gzmb", "Ncam1"),
      b_cell = c("Ms4a1", "Cd79a", "Cd79b", "Cd19"),
      plasma_cell = c("Mzb1", "Jchain", "Xbp1", "Igkc"),
      myeloid_monocyte = c("Lyz2", "S100a8", "S100a9", "Fcn1", "Fcgr3", "Ms4a7"),
      macrophage = c("C1qa", "C1qb", "Apoe", "Cd68"),
      dendritic_cell = c("Fcer1a", "Clec10a", "Cst3", "Siglech"),
      fibroblast = c("Col1a1", "Col1a2", "Dcn", "Lum"),
      endothelial = c("Pecam1", "Vwf", "Kdr", "Plvap"),
      platelet = c("Ppbp", "Pf4"),
      mast_cell = c("Tpsab1", "Tpsb2", "Kit", "Cpa3"),
      cycling = c("Mki67", "Top2a", "Stmn1", "Hmgb2"),
      stress_dissociation = c("Fos", "Jun", "Junb", "Dusp1", "Atf3", "Hspa1a")
    )
  }
}

#' Summarize which generic marker genes are present in an expression matrix
#'
#' @param mat Expression matrix or Seurat assay data with genes as rows
#' @param marker_panels Named list from generic_scrna_marker_panels()
#' @return Data frame with one row per marker gene
marker_presence_table <- function(mat, marker_panels) {
  genes <- rownames(mat)
  rows <- lapply(names(marker_panels), function(panel) {
    data.frame(
      panel = panel,
      gene = marker_panels[[panel]],
      present = marker_panels[[panel]] %in% genes,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

#' Select present marker genes for compact dot plots
#'
#' @param mat Expression matrix or Seurat assay data with genes as rows
#' @param marker_panels Named marker panel list
#' @param max_per_panel Maximum present genes retained per panel
#' @return Character vector of marker genes
present_marker_genes <- function(mat, marker_panels, max_per_panel = 3L) {
  presence <- marker_presence_table(mat, marker_panels)
  selected <- unlist(lapply(split(presence, presence$panel), function(tab) {
    head(tab$gene[tab$present], max_per_panel)
  }), use.names = FALSE)
  unique(selected)
}
