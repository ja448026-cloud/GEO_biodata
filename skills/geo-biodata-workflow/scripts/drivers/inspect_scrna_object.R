#!/usr/bin/env Rscript
#
# inspect_scrna_object.R
# Read-only inventory driver for scrna_author_object and scrna_raw_counts routes.
# Supports: Seurat RDS (.rds), H5AD (.h5ad via anndataR or zellkonverter)
# NEVER: dense-converts sparse matrices, reruns PCA/clustering/UMAP, modifies the object.
#
# Outputs:
#   tables/inventory.tsv        — structured inventory of slots/layers
#   tables/cell_metadata_fields.tsv — obs/meta.data column summary
#   tables/feature_metadata_fields.tsv — var/feature metadata summary
#   logs/inventory_summary.md   — human-readable audit report
#   workflow_status.tsv         — REVIEW_REQUIRED or MANIFEST_VALIDATED

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: inspect_scrna_object.R /path/to/run_manifest.yaml", call. = FALSE)
}

manifest_path <- args[[1L]]
if (!file.exists(manifest_path)) stop("Manifest file does not exist: ", manifest_path, call. = FALSE)

required <- c("yaml")
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
  normalizePath(file.path("skills", "geo-biodata-workflow", "scripts", "drivers", "inspect_scrna_object.R"), mustWork = TRUE)
}
driver_dir <- dirname(script_path)
script_dir <- normalizePath(file.path(driver_dir, ".."), winslash = "/", mustWork = TRUE)
manifest_dir <- dirname(normalizePath(manifest_path, winslash = "/", mustWork = TRUE))

resolve_manifest_path <- function(path) {
  if (!nzchar(path %||% "")) return("")
  if (grepl("^[A-Za-z]:[\\\\/]|^/", path)) return(path)
  file.path(manifest_dir, path)
}

manifest <- yaml::read_yaml(manifest_path)
route <- manifest$route %||% ""
if (!route %in% c("scrna_author_object", "scrna_raw_counts")) {
  stop("inspect_scrna_object.R supports routes: scrna_author_object, scrna_raw_counts", call. = FALSE)
}

input_path <- resolve_manifest_path(manifest$input$file)
if (!nzchar(input_path) || !file.exists(input_path)) {
  stop("Input file does not exist: ", manifest$input$file %||% "", call. = FALSE)
}

tables_dir <- file.path(manifest_dir, "tables")
logs_dir <- file.path(manifest_dir, "logs")
for (path in c(tables_dir, logs_dir)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# ── Format detection ─────────────────────────────────────────────────────────

ext <- tolower(tools::file_ext(input_path))
if (ext == "gz") {
  ext <- tolower(tools::file_ext(sub("\\.gz$", "", input_path)))
  if (ext == "gz") ext <- "rds"
}

inventory <- list()
inventory$object_format <- NA_character_
inventory$object_size_mb <- round(file.info(input_path)$size / 1e6, 2)

# ── Seurat RDS path ──────────────────────────────────────────────────────────

if (ext %in% c("rds", "rdata", "rda")) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required to read RDS objects. Install with: install.packages('Seurat')", call. = FALSE)
  }
  obj <- readRDS(input_path)
  class_info <- class(obj)

  if ("Seurat" %in% class_info) {
    inventory$object_format <- "Seurat"
    inventory$seurat_version <- paste(as.character(obj@version), collapse = ".")

    # Matrix orientation: Seurat stores genes x cells
    inventory$matrix_orientation <- "genes_x_cells"

    # Assays and layers
    assays_present <- names(obj@assays)
    inventory$assays <- paste(assays_present, collapse = "; ")

    assay_details <- list()
    for (assay_name in assays_present) {
      assay <- obj@assays[[assay_name]]
      layers <- if (inherits(assay, "Assay5")) {
        names(assay@layers)
      } else {
        c("counts", "data", "scale.data")
      }
      assay_details[[assay_name]] <- list(
        layers = layers,
        n_features = nrow(assay),
        sparsity_pct = if (inherits(assay, "Assay5") && "counts" %in% names(assay@layers)) {
          ct <- assay@layers[["counts"]]
          if (inherits(ct, "sparseMatrix")) round((1 - Matrix::nnzero(ct) / prod(dim(ct))) * 100, 1) else NA_real_
        } else NA_real_
      )
    }
    inventory$assay_details <- assay_details

    raw_count_layer <- character()
    normalized_layer <- character()
    scaled_layer <- character()
    for (aname in assays_present) {
      ad <- assay_details[[aname]]
      if (any(grepl("^counts$", ad$layers, ignore.case = TRUE))) {
        raw_count_layer <- c(raw_count_layer, paste0(aname, "/counts"))
      }
      if (any(grepl("^data$", ad$layers, ignore.case = TRUE))) {
        normalized_layer <- c(normalized_layer, paste0(aname, "/data"))
      }
      if (any(grepl("scale", ad$layers, ignore.case = TRUE))) {
        scaled_layer <- c(scaled_layer, paste0(aname, "/scale.data"))
      }
    }
    inventory$raw_count_layer <- paste(raw_count_layer, collapse = "; ")
    inventory$normalized_layer <- paste(normalized_layer, collapse = "; ")
    inventory$scaled_layer <- paste(scaled_layer, collapse = "; ")

    inventory$sparse_or_dense <- if (length(assay_details) > 0L && !is.na(assay_details[[1]]$sparsity_pct)) {
      if (assay_details[[1]]$sparsity_pct > 85) "sparse" else "dense"
    } else "unknown"

    # Cell metadata
    cell_meta <- obj@meta.data
    inventory$n_cells <- nrow(cell_meta)
    cell_fields <- data.frame(
      field = names(cell_meta),
      type = vapply(cell_meta, function(x) class(x)[[1L]], character(1L)),
      n_unique = vapply(cell_meta, function(x) length(unique(x)), integer(1L)),
      n_missing = vapply(cell_meta, function(x) sum(is.na(x)), integer(1L)),
      stringsAsFactors = FALSE
    )
    utils::write.table(cell_fields, file.path(tables_dir, "cell_metadata_fields.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    inventory$cell_metadata_fields <- nrow(cell_fields)

    # Feature metadata
    feature_meta <- NULL
    if (length(assays_present) > 0L) {
      assay1 <- obj@assays[[assays_present[[1L]]]]
      if (inherits(assay1, "Assay5")) {
        feature_meta <- assay1@meta.data
      } else {
        feature_meta <- assay1@meta.features
      }
    }
    if (!is.null(feature_meta) && ncol(feature_meta) > 0L) {
      feat_fields <- data.frame(
        field = names(feature_meta),
        type = vapply(feature_meta, function(x) class(x)[[1L]], character(1L)),
        n_unique = vapply(feature_meta, function(x) length(unique(x)), integer(1L)),
        stringsAsFactors = FALSE
      )
      utils::write.table(feat_fields, file.path(tables_dir, "feature_metadata_fields.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE, na = "")
      inventory$feature_metadata_fields <- nrow(feat_fields)
    } else {
      inventory$feature_metadata_fields <- 0L
    }

    # Embeddings
    embeddings <- names(obj@reductions)
    emb_details <- list()
    for (emb in embeddings) {
      emb_obj <- obj@reductions[[emb]]
      emb_details[[emb]] <- list(
        key = emb_obj@key,
        n_dimensions = ncol(emb_obj@cell.embeddings)
      )
    }
    inventory$embeddings <- paste(embeddings, collapse = "; ")
    inventory$embedding_details <- emb_details

    # Graphs
    graphs <- names(obj@graphs)
    graph_details <- list()
    for (gname in graphs) {
      graph_details[[gname]] <- list(
        n_cells = nrow(obj@graphs[[gname]]),
        nnz = if (inherits(obj@graphs[[gname]], "Matrix")) Matrix::nnzero(obj@graphs[[gname]]) else NA_integer_
      )
    }
    inventory$graphs <- paste(graphs, collapse = "; ")
    inventory$graph_details <- graph_details

    # Author labels — detect fields that look like cell-type annotations
    label_candidates <- grep("cell.?type|celltype|cluster|annotation|label|ident",
      names(cell_meta), value = TRUE, ignore.case = TRUE)
    label_info <- list()
    for (lc in label_candidates) {
      label_info[[lc]] <- list(
        n_levels = length(unique(cell_meta[[lc]])),
        levels_preview = paste(utils::head(unique(as.character(cell_meta[[lc]])), 10L), collapse = ", ")
      )
    }
    inventory$author_labels <- paste(label_candidates, collapse = "; ")
    inventory$author_label_details <- label_info

    # Conversion history (Seurat commands)
    commands <- names(obj@commands)
    inventory$conversion_history <- paste(commands, collapse = "; ")

    # Memory check — warn if large sparse matrix
    if (inventory$sparse_or_dense == "sparse") {
      inventory$dense_warning <- "Object contains sparse matrices. Do NOT call as.matrix() — it would materialize as dense and likely exhaust memory."
    } else {
      inventory$dense_warning <- ""
    }

  } else if ("SingleCellExperiment" %in% class_info) {
    inventory$object_format <- "SingleCellExperiment"
    inventory$matrix_orientation <- "genes_x_cells"
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("SingleCellExperiment package required.", call. = FALSE)
    }
    inventory$n_cells <- ncol(obj)
    n_assays <- length(SummarizedExperiment::assays(obj))
    inventory$assays <- paste(names(SummarizedExperiment::assays(obj)), collapse = "; ")
    inventory$sparse_or_dense <- if (inherits(SummarizedExperiment::assay(obj, 1L), "sparseMatrix")) "sparse" else "dense"
    cell_fields <- data.frame(
      field = names(SummarizedExperiment::colData(obj)),
      type = vapply(SummarizedExperiment::colData(obj), function(x) class(x)[[1L]], character(1L)),
      n_unique = vapply(SummarizedExperiment::colData(obj), function(x) length(unique(x)), integer(1L)),
      stringsAsFactors = FALSE
    )
    utils::write.table(cell_fields, file.path(tables_dir, "cell_metadata_fields.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    inventory$cell_metadata_fields <- nrow(cell_fields)
    red_dims <- names(SingleCellExperiment::reducedDims(obj))
    inventory$embeddings <- paste(red_dims, collapse = "; ")
    inventory$author_labels <- paste(grep("cell.?type|label|cluster",
      names(SummarizedExperiment::colData(obj)), value = TRUE, ignore.case = TRUE), collapse = "; ")
    inventory$conversion_history <- "not_recorded"
    inventory$raw_count_layer <- if ("counts" %in% names(SummarizedExperiment::assays(obj))) "counts" else "unknown"
    inventory$normalized_layer <- if ("logcounts" %in% names(SummarizedExperiment::assays(obj))) "logcounts" else ""
    inventory$scaled_layer <- ""
    inventory$feature_metadata_fields <- ncol(SummarizedExperiment::rowData(obj))
  } else {
    stop("Unsupported R object class: ", paste(class_info, collapse = ", "),
      ". Expected Seurat or SingleCellExperiment.", call. = FALSE)
  }

# ── H5AD path ────────────────────────────────────────────────────────────────

} else if (ext == "h5ad") {
  h5ad_method <- if (requireNamespace("anndataR", quietly = TRUE)) {
    "anndataR"
  } else if (requireNamespace("zellkonverter", quietly = TRUE)) {
    "zellkonverter"
  } else {
    stop("H5AD reading requires anndataR or zellkonverter. Install with: BiocManager::install('anndataR')", call. = FALSE)
  }

  if (h5ad_method == "anndataR") {
    obj <- anndataR::read_h5ad(input_path)
    inventory$object_format <- "AnnData_h5ad"
    inventory$matrix_orientation <- "cells_x_genes"
    inventory$n_cells <- nrow(obj)
    inventory$n_features <- ncol(obj)
    inventory$sparse_or_dense <- if (inherits(obj$X, "sparseMatrix") || inherits(obj$X, "Matrix")) "sparse" else "dense"

    if (inherits(obj$X, "sparseMatrix")) {
      inventory$sparsity_pct <- round((1 - Matrix::nnzero(obj$X) / prod(dim(obj$X))) * 100, 1)
      inventory$dense_warning <- "Object X matrix is sparse. Do NOT materialize as dense — it would likely exhaust memory."
    }

    inventory$raw_count_layer <- if ("counts" %in% names(obj$layers)) "layers/counts" else if (!is.null(obj$raw)) "raw.X" else ""
    inventory$normalized_layer <- if (is.null(obj$raw)) "X" else "X (after log1p)"
    inventory$scaled_layer <- if (any(grepl("scale", names(obj$layers), ignore.case = TRUE))) {
      paste(grep("scale", names(obj$layers), value = TRUE, ignore.case = TRUE), collapse = "; ")
    } else ""

    # Cell metadata
    obs_cols <- names(obj$obs)
    cell_fields <- data.frame(
      field = obs_cols,
      type = vapply(obj$obs, function(x) class(x)[[1L]], character(1L)),
      n_unique = vapply(obj$obs, function(x) length(unique(x)), integer(1L)),
      n_missing = vapply(obj$obs, function(x) sum(is.na(x)), integer(1L)),
      stringsAsFactors = FALSE
    )
    utils::write.table(cell_fields, file.path(tables_dir, "cell_metadata_fields.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    inventory$cell_metadata_fields <- nrow(cell_fields)

    # Feature metadata
    var_cols <- names(obj$var)
    if (length(var_cols) > 0L) {
      feat_fields <- data.frame(
        field = var_cols,
        type = vapply(obj$var, function(x) class(x)[[1L]], character(1L)),
        n_unique = vapply(obj$var, function(x) length(unique(x)), integer(1L)),
        stringsAsFactors = FALSE
      )
      utils::write.table(feat_fields, file.path(tables_dir, "feature_metadata_fields.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE, na = "")
      inventory$feature_metadata_fields <- nrow(feat_fields)
    } else {
      inventory$feature_metadata_fields <- 0L
    }

    # Embeddings
    obsm_keys <- names(obj$obsm)
    emb_details <- list()
    for (k in obsm_keys) {
      emb_details[[k]] <- list(n_dimensions = ncol(obj$obsm[[k]]))
    }
    inventory$embeddings <- paste(obsm_keys, collapse = "; ")
    inventory$embedding_details <- emb_details

    # Graphs
    obsp_keys <- names(obj$obsp)
    inventory$graphs <- paste(obsp_keys, collapse = "; ")

    # Author labels
    label_candidates <- grep("cell.?type|celltype|cluster|annotation|label|leiden|louvain",
      obs_cols, value = TRUE, ignore.case = TRUE)
    label_info <- list()
    for (lc in label_candidates) {
      vals <- obj$obs[[lc]]
      label_info[[lc]] <- list(
        n_levels = length(unique(vals)),
        levels_preview = paste(utils::head(unique(as.character(vals)), 10L), collapse = ", ")
      )
    }
    inventory$author_labels <- paste(label_candidates, collapse = "; ")
    inventory$author_label_details <- label_info

    inventory$conversion_history <- if ("log" %in% names(obj$uns)) "uns/log present" else "not_recorded"

  } else {
    # zellkonverter path: read as SingleCellExperiment
    obj <- zellkonverter::readH5AD(input_path)
    inventory$object_format <- "SingleCellExperiment_from_H5AD"
    inventory$matrix_orientation <- "genes_x_cells"
    inventory$n_cells <- ncol(obj)
    inventory$sparse_or_dense <- if (inherits(SummarizedExperiment::assay(obj, 1L), "sparseMatrix")) "sparse" else "dense"
    inventory$assays <- paste(names(SummarizedExperiment::assays(obj)), collapse = "; ")
    inventory$raw_count_layer <- if ("X" %in% names(SummarizedExperiment::assays(obj))) "X" else ""
    inventory$normalized_layer <- if ("logcounts" %in% names(SummarizedExperiment::assays(obj))) "logcounts" else ""
    inventory$scaled_layer <- ""
    cell_fields <- data.frame(
      field = names(SummarizedExperiment::colData(obj)),
      type = vapply(SummarizedExperiment::colData(obj), function(x) class(x)[[1L]], character(1L)),
      n_unique = vapply(SummarizedExperiment::colData(obj), function(x) length(unique(x)), integer(1L)),
      stringsAsFactors = FALSE
    )
    utils::write.table(cell_fields, file.path(tables_dir, "cell_metadata_fields.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    inventory$cell_metadata_fields <- nrow(cell_fields)
    inventory$feature_metadata_fields <- ncol(SummarizedExperiment::rowData(obj))
    red_dims <- names(SingleCellExperiment::reducedDims(obj))
    inventory$embeddings <- paste(red_dims, collapse = "; ")
    inventory$graphs <- ""
    inventory$author_labels <- paste(grep("cell.?type|label|cluster",
      names(SummarizedExperiment::colData(obj)), value = TRUE, ignore.case = TRUE), collapse = "; ")
    inventory$conversion_history <- "h5ad_roundtrip_via_zellkonverter"
  }

} else {
  stop("Unsupported file extension: ", ext, ". Expected .rds, .h5ad, or SingleCellExperiment RDS.", call. = FALSE)
}

# ── Write inventory table ────────────────────────────────────────────────────

inv_flat <- data.frame(
  field = c(
    "object_format", "object_size_mb", "n_cells", "n_features",
    "matrix_orientation", "sparse_or_dense",
    "raw_count_layer", "normalized_layer", "scaled_layer",
    "cell_metadata_fields", "feature_metadata_fields",
    "embeddings", "graphs", "author_labels", "conversion_history"
  ),
  value = c(
    inventory$object_format %||% "",
    as.character(inventory$object_size_mb %||% ""),
    as.character(inventory$n_cells %||% ""),
    as.character(inventory$n_features %||% ""),
    inventory$matrix_orientation %||% "",
    inventory$sparse_or_dense %||% "",
    inventory$raw_count_layer %||% "",
    inventory$normalized_layer %||% "",
    inventory$scaled_layer %||% "",
    as.character(inventory$cell_metadata_fields %||% ""),
    as.character(inventory$feature_metadata_fields %||% ""),
    inventory$embeddings %||% "",
    inventory$graphs %||% "",
    inventory$author_labels %||% "",
    inventory$conversion_history %||% ""
  ),
  stringsAsFactors = FALSE
)
utils::write.table(inv_flat, file.path(tables_dir, "inventory.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# ── Summary report ───────────────────────────────────────────────────────────

summary_lines <- c(
  paste("# scRNA Object Inventory"),
  "",
  paste("Generated:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "",
  paste("## Object"),
  sprintf("- Format: %s", inventory$object_format %||% "unknown"),
  sprintf("- File size: %.1f MB", inventory$object_size_mb %||% 0),
  sprintf("- Cells: %s", inventory$n_cells %||% "unknown"),
  sprintf("- Features: %s", inventory$n_features %||% "unknown"),
  sprintf("- Orientation: %s", inventory$matrix_orientation %||% "unknown"),
  sprintf("- Storage: %s", inventory$sparse_or_dense %||% "unknown"),
  "",
  "## Layers / Assays",
  sprintf("- Raw counts: %s", if (nzchar(inventory$raw_count_layer %||% "")) inventory$raw_count_layer else "NOT_FOUND"),
  sprintf("- Normalized: %s", if (nzchar(inventory$normalized_layer %||% "")) inventory$normalized_layer else (if (is.null(inventory$normalized_layer)) "NOT_FOUND" else "")),
  sprintf("- Scaled: %s", if (nzchar(inventory$scaled_layer %||% "")) inventory$scaled_layer else "none"),
  "",
  "## Metadata",
  sprintf("- Cell metadata fields: %d", inventory$cell_metadata_fields %||% 0L),
  sprintf("- Feature metadata fields: %d", inventory$feature_metadata_fields %||% 0L),
  "",
  "## Embeddings",
  sprintf("- Found: %s", if (nzchar(inventory$embeddings %||% "")) inventory$embeddings else "none"),
  "",
  "## Graphs",
  sprintf("- Found: %s", if (nzchar(inventory$graphs %||% "")) inventory$graphs else "none"),
  "",
  "## Author Labels",
  sprintf("- Candidate fields: %s", if (nzchar(inventory$author_labels %||% "")) inventory$author_labels else "NOT_FOUND"),
  "",
  "## Conversion History",
  sprintf("- %s", inventory$conversion_history %||% "not_recorded"),
  "",
  "## Warnings"
)

if (nzchar(inventory$dense_warning %||% "")) {
  summary_lines <- c(summary_lines, sprintf("- %s", inventory$dense_warning))
} else {
  summary_lines <- c(summary_lines, "- No critical warnings.")
}

writeLines(summary_lines, file.path(logs_dir, "inventory_summary.md"), useBytes = TRUE)

# ── Status ───────────────────────────────────────────────────────────────────

has_raw_counts <- nzchar(inventory$raw_count_layer %||% "")
has_author_labels <- nzchar(inventory$author_labels %||% "")
has_embeddings <- nzchar(inventory$embeddings %||% "")

flags <- character()
if (!has_raw_counts) flags <- c(flags, "No raw count layer found; raw-count QC route may not be possible.")
if (!has_author_labels) flags <- c(flags, "No author cell-type labels detected in metadata columns.")

# Status is REVIEW_REQUIRED until a human confirms the inventory
state <- "REVIEW_REQUIRED"
note <- paste(c(
  sprintf("Object inventory complete. Format: %s, %s cells, %s features.",
    inventory$object_format %||% "unknown", inventory$n_cells %||% "?",
    inventory$n_features %||% "?"),
  flags
), collapse = " ")

status <- data.frame(
  state = state,
  updated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  note = note,
  stringsAsFactors = FALSE
)
utils::write.table(status, file.path(manifest_dir, "workflow_status.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

session_lines <- utils::capture.output(utils::sessionInfo())
writeLines(session_lines, file.path(logs_dir, "sessionInfo_scrna_inventory.txt"), useBytes = TRUE)

cat(state, "\n", sep = "")
cat(normalizePath(manifest_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")
