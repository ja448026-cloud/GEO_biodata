#!/usr/bin/env Rscript
#
# inspect_scrna_object.R
# Read-only pre-analysis inventory driver for scrna_author_object and scrna_raw_counts routes.
# Supports: Seurat RDS (.rds), RData/.Rda (isolated load), H5AD (backed/HDF5 via anndataR)
# NEVER: dense-converts sparse matrices, reruns PCA/clustering/UMAP, modifies the object.
#
# Counts layer status: candidate (slot name found) → verified (sampled integer/nonneg check)
# Status: OBJECT_INVENTORY_COMPLETE, AGENT_ADJUDICATION_REQUIRED, USER_INPUT_REQUIRED, OBJECT_INTAKE_BLOCKED

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
  normalizePath(file.path("core", "R", "scrna", "inspect_object.R"), mustWork = TRUE)
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
inventory$count_layer_candidate <- ""
inventory$count_layer_verified <- FALSE
inventory$verification_method <- ""
inventory$sampled_integer_fraction <- NA_real_
inventory$nonnegative_fraction <- NA_real_
inventory$sparse_storage <- "unknown"
inventory$estimated_dense_memory_gb <- NA_real_
inventory$author_labels_found <- ""
inventory$donor_field_candidates <- ""
inventory$condition_field_candidates <- ""
inventory$route_recommendation <- ""
inventory$requires_agent_adjudication <- FALSE
inventory$requires_user_input <- FALSE

sample_matrix_values <- function(mat, sample_size = 100000L) {
  if (inherits(mat, "sparseMatrix")) {
    values <- mat@x
    if (length(values) == 0L) return(0)
    if (length(values) > sample_size) values <- sample(values, sample_size)
    return(as.numeric(values))
  }

  total_entries <- prod(dim(mat))
  if (!is.finite(total_entries) || total_entries < 1L) return(numeric())
  idx <- sample(total_entries, min(sample_size, total_entries))
  as.numeric(mat[idx])
}

verify_count_matrix <- function(mat, label) {
  sampled <- sample_matrix_values(mat)
  sampled <- sampled[is.finite(sampled)]
  if (length(sampled) == 0L) {
    return(data.frame(
      candidate = label, verified = FALSE, storage = "unknown",
      sampled_values = 0L, sampled_integer_fraction = NA_real_,
      nonnegative_fraction = NA_real_, stringsAsFactors = FALSE
    ))
  }
  data.frame(
    candidate = label,
    verified = sum(abs(sampled - round(sampled)) < 1e-8) / length(sampled) > 0.99 &&
      sum(sampled >= -1e-8) / length(sampled) > 0.99,
    storage = if (inherits(mat, "sparseMatrix")) "sparse" else "dense",
    sampled_values = length(sampled),
    sampled_integer_fraction = round(sum(abs(sampled - round(sampled)) < 1e-8) / length(sampled), 4),
    nonnegative_fraction = round(sum(sampled >= -1e-8) / length(sampled), 4),
    stringsAsFactors = FALSE
  )
}

get_seurat_counts_matrix <- function(assay) {
  if (inherits(assay, "Assay5")) {
    if ("counts" %in% names(assay@layers)) return(assay@layers[["counts"]])
    return(NULL)
  }
  if ("counts" %in% slotNames(assay)) return(slot(assay, "counts"))
  NULL
}

# ── Seurat RDS path ──────────────────────────────────────────────────────────
if (ext %in% c("rds")) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("SeuratObject package is required to read RDS. Install with: install.packages('SeuratObject')", call. = FALSE)
  }
  obj <- readRDS(input_path)
  class_info <- class(obj)

  if ("Seurat" %in% class_info) {
    inventory$object_format <- "Seurat_RDS"
    inventory$seurat_version <- paste(as.character(obj@version), collapse = ".")

    # Matrix orientation: Seurat stores genes x cells
    inventory$matrix_orientation <- "genes_x_cells"

    # Assays and layers
    assays_present <- names(obj@assays)
    inventory$assays <- paste(assays_present, collapse = "; ")

    assay_details <- list()
    for (assay_name in assays_present) {
      assay <- obj@assays[[assay_name]]
      layers <- if (inherits(assay, "Assay5")) names(assay@layers) else c("counts", "data", "scale.data")
      n_features <- nrow(assay)
      sparsity_pct <- NA_real_
      ct <- get_seurat_counts_matrix(assay)
      if (!is.null(ct)) {
        if (inherits(ct, "sparseMatrix")) {
          sparsity_pct <- round((1 - Matrix::nnzero(ct) / prod(dim(ct))) * 100, 1)
          inventory$sparse_storage <- "sparse"
          inventory$estimated_dense_memory_gb <- round(prod(dim(ct)) * 8 / 1e9, 2)
        }
      }
      assay_details[[assay_name]] <- list(layers = layers, n_features = n_features, sparsity_pct = sparsity_pct)
    }
    inventory$assay_details <- assay_details

    # Counts layer: candidate detection
    count_candidates <- character()
    for (aname in assays_present) {
      if (!is.null(get_seurat_counts_matrix(obj@assays[[aname]]))) {
        count_candidates <- c(count_candidates, paste0(aname, "/counts"))
      }
    }
    inventory$count_layer_candidate <- paste(count_candidates, collapse = "; ")

    # Sample verification of every candidate counts layer.
    if (length(count_candidates) > 0L && requireNamespace("Matrix", quietly = TRUE)) {
      count_checks <- do.call(rbind, lapply(count_candidates, function(candidate) {
        assay_name <- sub("/counts$", "", candidate)
        verify_count_matrix(get_seurat_counts_matrix(obj@assays[[assay_name]]), candidate)
      }))
      utils::write.table(count_checks, file.path(tables_dir, "count_layer_verification.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE, na = "")

      verified_idx <- which(count_checks$verified)
      if (length(verified_idx) > 0L) {
        best <- verified_idx[[1L]]
        inventory$count_layer_verified <- TRUE
        inventory$verification_method <- "sampled_nonzero_values_per_candidate"
        inventory$sampled_integer_fraction <- count_checks$sampled_integer_fraction[[best]]
        inventory$nonnegative_fraction <- count_checks$nonnegative_fraction[[best]]
        inventory$route_recommendation <- "scrna_raw_counts"
      } else {
        inventory$count_layer_verified <- FALSE
        inventory$verification_method <- "sampled_non_integer_per_candidate"
        inventory$sampled_integer_fraction <- max(count_checks$sampled_integer_fraction, na.rm = TRUE)
        inventory$nonnegative_fraction <- max(count_checks$nonnegative_fraction, na.rm = TRUE)
        inventory$route_recommendation <- "scrna_author_object"
      }
    }

    if (nzchar(inventory$count_layer_candidate) && !inventory$count_layer_verified) {
      inventory$requires_agent_adjudication <- TRUE
    }

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

    # Author labels detection
    label_candidates <- grep("cell.?type|celltype|cluster|annotation|label|ident",
      names(cell_meta), value = TRUE, ignore.case = TRUE)
    inventory$author_labels_found <- paste(label_candidates, collapse = "; ")

    # Donor/sample field candidates
    donor_candidates <- grep("donor|patient|sample_id|subject|individual|mouse_id",
      names(cell_meta), value = TRUE, ignore.case = TRUE)
    inventory$donor_field_candidates <- paste(donor_candidates, collapse = "; ")

    # Condition field candidates
    condition_candidates <- grep("condition|group|treatment|disease|status|timepoint|batch",
      names(cell_meta), value = TRUE, ignore.case = TRUE)
    inventory$condition_field_candidates <- paste(condition_candidates, collapse = "; ")

    # Embeddings
    embeddings <- names(obj@reductions)
    emb_details <- list()
    for (emb in embeddings) {
      emb_obj <- obj@reductions[[emb]]
      emb_details[[emb]] <- list(key = emb_obj@key, n_dimensions = ncol(emb_obj@cell.embeddings))
    }
    inventory$embeddings <- paste(embeddings, collapse = "; ")
    inventory$embedding_details <- emb_details

    # Graphs
    graphs <- names(obj@graphs)
    inventory$graphs <- paste(graphs, collapse = "; ")

    # Conversion history
    inventory$conversion_history <- paste(names(obj@commands), collapse = "; ")

    # Dense memory warning
    if (identical(inventory$sparse_storage, "sparse") && !is.na(inventory$estimated_dense_memory_gb)) {
      inventory$dense_warning <- sprintf(
        "Sparse matrix (~%.1f GB if dense). Do NOT call as.matrix() — it would exhaust memory.",
        inventory$estimated_dense_memory_gb)
    }

    inventory$feature_metadata_fields <- 0L
    if (length(assays_present) > 0L) {
      assay1 <- obj@assays[[assays_present[[1L]]]]
      if (inherits(assay1, "Assay5")) {
        feature_meta <- assay1@meta.data
      } else {
        feature_meta <- assay1@meta.features
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
      }
    }
    inventory$normalized_layer <- ""
    inventory$scaled_layer <- ""

  } else if ("SingleCellExperiment" %in% class_info) {
    inventory$object_format <- "SingleCellExperiment_RDS"
    inventory$matrix_orientation <- "genes_x_cells"
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
      stop("SingleCellExperiment package required.", call. = FALSE)
    }
    inventory$n_cells <- ncol(obj)
    inventory$sparse_storage <- if (inherits(SummarizedExperiment::assay(obj, 1L), "sparseMatrix")) "sparse" else "dense"

    if (inventory$sparse_storage == "sparse") {
      dims <- dim(SummarizedExperiment::assay(obj, 1L))
      inventory$estimated_dense_memory_gb <- round(prod(dims) * 8 / 1e9, 2)
    }

    assay_names <- names(SummarizedExperiment::assays(obj))
    inventory$count_layer_candidate <- if ("counts" %in% assay_names) "counts" else ""
    inventory$normalized_layer <- if ("logcounts" %in% assay_names) "logcounts" else ""
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
    inventory$embeddings <- paste(names(SingleCellExperiment::reducedDims(obj)), collapse = "; ")
    inventory$author_labels_found <- paste(grep("cell.?type|label|cluster",
      names(SummarizedExperiment::colData(obj)), value = TRUE, ignore.case = TRUE), collapse = "; ")
    inventory$conversion_history <- "SCE_direct_read"
    inventory$feature_metadata_fields <- ncol(SummarizedExperiment::rowData(obj))
  } else {
    stop("Unsupported R object class: ", paste(class_info, collapse = ", "), call. = FALSE)
  }

# ── RData / Rda → isolated environment load ──────────────────────────────────
} else if (ext %in% c("rdata", "rda")) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("SeuratObject package required. Install with: install.packages('SeuratObject')", call. = FALSE)
  }
  env <- new.env(parent = emptyenv())
  loaded_names <- load(input_path, envir = env)
  message("Loaded objects from RData: ", paste(loaded_names, collapse = ", "))

  # Find Seurat or SCE object in the environment
  obj <- NULL
  for (nm in loaded_names) {
    candidate <- env[[nm]]
    if (inherits(candidate, "Seurat") || inherits(candidate, "SingleCellExperiment")) {
      obj <- candidate
      inventory$conversion_history <- paste0("RData_isolated_load: ", nm)
      break
    }
  }
  if (is.null(obj)) {
    stop("No Seurat or SingleCellExperiment object found in RData file. Objects: ",
      paste(loaded_names, collapse = ", "), call. = FALSE)
  }

  # Same inventory logic as RDS path
  if (inherits(obj, "Seurat")) {
    inventory$object_format <- "Seurat_RData"
    inventory$matrix_orientation <- "genes_x_cells"
    inventory$n_cells <- ncol(obj)
    assays_present <- names(obj@assays)
    count_candidates <- character()
    for (aname in assays_present) {
      assay <- obj@assays[[aname]]
      ct <- get_seurat_counts_matrix(assay)
      if (!is.null(ct)) {
        count_candidates <- c(count_candidates, paste0(aname, "/counts"))
        inventory$sparse_storage <- if (inherits(ct, "sparseMatrix")) "sparse" else "dense"
        if (inherits(ct, "sparseMatrix")) {
          inventory$estimated_dense_memory_gb <- round(prod(dim(ct)) * 8 / 1e9, 2)
        }
      }
    }
    inventory$count_layer_candidate <- paste(count_candidates, collapse = "; ")
    if (length(count_candidates) > 0L && requireNamespace("Matrix", quietly = TRUE)) {
      count_checks <- do.call(rbind, lapply(count_candidates, function(candidate) {
        assay_name <- sub("/counts$", "", candidate)
        verify_count_matrix(get_seurat_counts_matrix(obj@assays[[assay_name]]), candidate)
      }))
      utils::write.table(count_checks, file.path(tables_dir, "count_layer_verification.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE, na = "")
      verified_idx <- which(count_checks$verified)
      inventory$count_layer_verified <- length(verified_idx) > 0L
      inventory$verification_method <- if (inventory$count_layer_verified) {
        "sampled_nonzero_values_per_candidate"
      } else {
        "sampled_non_integer_per_candidate"
      }
      if (length(verified_idx) > 0L) {
        best <- verified_idx[[1L]]
        inventory$sampled_integer_fraction <- count_checks$sampled_integer_fraction[[best]]
        inventory$nonnegative_fraction <- count_checks$nonnegative_fraction[[best]]
        inventory$route_recommendation <- "scrna_raw_counts"
      }
    }
    inventory$author_labels_found <- paste(grep("cell.?type|celltype|cluster|annotation|label|ident",
      names(obj@meta.data), value = TRUE, ignore.case = TRUE), collapse = "; ")
    inventory$embeddings <- paste(names(obj@reductions), collapse = "; ")
    inventory$graphs <- paste(names(obj@graphs), collapse = "; ")
    inventory$requires_agent_adjudication <- nzchar(inventory$count_layer_candidate)

    cell_meta <- obj@meta.data
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
  }

# ── H5AD path — default backed/HDF5 mode ─────────────────────────────────────
} else if (ext == "h5ad") {
  if (requireNamespace("anndataR", quietly = TRUE)) {
    message("Reading H5AD in backed/HDF5 mode via anndataR...")
    obj <- tryCatch(
      anndataR::read_h5ad(input_path, as = "HDF5AnnData", backed = TRUE, mode = "r"),
      error = function(e) {
        # Try alternate API
        tryCatch(
          anndataR::read_h5ad(input_path, backed = "r"),
          error = function(e2) {
            stop("Cannot open H5AD in backed mode. anndataR API may have changed. ",
              "Error: ", conditionMessage(e2), call. = FALSE)
          }
        )
      }
    )
    # Verify the object is actually backed
    obj_class <- class(obj)[[1L]]
    is_backed <- grepl("HDF5|Backed|backed", obj_class, ignore.case = TRUE)
    if (!is_backed) {
      status <- data.frame(
        inventory_state = "OBJECT_INTAKE_BLOCKED",
        updated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
        note = paste("BACKED_READER_CONTRACT_VIOLATION: anndataR returned", obj_class,
          "instead of a backed/HDF5AnnData object."),
        stringsAsFactors = FALSE
      )
      utils::write.table(status, file.path(manifest_dir, "workflow_status.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE, na = "")
      stop("BACKED_READER_CONTRACT_VIOLATION: anndataR returned ", obj_class,
        ". H5AD must be opened in backed mode. Check anndataR version.", call. = FALSE)
    }

    inventory$object_format <- obj_class
    inventory$h5ad_reader <- "anndataR"
    inventory$h5ad_backend <- "HDF5_backed"
    inventory$h5ad_mode <- "r"
    inventory$backed_verified <- TRUE
    inventory$materialization_policy <- "backed_read_only"
    inventory$reader_version <- as.character(utils::packageVersion("anndataR"))
    inventory$matrix_orientation <- "cells_x_genes"
    inventory$n_cells <- nrow(obj)
    inventory$n_features <- ncol(obj)

    # Check sparsity WITHOUT materializing
    if (inherits(obj$X, "sparseMatrix") || inherits(obj$X, "Matrix")) {
      inventory$sparse_storage <- "sparse"
      inventory$estimated_dense_memory_gb <- round(prod(dim(obj$X)) * 8 / 1e9, 2)
      inventory$dense_warning <- sprintf(
        "Sparse matrix (~%.1f GB if dense). Backed mode used — dense conversion blocked.",
        inventory$estimated_dense_memory_gb)
    } else {
      inventory$sparse_storage <- "dense"
    }

    # Counts layer candidates
    if ("counts" %in% names(obj$layers)) {
      inventory$count_layer_candidate <- "layers/counts"
    } else if (!is.null(obj$raw)) {
      inventory$count_layer_candidate <- "raw.X"
    }
    inventory$count_layer_verified <- FALSE
    inventory$verification_method <- "candidate_only"
    inventory$requires_agent_adjudication <- nzchar(inventory$count_layer_candidate)

    inventory$normalized_layer <- if (is.null(obj$raw)) "X" else "X (after log1p)"
    inventory$scaled_layer <- ""

    # Cell metadata
    cell_fields <- data.frame(
      field = names(obj$obs),
      type = vapply(obj$obs, function(x) class(x)[[1L]], character(1L)),
      n_unique = vapply(obj$obs, function(x) length(unique(x)), integer(1L)),
      n_missing = vapply(obj$obs, function(x) sum(is.na(x)), integer(1L)),
      stringsAsFactors = FALSE
    )
    utils::write.table(cell_fields, file.path(tables_dir, "cell_metadata_fields.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, na = "")
    inventory$cell_metadata_fields <- nrow(cell_fields)

    obs_cols <- names(obj$obs)
    inventory$author_labels_found <- paste(grep("cell.?type|celltype|cluster|annotation|label|leiden|louvain",
      obs_cols, value = TRUE, ignore.case = TRUE), collapse = "; ")
    inventory$donor_field_candidates <- paste(grep("donor|patient|sample_id|subject|individual",
      obs_cols, value = TRUE, ignore.case = TRUE), collapse = "; ")
    inventory$condition_field_candidates <- paste(grep("condition|group|treatment|disease|status|timepoint",
      obs_cols, value = TRUE, ignore.case = TRUE), collapse = "; ")
    inventory$embeddings <- paste(names(obj$obsm), collapse = "; ")
    inventory$graphs <- paste(names(obj$obsp), collapse = "; ")
    inventory$conversion_history <- if ("log" %in% names(obj$uns)) "uns/log present" else "not_recorded"
    inventory$feature_metadata_fields <- ncol(obj$var)

  } else if (requireNamespace("zellkonverter", quietly = TRUE)) {
    obj <- zellkonverter::readH5AD(input_path, use_hdf5 = TRUE, reader = "R")
    inventory$object_format <- "SingleCellExperiment_from_H5AD"
    inventory$h5ad_reader <- "zellkonverter"
    inventory$h5ad_backend <- "HDF5Array"
    inventory$h5ad_mode <- "r"
    inventory$matrix_orientation <- "genes_x_cells"
    inventory$n_cells <- ncol(obj)
    inventory$sparse_storage <- if (inherits(SummarizedExperiment::assay(obj, 1L), "sparseMatrix")) "sparse" else "dense"
    inventory$count_layer_candidate <- if ("X" %in% names(SummarizedExperiment::assays(obj))) "X" else ""
    inventory$normalized_layer <- if ("logcounts" %in% names(SummarizedExperiment::assays(obj))) "logcounts" else ""
    inventory$requires_agent_adjudication <- TRUE
    inventory$embeddings <- paste(names(SingleCellExperiment::reducedDims(obj)), collapse = "; ")
    inventory$author_labels_found <- paste(grep("cell.?type|label|cluster",
      names(SummarizedExperiment::colData(obj)), value = TRUE, ignore.case = TRUE), collapse = "; ")
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
  } else {
    stop("H5AD reading requires anndataR (preferred, backed mode) or zellkonverter. Install with: BiocManager::install('anndataR')", call. = FALSE)
  }
} else {
  stop("Unsupported file extension: ", ext, ". Expected .rds, .rdata, .rda, or .h5ad.", call. = FALSE)
}

# ── Write inventory table ────────────────────────────────────────────────────

inv_flat <- data.frame(
  field = c(
    "object_format", "object_size_mb", "n_cells", "n_features",
    "matrix_orientation", "sparse_storage",
    "count_layer_candidate", "count_layer_verified",
    "verification_method", "sampled_integer_fraction", "nonnegative_fraction",
    "estimated_dense_memory_gb",
    "cell_metadata_fields", "feature_metadata_fields",
    "embeddings", "graphs",
    "author_labels_found", "donor_field_candidates", "condition_field_candidates",
    "route_recommendation", "conversion_history",
    "requires_agent_adjudication", "requires_user_input",
    "h5ad_reader", "h5ad_backend", "h5ad_mode",
    "backed_verified", "materialization_policy", "reader_version"
  ),
  value = c(
    inventory$object_format %||% "",
    as.character(inventory$object_size_mb %||% ""),
    as.character(inventory$n_cells %||% ""),
    as.character(inventory$n_features %||% ""),
    inventory$matrix_orientation %||% "",
    inventory$sparse_storage %||% "unknown",
    inventory$count_layer_candidate %||% "",
    as.character(inventory$count_layer_verified %||% FALSE),
    inventory$verification_method %||% "",
    as.character(inventory$sampled_integer_fraction %||% NA),
    as.character(inventory$nonnegative_fraction %||% NA),
    as.character(inventory$estimated_dense_memory_gb %||% NA),
    as.character(inventory$cell_metadata_fields %||% ""),
    as.character(inventory$feature_metadata_fields %||% ""),
    inventory$embeddings %||% "",
    inventory$graphs %||% "",
    inventory$author_labels_found %||% "",
    inventory$donor_field_candidates %||% "",
    inventory$condition_field_candidates %||% "",
    inventory$route_recommendation %||% "",
    inventory$conversion_history %||% "",
    as.character(inventory$requires_agent_adjudication %||% FALSE),
    as.character(inventory$requires_user_input %||% FALSE),
    inventory$h5ad_reader %||% "",
    inventory$h5ad_backend %||% "",
    inventory$h5ad_mode %||% "",
    as.character(inventory$backed_verified %||% FALSE),
    inventory$materialization_policy %||% "",
    inventory$reader_version %||% ""
  ),
  stringsAsFactors = FALSE
)
utils::write.table(inv_flat, file.path(tables_dir, "inventory.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# ── Status determination ────────────────────────────────────────────────────

has_counts <- nzchar(inventory$count_layer_candidate %||% "")
counts_verified <- inventory$count_layer_verified %||% FALSE
has_labels <- nzchar(inventory$author_labels_found %||% "")
has_donor <- nzchar(inventory$donor_field_candidates %||% "")

# Determine status per the four-tier schema
if (inventory$object_format %||% "" == "" || identical(inventory$n_cells, 0L)) {
  obj_state <- "OBJECT_INTAKE_BLOCKED"
  obj_note <- "Could not read object or determine format."
} else if (has_counts && counts_verified && has_labels) {
  obj_state <- "OBJECT_INVENTORY_COMPLETE"
  obj_note <- sprintf("Object inventory complete. Format: %s, %s cells. Counts verified, labels found.",
    inventory$object_format, inventory$n_cells)
} else if (has_counts && !counts_verified) {
  obj_state <- "AGENT_ADJUDICATION_REQUIRED"
  obj_note <- sprintf("Counts layer candidate found (%s) but not verified as raw integer. Agent must adjudicate.",
    inventory$count_layer_candidate)
} else if (!has_counts && has_labels) {
  obj_state <- "USER_INPUT_REQUIRED"
  obj_note <- "No raw counts layer found. Author labels exist. User must confirm route."
} else {
  obj_state <- "AGENT_ADJUDICATION_REQUIRED"
  obj_note <- sprintf("Object format '%s' recognized but routing requires agent review.", inventory$object_format %||% "unknown")
}

status <- data.frame(
  inventory_state = obj_state,
  updated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  note = obj_note,
  stringsAsFactors = FALSE
)
utils::write.table(status, file.path(manifest_dir, "workflow_status.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# Summary
summary_lines <- c(
  "# scRNA Object Inventory",
  "",
  paste("Generated:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("State:", obj_state),
  "",
  sprintf("Format: %s | Size: %.1f MB | Cells: %s | Features: %s",
    inventory$object_format %||% "?", inventory$object_size_mb %||% 0,
    inventory$n_cells %||% "?", inventory$n_features %||% "?"),
  sprintf("Storage: %s | Orientation: %s",
    inventory$sparse_storage %||% "?", inventory$matrix_orientation %||% "?"),
  sprintf("Counts layer candidate: %s | Verified: %s",
    if (nzchar(inventory$count_layer_candidate %||% "")) inventory$count_layer_candidate else "NOT_FOUND",
    as.character(counts_verified)),
  sprintf("Labels: %s | Donor fields: %s",
    if (nzchar(inventory$author_labels_found %||% "")) inventory$author_labels_found else "NOT_FOUND",
    if (nzchar(inventory$donor_field_candidates %||% "")) inventory$donor_field_candidates else "NOT_FOUND"),
  sprintf("Route recommendation: %s",
    if (nzchar(inventory$route_recommendation %||% "")) inventory$route_recommendation else "review_required")
)

if (nzchar(inventory$dense_warning %||% "")) {
  summary_lines <- c(summary_lines, "", paste("WARNING:", inventory$dense_warning))
}

writeLines(summary_lines, file.path(logs_dir, "inventory_summary.md"), useBytes = TRUE)

session_lines <- utils::capture.output(utils::sessionInfo())
writeLines(session_lines, file.path(logs_dir, "sessionInfo_scrna_inventory.txt"), useBytes = TRUE)

cat(obj_state, "\n", sep = "")
cat(normalizePath(manifest_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")
