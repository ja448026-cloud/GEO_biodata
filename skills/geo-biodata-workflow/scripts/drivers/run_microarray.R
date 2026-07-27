#!/usr/bin/env Rscript
#
# run_microarray.R
# Manifest-driven limma DE driver for microarray_series_matrix route.
# Handles: Series Matrix files, ExpressionSet (RDS), tabular intensity matrices.
# Performs: platform annotation, probe-to-gene mapping, duplicate probe aggregation,
#           then limma DE via bulk_limma_common.R.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: run_microarray.R /path/to/run_manifest.yaml", call. = FALSE)
}

manifest_path <- args[[1L]]
if (!file.exists(manifest_path)) stop("Manifest file does not exist: ", manifest_path, call. = FALSE)

required <- c("yaml", "limma", "ggplot2", "Biobase")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages({ library(limma) })

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(file.path("skills", "geo-biodata-workflow", "scripts", "drivers", "run_microarray.R"), mustWork = TRUE)
}
driver_dir <- dirname(script_path)
script_dir <- normalizePath(file.path(driver_dir, ".."), winslash = "/", mustWork = TRUE)
manifest_dir <- dirname(normalizePath(manifest_path, winslash = "/", mustWork = TRUE))

# Source shared limma functions
limma_common <- file.path(script_dir, "bulk_limma_common.R")
if (!file.exists(limma_common)) stop("Could not find bulk_limma_common.R.", call. = FALSE)
source(limma_common)

# Validate manifest
validator <- file.path(script_dir, "validate_manifest.R")
if (!file.exists(validator)) stop("Could not find validate_manifest.R.", call. = FALSE)
validation_output <- system2("Rscript", c(shQuote(validator), shQuote(manifest_path)),
  stdout = TRUE, stderr = TRUE)
validation_status <- attr(validation_output, "status") %||% 0L
if (!identical(as.integer(validation_status), 0L)) {
  cat(paste(validation_output, collapse = "\n"), "\n")
  stop("Manifest validation failed; microarray analysis did not run.", call. = FALSE)
}

manifest <- yaml::read_yaml(manifest_path)
if (!identical(manifest$route, "microarray_series_matrix")) {
  stop("run_microarray.R only supports route: microarray_series_matrix", call. = FALSE)
}

allowed_inputs <- c("expression_set", "series_matrix", "microarray_intensity")
input_type <- manifest$input$input_type %||% ""
if (!input_type %in% allowed_inputs) {
  stop("microarray_series_matrix route requires input_type in: ", paste(allowed_inputs, collapse = ", "),
    ", got: ", input_type, call. = FALSE)
}

resolve_manifest_path <- function(path) {
  if (!nzchar(path %||% "")) return("")
  if (grepl("^[A-Za-z]:[\\\\/]|^/", path)) return(path)
  file.path(manifest_dir, path)
}

input_path <- resolve_manifest_path(manifest$input$file)
sample_path <- resolve_manifest_path(manifest$sample_mapping$file)
sample_id_col <- manifest$sample_mapping$sample_id_column
contrast <- manifest$design$contrast
contrast_factor <- contrast$factor
design_formula <- stats::as.formula(manifest$design$formula)

# Optional: platform annotation path in manifest
platform_id <- manifest$input$platform_id %||% ""
probe_map_path <- if (!is.null(manifest$input$probe_map_file) && nzchar(manifest$input$probe_map_file %||% "")) {
  resolve_manifest_path(manifest$input$probe_map_file)
} else {
  ""
}

tables_dir <- file.path(manifest_dir, "tables")
figures_dir <- file.path(manifest_dir, "figures")
logs_dir <- file.path(manifest_dir, "logs")
for (path in c(tables_dir, figures_dir, logs_dir)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# ── Input parsing ────────────────────────────────────────────────────────────

parse_series_matrix <- function(path) {
  # Parse a GEO Series Matrix file: header lines starting with ! then data
  con <- file(path, "r")
  on.exit(close(con))
  header <- list()
  lines <- readLines(con, warn = FALSE)
  data_start <- which(grepl("^!series_matrix_table_begin", lines))
  data_end <- which(grepl("^!series_matrix_table_end", lines))
  if (length(data_start) == 0L) stop("Series matrix missing '!series_matrix_table_begin' marker.", call. = FALSE)
  for (i in seq_len(data_start - 1L)) {
    line <- lines[i]
    if (grepl("^!", line)) {
      parts <- strsplit(sub("^!", "", line), "\t")[[1L]]
      key <- trimws(parts[[1L]])
      vals <- if (length(parts) > 1L) trimws(gsub('"', '', parts[-1L])) else ""
      header[[key]] <- vals
    }
  }
  data_lines <- lines[(data_start + 1L):(data_end - 1L)]
  if (length(data_lines) == 0L) stop("Series matrix data section is empty.", call. = FALSE)
  tf <- textConnection(data_lines)
  mat <- utils::read.delim(tf, check.names = FALSE, stringsAsFactors = FALSE)
  close(tf)
  feature_ids <- as.character(mat[[1L]])
  if (any(duplicated(feature_ids))) {
    warning("Duplicate probe IDs in series matrix; keeping first occurrence.", call. = FALSE)
    mat <- mat[!duplicated(feature_ids), ]
    feature_ids <- feature_ids[!duplicated(feature_ids)]
  }
  rownames(mat) <- feature_ids
  mat <- mat[, -1L, drop = FALSE]
  # Extract sample data processing notes
  processing_notes <- header[["Sample_data_processing"]] %||% "not_provided"
  # Extract platform IDs
  platform_ids <- header[["Series_platform_id"]] %||% character()
  storage.mode(mat) <- "numeric"
  list(
    matrix = as.matrix(mat),
    header = header,
    processing_notes = processing_notes,
    platform_ids = platform_ids
  )
}

parse_expression_set <- function(path) {
  if (!requireNamespace("Biobase", quietly = TRUE)) stop("Biobase is required for ExpressionSet input.", call. = FALSE)
  eset <- readRDS(path)
  if (!inherits(eset, "ExpressionSet")) stop("File is not an ExpressionSet object.", call. = FALSE)
  mat <- Biobase::exprs(eset)
  pdata <- Biobase::pData(eset)
  fdata <- Biobase::fData(eset)
  list(
    matrix = mat,
    phenotype_data = pdata,
    feature_data = fdata,
    processing_notes = if ("data_processing" %in% names(Biobase::experimentData(eset)@other)) {
      Biobase::experimentData(eset)@other[["data_processing"]]
    } else {
      "not_provided"
    },
    platform_ids = as.character(unique(Biobase::annotation(eset)))
  )
}

# Load input
input_data <- switch(input_type,
  series_matrix = parse_series_matrix(input_path),
  expression_set = parse_expression_set(input_path),
  microarray_intensity = list(matrix = read_matrix_input(input_path), processing_notes = "not_provided", platform_ids = character()),
  stop("Unsupported input_type: ", input_type, call. = FALSE)
)

mat <- input_data$matrix
processing_notes <- input_data$processing_notes %||% "not_provided"
platform_ids <- input_data$platform_ids %||% character()

# Log processing notes for audit
if (length(processing_notes) > 0L) {
  writeLines(c("Sample data processing notes from input:", processing_notes),
    file.path(logs_dir, "processing_notes.txt"), useBytes = TRUE)
}

message(sprintf("Input matrix: %d probes x %d samples", nrow(mat), ncol(mat)))
if (length(platform_ids) > 0L) {
  message("Platform(s): ", paste(platform_ids, collapse = ", "))
}

# ── Platform annotation and probe-to-gene mapping ────────────────────────────

probe_to_gene <- function(mat, platform, probe_map_file = "") {
  # Returns gene-level matrix. Priority: user file > Bioc annotation package > fail
  probe_gene_map <- NULL

  if (nzchar(probe_map_file) && file.exists(probe_map_file)) {
    probe_gene_map <- utils::read.delim(probe_map_file, stringsAsFactors = FALSE, check.names = FALSE)
    if (ncol(probe_gene_map) < 2L) stop("Probe map file must have at least 2 columns (probe_id, gene_symbol).", call. = FALSE)
    colnames(probe_gene_map)[1:2] <- c("probe_id", "gene_symbol")
  } else if (nzchar(platform) && length(platform) > 0L) {
    platform_pkg <- paste0(platform, ".db")
    if (requireNamespace(platform_pkg, quietly = TRUE)) {
      pkg_env <- asNamespace(platform_pkg)
      probe_gene_map <- tryCatch({
        AnnotationDbi::select(
          getExportedValue(platform_pkg, platform_pkg),
          keys = rownames(mat),
          columns = "SYMBOL",
          keytype = "PROBEID"
        )
      }, error = function(e) NULL)
      if (!is.null(probe_gene_map)) {
        colnames(probe_gene_map)[1:2] <- c("probe_id", "gene_symbol")
      }
    }
  }

  if (is.null(probe_gene_map) || nrow(probe_gene_map) == 0L) {
    warning("No probe-to-gene annotation available. Using probe IDs as gene-level features. ",
      "Set platform_id in manifest or provide probe_map_file.", call. = FALSE)
    return(mat)
  }

  probe_gene_map <- probe_gene_map[!is.na(probe_gene_map$gene_symbol) &
    nzchar(probe_gene_map$gene_symbol) &
    probe_gene_map$gene_symbol != "NA", ]

  # Match probes
  common_probes <- intersect(rownames(mat), probe_gene_map$probe_id)
  if (length(common_probes) == 0L) {
    warning("No probe IDs matched. Probe IDs in matrix might differ from annotation keys. Using probe IDs as features.", call. = FALSE)
    return(mat)
  }

  message(sprintf("Mapped %d / %d probes to %d unique genes",
    length(common_probes), nrow(mat),
    length(unique(probe_gene_map$gene_symbol[probe_gene_map$probe_id %in% common_probes]))))

  mat_mapped <- mat[common_probes, , drop = FALSE]
  gene_symbols <- probe_gene_map$gene_symbol[match(common_probes, probe_gene_map$probe_id)]

  # Aggregate multi-probe genes: take the probe with maximum variance (most informative)
  dup_genes <- unique(gene_symbols[duplicated(gene_symbols)])
  if (length(dup_genes) > 0L) {
    message(sprintf("Aggregating %d genes with multiple probes (max-variance method)", length(dup_genes)))
  }

  probe_var <- apply(mat_mapped, 1L, stats::var, na.rm = TRUE)
  gene_list <- split(seq_len(nrow(mat_mapped)), gene_symbols)
  gene_mat <- do.call(rbind, lapply(gene_list, function(idx) {
    if (length(idx) == 1L) return(mat_mapped[idx, , drop = FALSE])
    best <- idx[which.max(probe_var[idx])]
    mat_mapped[best, , drop = FALSE]
  }))
  rownames(gene_mat) <- names(gene_list)

  # Record probe mapping table
  probe_table <- data.frame(
    probe_id = common_probes,
    gene_symbol = gene_symbols,
    probe_variance = probe_var,
    kept_in_aggregation = !duplicated(gene_symbols) & !gene_symbols %in% dup_genes |
      seq_along(gene_symbols) %in% unlist(lapply(gene_list, function(idx) idx[which.max(probe_var[idx])])),
    stringsAsFactors = FALSE
  )
  utils::write.table(probe_table, file.path(tables_dir, "probe_mapping.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  gene_mat
}

platform_to_use <- if (nzchar(platform_id)) platform_id else if (length(platform_ids) > 0L) platform_ids[[1L]] else ""
mat <- probe_to_gene(mat, platform_to_use, probe_map_path)

# ── Sample alignment ─────────────────────────────────────────────────────────

sample_map <- read_sample_mapping(sample_path, sample_id_col, contrast_factor)
aligned <- align_samples(mat, sample_map, sample_id_col)
mat <- aligned$matrix
sample_map <- aligned$metadata

utils::write.table(sample_map, file.path(tables_dir, "sample_mapping_used.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# ── Low-expression filter ────────────────────────────────────────────────────

mat_filtered <- filter_low_expression(mat)
message(sprintf("Genes after low-expression filter: %d / %d", nrow(mat_filtered), nrow(mat)))

# ── limma DE ─────────────────────────────────────────────────────────────────

fit_result <- run_limma_de(mat_filtered, sample_map, design_formula, contrast)
result_df <- fit_result$result

outputs <- write_limma_outputs(result_df, fit_result$ebayes_fit, contrast_factor,
  contrast$numerator, contrast$denominator,
  sample_map, mat_filtered, tables_dir, figures_dir, logs_dir)

# ── QC checks ────────────────────────────────────────────────────────────────

qc <- run_limma_qc(fit_result, mat_filtered, sample_map, contrast_factor)
utils::write.table(qc$qc_table, file.path(tables_dir, "qc_checks.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# ── Status ───────────────────────────────────────────────────────────────────

critical_outputs <- c(
  file.path(tables_dir, "sample_mapping_used.tsv"),
  file.path(tables_dir, "design_matrix.tsv"),
  file.path(tables_dir, "library_sizes.tsv"),
  file.path(tables_dir, "pca_coordinates.tsv"),
  file.path(tables_dir, "qc_checks.tsv"),
  outputs["de_path"],
  file.path(figures_dir, "bulk_library_sizes.pdf"),
  file.path(figures_dir, "bulk_pca.pdf"),
  file.path(figures_dir, paste0("bulk_pvalue_histogram_", outputs["contrast_name"], ".pdf")),
  file.path(figures_dir, paste0("bulk_meanvar_", outputs["contrast_name"], ".pdf"))
)

n_total <- nrow(result_df)
n_de <- sum(result_df$adj.P.Val < 0.05, na.rm = TRUE)
status_out <- determine_limma_status(critical_outputs, qc$qc_flags, n_total, n_de)

status <- data.frame(
  state = status_out$state,
  updated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  note = status_out$note,
  stringsAsFactors = FALSE
)
utils::write.table(status, file.path(manifest_dir, "workflow_status.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

session_lines <- utils::capture.output(utils::sessionInfo())
writeLines(session_lines, file.path(logs_dir, "sessionInfo_microarray.txt"), useBytes = TRUE)

cat(status_out$state, "\n", sep = "")
cat(normalizePath(manifest_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")
