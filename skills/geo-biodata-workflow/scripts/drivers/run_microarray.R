#!/usr/bin/env Rscript
#
# run_microarray.R
# Manifest-driven limma DE driver for microarray_series_matrix route.
# Handles: Series Matrix files (.txt/.txt.gz), ExpressionSet (RDS), tabular intensity matrices.
# Uses knowledge/platform_registry.tsv for annotation lookup.
# Supports: probe-level DE (no mapping), gene-level DE (with mapping), coverage gates.

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

# Find repo root for platform registry
find_repo_root <- function(dirs) {
  for (start_dir in dirs) {
    current <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
    for (i in seq_len(10L)) {
      if (file.exists(file.path(current, "knowledge", "platform_registry.tsv"))) return(current)
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  stop("Could not find knowledge/platform_registry.tsv from: ", paste(dirs, collapse = ", "), call. = FALSE)
}
repo_root <- find_repo_root(c(script_dir, manifest_dir, getwd()))

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

# Platform: from manifest or auto-detect
platform_id <- manifest$input$platform_id %||% ""
probe_map_path <- if (!is.null(manifest$input$probe_map_file) &&
  nzchar(manifest$input$probe_map_file %||% "")) {
  resolve_manifest_path(manifest$input$probe_map_file)
} else ""

tables_dir <- file.path(manifest_dir, "tables")
figures_dir <- file.path(manifest_dir, "figures")
logs_dir <- file.path(manifest_dir, "logs")
for (path in c(tables_dir, figures_dir, logs_dir)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

fallback_events <- init_fallback_events()

# ── Platform registry lookup ──────────────────────────────────────────────────

load_platform_registry <- function(repo_root) {
  reg_path <- file.path(repo_root, "knowledge", "platform_registry.tsv")
  if (!file.exists(reg_path)) {
    stop("Platform registry not found: ", reg_path, call. = FALSE)
  }
  utils::read.delim(reg_path, stringsAsFactors = FALSE, check.names = FALSE)
}

resolve_platform <- function(gpl_id, registry) {
  # Only use registry lookup — no paste0(GPL, ".db") guessing
  if (!nzchar(gpl_id)) return(NULL)

  entry <- registry[registry$gpl_id == gpl_id, , drop = FALSE]
  if (nrow(entry) == 0L) {
    warning("GPL ", gpl_id, " not found in platform_registry.tsv. ",
      "Add it before using this platform for gene-level analysis.", call. = FALSE)
    return(NULL)
  }
  as.list(entry[1L, ])
}

platform_registry <- load_platform_registry(repo_root)
platform_info <- NULL

# Detect platforms from input
detected_platforms <- platform_id %||% ""

if (length(detected_platforms) > 1L) {
  stop("Multi-platform GSE detected: ", paste(detected_platforms, collapse = ", "),
    ". Multi-platform analysis requires per-platform splitting. ",
    "Process each platform separately with its own manifest.",
    call. = FALSE)
}

if (nzchar(detected_platforms)) {
  platform_info <- resolve_platform(detected_platforms, platform_registry)
}

# Write platform resolution
if (!is.null(platform_info)) {
  utils::write.table(
    as.data.frame(platform_info, stringsAsFactors = FALSE),
    file.path(tables_dir, "platform_resolution.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )
}

# ── Input parsing ────────────────────────────────────────────────────────────

parse_series_matrix <- function(path) {
  # Support both .txt and .txt.gz
  con <- if (grepl("\\.gz$", path)) gzfile(path, "r") else file(path, "r")
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
  storage.mode(mat) <- "numeric"
  platform_from_header <- header[["Series_platform_id"]] %||% character()
  list(
    matrix = as.matrix(mat),
    header = header,
    processing_notes = paste(header[["Sample_data_processing"]] %||% "not_provided", collapse = "; "),
    platform_ids = as.character(platform_from_header)
  )
}

parse_expression_set <- function(path) {
  eset <- readRDS(path)
  if (!inherits(eset, "ExpressionSet")) stop("File is not an ExpressionSet.", call. = FALSE)
  mat <- Biobase::exprs(eset)
  list(
    matrix = mat,
    processing_notes = "not_provided",
    platform_ids = as.character(unique(Biobase::annotation(eset)))
  )
}

input_data <- switch(input_type,
  series_matrix = parse_series_matrix(input_path),
  expression_set = parse_expression_set(input_path),
  microarray_intensity = list(matrix = read_matrix_input(input_path),
    processing_notes = "not_provided", platform_ids = character()),
  stop("Unsupported input_type: ", input_type, call. = FALSE)
)

mat <- input_data$matrix
processing_notes <- input_data$processing_notes %||% "not_provided"

# Detect platforms from input if not explicitly set
detected_from_input <- input_data$platform_ids %||% character()
if (!nzchar(detected_platforms) && length(detected_from_input) > 0L) {
  if (length(detected_from_input) > 1L) {
    stop("Multi-platform input: ", paste(detected_from_input, collapse = ", "),
      ". Split by platform before analysis.", call. = FALSE)
  }
  detected_platforms <- detected_from_input[[1L]]
  if (nzchar(detected_platforms) && is.null(platform_info)) {
    platform_info <- resolve_platform(detected_platforms, platform_registry)
  }
}

writeLines(c("Sample data processing notes from input:", processing_notes),
  file.path(logs_dir, "processing_notes.txt"), useBytes = TRUE)

message(sprintf("Input matrix: %d probes x %d samples", nrow(mat), ncol(mat)))

# ── Probe-to-gene mapping ────────────────────────────────────────────────────

probe_to_gene <- function(mat, platform_info, probe_map_file, fallback_events) {
  probe_gene_map <- NULL
  mapping_method <- "none"
  annotation_pkg <- platform_info$annotation_package %||% ""
  mapping_status_info <- platform_info$mapping_status %||% "unknown"

  if (identical(mapping_status_info, "no_mapping_needed")) {
    return(list(matrix = mat, mapping_method = "none_needed",
      coverage_pct = 100, warnings = character()))
  }

  if (!nzchar(annotation_pkg) && !nzchar(probe_map_file)) {
    fallback_events <- log_fallback(fallback_events,
      "probe_to_gene", "probe_id_used",
      "platform annotation missing",
      "gene_symbol_mapping", "probe_id_as_feature",
      "DE results use probe IDs; gene-level interpretation is not supported.",
      requires_review = TRUE)
    return(list(matrix = mat, mapping_method = "probe_id_only",
      coverage_pct = 0, warnings = "No gene-level annotation available. Results are probe-level.",
      fallback_events = fallback_events))
  }

  if (nzchar(probe_map_file) && file.exists(probe_map_file)) {
    probe_gene_map <- utils::read.delim(probe_map_file, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(probe_gene_map)[1:2] <- c("probe_id", "gene_symbol")
    mapping_method <- "user_provided_map"
  } else if (nzchar(annotation_pkg)) {
    if (!requireNamespace(annotation_pkg, quietly = TRUE)) {
      fallback_events <- log_fallback(fallback_events,
        "probe_to_gene", "annotation_pkg_missing",
        paste("package", annotation_pkg, "not installed"),
        "Bioconductor annotation package", "probe_id_as_feature",
        paste("Annotation package", annotation_pkg, "not available. Results are probe-level."),
        requires_review = TRUE)
      return(list(matrix = mat, mapping_method = "probe_id_only",
        coverage_pct = 0, warnings = paste("Package", annotation_pkg, "not installed."),
        fallback_events = fallback_events))
    }
    probe_gene_map <- AnnotationDbi::select(
      getExportedValue(annotation_pkg, annotation_pkg),
      keys = rownames(mat),
      columns = platform_info$gene_column %||% "SYMBOL",
      keytype = platform_info$probe_keytype %||% "PROBEID"
    )
    colnames(probe_gene_map)[1:2] <- c("probe_id", "gene_symbol")
    mapping_method <- paste0("bioc:", annotation_pkg)
  }

  probe_gene_map <- probe_gene_map[!is.na(probe_gene_map$gene_symbol) &
    nzchar(probe_gene_map$gene_symbol) & probe_gene_map$gene_symbol != "NA", ]

  common_probes <- intersect(rownames(mat), probe_gene_map$probe_id)
  coverage_pct <- round(length(common_probes) / nrow(mat) * 100, 1)

  if (coverage_pct < 50) {
    fallback_events <- log_fallback(fallback_events,
      "probe_to_gene", "low_coverage",
      paste("only", coverage_pct, "% probes mapped to genes"),
      "gene_symbol_mapping", "probe_id_as_feature",
      paste("Gene mapping coverage too low (", coverage_pct, "%)."),
      requires_review = TRUE)
    return(list(matrix = mat, mapping_method = "probe_id_only",
      coverage_pct = coverage_pct,
      warnings = paste("Low mapping coverage:", coverage_pct, "%."),
      fallback_events = fallback_events))
  }

  mat_mapped <- mat[common_probes, , drop = FALSE]
  gene_symbols <- probe_gene_map$gene_symbol[match(common_probes, probe_gene_map$probe_id)]

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

  # Save probe mapping
  probe_table <- data.frame(
    probe_id = common_probes,
    gene_symbol = gene_symbols,
    probe_variance = probe_var,
    kept_in_aggregation = seq_along(gene_symbols) %in%
      unlist(lapply(gene_list, function(idx) idx[which.max(probe_var[idx])])),
    stringsAsFactors = FALSE
  )
  utils::write.table(probe_table, file.path(tables_dir, "probe_mapping.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "")

  # Write mapping coverage
  utils::write.table(
    data.frame(total_probes = nrow(mat), mapped_probes = length(common_probes),
      mapped_genes = nrow(gene_mat), coverage_pct = coverage_pct,
      mapping_method = mapping_method,
      stringsAsFactors = FALSE),
    file.path(tables_dir, "mapping_coverage.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = ""
  )

  list(matrix = gene_mat, mapping_method = mapping_method,
    coverage_pct = coverage_pct, warnings = character(),
    fallback_events = fallback_events)
}

mapping_result <- probe_to_gene(mat, platform_info, probe_map_path, fallback_events)
mat <- mapping_result$matrix
fallback_events <- mapping_result$fallback_events %||% fallback_events

# Determine mapping level
gene_level_de <- !identical(mapping_result$mapping_method, "probe_id_only")
mapping_status_label <- if (gene_level_de) "GENE_LEVEL_ANALYSIS_COMPLETE" else "PROBE_LEVEL_ANALYSIS_COMPLETE"
if (mapping_result$coverage_pct < 50) {
  mapping_status_label <- "GENE_MAPPING_REVIEW_REQUIRED"
}

# ── Sample alignment ─────────────────────────────────────────────────────────

sample_map <- read_sample_mapping(sample_path, sample_id_col, contrast_factor)
aligned <- align_samples(mat, sample_map, sample_id_col)
mat <- aligned$matrix
sample_map <- aligned$metadata

utils::write.table(sample_map, file.path(tables_dir, "sample_mapping_used.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

# ── Low-expression filter ────────────────────────────────────────────────────

mat_filtered <- filter_low_expression(mat)
message(sprintf("Features after filter: %d / %d", nrow(mat_filtered), nrow(mat)))

# ── limma DE ─────────────────────────────────────────────────────────────────

fit_result <- run_limma_de(mat_filtered, sample_map, design_formula, contrast)
result_df <- fit_result$result

outputs <- write_limma_outputs(result_df, fit_result$ebayes_fit,
  fit_result$design, fit_result$contrast_matrix,
  contrast_factor, contrast$numerator, contrast$denominator,
  sample_map, mat_filtered, tables_dir, figures_dir, logs_dir)

# ── QC checks ────────────────────────────────────────────────────────────────

qc <- run_limma_qc(fit_result, mat_filtered, sample_map, contrast_factor)
utils::write.table(qc$qc_table, file.path(tables_dir, "qc_checks.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

if (nrow(fallback_events) > 0L) {
  utils::write.table(fallback_events, file.path(tables_dir, "fallback_events.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

# ── Status ───────────────────────────────────────────────────────────────────

critical_outputs <- c(
  file.path(tables_dir, "sample_mapping_used.tsv"),
  file.path(tables_dir, "design_matrix_used.tsv"),
  file.path(tables_dir, "contrast_matrix_used.tsv"),
  file.path(tables_dir, "factor_levels_used.tsv"),
  file.path(tables_dir, "library_sizes.tsv"),
  file.path(tables_dir, "pca_coordinates.tsv"),
  file.path(tables_dir, "qc_checks.tsv"),
  file.path(tables_dir, "platform_resolution.tsv"),
  file.path(tables_dir, "mapping_coverage.tsv"),
  outputs["de_path"],
  file.path(figures_dir, "bulk_library_sizes.pdf"),
  file.path(figures_dir, "bulk_pca.pdf"),
  file.path(figures_dir, paste0("bulk_pvalue_histogram_", outputs["contrast_name"], ".pdf")),
  file.path(figures_dir, paste0("bulk_meanvar_", outputs["contrast_name"], ".pdf"))
)

status_out <- determine_limma_status(critical_outputs, qc, fallback_events)
status_out$note <- paste(status_out$note, "| mapping:", mapping_status_label)

status <- data.frame(
  execution_state = status_out$execution_state,
  technical_qc = status_out$technical_qc,
  result_signal = status_out$result_signal,
  mapping_status = mapping_status_label,
  updated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  note = status_out$note,
  stringsAsFactors = FALSE
)
utils::write.table(status, file.path(manifest_dir, "workflow_status.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

session_lines <- utils::capture.output(utils::sessionInfo())
writeLines(session_lines, file.path(logs_dir, "sessionInfo_microarray.txt"), useBytes = TRUE)

cat(status_out$execution_state, "\n", sep = "")
cat(normalizePath(manifest_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")
