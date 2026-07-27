#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: validate_manifest.R /path/to/run_manifest.yaml", call. = FALSE)
}

manifest_path <- args[[1L]]
if (!file.exists(manifest_path)) {
  stop("Manifest file does not exist: ", manifest_path, call. = FALSE)
}
if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Missing required package: yaml", call. = FALSE)
}

manifest <- yaml::read_yaml(manifest_path)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

find_schema <- function(starts) {
  for (start in starts) {
    current <- normalizePath(start, winslash = "/", mustWork = FALSE)
    for (i in seq_len(8L)) {
      candidate <- file.path(current, "schemas", "run_manifest.schema.yaml")
      if (file.exists(candidate)) return(candidate)
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  NA_character_
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)) else getwd()
schema_path <- find_schema(c(getwd(), dirname(normalizePath(manifest_path, mustWork = TRUE)), script_dir))
if (is.na(schema_path)) {
  stop("Could not find schemas/run_manifest.schema.yaml from manifest, script, or working directory.", call. = FALSE)
}
schema <- yaml::read_yaml(schema_path)

errors <- character()
warnings <- character()
add_error <- function(...) errors <<- c(errors, paste0(...))
add_warning <- function(...) warnings <<- c(warnings, paste0(...))

is_true <- function(x) isTRUE(x) || identical(tolower(as.character(x)), "true")
has_field <- function(x, field) !is.null(x[[field]]) && length(x[[field]]) > 0L && !identical(x[[field]], "")

for (field in schema$required_top_level) {
  if (!has_field(manifest, field)) add_error("Missing top-level field: ", field)
}

accession <- manifest$accession %||% ""
if (!grepl("^GSE[0-9]+$", toupper(accession))) {
  add_error("accession must be a GEO Series accession like GSE123456.")
}

route <- manifest$route %||% ""
allowed_routes <- names(schema$allowed_routes)
if (!route %in% allowed_routes) {
  add_error("route must be one of: ", paste(allowed_routes, collapse = ", "))
} else {
  input_type <- manifest$input$input_type %||% ""
  allowed_input_types <- schema$allowed_routes[[route]]$input_types
  if (!input_type %in% allowed_input_types) {
    add_error(
      "input.input_type '", input_type, "' is not compatible with route '", route,
      "'. Allowed: ", paste(allowed_input_types, collapse = ", ")
    )
  }
}

for (field in schema$required_input_fields) {
  if (!has_field(manifest$input, field)) add_error("Missing input field: input.", field)
}
for (field in schema$required_sample_mapping_fields) {
  if (!has_field(manifest$sample_mapping, field)) add_error("Missing sample_mapping field: sample_mapping.", field)
}
for (field in schema$required_design_fields) {
  if (!has_field(manifest$design, field)) add_error("Missing design field: design.", field)
}
for (field in schema$required_contrast_fields) {
  if (!has_field(manifest$design$contrast, field)) add_error("Missing contrast field: design.contrast.", field)
}
for (field in schema$required_review_flags) {
  if (!has_field(manifest$review, field) || !is_true(manifest$review[[field]])) {
    add_error("Review gate is not confirmed: review.", field)
  }
}

manifest_dir <- dirname(normalizePath(manifest_path, winslash = "/", mustWork = TRUE))
resolve_manifest_path <- function(path) {
  if (grepl("^[A-Za-z]:[\\\\/]|^/", path)) return(path)
  file.path(manifest_dir, path)
}

input_file <- manifest$input$file %||% ""
sample_file <- manifest$sample_mapping$file %||% ""
if (nzchar(input_file) && !file.exists(resolve_manifest_path(input_file))) {
  add_warning("Input file not found relative to manifest: ", input_file)
}
sample_path <- if (nzchar(sample_file)) resolve_manifest_path(sample_file) else ""
sample_map <- NULL
if (nzchar(sample_path) && file.exists(sample_path)) {
  sample_map <- utils::read.delim(sample_path, stringsAsFactors = FALSE, check.names = FALSE)
  sample_id_col <- manifest$sample_mapping$sample_id_column %||% ""
  group_col <- manifest$sample_mapping$group_column %||% ""
  design_formula <- manifest$design$formula %||% ""
  contrast <- manifest$design$contrast

  for (col in c(sample_id_col, group_col)) {
    if (!col %in% names(sample_map)) add_error("Sample mapping column not found: ", col)
  }
  if (sample_id_col %in% names(sample_map)) {
    ids <- sample_map[[sample_id_col]]
    if (any(is.na(ids) | !nzchar(as.character(ids)))) add_error("Sample IDs contain missing or empty values.")
    if (any(duplicated(ids))) add_error("Sample IDs are duplicated.")
  }
  formula_vars <- all.vars(stats::as.formula(design_formula))
  missing_formula_vars <- setdiff(formula_vars, names(sample_map))
  if (length(missing_formula_vars) > 0L) {
    add_error("Design formula fields missing from sample mapping: ", paste(missing_formula_vars, collapse = ", "))
  }
  if (length(missing_formula_vars) == 0L && length(formula_vars) > 0L) {
    design_matrix <- tryCatch(
      stats::model.matrix(stats::as.formula(design_formula), data = sample_map),
      error = function(e) {
        add_error("Design matrix could not be built: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(design_matrix) && qr(design_matrix)$rank < ncol(design_matrix)) {
      add_error("Design matrix is not full rank.")
    }
  }
  if (!is.null(contrast)) {
    contrast_factor <- contrast$factor %||% group_col
    if (!contrast_factor %in% names(sample_map)) {
      add_error("Contrast factor missing from sample mapping: ", contrast_factor)
    } else {
      levels_present <- unique(as.character(sample_map[[contrast_factor]]))
      if (!contrast$numerator %in% levels_present) add_error("Contrast numerator not present: ", contrast$numerator)
      if (!contrast$denominator %in% levels_present) add_error("Contrast denominator not present: ", contrast$denominator)
    }
  }
} else if (nzchar(sample_file)) {
  add_warning("Sample mapping file not found relative to manifest: ", sample_file)
}

out_dir <- dirname(manifest_path)
status <- if (length(errors) == 0L) "MANIFEST_VALIDATED" else "MANIFEST_INVALID"
validation <- data.frame(
  status = c(rep("ERROR", length(errors)), rep("WARNING", length(warnings))),
  message = c(errors, warnings),
  stringsAsFactors = FALSE
)
if (nrow(validation) == 0L) {
  validation <- data.frame(status = "OK", message = "Manifest validation passed.", stringsAsFactors = FALSE)
}
utils::write.table(
  validation,
  file.path(out_dir, "manifest_validation.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

workflow_status <- data.frame(
  state = status,
  updated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  note = if (length(errors) == 0L) "Manifest structure and review gates passed." else paste(errors, collapse = " | "),
  stringsAsFactors = FALSE
)
utils::write.table(
  workflow_status,
  file.path(out_dir, "workflow_status.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

cat(status, "\n", sep = "")
cat(normalizePath(file.path(out_dir, "manifest_validation.tsv"), winslash = "/", mustWork = TRUE), "\n")
if (length(errors) > 0L) quit(status = 1L)
