#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L || length(args) > 2L) {
  stop("Usage: generate_manifest_draft.R <intake_handoff.yaml> [output-yaml]", call. = FALSE)
}
if (!requireNamespace("yaml", quietly = TRUE)) stop("Missing required package: yaml", call. = FALSE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
handoff_path <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
handoff <- yaml::read_yaml(handoff_path)
run_dir <- handoff$run_dir
out_path <- if (length(args) == 2L) args[[2L]] else file.path(run_dir, "run_manifest.draft.yaml")

route <- handoff$selected_input$route
input_type <- handoff$selected_input$input_type
if (identical(input_type, "tpm_fpkm_cpm_or_unlogged_normalized")) input_type <- "normalized_expression"
if (identical(input_type, "log_or_transformed_expression")) input_type <- "log_normalized"
if (identical(input_type, "raw_integer_counts")) input_type <- "raw_integer_counts"

manifest <- list(
  accession = handoff$accession,
  route = route,
  input = list(
    file = handoff$selected_input$path,
    input_type = input_type,
    location_type = "local",
    species = handoff$selected_input$species %||% "",
    scale = handoff$scale
  ),
  sample_mapping = list(
    file = handoff$sample_mapping$file %||% "",
    biological_unit = handoff$sample_mapping$biological_unit %||% "",
    sample_id_column = handoff$sample_mapping$sample_id_column %||% "",
    group_column = handoff$sample_mapping$group_column %||% ""
  ),
  design = list(
    formula = handoff$design_readiness$suggested_formula %||% "",
    contrast = handoff$design_readiness$suggested_contrast
  ),
  review = list(
    mode = "manual",
    input_type_confirmed = FALSE,
    sample_mapping_confirmed = FALSE,
    design_confirmed = FALSE,
    contrast_confirmed = FALSE
  )
)

yaml::write_yaml(manifest, out_path)
cat("MANIFEST_DRAFT_WRITTEN\n", normalizePath(out_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
