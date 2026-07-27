#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
this_file <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(file.path("validation", "run_smoke_checks.R"), mustWork = TRUE)
}
repo_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = TRUE)
script_dir <- file.path(repo_root, "skills", "geo-biodata-workflow", "scripts")

fail <- function(message) {
  stop(message, call. = FALSE)
}

expect_status <- function(label, command, args, expected_status) {
  status <- suppressWarnings(system2(command, shQuote(args), stdout = TRUE, stderr = TRUE))
  exit_status <- attr(status, "status") %||% 0L
  if (!identical(as.integer(exit_status), as.integer(expected_status))) {
    cat(paste(status, collapse = "\n"), "\n")
    fail(sprintf("%s: expected status %d, got %d", label, expected_status, exit_status))
  }
  invisible(status)
}

cat("== Parse R scripts ==\n")
r_files <- list.files(script_dir, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
for (path in sort(r_files)) {
  parse(file = path)
  cat("PARSE_OK\t", basename(path), "\n", sep = "")
}

cat("== Download guard checks ==\n")
scratch <- file.path(tempdir(), paste0("geo_biodata_smoke_", Sys.getpid()))
dir.create(file.path(scratch, "resources"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(scratch, "raw"), recursive = TRUE, showWarnings = FALSE)

supplement_index <- data.frame(
  supplement_url = c("not_a_url", "not_a_url"),
  file_name = c("GSE000000_counts.tsv.gz", "GSE000000_notes.txt"),
  size = c(100, 20),
  stringsAsFactors = FALSE
)
utils::write.table(
  supplement_index,
  file.path(scratch, "resources", "supplement_index.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

expect_status(
  "Broad download-plan regex refusal",
  "Rscript",
  c(file.path(script_dir, "generate_download_plan.R"), scratch, ".*"),
  1L
)

expect_status(
  "Specific download-plan generation",
  "Rscript",
  c(file.path(script_dir, "generate_download_plan.R"), scratch, "counts", "count matrix candidate"),
  0L
)

plan_path <- file.path(scratch, "plans", "download_plan.tsv")
plan <- utils::read.delim(plan_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!any(tolower(plan$selected) %in% c("true", "1", "yes"))) {
  fail("Generated plan did not select the expected counts file.")
}

expect_status(
  "Unreviewed plan download refusal",
  "Rscript",
  c(file.path(script_dir, "download_geo_supp.R"), plan_path, file.path(scratch, "raw")),
  1L
)

selected_mask <- tolower(as.character(plan$selected)) %in% c("true", "1", "yes")
plan$reviewed[selected_mask] <- TRUE
utils::write.table(plan, plan_path, sep = "\t", quote = FALSE, row.names = FALSE)
expect_status(
  "Invalid URL download refusal",
  "Rscript",
  c(file.path(script_dir, "download_geo_supp.R"), plan_path, file.path(scratch, "raw")),
  1L
)

expect_status(
  "Legacy broad regex refusal",
  "Rscript",
  c(file.path(script_dir, "download_geo_supp.R"), "GSE000000", file.path(scratch, "raw"), ".*"),
  1L
)

cat("== Manifest validation checks ==\n")
manifest_fixture <- file.path(repo_root, "validation", "fixtures", "manifest_valid", "run_manifest.yaml")
expect_status(
  "Valid manifest",
  "Rscript",
  c(file.path(script_dir, "validate_manifest.R"), manifest_fixture),
  0L
)

expect_status(
  "Metadata-only manifest",
  "Rscript",
  c(file.path(script_dir, "validate_manifest.R"), file.path(repo_root, "validation", "fixtures", "manifest_metadata_only", "run_manifest.yaml")),
  0L
)

expect_status(
  "scRNA author-object manifest",
  "Rscript",
  c(file.path(script_dir, "validate_manifest.R"), file.path(repo_root, "validation", "fixtures", "manifest_scrna_author_object", "run_manifest.yaml")),
  0L
)

expect_status(
  "Normalized bulk manifest",
  "Rscript",
  c(file.path(script_dir, "validate_manifest.R"), file.path(repo_root, "validation", "fixtures", "manifest_bulk_normalized", "run_manifest.yaml")),
  0L
)

expect_status(
  "Microarray manifest",
  "Rscript",
  c(file.path(script_dir, "validate_manifest.R"), file.path(repo_root, "validation", "fixtures", "manifest_microarray", "run_manifest.yaml")),
  0L
)

invalid_dir <- file.path(scratch, "manifest_invalid")
dir.create(invalid_dir, recursive = TRUE, showWarnings = FALSE)
invalid_manifest <- file.path(invalid_dir, "run_manifest.yaml")
valid_text <- readLines(manifest_fixture, warn = FALSE)
invalid_text <- sub("input_type_confirmed: true", "input_type_confirmed: false", valid_text, fixed = TRUE)
writeLines(invalid_text, invalid_manifest, useBytes = TRUE)
dir.create(file.path(invalid_dir, "resources"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(invalid_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
invisible(file.copy(
  file.path(dirname(manifest_fixture), "resources", "sample_mapping_reviewed.tsv"),
  file.path(invalid_dir, "resources", "sample_mapping_reviewed.tsv"),
  overwrite = TRUE
))
invisible(file.copy(
  file.path(dirname(manifest_fixture), "raw", "counts.tsv"),
  file.path(invalid_dir, "raw", "counts.tsv"),
  overwrite = TRUE
))
expect_status(
  "Unconfirmed manifest review gate refusal",
  "Rscript",
  c(file.path(script_dir, "validate_manifest.R"), invalid_manifest),
  1L
)

missing_input_dir <- file.path(scratch, "manifest_missing_input")
dir.create(missing_input_dir, recursive = TRUE, showWarnings = FALSE)
missing_input_manifest <- file.path(missing_input_dir, "run_manifest.yaml")
missing_input_text <- sub("file: raw/counts.tsv", "file: raw/missing_counts.tsv", valid_text, fixed = TRUE)
writeLines(missing_input_text, missing_input_manifest, useBytes = TRUE)
dir.create(file.path(missing_input_dir, "resources"), recursive = TRUE, showWarnings = FALSE)
invisible(file.copy(
  file.path(dirname(manifest_fixture), "resources", "sample_mapping_reviewed.tsv"),
  file.path(missing_input_dir, "resources", "sample_mapping_reviewed.tsv"),
  overwrite = TRUE
))
expect_status(
  "Missing local input refusal",
  "Rscript",
  c(file.path(script_dir, "validate_manifest.R"), missing_input_manifest),
  1L
)

cat("== Dependency bootstrap plan ==\n")
expect_status(
  "Core dependency plan",
  "Rscript",
  c(file.path(script_dir, "bootstrap_environment.R"), "--profile", "core", "--plan"),
  0L
)

cat("== Bulk driver optional check ==\n")
if (all(vapply(c("DESeq2", "ggplot2", "SummarizedExperiment"), requireNamespace, logical(1), quietly = TRUE))) {
  bulk_dir <- file.path(scratch, "bulk_driver")
  dir.create(file.path(bulk_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(bulk_dir, "resources"), recursive = TRUE, showWarnings = FALSE)
  invisible(file.copy(manifest_fixture, file.path(bulk_dir, "run_manifest.yaml"), overwrite = TRUE))
  invisible(file.copy(file.path(dirname(manifest_fixture), "raw", "counts.tsv"), file.path(bulk_dir, "raw", "counts.tsv"), overwrite = TRUE))
  invisible(file.copy(file.path(dirname(manifest_fixture), "resources", "sample_mapping_reviewed.tsv"), file.path(bulk_dir, "resources", "sample_mapping_reviewed.tsv"), overwrite = TRUE))
  expect_status(
    "Bulk raw-count driver",
    "Rscript",
    c(file.path(script_dir, "drivers", "run_bulk_counts.R"), file.path(bulk_dir, "run_manifest.yaml")),
    0L
  )
  status <- utils::read.delim(file.path(bulk_dir, "workflow_status.tsv"), stringsAsFactors = FALSE)
  valid_completion <- c("BASIC_ANALYSIS_COMPLETE", "QC_REVIEW_REQUIRED")
  if (!status$state[[1L]] %in% valid_completion) {
    fail(sprintf("Bulk driver produced unexpected state: %s (expected BASIC_ANALYSIS_COMPLETE or QC_REVIEW_REQUIRED)",
          status$state[[1L]]))
  }
  qc_table <- utils::read.delim(file.path(bulk_dir, "tables", "qc_checks.tsv"), stringsAsFactors = FALSE)
  if (nrow(qc_table) == 0L) fail("Bulk driver did not produce QC checks table.")
} else {
  cat("SKIP\tBulk driver optional check requires DESeq2, ggplot2, and SummarizedExperiment.\n")
}

cat("== Normalized bulk driver optional check ==\n")
if (all(vapply(c("limma", "ggplot2", "yaml"), requireNamespace, logical(1), quietly = TRUE))) {
  norm_dir <- file.path(scratch, "bulk_normalized_driver")
  norm_fixture <- file.path(repo_root, "validation", "fixtures", "manifest_bulk_normalized")
  dir.create(file.path(norm_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(norm_dir, "resources"), recursive = TRUE, showWarnings = FALSE)
  invisible(file.copy(file.path(norm_fixture, "run_manifest.yaml"), file.path(norm_dir, "run_manifest.yaml"), overwrite = TRUE))
  invisible(file.copy(file.path(norm_fixture, "raw", "counts.tsv"), file.path(norm_dir, "raw", "counts.tsv"), overwrite = TRUE))
  invisible(file.copy(file.path(norm_fixture, "resources", "sample_mapping_reviewed.tsv"), file.path(norm_dir, "resources", "sample_mapping_reviewed.tsv"), overwrite = TRUE))
  expect_status(
    "Bulk normalized driver",
    "Rscript",
    c(file.path(script_dir, "drivers", "run_bulk_normalized.R"), file.path(norm_dir, "run_manifest.yaml")),
    0L
  )
  norm_status <- utils::read.delim(file.path(norm_dir, "workflow_status.tsv"), stringsAsFactors = FALSE)
  if (!norm_status$state[[1L]] %in% c("BASIC_ANALYSIS_COMPLETE", "QC_REVIEW_REQUIRED")) {
    fail(sprintf("Bulk normalized driver produced unexpected state: %s", norm_status$state[[1L]]))
  }
} else {
  cat("SKIP\tBulk normalized driver optional check requires limma, ggplot2, and yaml.\n")
}

cat("== Microarray driver optional check ==\n")
if (all(vapply(c("limma", "ggplot2", "yaml", "Biobase"), requireNamespace, logical(1), quietly = TRUE))) {
  ma_dir <- file.path(scratch, "microarray_driver")
  ma_fixture <- file.path(repo_root, "validation", "fixtures", "manifest_microarray")
  dir.create(file.path(ma_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(ma_dir, "resources"), recursive = TRUE, showWarnings = FALSE)
  invisible(file.copy(file.path(ma_fixture, "run_manifest.yaml"), file.path(ma_dir, "run_manifest.yaml"), overwrite = TRUE))
  invisible(file.copy(file.path(ma_fixture, "raw", "array_matrix.tsv"), file.path(ma_dir, "raw", "array_matrix.tsv"), overwrite = TRUE))
  invisible(file.copy(file.path(ma_fixture, "resources", "sample_mapping_reviewed.tsv"), file.path(ma_dir, "resources", "sample_mapping_reviewed.tsv"), overwrite = TRUE))
  expect_status(
    "Microarray driver",
    "Rscript",
    c(file.path(script_dir, "drivers", "run_microarray.R"), file.path(ma_dir, "run_manifest.yaml")),
    0L
  )
  ma_status <- utils::read.delim(file.path(ma_dir, "workflow_status.tsv"), stringsAsFactors = FALSE)
  if (!ma_status$state[[1L]] %in% c("BASIC_ANALYSIS_COMPLETE", "QC_REVIEW_REQUIRED")) {
    fail(sprintf("Microarray driver produced unexpected state: %s", ma_status$state[[1L]]))
  }
} else {
  cat("SKIP\tMicroarray driver optional check requires limma, ggplot2, yaml, and Biobase.\n")
}

cat("== Dangerous-pattern checks ==\n")
script_text <- unlist(lapply(list.files(script_dir, pattern = "\\.R$", full.names = TRUE, recursive = TRUE), readLines, warn = FALSE))
dangerous_patterns <- c(
  "max\\(vals.*<\\s*50",
  "as\\.matrix\\(counts_raw\\)"
)
for (pattern in dangerous_patterns) {
  if (any(grepl(pattern, script_text))) {
    fail(paste("Dangerous legacy pattern is still present:", pattern))
  }
}

cat("== Public tree path scan ==\n")
all_files <- list.files(repo_root, all.files = TRUE, recursive = TRUE, full.names = TRUE, no.. = TRUE)
text_files <- all_files[file.info(all_files)$isdir == FALSE]
bad_hits <- character()
private_user_pattern <- paste0("Users", "\\\\", "11495")
private_project_pattern <- paste0("SE", "_", "Program")
private_eval_pattern <- paste0("geo_biodata", "_", "eval")
windows_abs_path_pattern <- "[A-Za-z]:[\\\\/](Users|Program Files|ProgramData|Windows|Temp)"
for (path in text_files) {
  rel <- sub(paste0("^", gsub("\\\\", "/", normalizePath(repo_root, winslash = "/")), "/?"), "", normalizePath(path, winslash = "/"))
  if (grepl("^\\.git/", rel)) next
  content <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
  if (any(grepl(windows_abs_path_pattern, content) | grepl(private_user_pattern, content) | grepl(private_project_pattern, content) | grepl(private_eval_pattern, content))) {
    bad_hits <- c(bad_hits, rel)
  }
}
if (length(bad_hits) > 0L) {
  fail(paste("Public tree contains local/private strings:", paste(unique(bad_hits), collapse = ", ")))
}

cat("SMOKE_CHECKS_PASS\n")
