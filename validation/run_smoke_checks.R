#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

args <- commandArgs(trailingOnly = TRUE)
quick_mode <- "--quick" %in% args

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
this_file <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(file.path("validation", "run_smoke_checks.R"), mustWork = TRUE)
}
repo_root <- normalizePath(file.path(dirname(this_file), ".."), mustWork = TRUE)
script_dir <- file.path(repo_root, "skills", "geo-biodata-workflow", "scripts")
core_dir <- file.path(repo_root, "core", "R")

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
script_roots <- c(script_dir, core_dir)
r_files <- unlist(lapply(script_roots, list.files, pattern = "\\.R$", full.names = TRUE, recursive = TRUE),
  use.names = FALSE)
for (path in sort(r_files)) {
  parse(file = path)
  rel_path <- sub(paste0("^", gsub("\\\\", "/", normalizePath(repo_root, winslash = "/")), "/?"),
    "", normalizePath(path, winslash = "/"))
  cat("PARSE_OK\t", rel_path, "\n", sep = "")
}

cat("== Skill frontmatter checks ==\n")
skill_files <- list.files(file.path(repo_root, "skills"), pattern = "^SKILL\\.md$",
  full.names = TRUE, recursive = TRUE)
for (path in sort(skill_files)) {
  lines <- readLines(path, warn = FALSE)
  if (length(lines) < 4L || !identical(lines[[1L]], "---")) {
    fail(sprintf("Skill missing YAML frontmatter: %s", path))
  }
  closing <- which(lines[-1L] == "---")
  if (length(closing) == 0L) {
    fail(sprintf("Skill frontmatter is not closed: %s", path))
  }
  header <- lines[2L:closing[[1L]]]
  parsed_header <- tryCatch(yaml::yaml.load(paste(header, collapse = "\n")),
    error = function(e) e)
  if (inherits(parsed_header, "error")) {
    fail(sprintf("Skill frontmatter is invalid YAML: %s", path))
  }
  if (!any(grepl("^name:\\s*[a-z0-9-]+\\s*$", header))) {
    fail(sprintf("Skill frontmatter missing valid name: %s", path))
  }
  if (is.null(parsed_header$description) || !nzchar(parsed_header$description)) {
    fail(sprintf("Skill frontmatter missing description: %s", path))
  }
  cat("SKILL_OK\t", basename(dirname(path)), "\n", sep = "")
}

cat("== Figure playbook checks ==\n")
figure_skill_dir <- file.path(repo_root, "skills", "geo-biodata-figure")
required_figure_files <- file.path(figure_skill_dir, c(
  "SKILL.md",
  "references/principles-and-qa.md",
  "references/figure-decision-matrix.md",
  "references/omics-plot-recipes.md",
  "references/manuscript-figure-planning.md",
  "templates/figure_plan.yaml",
  "templates/figure_spec.yaml",
  "templates/figure_qa.tsv",
  "templates/figure_caption.md"
))
required_figure_files <- c(required_figure_files, file.path(repo_root, "knowledge", "visualization_source_registry.yaml"))
missing_figure_files <- required_figure_files[!file.exists(required_figure_files)]
if (length(missing_figure_files) > 0L) {
  fail(paste("Figure playbook files missing:", paste(missing_figure_files, collapse = ", ")))
}
figure_skill_text <- paste(readLines(file.path(figure_skill_dir, "SKILL.md"), warn = FALSE), collapse = "\n")
if (!grepl("skill_mode`: playbook", figure_skill_text, fixed = TRUE) ||
    !grepl("external_skill_dependency`: none", figure_skill_text, fixed = TRUE)) {
  fail("Figure skill must declare playbook mode and no external skill dependency.")
}
invisible(yaml::yaml.load_file(file.path(repo_root, "knowledge", "visualization_source_registry.yaml")))
blocked_skill_names <- c(
  paste0("figure", "-", "planner"),
  paste0("nature", "-", "figure"),
  paste0("omics", "-", "figure", "-", "qa")
)
repo_text_files <- list.files(repo_root, all.files = TRUE, recursive = TRUE, full.names = TRUE, no.. = TRUE)
repo_text_files <- repo_text_files[file.info(repo_text_files)$isdir == FALSE]
blocked_hits <- character()
for (path in repo_text_files) {
  rel <- sub(paste0("^", gsub("\\\\", "/", normalizePath(repo_root, winslash = "/")), "/?"), "", normalizePath(path, winslash = "/"))
  if (grepl("^\\.git/", rel)) next
  content <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
  if (any(vapply(blocked_skill_names, function(pattern) any(grepl(pattern, content, fixed = TRUE)), logical(1)))) {
    blocked_hits <- c(blocked_hits, rel)
  }
}
if (length(blocked_hits) > 0L) {
  fail(paste("Figure playbook still names removed external skill dependencies:", paste(unique(blocked_hits), collapse = ", ")))
}
cat("FIGURE_PLAYBOOK_OK\n")

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

legacy_url_dir <- file.path(scratch, "legacy_url_index")
dir.create(file.path(legacy_url_dir, "resources"), recursive = TRUE, showWarnings = FALSE)
legacy_url <- "https://example.org/GSE000000_counts.tsv.gz"
utils::write.table(
  data.frame(
    fname = "GSE000000_counts.tsv.gz",
    url = legacy_url,
    stringsAsFactors = FALSE
  ),
  file.path(legacy_url_dir, "resources", "supplement_index.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
expect_status(
  "Download-plan URL-column fallback",
  "Rscript",
  c(file.path(script_dir, "generate_download_plan.R"), legacy_url_dir, "counts", "URL fallback candidate"),
  0L
)
legacy_plan <- utils::read.delim(file.path(legacy_url_dir, "plans", "download_plan.tsv"),
  stringsAsFactors = FALSE, check.names = FALSE)
if (!identical(legacy_plan$supplement_url[[1L]], legacy_url)) {
  fail("Download-plan URL-column fallback did not populate supplement_url.")
}

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
  "Manifest dependency plan",
  "Rscript",
  c(file.path(core_dir, "bootstrap_environment.R"), "--profile", "manifest", "--plan"),
  0L
)

expect_status(
  "Legacy core dependency plan",
  "Rscript",
  c(file.path(script_dir, "bootstrap_environment.R"), "--profile", "core", "--plan"),
  0L
)

cat("== Core wrapper checks ==\n")
expect_status(
  "Core manifest validation wrapper",
  "Rscript",
  c(file.path(core_dir, "validate_manifest.R"), manifest_fixture),
  0L
)

cat("== Bulk driver optional check ==\n")
if (quick_mode) {
  cat("SKIP\tBulk driver optional check skipped in quick mode.\n")
} else if (all(vapply(c("DESeq2", "ggplot2", "SummarizedExperiment"), requireNamespace, logical(1), quietly = TRUE))) {
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
  valid_execution <- c("EXECUTION_COMPLETE", "BASIC_ANALYSIS_COMPLETE", "QC_REVIEW_REQUIRED")
  exec_state <- if ("execution_state" %in% names(status)) status$execution_state[[1L]] else status$state[[1L]]
  if (!exec_state %in% valid_execution) {
    fail(sprintf("Bulk driver produced unexpected execution_state: %s", exec_state))
  }
  if ("technical_qc" %in% names(status) && status$technical_qc[[1L]] == "REVIEW_REQUIRED") {
    cat("NOTE: Bulk driver technical_qc=REVIEW_REQUIRED (expected with small fixture).\n")
  }
  qc_table <- utils::read.delim(file.path(bulk_dir, "tables", "qc_checks.tsv"), stringsAsFactors = FALSE)
  if (nrow(qc_table) == 0L) fail("Bulk driver did not produce QC checks table.")
} else {
  cat("SKIP\tBulk driver optional check requires DESeq2, ggplot2, and SummarizedExperiment.\n")
}

cat("== Normalized bulk driver optional check ==\n")
if (quick_mode) {
  cat("SKIP\tBulk normalized driver optional check skipped in quick mode.\n")
} else if (all(vapply(c("limma", "ggplot2", "yaml"), requireNamespace, logical(1), quietly = TRUE))) {
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
  fallback_path <- file.path(norm_dir, "tables", "fallback_events.tsv")
  if (!file.exists(fallback_path) || file.info(fallback_path)$size <= 0L) {
    fail("Bulk normalized driver did not write fallback_events.tsv contract file.")
  }
  cor_path <- file.path(norm_dir, "tables", "sample_correlation.tsv")
  cor_fig <- file.path(norm_dir, "figures", "sample_correlation_heatmap.pdf")
  if (!file.exists(cor_path) || file.info(cor_path)$size <= 0L ||
      !file.exists(cor_fig) || file.info(cor_fig)$size <= 0L) {
    fail("Bulk normalized driver did not write sample-correlation table and heatmap.")
  }
  norm_status <- utils::read.delim(file.path(norm_dir, "workflow_status.tsv"), stringsAsFactors = FALSE)
  valid_exec <- c("EXECUTION_COMPLETE", "BASIC_ANALYSIS_COMPLETE", "QC_REVIEW_REQUIRED")
  exec_s <- if ("execution_state" %in% names(norm_status)) norm_status$execution_state[[1L]] else norm_status$state[[1L]]
  if (!exec_s %in% valid_exec) {
    fail(sprintf("Bulk normalized driver produced unexpected state: %s", exec_s))
  }
} else {
  cat("SKIP\tBulk normalized driver optional check requires limma, ggplot2, and yaml.\n")
}

cat("== Microarray driver optional check ==\n")
if (quick_mode) {
  cat("SKIP\tMicroarray driver optional check skipped in quick mode.\n")
} else if (all(vapply(c("limma", "ggplot2", "yaml", "Biobase"), requireNamespace, logical(1), quietly = TRUE))) {
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
  valid_ma <- c("EXECUTION_COMPLETE", "BASIC_ANALYSIS_COMPLETE", "QC_REVIEW_REQUIRED", "EXECUTION_FAILED")
  ma_s <- if ("execution_state" %in% names(ma_status)) ma_status$execution_state[[1L]] else ma_status$state[[1L]]
  if (!ma_s %in% valid_ma) {
    fail(sprintf("Microarray driver produced unexpected state: %s", ma_s))
  }
} else {
  cat("SKIP\tMicroarray driver optional check requires limma, ggplot2, yaml, and Biobase.\n")
}

cat("== Negative fixture checks ==\n")

# raw-like mislabeled: expected fail (INPUT_SCALE_CONFLICT)
raw_mislabel_dir <- file.path(scratch, "neg_raw_mislabel")
dir.create(file.path(raw_mislabel_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(raw_mislabel_dir, "resources"), recursive = TRUE, showWarnings = FALSE)
neg_fixture_base <- file.path(repo_root, "validation", "fixtures", "negative")
invisible(file.copy(file.path(neg_fixture_base, "raw_mislabeled_normalized", "run_manifest.yaml"),
  file.path(raw_mislabel_dir, "run_manifest.yaml"), overwrite = TRUE))
invisible(file.copy(file.path(neg_fixture_base, "raw_mislabeled_normalized", "raw", "counts.tsv"),
  file.path(raw_mislabel_dir, "raw", "counts.tsv"), overwrite = TRUE))
invisible(file.copy(file.path(neg_fixture_base, "raw_mislabeled_normalized", "resources", "sample_mapping_reviewed.tsv"),
  file.path(raw_mislabel_dir, "resources", "sample_mapping_reviewed.tsv"), overwrite = TRUE))
if (quick_mode) {
  cat("SKIP\tRaw-count mislabeled driver check skipped in quick mode.\n")
} else {
  expect_status(
    "Raw-count mislabeled as log_normalized failure",
    "Rscript",
    c(file.path(script_dir, "drivers", "run_bulk_normalized.R"), file.path(raw_mislabel_dir, "run_manifest.yaml")),
    1L
  )
}

# missing scale + normalized_expression + DE: expected MANIFEST_INVALID
missing_scale_dir <- file.path(scratch, "neg_missing_scale")
dir.create(file.path(missing_scale_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(missing_scale_dir, "resources"), recursive = TRUE, showWarnings = FALSE)
invisible(file.copy(file.path(neg_fixture_base, "missing_scale_contract", "run_manifest.yaml"),
  file.path(missing_scale_dir, "run_manifest.yaml"), overwrite = TRUE))
invisible(file.copy(file.path(neg_fixture_base, "missing_scale_contract", "raw", "counts.tsv"),
  file.path(missing_scale_dir, "raw", "counts.tsv"), overwrite = TRUE))
invisible(file.copy(file.path(neg_fixture_base, "missing_scale_contract", "resources", "sample_mapping_reviewed.tsv"),
  file.path(missing_scale_dir, "resources", "sample_mapping_reviewed.tsv"), overwrite = TRUE))
expect_status(
  "Missing scale contract refusal",
  "Rscript",
  c(file.path(script_dir, "validate_manifest.R"), file.path(missing_scale_dir, "run_manifest.yaml")),
  1L
)

# TPM untransformed DE
tpm_dir <- file.path(scratch, "neg_tpm")
dir.create(file.path(tpm_dir, "raw"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(tpm_dir, "resources"), recursive = TRUE, showWarnings = FALSE)
invisible(file.copy(file.path(neg_fixture_base, "tpm_untransformed_de", "run_manifest.yaml"),
  file.path(tpm_dir, "run_manifest.yaml"), overwrite = TRUE))
invisible(file.copy(file.path(neg_fixture_base, "tpm_untransformed_de", "raw", "counts.tsv"),
  file.path(tpm_dir, "raw", "counts.tsv"), overwrite = TRUE))
invisible(file.copy(file.path(neg_fixture_base, "tpm_untransformed_de", "resources", "sample_mapping_reviewed.tsv"),
  file.path(tpm_dir, "resources", "sample_mapping_reviewed.tsv"), overwrite = TRUE))
expect_status(
  "TPM untransformed DE refusal",
  "Rscript",
  c(file.path(script_dir, "validate_manifest.R"), file.path(tpm_dir, "run_manifest.yaml")),
  1L
)

cat("== Dangerous-pattern checks ==\n")
script_text <- unlist(lapply(list.files(script_dir, pattern = "\\.R$", full.names = TRUE, recursive = TRUE), readLines, warn = FALSE))
script_text <- c(
  script_text,
  unlist(lapply(list.files(core_dir, pattern = "\\.R$", full.names = TRUE, recursive = TRUE), readLines, warn = FALSE))
)
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
