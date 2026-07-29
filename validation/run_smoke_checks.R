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
core_dir <- file.path(repo_root, "core", "R")
legacy_script_dir <- file.path(repo_root, "skills", "geo-biodata-workflow", "scripts")

fail <- function(message) stop(message, call. = FALSE)

expect_status <- function(label, command, args, expected_status) {
  status <- suppressWarnings(system2(command, shQuote(args), stdout = TRUE, stderr = TRUE))
  exit_status <- attr(status, "status") %||% 0L
  if (!identical(as.integer(exit_status), as.integer(expected_status))) {
    cat(paste(status, collapse = "\n"), "\n")
    fail(sprintf("%s: expected status %d, got %d", label, expected_status, exit_status))
  }
  invisible(status)
}

copy_fixture <- function(name, target) {
  source <- file.path(repo_root, "validation", "fixtures", name)
  files <- list.files(source, all.files = TRUE, recursive = TRUE, full.names = TRUE, no.. = TRUE)
  for (from in files) {
    rel <- substring(from, nchar(source) + 2L)
    to <- file.path(target, rel)
    dir.create(dirname(to), recursive = TRUE, showWarnings = FALSE)
    file.copy(from, to, overwrite = TRUE)
  }
  file.path(target, "run_manifest.yaml")
}

cat("== Parse R scripts ==\n")
r_files <- unlist(lapply(c(core_dir, legacy_script_dir), list.files,
  pattern = "\\.R$", full.names = TRUE, recursive = TRUE), use.names = FALSE)
for (path in sort(r_files)) {
  parse(file = path)
  rel <- sub(paste0("^", gsub("\\\\", "/", normalizePath(repo_root, winslash = "/")), "/?"),
    "", normalizePath(path, winslash = "/"))
  cat("PARSE_OK\t", rel, "\n", sep = "")
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
  if (length(closing) == 0L) fail(sprintf("Skill frontmatter is not closed: %s", path))
  header <- lines[2L:closing[[1L]]]
  parsed_header <- tryCatch(yaml::yaml.load(paste(header, collapse = "\n")), error = function(e) e)
  if (inherits(parsed_header, "error")) fail(sprintf("Skill frontmatter is invalid YAML: %s", path))
  if (!any(grepl("^name:\\s*[a-z0-9-]+\\s*$", header))) {
    fail(sprintf("Skill frontmatter missing valid name: %s", path))
  }
  if (is.null(parsed_header$description) || !nzchar(parsed_header$description)) {
    fail(sprintf("Skill frontmatter missing description: %s", path))
  }
  cat("SKILL_OK\t", basename(dirname(path)), "\n", sep = "")
}

cat("== Legacy shim size checks ==\n")
legacy_files <- list.files(legacy_script_dir, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
for (path in sort(legacy_files)) {
  size <- file.info(path)$size
  rel <- sub(paste0("^", gsub("\\\\", "/", normalizePath(repo_root, winslash = "/")), "/?"),
    "", normalizePath(path, winslash = "/"))
  if (size >= 600L) fail(sprintf("Legacy shim is >=600 bytes: %s (%d)", rel, size))
  cat("SHIM_OK\t", rel, "\t", size, "\n", sep = "")
}

scratch <- file.path(tempdir(), paste0("geo_biodata_smoke_", Sys.getpid()))
dir.create(scratch, recursive = TRUE, showWarnings = FALSE)

cat("== Manifest validation checks ==\n")
manifest_cases <- c(
  manifest_valid = 0L,
  manifest_metadata_only = 0L,
  manifest_scrna_author_object = 0L,
  manifest_bulk_normalized = 0L,
  manifest_microarray = 0L
)
for (name in names(manifest_cases)) {
  manifest <- copy_fixture(name, file.path(scratch, name))
  expect_status(paste("Manifest", name), "Rscript",
    c(file.path(core_dir, "validate_manifest.R"), manifest), manifest_cases[[name]])
}

required_field_refusals <- list(
  missing_accession = "^accession:",
  missing_input_species = "^  species:",
  missing_sample_mapping_file = "^  file: resources/sample_mapping_reviewed.tsv",
  missing_design_formula = "^  formula:",
  missing_contrast_numerator = "^    numerator:"
)
for (name in names(required_field_refusals)) {
  manifest <- copy_fixture("manifest_valid", file.path(scratch, name))
  lines <- readLines(manifest, warn = FALSE)
  writeLines(lines[!grepl(required_field_refusals[[name]], lines)], manifest, useBytes = TRUE)
  expect_status(paste("Required field refusal", name), "Rscript",
    c(file.path(core_dir, "validate_manifest.R"), manifest), 1L)
}

cat("== Enrichment quick checks ==\n")
if (!requireNamespace("fgsea", quietly = TRUE)) fail("fgsea is required for enrichment smoke checks.")
gmt_dir <- file.path(scratch, "gmt")
dir.create(gmt_dir, recursive = TRUE, showWarnings = FALSE)
gmt <- file.path(gmt_dir, "h.all.v1.Hs.symbols.gmt")
writeLines(c(
  paste(c("PATHWAY_A", "desc", sprintf("GENE%03d", 1:15)), collapse = "\t"),
  paste(c("PATHWAY_B", "desc", sprintf("GENE%03d", 10:30)), collapse = "\t")
), gmt, useBytes = TRUE)
query <- file.path(scratch, "ora_query.tsv")
utils::write.table(data.frame(gene_symbol = sprintf("GENE%03d", 1:6)),
  query, sep = "\t", quote = FALSE, row.names = FALSE)
universe <- file.path(scratch, "ora_universe.txt")
writeLines(sprintf("GENE%03d", 1:120), universe, useBytes = TRUE)
expect_status("ORA positive universe contract", "Rscript",
  c(file.path(core_dir, "enrichment", "run_ora_enrichment.R"),
    "--gene-list", query, "--gene-column", "gene_symbol",
    "--universe", universe, "--gmt-dir", gmt_dir,
    "--out-dir", file.path(scratch, "ora_positive"),
    "--collections", "h.all"), 0L)

rank <- file.path(scratch, "gsea_rank.tsv")
utils::write.table(data.frame(gene_symbol = sprintf("GENE%03d", 1:80), t = seq(4, -4, length.out = 80)),
  rank, sep = "\t", quote = FALSE, row.names = FALSE)
expect_status("GSEA ID-overlap gate positive", "Rscript",
  c(file.path(core_dir, "enrichment", "run_preranked_gsea.R"),
    "--rank-table", rank, "--gene-column", "gene_symbol", "--rank-column", "t",
    "--gmt", gmt, "--out-dir", file.path(scratch, "gsea_positive"),
    "--collection", "h.all"), 0L)

cat("== Negative fixture checks ==\n")
negative_base <- file.path(scratch, "negative")
missing_scale <- copy_fixture(file.path("negative", "missing_scale_contract"),
  file.path(negative_base, "missing_scale_contract"))
tpm_untransformed <- copy_fixture(file.path("negative", "tpm_untransformed_de"),
  file.path(negative_base, "tpm_untransformed_de"))
raw_mislabeled <- copy_fixture(file.path("negative", "raw_mislabeled_normalized"),
  file.path(negative_base, "raw_mislabeled_normalized"))
expect_status("Missing scale contract refusal", "Rscript",
  c(file.path(core_dir, "validate_manifest.R"), missing_scale), 1L)
expect_status("TPM untransformed DE refusal", "Rscript",
  c(file.path(core_dir, "validate_manifest.R"), tpm_untransformed), 1L)
expect_status("Raw-count mislabeled normalized refusal", "Rscript",
  c(file.path(core_dir, "drivers", "run_bulk_normalized.R"), raw_mislabeled), 1L)

if (quick_mode) cat("QUICK_MODE_COMPLETE\n")
cat("SMOKE_CHECKS_PASS\n")
