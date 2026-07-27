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
r_files <- list.files(script_dir, pattern = "\\.R$", full.names = TRUE)
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
