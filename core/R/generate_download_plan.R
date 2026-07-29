#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop(
    "Usage: generate_download_plan.R /path/to/run-directory 'reviewed-file-regex' [selection-reason]",
    call. = FALSE
  )
}

run_dir <- args[[1L]]
reviewed_regex <- args[[2L]]
selection_reason <- if (length(args) == 3L) args[[3L]] else "selected by reviewed filename regex"

is_broad_regex <- function(pattern) {
  normalized <- gsub("\\s+", "", tolower(pattern))
  normalized %in% c("", ".*", ".+", "*", "all", "^.*$", "^.+$") ||
    grepl("^\\.?\\*\\??$", normalized)
}

if (is_broad_regex(reviewed_regex)) {
  stop("Refusing to generate a broad download plan. Use a route-specific filename regex.", call. = FALSE)
}

supplement_path <- file.path(run_dir, "resources", "supplement_index.tsv")
if (!file.exists(supplement_path)) {
  stop("Missing supplement index: ", supplement_path, call. = FALSE)
}

supplements <- utils::read.delim(
  supplement_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
if (nrow(supplements) < 1L) {
  stop("Supplement index is empty.", call. = FALSE)
}
if ("status" %in% names(supplements)) {
  status_values <- supplements$status[!is.na(supplements$status) & nzchar(supplements$status)]
}
if (exists("status_values") && length(status_values) > 0L) {
  stop(
    "Supplement index is not downloadable; status column reports: ",
    paste(unique(status_values), collapse = ", "),
    call. = FALSE
  )
}
if (!"supplement_url" %in% names(supplements)) {
  if ("url" %in% names(supplements)) {
    supplements$supplement_url <- supplements$url
  } else {
    stop("Supplement index must include supplement_url. Re-run discover_geo.R with the updated workflow.", call. = FALSE)
  }
}
missing_url <- is.na(supplements$supplement_url) | !nzchar(supplements$supplement_url)
if (any(missing_url) && "url" %in% names(supplements)) {
  url_values <- as.character(supplements$url)
  url_values[is.na(url_values)] <- ""
  usable_url <- grepl("^(https?|ftp)://", url_values, ignore.case = TRUE)
  supplements$supplement_url[missing_url & usable_url] <- url_values[missing_url & usable_url]
}

file_name <- if ("file_name" %in% names(supplements)) {
  supplements$file_name
} else {
  basename(supplements$supplement_url)
}
search_text <- paste(file_name, supplements$supplement_url)
selected <- grepl(reviewed_regex, search_text, ignore.case = TRUE, perl = TRUE)
if (!any(selected)) {
  stop("The reviewed regex matched no supplementary files.", call. = FALSE)
}

accession <- basename(normalizePath(run_dir, winslash = "/", mustWork = FALSE))
series_meta <- file.path(run_dir, "resources", "series_metadata.tsv")
if (file.exists(series_meta)) {
  metadata <- utils::read.delim(series_meta, stringsAsFactors = FALSE)
  if (all(c("field", "value") %in% names(metadata))) {
    accession_row <- metadata$value[metadata$field %in% c("accession", "gse")]
    if (length(accession_row) > 0L && nzchar(accession_row[[1L]])) {
      accession <- accession_row[[1L]]
    }
  }
}

plans_dir <- file.path(run_dir, "plans")
if (!dir.exists(plans_dir)) dir.create(plans_dir, recursive = TRUE)

plan <- data.frame(
  accession = accession,
  source_accession = if ("source_accession" %in% names(supplements)) supplements$source_accession else accession,
  source_scope = if ("source_scope" %in% names(supplements)) supplements$source_scope else "series",
  selected = selected,
  reviewed = FALSE,
  file_name = file_name,
  supplement_url = supplements$supplement_url,
  size_bytes = if ("size" %in% names(supplements)) supplements$size else NA_real_,
  selection_reason = ifelse(selected, selection_reason, "not selected"),
  route_hint = "",
  review_note = ifelse(selected, "Set reviewed=TRUE only after checking route relevance.", ""),
  stringsAsFactors = FALSE
)

out_path <- file.path(plans_dir, "download_plan.tsv")
utils::write.table(
  plan,
  out_path,
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

cat("DOWNLOAD_PLAN_CREATED\n")
cat(normalizePath(out_path, winslash = "/", mustWork = TRUE), "\n")
