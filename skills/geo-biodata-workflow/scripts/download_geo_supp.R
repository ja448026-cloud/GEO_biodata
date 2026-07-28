#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L, 3L)) {
  stop(
    paste(
      "Usage:",
      "download_geo_supp.R /path/to/download_plan.tsv /path/to/raw-directory",
      "or legacy: download_geo_supp.R GSE000000 /path/to/raw-directory 'reviewed-regex'",
      sep = "\n  "
    ),
    call. = FALSE
  )
}

truthy <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "yes", "y", "1")
}

is_broad_regex <- function(pattern) {
  normalized <- gsub("\\s+", "", tolower(pattern))
  normalized %in% c("", ".*", ".+", "*", "all", "^.*$", "^.+$") ||
    grepl("^\\.?\\*\\??$", normalized)
}

max_total_gb <- as.numeric(Sys.getenv("GEO_BIODATA_MAX_TOTAL_GB", unset = "20"))
max_single_gb <- as.numeric(Sys.getenv("GEO_BIODATA_MAX_SINGLE_GB", unset = "10"))
max_total_bytes <- max_total_gb * 1024^3
max_single_bytes <- max_single_gb * 1024^3
download_timeout_sec <- as.numeric(Sys.getenv("GEO_BIODATA_DOWNLOAD_TIMEOUT_SEC", unset = "300"))
download_retries <- as.integer(Sys.getenv("GEO_BIODATA_DOWNLOAD_RETRIES", unset = "3"))
if (!is.finite(download_timeout_sec) || download_timeout_sec < 1) download_timeout_sec <- 300
if (is.na(download_retries) || download_retries < 1L) download_retries <- 3L
old_timeout <- getOption("timeout")
options(timeout = max(old_timeout, download_timeout_sec))
on.exit(options(timeout = old_timeout), add = TRUE)

required <- c("digest")
if (length(args) == 3L) {
  required <- c(required, "GEOquery")
}
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

download_with_retries <- function(url, dest) {
  last_error <- ""
  for (attempt in seq_len(download_retries)) {
    if (file.exists(dest)) unlink(dest)
    status <- tryCatch(
      utils::download.file(url, destfile = dest, mode = "wb", quiet = FALSE),
      error = function(e) {
        last_error <<- conditionMessage(e)
        1L
      }
    )
    if (identical(status, 0L) && file.exists(dest) && file.info(dest)$size > 0L) {
      return(dest)
    }
    if (attempt < download_retries) {
      wait_sec <- min(30, 2^attempt)
      message("Download attempt ", attempt, " failed; retrying in ", wait_sec, " seconds.")
      Sys.sleep(wait_sec)
    }
  }
  if (file.exists(dest)) unlink(dest)
  stop("Download failed after ", download_retries, " attempt(s) for ", url,
    if (nzchar(last_error)) paste0(": ", last_error) else "", call. = FALSE)
}

if (length(args) == 2L) {
  plan_path <- args[[1L]]
  raw_dir <- args[[2L]]
  if (!file.exists(plan_path)) stop("Missing download plan: ", plan_path, call. = FALSE)

  plan <- utils::read.delim(plan_path, stringsAsFactors = FALSE, check.names = FALSE)
  required_cols <- c("accession", "selected", "reviewed", "file_name", "supplement_url")
  missing_cols <- setdiff(required_cols, names(plan))
  if (length(missing_cols) > 0L) {
    stop("Download plan missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  selected <- plan[truthy(plan$selected), , drop = FALSE]
  if (nrow(selected) < 1L) {
    stop("Download plan has no selected files.", call. = FALSE)
  }
  if (!all(truthy(selected$reviewed))) {
    stop("All selected files must have reviewed=TRUE before download.", call. = FALSE)
  }
  if (any(!nzchar(selected$supplement_url))) {
    stop("Selected rows must include supplement_url.", call. = FALSE)
  }
  if (any(!grepl("^(https?|ftp)://", selected$supplement_url, ignore.case = TRUE))) {
    stop("Selected supplement URLs must start with http://, https://, or ftp://.", call. = FALSE)
  }
  if ("size_bytes" %in% names(selected)) {
    sizes <- suppressWarnings(as.numeric(selected$size_bytes))
    known_sizes <- sizes[is.finite(sizes)]
    if (length(known_sizes) > 0L) {
      if (any(known_sizes > max_single_bytes)) {
        stop("At least one selected file exceeds GEO_BIODATA_MAX_SINGLE_GB.", call. = FALSE)
      }
      if (sum(known_sizes) > max_total_bytes) {
        stop("Selected files exceed GEO_BIODATA_MAX_TOTAL_GB.", call. = FALSE)
      }
    }
  }

  if (!dir.exists(raw_dir)) dir.create(raw_dir, recursive = TRUE)
  paths <- character(nrow(selected))
  for (i in seq_len(nrow(selected))) {
    dest <- file.path(raw_dir, basename(selected$file_name[[i]]))
    download_with_retries(selected$supplement_url[[i]], dest)
    paths[[i]] <- dest
  }
  accession <- selected$accession
} else {
  accession <- toupper(args[[1L]])
  raw_dir <- args[[2L]]
  filter_regex <- args[[3L]]
  if (!grepl("^GSE[0-9]+$", accession)) {
    stop("The accession must match GSE followed by digits.", call. = FALSE)
  }
  if (is_broad_regex(filter_regex)) {
    stop("Provide a route-specific filename regex; refusing an unfiltered bulk download.", call. = FALSE)
  }
  cat("Warning: legacy regex download mode is retained for compatibility. Prefer generate_download_plan.R followed by plan review.\n")
  if (!dir.exists(raw_dir)) dir.create(raw_dir, recursive = TRUE)

  listed <- GEOquery::getGEOSuppFiles(
    accession,
    makeDirectory = FALSE,
    baseDir = raw_dir,
    fetch_files = FALSE,
    filter_regex = filter_regex
  )
  if (nrow(listed) < 1L) {
    stop("The reviewed regex matched no GEO supplementary files.", call. = FALSE)
  }

  downloaded <- GEOquery::getGEOSuppFiles(
    accession,
    makeDirectory = FALSE,
    baseDir = raw_dir,
    fetch_files = TRUE,
    filter_regex = filter_regex
  )
  paths <- rownames(downloaded)
  paths <- paths[file.exists(paths)]
  if (length(paths) < 1L) stop("No downloaded file exists after transfer.", call. = FALSE)
}

info <- file.info(paths)
manifest <- data.frame(
  accession = accession,
  file = basename(paths),
  bytes = info$size,
  sha256 = vapply(
    paths,
    digest::digest,
    character(1),
    algo = "sha256",
    file = TRUE,
    serialize = FALSE
  ),
  downloaded_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  stringsAsFactors = FALSE
)
utils::write.table(
  manifest,
  file.path(raw_dir, "download_manifest.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

if (any(manifest$bytes <= 0L)) stop("At least one downloaded file is empty.", call. = FALSE)
cat("DOWNLOAD_COMPLETE\n")
cat(normalizePath(file.path(raw_dir, "download_manifest.tsv"), winslash = "/", mustWork = TRUE), "\n")
