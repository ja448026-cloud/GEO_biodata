#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: download_geo_supp.R GSE000000 /path/to/raw-directory 'reviewed-regex'",
    call. = FALSE
  )
}

accession <- toupper(args[[1L]])
raw_dir <- args[[2L]]
filter_regex <- args[[3L]]
if (!grepl("^GSE[0-9]+$", accession)) {
  stop("The accession must match GSE followed by digits.", call. = FALSE)
}
if (!nzchar(filter_regex) || identical(filter_regex, ".*")) {
  stop("Provide a reviewed filename regex; refusing an unfiltered bulk download.", call. = FALSE)
}

required <- c("GEOquery", "digest")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}
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
