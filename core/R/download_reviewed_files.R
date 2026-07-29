#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(2L, 3L)) {
  stop(
    paste(
      "Usage:",
      "download_reviewed_files.R /path/to/download_plan.tsv /path/to/raw-directory",
      "or legacy: download_reviewed_files.R GSE000000 /path/to/raw-directory 'reviewed-regex'",
      sep = "\n  "
    ),
    call. = FALSE
  )
}

truthy <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "yes", "y", "1")
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

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
download_methods <- tolower(Sys.getenv("GEO_BIODATA_DOWNLOAD_METHODS", unset = "geoquery,aria2c,curl,download.file"))
backup_timeout_sec <- as.numeric(Sys.getenv("GEO_BIODATA_BACKUP_TIMEOUT_SEC", unset = "1200"))
if (!is.finite(download_timeout_sec) || download_timeout_sec < 1) download_timeout_sec <- 300
if (is.na(download_retries) || download_retries < 1L) download_retries <- 3L
if (!is.finite(backup_timeout_sec) || backup_timeout_sec < 1) backup_timeout_sec <- 1200
download_methods <- trimws(strsplit(download_methods, ",", fixed = TRUE)[[1L]])
download_methods <- download_methods[nzchar(download_methods) & download_methods != "none"]
download_methods <- unique(download_methods)
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

download_ok <- function(dest) {
  file.exists(dest) && file.info(dest)$size > 0L
}

escape_regex <- function(x) {
  gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)
}

geo_accession_from_url <- function(url, fallback_accession = NA_character_) {
  fallback_accession <- toupper(fallback_accession %||% "")
  if (grepl("^(GSE|GSM)[0-9]+$", fallback_accession)) return(fallback_accession)
  gsm <- regmatches(url, regexpr("GSM[0-9]+", url, ignore.case = TRUE))
  if (length(gsm) > 0L && nzchar(gsm)) return(toupper(gsm))
  gse <- regmatches(url, regexpr("GSE[0-9]+", url, ignore.case = TRUE))
  if (length(gse) > 0L && nzchar(gse)) return(toupper(gse))
  NA_character_
}

download_with_geoquery <- function(url, dest, accession = NA_character_) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) return(FALSE)
  geo_acc <- geo_accession_from_url(url, accession)
  if (is.na(geo_acc) || !grepl("^(GSE|GSM)[0-9]+$", geo_acc)) return(FALSE)

  if (download_ok(dest)) unlink(dest)
  target_name <- basename(dest)
  target_pattern <- paste0("^", escape_regex(target_name), "$")
  raw_dir <- dirname(dest)
  downloaded <- tryCatch(
    GEOquery::getGEOSuppFiles(
      geo_acc,
      makeDirectory = FALSE,
      baseDir = raw_dir,
      fetch_files = TRUE,
      filter_regex = target_pattern
    ),
    error = function(e) NULL
  )
  if (is.null(downloaded)) return(FALSE)

  candidates <- unique(c(dest, rownames(downloaded), file.path(raw_dir, target_name)))
  candidates <- candidates[file.exists(candidates) & basename(candidates) == target_name]
  candidates <- candidates[file.info(candidates)$size > 0L]
  if (length(candidates) < 1L) return(FALSE)

  if (!identical(normalizePath(candidates[[1L]], winslash = "/", mustWork = TRUE),
                 normalizePath(dest, winslash = "/", mustWork = FALSE))) {
    file.copy(candidates[[1L]], dest, overwrite = TRUE)
  }
  download_ok(dest)
}

download_with_aria2c <- function(url, dest) {
  aria2c <- Sys.which("aria2c")
  if (!nzchar(aria2c)) return(FALSE)
  if (file.exists(dest)) unlink(dest)
  status <- tryCatch(
    system2(aria2c, c(
      "--allow-overwrite=true",
      "--auto-file-renaming=false",
      "--continue=true",
      "--max-connection-per-server=4",
      "--split=4",
      "--max-tries", as.character(download_retries),
      "--timeout", as.character(min(backup_timeout_sec, 1200)),
      "--dir", dirname(dest),
      "--out", basename(dest),
      url
    )),
    error = function(e) 1L
  )
  identical(as.integer(status), 0L) && download_ok(dest)
}

download_with_curl <- function(url, dest) {
  curl <- Sys.which("curl")
  if (!nzchar(curl)) return(FALSE)
  if (file.exists(dest)) unlink(dest)
  status <- tryCatch(
    system2(curl, c(
      "-L",
      "--retry", as.character(download_retries),
      "--max-time", as.character(backup_timeout_sec),
      "-o", dest,
      url
    )),
    error = function(e) 1L
  )
  identical(as.integer(status), 0L) && download_ok(dest)
}

download_with_download_file <- function(url, dest) {
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
      return(TRUE)
    }
    if (attempt < download_retries) {
      wait_sec <- min(30, 2^attempt)
      message("Download attempt ", attempt, " failed; retrying in ", wait_sec, " seconds.")
      Sys.sleep(wait_sec)
    }
  }
  if (file.exists(dest)) unlink(dest)
  attr(FALSE, "last_error") <- last_error
  FALSE
}

download_one <- function(url, dest, accession = NA_character_) {
  errors <- character()
  for (method in download_methods) {
    ok <- switch(method,
      "geoquery" = download_with_geoquery(url, dest, accession),
      "aria2" = download_with_aria2c(url, dest),
      "aria2c" = download_with_aria2c(url, dest),
      "curl" = download_with_curl(url, dest),
      "download.file" = download_with_download_file(url, dest),
      "download_file" = download_with_download_file(url, dest),
      FALSE
    )
    if (isTRUE(ok) && download_ok(dest)) {
      message("Download OK via ", method, ": ", basename(dest))
      return(list(path = dest, method = method))
    }
    errors <- c(errors, method)
  }
  if (file.exists(dest)) unlink(dest)
  stop("All download methods failed for ", url, ". Tried: ", paste(errors, collapse = ", "), call. = FALSE)
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
  planned_names <- basename(as.character(selected$file_name))
  if (any(!nzchar(planned_names)) || any(is.na(planned_names))) {
    stop("Selected rows must include non-empty local file_name values.", call. = FALSE)
  }
  duplicate_names <- unique(planned_names[duplicated(tolower(planned_names))])
  if (length(duplicate_names) > 0L) {
    stop("Selected download plan would overwrite local files. Regenerate the plan or assign unique file_name values: ",
      paste(duplicate_names, collapse = ", "), call. = FALSE)
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
  methods_used <- character(nrow(selected))
  for (i in seq_len(nrow(selected))) {
    dest <- file.path(raw_dir, basename(selected$file_name[[i]]))
    source_i <- if ("source_accession" %in% names(selected)) as.character(selected$source_accession[[i]]) else ""
    if (is.na(source_i)) source_i <- ""
    accession_i <- if (nzchar(source_i)) {
      selected$source_accession[[i]]
    } else {
      selected$accession[[i]]
    }
    downloaded <- download_one(
      selected$supplement_url[[i]],
      dest,
      accession = accession_i
    )
    paths[[i]] <- downloaded$path
    methods_used[[i]] <- downloaded$method
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
  methods_used <- rep("geoquery", length(paths))
}

info <- file.info(paths)
manifest <- data.frame(
  accession = accession,
  file = basename(paths),
  bytes = info$size,
  download_method = methods_used,
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
