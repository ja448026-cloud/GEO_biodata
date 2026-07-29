#!/usr/bin/env Rscript
#
# generate_dataset_report.R
# Post-discovery dataset report: reads resources/ directory and produces a rich markdown
# report with dataset overview, publication details (DOI, abstract), and route summary.
#
# Usage:
#   Rscript generate_dataset_report.R runs/GSE000000
# Output:
#   runs/GSE000000/dataset_report.md

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: generate_dataset_report.R /path/to/run-directory", call. = FALSE)
}

run_dir <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
resources_dir <- file.path(run_dir, "resources")

meta_path <- file.path(resources_dir, "series_metadata.tsv")
pub_path  <- file.path(resources_dir, "publication_links.tsv")
route_path <- file.path(resources_dir, "routing_hint.tsv")
supp_path  <- file.path(resources_dir, "supplement_index.tsv")
samp_path  <- file.path(resources_dir, "sample_index.tsv")

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

# Required packages
for (pkg in c("httr2", "jsonlite")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Note:", pkg, "not available; DOI/abstract lookup will be skipped.\n")
  }
}

read_if <- function(path) {
  if (file.exists(path)) utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE) else NULL
}

has_quant_headers <- function(run_dir) {
  raw_dir <- file.path(run_dir, "raw")
  if (!dir.exists(raw_dir)) return(FALSE)
  paths <- list.files(raw_dir, pattern = "\\.(txt|tsv|csv)(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
  if (!length(paths)) return(FALSE)
  for (path in head(sort(paths), 5L)) {
    tab <- tryCatch(utils::read.delim(path, nrows = 5L, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL)
    if (is.null(tab)) next
    cols <- tolower(names(tab))
    if (any(cols %in% c("expected_count", "tpm", "fpkm"))) return(TRUE)
  }
  FALSE
}

meta   <- read_if(meta_path)
pub    <- read_if(pub_path)
route  <- read_if(route_path)
supp   <- read_if(supp_path)
samp   <- read_if(samp_path)

accession <- if (!is.null(route) && "accession" %in% names(route)) route$accession[[1L]] else basename(run_dir)
route_recommendation <- if (!is.null(route)) route$recommended_route[[1L]] %||% "unknown" else "unknown"
review_required <- if (!is.null(route)) isTRUE(route$review_required[[1L]]) else FALSE
if (has_quant_headers(run_dir) && identical(route_recommendation, "metadata_only")) {
  route_recommendation <- "bulk_quant_review_required"
  review_required <- TRUE
}

# Helper: extract a metadata field value
meta_val <- function(field_name) {
  if (is.null(meta)) return("unknown")
  row <- meta[meta$field == field_name, ]
  if (nrow(row) > 0L) row$value[[1L]] else "unknown"
}

query_europe_pmc <- function(pmid) {
  empty <- list(doi = NULL, pmcid = NULL, abstract = NULL, lookup_status = "not attempted")
  if (is.na(pmid) || !nzchar(pmid) || pmid == "not available") return(empty)
  if (!requireNamespace("httr2", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) {
    empty$lookup_status <- "skipped: httr2/jsonlite unavailable"
    return(empty)
  }
  tryCatch({
    resp <- httr2::request("https://www.ebi.ac.uk/europepmc/webservices/rest/search") |>
      httr2::req_url_query(
        query = paste0("EXT_ID:", pmid),
        resultType = "core",
        format = "json"
      ) |>
      httr2::req_user_agent("geo-biodata/0.1") |>
      httr2::req_retry(max_tries = 2L) |>
      httr2::req_perform()
    parsed <- jsonlite::fromJSON(httr2::resp_body_string(resp), simplifyVector = FALSE)
    results <- parsed$resultList$result
    if (length(results) < 1L) {
      return(list(doi = NULL, pmcid = NULL, abstract = NULL, lookup_status = "not found"))
    }
    r <- results[[1L]]
    list(
      doi = r$doi %||% NULL,
      pmcid = r$pmcid %||% NULL,
      abstract = r$abstractText %||% NULL,
      lookup_status = "ok"
    )
  }, error = function(e) {
    list(doi = NULL, pmcid = NULL, abstract = NULL,
      lookup_status = paste("failed:", conditionMessage(e)))
  })
}

# ── Fetch DOI + abstract from Europe PMC if not already in pub table ──────────

doi_str     <- "not available"
abstract_str <- NULL
pmid_str    <- "not available"
pmc_str     <- "not available"

if (!is.null(pub) && nrow(pub) > 0L) {
  pmid_str <- as.character(pub$identifier[[1L]])
  if (is.na(pmid_str) || !nzchar(pmid_str)) pmid_str <- "not available"

  # Check if pub already has doi/abstract columns from enhanced discover_geo
  if ("doi" %in% names(pub) && !is.na(pub$doi[[1L]]) && nzchar(pub$doi[[1L]])) {
    doi_str <- pub$doi[[1L]]
  }
  if ("abstract" %in% names(pub) && !is.na(pub$abstract[[1L]]) && nzchar(pub$abstract[[1L]])) {
    abstract_str <- pub$abstract[[1L]]
  }
  if ("pmcid" %in% names(pub) && !is.na(pub$pmcid[[1L]]) && nzchar(pub$pmcid[[1L]])) {
    pmc_str <- pub$pmcid[[1L]]
  }

  # Fetch from Europe PMC if missing
  if (doi_str == "not available" || is.null(abstract_str)) {
    epm <- query_europe_pmc(pmid_str)
    if (doi_str == "not available") {
      doi_str <- epm$doi %||% "not available"
    }
    if (is.null(abstract_str)) {
      abstract_str <- epm$abstract %||% NULL
    }
    if (pmc_str == "not available") {
      pmc_str <- epm$pmcid %||% "not available"
    }
    epm_lookup_status <- epm$lookup_status
  }
}
if (!exists("epm_lookup_status")) epm_lookup_status <- "not needed"

# ── Assemble report ──────────────────────────────────────────────────────────

now_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)

lines <- c(
  paste0("# ", accession, " — GEO Dataset Report"),
  "",
  paste("**Generated**:", now_utc),
  "",
  "## 1. Dataset Overview",
  "",
  paste("| Field | Value |"),
  paste("|-------|-------|"),
  paste("| **Title** |", meta_val("title"), "|"),
  paste("| **Organism** |", meta_val("taxon"), "|"),
  paste("| **GEO Status** |", meta_val("status"), "|"),
  paste("| **Experiment Type** |", meta_val("gdstype"), "|"),
  paste("| **Submission Date** |", meta_val("submission_date"), "|"),
  paste("| **Last Update** |", meta_val("last_update_date"), "|")
)

if (!is.null(samp)) {
  lines <- c(lines, paste("| **Samples** |", nrow(samp), "|"))
}
if (!is.null(route)) {
  lines <- c(lines,
    paste("| **Assay Type** |", route$assay_type[[1L]] %||% "unknown", "|"),
    paste("| **Recommended Route** |", route_recommendation, "|")
  )
}
if (!is.null(supp) && nrow(supp) > 0L && "file_name" %in% names(supp)) {
  n_supp <- sum(nzchar(supp$file_name) & !is.na(supp$file_name))
  lines <- c(lines, paste("| **Supplementary Files** |", n_supp, "|"))
}

# ── Publication section ──────────────────────────────────────────────────────

lines <- c(lines,
  "",
  "## 2. Publication",
  "",
  paste("| Field | Value |"),
  paste("|-------|-------|"),
  paste("| **PMID** |", pmid_str, "|"),
  paste("| **DOI** |", doi_str, "|"),
  paste("| **DOI/Abstract Lookup** |", epm_lookup_status, "|")
)

if (pmid_str != "not available") {
  lines <- c(lines, paste("| **PubMed** | https://pubmed.ncbi.nlm.nih.gov/", pmid_str, "/ |"))
}
if (pmc_str != "not available") {
  lines <- c(lines, paste("| **Full Text (PMC)** | https://www.ncbi.nlm.nih.gov/pmc/articles/", pmc_str, "/ |"))
}

# ── Abstract ─────────────────────────────────────────────────────────────────

lines <- c(lines, "", "## 3. GEO Summary", "")
geo_summary <- meta_val("summary")
if (nzchar(geo_summary) && geo_summary != "unknown") {
  lines <- c(lines, geo_summary)
} else {
  lines <- c(lines, "_No GEO summary available._")
}

if (!is.null(abstract_str) && nzchar(abstract_str)) {
  lines <- c(lines, "", "## 4. Publication Abstract", "", abstract_str)
}

# ── Data & Routes ────────────────────────────────────────────────────────────

lines <- c(lines, "", "## 5. Data Files & Routes")

if (!is.null(supp) && nrow(supp) > 0L && "file_name" %in% names(supp)) {
  lines <- c(lines, "", "### Supplementary Files", "")
  for (i in seq_len(min(nrow(supp), 20L))) {
    fname <- supp$file_name[[i]]
    fsize <- if ("size_mb" %in% names(supp) && !is.na(supp$size_mb[[i]])) {
      paste0(" (", supp$size_mb[[i]], " Mb)")
    } else if ("size" %in% names(supp) && !is.na(supp$size[[i]])) {
      paste0(" (", round(supp$size[[i]] / 1e6, 1), " Mb)")
    } else { "" }
    lines <- c(lines, paste("- ", fname, fsize))
  }
} else {
  lines <- c(lines, "", "_No supplementary files listed._")
}

if (!is.null(route)) {
  lines <- c(lines,
    "",
    paste("**Route recommendation**: `", route_recommendation, "`"),
    paste("  - Assay type:", route$assay_type[[1L]] %||% "unknown"),
    paste("  - Review required:", if (isTRUE(review_required)) "yes" else "no")
  )
}

# ── Footer ───────────────────────────────────────────────────────────────────

lines <- c(lines,
  "",
  "---",
  "",
  paste("Report generated by `generate_dataset_report.R` at", now_utc),
  paste("Run directory:", run_dir),
  ""
)

writeLines(lines, file.path(run_dir, "dataset_report.md"), useBytes = TRUE)

cat("DATASET_REPORT_GENERATED\n")
cat(normalizePath(file.path(run_dir, "dataset_report.md"), winslash = "/", mustWork = TRUE), "\n")
