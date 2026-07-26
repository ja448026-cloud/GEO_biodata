#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: collect_publication_supplements.R /path/to/run-dir [pmid|doi]", call. = FALSE)
}

run_dir <- args[[1L]]
identifier <- args[[2L]]

required <- c("httr2", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

resources_dir <- file.path(run_dir, "resources")
if (!dir.exists(resources_dir)) dir.create(resources_dir, recursive = TRUE)

request_json <- function(url, query) {
  response <- httr2::request(url) |>
    httr2::req_url_query(!!!query) |>
    httr2::req_user_agent("geo-biodata-workflow/0.1") |>
    httr2::req_retry(max_tries = 3L) |>
    httr2::req_perform()
  jsonlite::fromJSON(httr2::resp_body_string(response), simplifyVector = FALSE)
}

is_pmid <- grepl("^[0-9]+$", identifier)
is_doi <- grepl("^10\\.", identifier)
is_pmcid <- grepl("^PMC[0-9]+$", identifier, ignore.case = TRUE)

if (!is_pmid && !is_doi && !is_pmcid) {
  stop("Identifier must be a PMID (digits), DOI (10.xxx), or PMCID (PMCxxx).", call. = FALSE)
}

query_field <- if (is_pmid) {
  paste0("EXT_ID:", identifier)
} else if (is_pmcid) {
  paste0("PMCID:", identifier)
} else {
  paste0("DOI:", gsub("/", "%2F", identifier))
}

cat(sprintf("Querying Europe PMC for %s...\n", identifier))

search_result <- request_json(
  "https://www.ebi.ac.uk/europepmc/webservices/rest/search",
  list(query = query_field, resultType = "core", format = "json", pageSize = "1")
)

results <- search_result$resultList$result
if (length(results) < 1L) {
  cat("No Europe PMC record found for", identifier, "\n")
  empty_df <- data.frame(
    identifier = character(), title = character(),
    has_supplementary = logical(), supplementary_count = integer(),
    stringsAsFactors = FALSE
  )
  utils::write.table(
    empty_df,
    file.path(resources_dir, "publication_supplements.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  quit(status = 0)
}

article <- results[[1L]]

pub_meta <- data.frame(
  identifier = identifier,
  title = if (is.null(article$title)) NA_character_ else article$title,
  journal = if (is.null(article$journalTitle)) NA_character_ else article$journalTitle,
  pub_year = if (is.null(article$pubYear)) NA_character_ else article$pubYear,
  doi = if (is.null(article$doi)) NA_character_ else article$doi,
  pmid = if (is.null(article$pmid)) NA_character_ else article$pmid,
  pmcid = if (is.null(article$pmcid)) NA_character_ else article$pmcid,
  is_open_access = if (is.null(article$isOpenAccess)) NA else article$isOpenAccess,
  stringsAsFactors = FALSE
)

utils::write.table(
  pub_meta,
  file.path(resources_dir, "publication_meta.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

# Collect supplementary file links when PMCID is available
pmcid <- article$pmcid
supp_data <- if (!is.null(pmcid) && nzchar(pmcid)) {
  cat(sprintf("Fetching supplementary material list for %s...\n", pmcid))

  supp_result <- tryCatch({
    request_json(
      "https://www.ebi.ac.uk/europepmc/webservices/rest/pmc",
      list(resultType = "supplementary", format = "json", id = pmcid)
    )
  }, error = function(e) {
    cat("Could not fetch supplementary list:", conditionMessage(e), "\n")
    NULL
  })

  supp_files <- supp_result$supplementaryFiles
  if (length(supp_files) > 0L) {
    tab <- do.call(rbind, lapply(supp_files, function(f) {
      data.frame(
        pmcid = pmcid,
        file_name = if (is.null(f$fileName)) NA_character_ else f$fileName,
        file_url = if (is.null(f$url)) NA_character_ else f$url,
        description = if (is.null(f$description)) NA_character_ else f$description,
        content_type = if (is.null(f$contentType)) NA_character_ else f$contentType,
        stringsAsFactors = FALSE
      )
    }))
  } else {
    data.frame(
      pmcid = pmcid,
      file_name = character(), file_url = character(),
      description = character(), content_type = character(),
      stringsAsFactors = FALSE
    )
  }
} else {
  cat("No PMCID available; cannot list supplementary files from publisher.\n")
  data.frame(
    pmcid = character(), file_name = character(), file_url = character(),
    description = character(), content_type = character(),
    stringsAsFactors = FALSE
  )
}

utils::write.table(
  supp_data,
  file.path(resources_dir, "publication_supplements.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

cat(sprintf(
  "Publication: %s\nOpen access: %s\nSupplementary files listed: %d\n",
  pub_meta$title,
  if (isTRUE(pub_meta$is_open_access)) "yes" else "no",
  nrow(supp_data)
))
