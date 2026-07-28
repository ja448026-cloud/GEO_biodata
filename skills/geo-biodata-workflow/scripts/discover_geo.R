#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop("Usage: discover_geo.R GSE000000 /path/to/run-directory [--no-characteristics]", call. = FALSE)
}

accession <- toupper(args[[1L]])
run_dir <- args[[2L]]
skip_characteristics <- length(args) == 3L && identical(args[[3L]], "--no-characteristics")
if (!grepl("^GSE[0-9]+$", accession)) {
  stop("The accession must match GSE followed by digits.", call. = FALSE)
}

required <- c("GEOquery", "httr2", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0L) {
  stop("Missing required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

ncbi_email <- Sys.getenv("ENTREZ_EMAIL", unset = "")
if (!nzchar(ncbi_email)) {
  ncbi_email <- getOption("Entrez.email", default = "")
}
if (!nzchar(ncbi_email)) {
  cat("Warning: ENTREZ_EMAIL or options('Entrez.email') is not set; set one for repeated NCBI requests.\n")
}

resources_dir <- file.path(run_dir, "resources")
for (path in c(
  run_dir, resources_dir, file.path(run_dir, "raw"), file.path(run_dir, "derived"),
  file.path(run_dir, "tables"), file.path(run_dir, "figures"),
  file.path(run_dir, "logs"), file.path(run_dir, "scripts")
)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
}

workflow_events <- data.frame(
  stage = character(),
  severity = character(),
  status = character(),
  message = character(),
  stringsAsFactors = FALSE
)

add_event <- function(stage, severity, status, message) {
  workflow_events <<- rbind(
    workflow_events,
    data.frame(
      stage = stage,
      severity = severity,
      status = status,
      message = message,
      stringsAsFactors = FALSE
    )
  )
}

request_json <- function(url, query) {
  request <- httr2::request(url) |>
    httr2::req_url_query(!!!query) |>
    httr2::req_user_agent("geo-biodata-workflow/0.1")
  if (nzchar(ncbi_email)) {
    request <- request |>
      httr2::req_headers("NCBI-Email" = ncbi_email)
  }
  response <- request |>
    httr2::req_retry(max_tries = 3L) |>
    httr2::req_perform()
  jsonlite::fromJSON(httr2::resp_body_string(response), simplifyVector = FALSE)
}

collapse_value <- function(value) {
  value <- unlist(value, recursive = TRUE, use.names = FALSE)
  value <- value[!is.na(value)]
  paste(as.character(value), collapse = "; ")
}

# ── Series-level metadata via NCBI E-utilities ───────────────────────────

search <- request_json(
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
  list(db = "gds", term = paste0(accession, "[ACCN]"), retmode = "json")
)
uids <- unlist(search$esearchresult$idlist, use.names = FALSE)
if (length(uids) < 1L) stop("NCBI GEO did not return a record for ", accession, call. = FALSE)

summary_result <- request_json(
  "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi",
  list(db = "gds", id = uids[[1L]], retmode = "json")
)
record <- summary_result$result[[uids[[1L]]]]

metadata_names <- setdiff(names(record), "uid")
metadata <- data.frame(
  field = metadata_names,
  value = vapply(record[metadata_names], collapse_value, character(1)),
  stringsAsFactors = FALSE
)
utils::write.table(
  metadata,
  file.path(resources_dir, "series_metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

# ── Sample index from esummary ───────────────────────────────────────────

sample_records <- record$samples
samples <- if (length(sample_records) > 0L) {
  data.frame(
    gsm = vapply(sample_records, function(x) collapse_value(x$accession), character(1)),
    title = vapply(sample_records, function(x) collapse_value(x$title), character(1)),
    stringsAsFactors = FALSE
  )
} else {
  data.frame(gsm = character(), title = character())
}
utils::write.table(
  samples,
  file.path(resources_dir, "sample_index.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

# ── Sample characteristics via GEOquery SOFT parsing ─────────────────────

characteristics <- data.frame(
  gsm = character(), characteristic = character(), value = character(),
  stringsAsFactors = FALSE
)

if (!skip_characteristics && length(sample_records) > 0L) {
  cat("Extracting sample characteristics from GEO SOFT file...\n")
  characteristics_status <- "OK"
  characteristics_message <- ""
  characteristics <- tryCatch({
    gse_soft <- GEOquery::getGEO(
      accession,
      GSEMatrix = FALSE,
      getGPL = FALSE,
      AnnotGPL = FALSE
    )
    gsm_list <- GEOquery::GSMList(gse_soft)
    all_chars <- list()
    for (gsm_name in names(gsm_list)) {
      gsm_meta <- GEOquery::Meta(gsm_list[[gsm_name]])
      char_keys <- grep("^characteristics_ch", names(gsm_meta), value = TRUE)
      if (length(char_keys) == 0L) next
      for (key in char_keys) {
        raw_values <- gsm_meta[[key]]
        for (raw in raw_values) {
          parts <- strsplit(raw, ":\\s*", perl = TRUE)[[1L]]
          if (length(parts) >= 2L) {
            all_chars[[length(all_chars) + 1L]] <- data.frame(
              gsm = gsm_name,
              characteristic = trimws(parts[[1L]]),
              value = trimws(paste(parts[-1L], collapse = ": ")),
              stringsAsFactors = FALSE
            )
          } else if (nzchar(trimws(raw))) {
            all_chars[[length(all_chars) + 1L]] <- data.frame(
              gsm = gsm_name,
              characteristic = key,
              value = trimws(raw),
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }
    if (length(all_chars) > 0L) {
      do.call(rbind, all_chars)
    } else {
      data.frame(
        gsm = character(), characteristic = character(), value = character(),
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) {
    characteristics_status <<- "FAILED"
    characteristics_message <<- conditionMessage(e)
    cat("Warning: Could not extract sample characteristics:", characteristics_message, "\n")
    add_event(
      "sample_characteristics",
      "error",
      "FAILED",
      characteristics_message
    )
    data.frame(
      gsm = character(), characteristic = character(), value = character(),
      stringsAsFactors = FALSE
    )
  })
} else if (skip_characteristics && length(sample_records) > 0L) {
  characteristics_status <- "SKIPPED"
  characteristics_message <- "SOFT sample-characteristics extraction was skipped by --no-characteristics."
  add_event(
    "sample_characteristics",
    "warning",
    "SKIPPED",
    characteristics_message
  )
} else {
  characteristics_status <- "NO_SAMPLES"
  characteristics_message <- "No sample records were available for characteristic extraction."
  add_event(
    "sample_characteristics",
    "error",
    "NO_SAMPLES",
    characteristics_message
  )
}

if (identical(characteristics_status, "OK") && length(sample_records) > 0L && nrow(characteristics) == 0L) {
  characteristics_status <- "EMPTY"
  characteristics_message <- "No structured characteristics were extracted from GEO SOFT."
  add_event(
    "sample_characteristics",
    "warning",
    "EMPTY",
    characteristics_message
  )
}

if (nrow(characteristics) > 0L) {
  char_wide <- reshape(
    characteristics,
    idvar = "gsm",
    timevar = "characteristic",
    direction = "wide"
  )
  names(char_wide) <- gsub("^value\\.", "", names(char_wide))
} else {
  char_wide <- data.frame(
    gsm = if (nrow(samples) > 0L) samples$gsm else character(),
    note = rep("No structured characteristics extracted", nrow(samples)),
    stringsAsFactors = FALSE
  )
  if (nrow(char_wide) == 0L) {
    char_wide <- data.frame(gsm = character(), note = character(), stringsAsFactors = FALSE)
  }
}

utils::write.table(
  char_wide,
  file.path(resources_dir, "sample_characteristics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

# ── Geo-series-type summary for routing ──────────────────────────────────

detected_fields <- names(record)
assay_hint <- NA_character_
if ("gdstype" %in% detected_fields) {
  assay_hint <- collapse_value(record$gdstype)
}

series_relation <- if (!is.null(record$relations)) collapse_value(record$relations) else ""
superseries_warning <- grepl("super.?series|sub.?series", series_relation, ignore.case = TRUE)

platform_ids <- if (!is.null(record$gpl)) {
  gsub("[^0-9]", "", unlist(record$gpl, recursive = TRUE, use.names = FALSE))
} else {
  character()
}

cat(sprintf(
  "Accession: %s\nAssay hint: %s\nPlatform(s): %s\nSamples: %d\n",
  accession,
  if (is.na(assay_hint)) "unknown" else assay_hint,
  if (length(platform_ids) > 0L) paste(platform_ids, collapse = ", ") else "unknown",
  nrow(samples)
))
if (superseries_warning) {
  cat("Series relation warning: this record appears to reference a SuperSeries/SubSeries relationship; review related GSE accessions before analysis.\n")
}

# ── Supplementary file index ─────────────────────────────────────────────

cat("Listing GEO supplementary files...\n")
supplements_status <- "OK"
supplements_message <- ""
supplements <- tryCatch(
  GEOquery::getGEOSuppFiles(
    accession,
    makeDirectory = FALSE,
    baseDir = resources_dir,
    fetch_files = FALSE
  ),
  error = function(error) {
    supplements_status <<- "FAILED"
    supplements_message <<- conditionMessage(error)
    add_event(
      "supplements",
      "error",
      "FAILED",
      supplements_message
    )
    data.frame(
      status = "SUPPLEMENT_QUERY_FAILED",
      message = supplements_message,
      stringsAsFactors = FALSE
    )
  }
)
supplement_rows <- rownames(supplements)
supplements <- data.frame(supplements, check.names = FALSE)
if (!identical(supplements_status, "FAILED") && (ncol(supplements) == 0L || nrow(supplements) == 0L)) {
  supplements_status <- "NO_SUPPLEMENT_FILES"
  supplements_message <- "GEOquery returned no series-level supplementary files."
  add_event(
    "supplements",
    "info",
    "NO_SUPPLEMENT_FILES",
    supplements_message
  )
  supplements <- data.frame(
    status = "NO_SUPPLEMENT_FILES",
    message = supplements_message,
    stringsAsFactors = FALSE
  )
}
if (identical(supplements_status, "OK") && nrow(supplements) > 0L) {
  row_is_url <- grepl("^(https?|ftp)://", supplement_rows, ignore.case = TRUE)
  url_col <- if ("url" %in% names(supplements)) as.character(supplements$url) else rep("", nrow(supplements))
  url_col[is.na(url_col)] <- ""
  col_is_url <- grepl("^(https?|ftp)://", url_col, ignore.case = TRUE)
  supplements$supplement_url <- ifelse(row_is_url, supplement_rows, ifelse(col_is_url, url_col, NA_character_))
  supplements$file_name <- ifelse(
    row_is_url,
    basename(supplement_rows),
    if ("fname" %in% names(supplements)) basename(supplements$fname) else supplement_rows
  )
}
if ("size" %in% names(supplements)) {
  supplements$size_mb <- round(supplements$size / 1e6, 2)
}
utils::write.table(
  supplements,
  file.path(resources_dir, "supplement_index.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

# ── Publication links ────────────────────────────────────────────────────

pubmed_ids <- unique(unlist(record$pubmedids, recursive = TRUE, use.names = FALSE))
pubmed_ids <- pubmed_ids[nzchar(pubmed_ids)]
publications <- if (length(pubmed_ids) > 0L) {
  pmc_ids <- vapply(pubmed_ids, function(pmid) {
    tryCatch({
      resp <- httr2::request(
        "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
      ) |>
        httr2::req_url_query(
          query = paste0("EXT_ID:", pmid),
          resultType = "core",
          format = "json"
        ) |>
        httr2::req_user_agent("geo-biodata-workflow/0.1") |>
        httr2::req_retry(max_tries = 2L) |>
        httr2::req_perform()
      result <- jsonlite::fromJSON(
        httr2::resp_body_string(resp),
        simplifyVector = FALSE
      )$resultList$result
      if (length(result) > 0L) {
        pmcid <- result[[1L]]$pmcid
        if (is.null(pmcid)) NA_character_ else pmcid
      } else {
        NA_character_
      }
    }, error = function(e) NA_character_)
  }, character(1), USE.NAMES = FALSE)

  data.frame(
    identifier_type = "PMID",
    identifier = pubmed_ids,
    pmcid = pmc_ids,
    url = paste0("https://pubmed.ncbi.nlm.nih.gov/", pubmed_ids, "/"),
    open_access_url = ifelse(
      is.na(pmc_ids), NA_character_,
      paste0("https://www.ncbi.nlm.nih.gov/pmc/articles/", pmc_ids, "/")
    ),
    stringsAsFactors = FALSE
  )
} else {
  data.frame(
    identifier_type = character(), identifier = character(),
    pmcid = character(), url = character(),
    open_access_url = character(), stringsAsFactors = FALSE
  )
}
utils::write.table(
  publications,
  file.path(resources_dir, "publication_links.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

# ── SRA/BioProject links for raw-data handoff ────────────────────────────

bioproject <- if (!is.null(record$bioproject)) collapse_value(record$bioproject) else ""
sra_links <- data.frame(
  accession = accession,
  bioproject = bioproject,
  sra_url = if (nzchar(bioproject)) {
    paste0("https://www.ncbi.nlm.nih.gov/sra/?term=", bioproject)
  } else {
    "not_available"
  },
  fetchngs_hint = if (nzchar(bioproject)) {
    paste0(
      "Write ", bioproject, " to a one-line accessions file, then run: ",
      "nf-core fetchngs --input accessions.txt --outdir raw_fastq/"
    )
  } else {
    "not_available"
  },
  stringsAsFactors = FALSE
)
utils::write.table(
  sra_links,
  file.path(resources_dir, "sra_links.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

# ── Routing hint ─────────────────────────────────────────────────────────

assay_str <- if (is.na(assay_hint)) "" else tolower(assay_hint)
infer_primary_route <- function(assay_text) {
  if (grepl("single.cell|scrna|10x|droplet|smart-seq", assay_text, ignore.case = TRUE)) {
    return("scrna_raw_counts")
  }
  if (grepl("array", assay_text, ignore.case = TRUE)) {
    return("microarray_series_matrix")
  }
  if (grepl("high.throughput|rna.seq|rna-seq|sequencing", assay_text, ignore.case = TRUE)) {
    return("metadata_only")
  }
  "metadata_only"
}

routing <- data.frame(
  accession = accession,
  assay_type = if (is.na(assay_hint)) "unknown" else assay_hint,
  series_relation = series_relation,
  superseries_or_subseries = superseries_warning,
  n_samples = nrow(samples),
  n_platforms = length(platform_ids),
  has_supplements = identical(supplements_status, "OK") && nrow(supplements) > 0L,
  sample_characteristics_status = characteristics_status,
  supplement_status = supplements_status,
  likely_scRNA = grepl(
    "single.cell|scrna|scRNA|10x|droplet|smart-seq",
    assay_str, ignore.case = TRUE
  ),
  likely_bulk = grepl(
    "expression.profiling.by.array|expression.profiling.by.high.throughput",
    assay_str, ignore.case = TRUE
  ),
  recommended_route = infer_primary_route(assay_str),
  stringsAsFactors = FALSE
)
utils::write.table(
  routing,
  file.path(resources_dir, "routing_hint.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

# ── Routing evidence ─────────────────────────────────────────────────────

evidence_rows <- list()
add_evidence <- function(candidate_route, source, value, weight, supports, conflicts, confidence, review_required, note) {
  evidence_rows[[length(evidence_rows) + 1L]] <<- data.frame(
    candidate_route = candidate_route,
    evidence_source = source,
    evidence_value = value,
    evidence_weight = weight,
    supports = supports,
    conflicts = conflicts,
    confidence = confidence,
    review_required = review_required,
    note = note,
    stringsAsFactors = FALSE
  )
}

if (nzchar(assay_str)) {
  add_evidence(
    routing$recommended_route,
    "gds_assay_type",
    routing$assay_type,
    0.35,
    TRUE,
    FALSE,
    if (routing$recommended_route == "metadata_only") 0.3 else 0.65,
    routing$recommended_route == "metadata_only",
    "GDS assay type is useful for broad routing but does not prove input scale or sample mapping."
  )
  if (grepl("high.throughput|rna.seq|rna-seq|sequencing", assay_str, ignore.case = TRUE)) {
    add_evidence(
      "bulk_raw_counts",
      "gds_assay_type",
      routing$assay_type,
      0.2,
      TRUE,
      FALSE,
      0.35,
      TRUE,
      "Sequencing assay type can support a raw-count route only after content inspection confirms integer counts."
    )
    add_evidence(
      "bulk_normalized",
      "gds_assay_type",
      routing$assay_type,
      0.15,
      TRUE,
      FALSE,
      0.3,
      TRUE,
      "Sequencing assay type can also produce normalized author matrices; scale remains unresolved."
    )
  }
}
if (identical(characteristics_status, "OK") && nrow(characteristics) > 0L) {
  add_evidence(
    routing$recommended_route,
    "sample_characteristics",
    paste0(nrow(characteristics), " characteristic entries"),
    0.25,
    TRUE,
    FALSE,
    0.55,
    TRUE,
    "Sample characteristics exist, but biological groups and units still require review."
  )
} else {
  add_evidence(
    "metadata_only",
    "sample_characteristics",
    characteristics_status,
    0.3,
    FALSE,
    TRUE,
    0.8,
    TRUE,
    "Missing, skipped, or empty sample characteristics prevent fully automatic analysis routing."
  )
}
if (identical(supplements_status, "OK") && nrow(supplements) > 0L) {
  supplement_text <- paste(tolower(supplements$file_name), collapse = " ")
  if (grepl("barcodes|features|matrix\\.mtx|filtered_feature_bc_matrix|10x", supplement_text)) {
    add_evidence("scrna_raw_counts", "supplement_filenames", supplement_text, 0.45, TRUE, FALSE, 0.75, TRUE, "10x-style filenames suggest scRNA counts; object content still needs inspection.")
  }
  if (grepl("h5ad|seurat|rds|rdata", supplement_text)) {
    add_evidence("scrna_author_object", "supplement_filenames", supplement_text, 0.4, TRUE, FALSE, 0.65, TRUE, "Author object filenames require inspection of raw-count layers and metadata.")
  }
  if (grepl("count|counts|matrix|expression|series", supplement_text)) {
    add_evidence("bulk_raw_counts", "supplement_filenames", supplement_text, 0.25, TRUE, FALSE, 0.45, TRUE, "Counts-like filenames require content inspection before raw-count DE.")
    add_evidence("bulk_normalized", "supplement_filenames", supplement_text, 0.2, TRUE, FALSE, 0.4, TRUE, "Expression-like filenames may indicate normalized values rather than raw counts.")
  }
  if (grepl("series_matrix|matrix\\.txt|gse.*series", supplement_text)) {
    add_evidence("microarray_series_matrix", "supplement_filenames", supplement_text, 0.25, TRUE, FALSE, 0.45, TRUE, "Series matrix filenames require platform and normalization review.")
  }
} else {
  add_evidence("metadata_only", "supplements", supplements_status, 0.2, FALSE, FALSE, 0.4, TRUE, "No downloadable supplement evidence was available from GEO series-level listing.")
}
if (superseries_warning) {
  add_evidence(
    "metadata_only",
    "series_relation",
    series_relation,
    0.6,
    FALSE,
    TRUE,
    0.9,
    TRUE,
    "SuperSeries/SubSeries records require choosing the analysis unit before download or analysis."
  )
}
if (length(evidence_rows) == 0L) {
  add_evidence("metadata_only", "none", "no routing evidence generated", 0, FALSE, FALSE, 0, TRUE, "Review GEO page manually.")
}
routing_evidence <- do.call(rbind, evidence_rows)
utils::write.table(
  routing_evidence,
  file.path(resources_dir, "routing_evidence.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

global_route_review_required <- (
  skip_characteristics ||
    !identical(characteristics_status, "OK") ||
    superseries_warning ||
    !identical(supplements_status, "OK")
)

route_candidates <- do.call(
  rbind,
  lapply(split(routing_evidence, routing_evidence$candidate_route), function(df) {
    support_score <- sum(df$evidence_weight[df$supports], na.rm = TRUE)
    conflict_score <- sum(df$evidence_weight[df$conflicts], na.rm = TRUE)
    confidence <- max(df$confidence, na.rm = TRUE)
    review_required <- global_route_review_required || any(df$review_required) || conflict_score > 0 || df$candidate_route[[1L]] == "metadata_only"
    confidence_label <- if (confidence >= 0.75 && !review_required) {
      "high"
    } else if (confidence >= 0.5) {
      "medium"
    } else {
      "low"
    }
    data.frame(
      route = df$candidate_route[[1L]],
      support_score = round(support_score, 3),
      conflict_score = round(conflict_score, 3),
      confidence = round(confidence, 3),
      confidence_label = confidence_label,
      review_required = review_required,
      evidence_sources = paste(unique(df$evidence_source), collapse = ";"),
      note = paste(unique(df$note), collapse = " | "),
      stringsAsFactors = FALSE
    )
  })
)
route_candidates <- route_candidates[order(-route_candidates$support_score, route_candidates$conflict_score, route_candidates$route), ]
route_candidates$decision <- "secondary"
if (nrow(route_candidates) > 0L) {
  if (all(route_candidates$review_required)) {
    route_candidates$decision[[1L]] <- "review_required"
  } else {
    route_candidates$decision[[1L]] <- "selected"
  }
}
utils::write.table(
  route_candidates,
  file.path(resources_dir, "route_candidates.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

selected_candidate <- if (nrow(route_candidates) > 0L) route_candidates$route[[1L]] else "metadata_only"
analysis_decisions <- data.frame(
  decision_id = paste0(accession, "_route_candidate_001"),
  decision_type = "route_candidate",
  selected_value = selected_candidate,
  alternative_values = paste(setdiff(route_candidates$route, selected_candidate), collapse = ";"),
  evidence_sources = paste(unique(routing_evidence$evidence_source), collapse = ";"),
  source_tier = "public_metadata",
  supporting_text = paste(
    "Assay hint:", if (is.na(assay_hint)) "unknown" else assay_hint,
    "| characteristics:", characteristics_status,
    "| supplements:", supplements_status
  ),
  conflicting_evidence = if (any(routing_evidence$conflicts)) {
    paste(unique(routing_evidence$note[routing_evidence$conflicts]), collapse = " | ")
  } else {
    ""
  },
  confidence = if (nrow(route_candidates) > 0L) route_candidates$confidence[[1L]] else 0,
  decision_rule = "superseries_selection_v0.1",
  agent_model = "rules_only",
  conflict_status = if (any(routing_evidence$conflicts)) "unresolved" else "none",
  evidence_sources = paste(unique(routing_evidence$evidence_source), collapse = ";"),
  retrieved_at = now_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE),
  decided_at = now_utc,
  requires_user_input = if (nrow(route_candidates) > 0L) route_candidates$review_required[[1L]] else TRUE,
  stringsAsFactors = FALSE
)
utils::write.table(
  analysis_decisions,
  file.path(resources_dir, "analysis_decisions.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

# ── Status ───────────────────────────────────────────────────────────────

now_utc <- format(Sys.time(), tz = "UTC", usetz = TRUE)
n_chars <- nrow(characteristics)
n_char_fields <- if (nrow(char_wide) > 0L) ncol(char_wide) - 1L else 0L

if (superseries_warning) {
  add_event(
    "routing",
    "warning",
    "REVIEW_REQUIRED",
    "SuperSeries/SubSeries relation requires manual review before choosing the analysis unit."
  )
}

status_state <- "RESOURCE_INVENTORY_COMPLETE"
if (nrow(samples) < 1L) {
  status_state <- "BLOCKED_METADATA"
} else if (any(workflow_events$severity == "error")) {
  status_state <- "DISCOVERY_PARTIAL"
} else if (
  skip_characteristics ||
  identical(characteristics_status, "EMPTY") ||
  superseries_warning
) {
  status_state <- "REVIEW_REQUIRED"
}

if (nrow(workflow_events) == 0L) {
  workflow_events <- data.frame(
    stage = "discovery",
    severity = "info",
    status = "OK",
    message = "No discovery warnings or errors were recorded.",
    stringsAsFactors = FALSE
  )
}
workflow_events$updated_at_utc <- now_utc
utils::write.table(
  workflow_events,
  file.path(run_dir, "workflow_events.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)

status <- data.frame(
  accession = accession,
  state = status_state,
  updated_at_utc = now_utc,
  note = paste0(
    "Review metadata and supplement_index.tsv before downloading. ",
    "Extracted characteristics for ", n_chars, " entries across ", n_char_fields, " fields.",
    " Characteristic status: ", characteristics_status, ".",
    " Supplement status: ", supplements_status, ".",
    if (superseries_warning) " SuperSeries/SubSeries relation requires manual review." else ""
  ),
  stringsAsFactors = FALSE
)
utils::write.table(
  status,
  file.path(run_dir, "workflow_status.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

summary_lines <- c(
  paste("#", accession, "resource inventory"),
  "",
  paste("Generated:", now_utc),
  "",
  paste0("Status: `", status_state, "`"),
  "",
  sprintf("Assay: %s | Samples: %d | Platforms: %s",
          if (is.na(assay_hint)) "unknown" else assay_hint,
          nrow(samples),
          if (length(platform_ids) > 0L) paste(platform_ids, collapse = ", ") else "unknown"),
  "",
  paste("Recommended route:", routing$recommended_route),
  "",
  paste("Sample characteristics status:", characteristics_status),
  paste("Supplement status:", supplements_status),
  "",
  "Review the following before selecting an input route:",
  "",
  "- `resources/series_metadata.tsv` — GEO series-level metadata",
  "- `resources/sample_index.tsv` — GSM accessions and titles",
  "- `resources/sample_characteristics.tsv` — extracted sample-level characteristics",
  "- `resources/supplement_index.tsv` — GEO supplementary file listing",
  "- `resources/publication_links.tsv` — publication identifiers and open-access URLs",
  "- `resources/routing_hint.tsv` — automatically inferred assay type and route",
  "- `resources/routing_evidence.tsv` — evidence rows behind route candidates",
  "- `resources/route_candidates.tsv` — scored route candidates using the unified route ontology",
  "- `resources/analysis_decisions.tsv` — deterministic decision record for agent/user adjudication",
  "- `resources/sra_links.tsv` — SRA/BioProject links for raw-data handoff",
  "- `workflow_events.tsv` — warnings and partial-failure events recorded during discovery",
  "",
  if (superseries_warning) {
    "SuperSeries/SubSeries warning: review related accessions before choosing the analysis unit."
  } else {
    NULL
  },
  "If sample characteristics were not extractable, review the GEO page directly.",
  "Use `--no-characteristics` to skip the SOFT download for very large series."
)
writeLines(summary_lines, file.path(run_dir, "summary.md"), useBytes = TRUE)

cat(status_state, "\n", sep = "")
cat(normalizePath(run_dir, winslash = "/", mustWork = TRUE), "\n")
