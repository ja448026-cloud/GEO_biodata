#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(1L, 2L)) {
  stop("Usage: merge_quant_tables.R /path/to/run-or-raw-dir [GSE000000]", call. = FALSE)
}

root <- normalizePath(args[[1L]], winslash = "/", mustWork = TRUE)
run_dir <- if (basename(root) == "raw") dirname(root) else root
raw_dir <- if (basename(root) == "raw") root else file.path(root, "raw")
if (!dir.exists(raw_dir)) stop("Raw directory does not exist: ", raw_dir, call. = FALSE)
accession <- if (length(args) == 2L) toupper(args[[2L]]) else basename(run_dir)

tables_dir <- file.path(run_dir, "tables")
derived_dir <- file.path(run_dir, "derived")
for (path in c(tables_dir, derived_dir)) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

read_quant <- function(path) {
  tryCatch(utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL)
}

sample_id_from <- function(path) {
  hit <- regmatches(basename(path), regexpr("GSM[0-9]+", basename(path), ignore.case = TRUE))
  if (length(hit) && nzchar(hit)) toupper(hit) else tools::file_path_sans_ext(basename(path))
}

find_col <- function(tab, patterns) {
  lname <- tolower(names(tab))
  for (pattern in patterns) {
    idx <- which(grepl(pattern, lname))
    if (length(idx)) return(names(tab)[idx[[1L]]])
  }
  NA_character_
}

merge_scale <- function(items, scale) {
  rows <- lapply(items, function(item) {
    vals <- item$tab[, c(item$feature_col, item[[scale]]), drop = FALSE]
    names(vals) <- c("feature_id", item$sample_id)
    vals
  })
  Reduce(function(x, y) merge(x, y, by = "feature_id", all = TRUE), rows)
}

write_matrix <- function(items, scale, suffix) {
  usable <- items[vapply(items, function(x) !is.na(x[[scale]]), logical(1))]
  if (!length(usable)) return("")
  mat <- merge_scale(usable, scale)
  out <- file.path(derived_dir, paste0(accession, "_", suffix, "_matrix.tsv"))
  utils::write.table(mat, out, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  file.path("derived", basename(out))
}

input_files <- list.files(raw_dir, pattern = "\\.(txt|tsv|csv)(\\.gz)?$", full.names = TRUE, ignore.case = TRUE)
items <- list()
metrics <- list()

for (path in sort(input_files)) {
  tab <- read_quant(path)
  if (is.null(tab) || ncol(tab) < 2L) next
  feature_col <- names(tab)[[1L]]
  expected_col <- find_col(tab, c("^expected_count$", "expected.?count"))
  tpm_col <- find_col(tab, c("^tpm$"))
  fpkm_col <- find_col(tab, c("^fpkm$"))
  if (all(is.na(c(expected_col, tpm_col, fpkm_col)))) next

  numeric_col <- function(col) {
    if (is.na(col)) return(numeric())
    suppressWarnings(as.numeric(tab[[col]]))
  }
  exp_vals <- numeric_col(expected_col)
  tpm_vals <- numeric_col(tpm_col)
  fpkm_vals <- numeric_col(fpkm_col)
  sample_id <- sample_id_from(path)
  items[[length(items) + 1L]] <- list(
    sample_id = sample_id,
    tab = tab,
    feature_col = feature_col,
    expected_count = expected_col,
    tpm = tpm_col,
    fpkm = fpkm_col
  )
  integer_frac <- if (length(exp_vals)) mean(abs(exp_vals - round(exp_vals)) < 1e-8, na.rm = TRUE) else NA_real_
  metrics[[length(metrics) + 1L]] <- data.frame(
    sample_id = sample_id,
    file_name = basename(path),
    n_features = nrow(tab),
    expected_count_sum = if (length(exp_vals)) sum(exp_vals, na.rm = TRUE) else NA_real_,
    tpm_sum = if (length(tpm_vals)) sum(tpm_vals, na.rm = TRUE) else NA_real_,
    fpkm_sum = if (length(fpkm_vals)) sum(fpkm_vals, na.rm = TRUE) else NA_real_,
    expected_count_integer_like_fraction = integer_frac,
    tpm_sum_flag = if (length(tpm_vals) && abs(sum(tpm_vals, na.rm = TRUE) - 1e6) > 1e4) "review" else "ok",
    expected_count_integer_flag = if (is.finite(integer_frac) && integer_frac < 0.9) "review" else "ok",
    stringsAsFactors = FALSE
  )
}

if (!length(items)) stop("No expected_count, TPM, or FPKM sample quantification tables found.", call. = FALSE)

metrics_df <- do.call(rbind, metrics)
utils::write.table(metrics_df, file.path(tables_dir, "sample_quant_metrics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "")

expected_path <- write_matrix(items, "expected_count", "expected_count")
tpm_path <- write_matrix(items, "tpm", "tpm")
fpkm_path <- write_matrix(items, "fpkm", "fpkm")

draft <- c(
  paste0("accession: ", accession),
  "route: bulk_normalized",
  "input:",
  paste0("  file: ", if (nzchar(tpm_path)) tpm_path else expected_path),
  "  input_type: tpm",
  "  scale:",
  "    transformed: false",
  "    evidence_source: tables/sample_quant_metrics.tsv",
  "    evidence_note: Sample quantification tables include expected_count, TPM, and/or FPKM; review before DE.",
  "analysis:",
  "  intent: eda_only",
  "review:",
  "  mode: manual",
  "  reviewed: false",
  "outputs:",
  paste0("  expected_count_matrix: ", expected_path),
  paste0("  tpm_matrix: ", tpm_path),
  paste0("  fpkm_matrix: ", fpkm_path)
)
writeLines(draft, file.path(run_dir, "run_manifest.quant.draft.yaml"), useBytes = TRUE)

cat("QUANT_TABLES_MERGED\n")
cat(normalizePath(file.path(tables_dir, "sample_quant_metrics.tsv"), winslash = "/", mustWork = TRUE), "\n")
