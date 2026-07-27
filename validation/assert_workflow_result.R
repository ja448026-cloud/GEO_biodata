#!/usr/bin/env Rscript
#
# assert_workflow_result.R
# Result-level CI assertion for geo_biodata_workflow drivers.
# Usage: Rscript assert_workflow_result.R <run_dir> <expected_execution_state>
#   [--expect-de] [--expect-eda] [--expect-blocked]
#   [--check-logfc-direction] [--expect-mapping <status>]
#   [--expect-contract <state>] [--expect-no-de-output]

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: assert_workflow_result.R <run_dir> <expected_execution_state> [options]", call. = FALSE)
}

run_dir <- args[[1L]]
expected_execution <- args[[2L]]
opts <- args[-(1:2)]

fail <- function(msg) { stop(msg, call. = FALSE) }

has_opt <- function(flag) any(opts == flag)
get_opt_val <- function(flag) {
  idx <- which(opts == flag)
  if (length(idx) == 0L || idx == length(opts)) NA_character_
  else opts[[idx + 1L]]
}

# ── Check workflow_status.tsv ────────────────────────────────────────────────

status_path <- file.path(run_dir, "workflow_status.tsv")
if (!file.exists(status_path)) fail("workflow_status.tsv does not exist.")
status <- utils::read.delim(status_path, stringsAsFactors = FALSE)

# Execution state (supports execution_state, state, or inventory_state)
exec_state <- if ("execution_state" %in% names(status)) {
  status$execution_state[[1L]]
} else if ("state" %in% names(status)) {
  status$state[[1L]]
} else if ("inventory_state" %in% names(status)) {
  status$inventory_state[[1L]]
} else {
  fail("workflow_status.tsv has no execution_state, state, or inventory_state column.")
}

if (exec_state != expected_execution) {
  fail(sprintf("Expected execution_state='%s', got '%s'.", expected_execution, exec_state))
}
cat(sprintf("PASS execution_state=%s\n", exec_state))

# ── Optional: contract_state ─────────────────────────────────────────────────

expected_contract <- get_opt_val("--expect-contract")
if (!is.na(expected_contract)) {
  if (!"contract_state" %in% names(status)) {
    fail("--expect-contract specified but workflow_status.tsv has no contract_state column.")
  }
  cs <- status$contract_state[[1L]]
  if (cs != expected_contract) {
    fail(sprintf("Expected contract_state='%s', got '%s'.", expected_contract, cs))
  }
  cat(sprintf("PASS contract_state=%s\n", cs))
}

# ── Optional: technical QC ──────────────────────────────────────────────────

expected_tqc <- get_opt_val("--expect-technical-qc")
if (!is.na(expected_tqc)) {
  if (!"technical_qc" %in% names(status)) {
    fail("--expect-technical-qc specified but workflow_status.tsv has no technical_qc column.")
  }
  tqc <- status$technical_qc[[1L]]
  if (tqc != expected_tqc) {
    fail(sprintf("Expected technical_qc='%s', got '%s'.", expected_tqc, tqc))
  }
  cat(sprintf("PASS technical_qc=%s\n", tqc))
}

# ── Optional: result_signal ─────────────────────────────────────────────────

expected_rs <- get_opt_val("--expect-result-signal")
if (!is.na(expected_rs)) {
  if (!"result_signal" %in% names(status)) {
    fail("--expect-result-signal specified but workflow_status.tsv has no result_signal column.")
  }
  rs <- status$result_signal[[1L]]
  if (rs != expected_rs) {
    fail(sprintf("Expected result_signal='%s', got '%s'.", expected_rs, rs))
  }
  cat(sprintf("PASS result_signal=%s\n", rs))
}

# ── Optional: mapping_status ─────────────────────────────────────────────────

expected_map <- get_opt_val("--expect-mapping")
if (!is.na(expected_map)) {
  if (!"mapping_status" %in% names(status)) {
    fail("--expect-mapping specified but workflow_status.tsv has no mapping_status column.")
  }
  ms <- status$mapping_status[[1L]]
  if (ms != expected_map) {
    fail(sprintf("Expected mapping_status='%s', got '%s'.", expected_map, ms))
  }
  cat(sprintf("PASS mapping_status=%s\n", ms))
}

# ── Optional: inventory_state ────────────────────────────────────────────────

expected_inv <- get_opt_val("--expect-inventory")
if (!is.na(expected_inv)) {
  if (!"inventory_state" %in% names(status)) {
    fail("--expect-inventory specified but workflow_status.tsv has no inventory_state column.")
  }
  ist <- status$inventory_state[[1L]]
  if (ist != expected_inv) {
    fail(sprintf("Expected inventory_state='%s', got '%s'.", expected_inv, ist))
  }
  cat(sprintf("PASS inventory_state=%s\n", ist))
}

# ── EDA-only: must not have DE outputs ───────────────────────────────────────

if (has_opt("--expect-no-de-output")) {
  de_files <- list.files(file.path(run_dir, "tables"), pattern = "^de_results_", full.names = TRUE)
  if (length(de_files) > 0L) {
    fail(sprintf("EDA-only should not generate DE tables, found: %s", paste(basename(de_files), collapse = ", ")))
  }
  pval_files <- list.files(file.path(run_dir, "figures"), pattern = "bulk_pvalue_histogram|bulk_meanvar", full.names = TRUE)
  if (length(pval_files) > 0L) {
    fail(sprintf("EDA-only should not generate p-value/meanvar plots, found: %s", paste(basename(pval_files), collapse = ", ")))
  }
  cat("PASS no DE outputs in EDA-only mode\n")
}

# ── Blocked: must not have DE outputs ────────────────────────────────────────

if (has_opt("--expect-blocked")) {
  de_files <- list.files(file.path(run_dir, "tables"), pattern = "^de_results_", full.names = TRUE)
  if (length(de_files) > 0L) {
    fail(sprintf("BLOCKED should not generate DE tables, found: %s", paste(basename(de_files), collapse = ", ")))
  }
  cat("PASS blocked state: no DE outputs\n")
}

# ── Optional: check logFC direction ──────────────────────────────────────────

if (has_opt("--check-logfc")) {
  gene_name <- get_opt_val("--gene")
  expected_sign <- get_opt_val("--sign")
  if (is.na(gene_name) || is.na(expected_sign)) {
    fail("--check-logfc requires --gene <name> and --sign <positive|negative>")
  }
  de_files <- list.files(file.path(run_dir, "tables"), pattern = "^de_results_", full.names = TRUE)
  if (length(de_files) == 0L) fail("No DE results found for logFC check.")
  de <- utils::read.delim(de_files[[1L]], stringsAsFactors = FALSE)
  id_col <- if ("feature_id" %in% names(de)) "feature_id" else names(de)[[1L]]
  logfc_col <- if ("logFC" %in% names(de)) "logFC" else names(de)[grep("logFC|log2FoldChange", names(de))][[1L]]

  if (!gene_name %in% de[[id_col]]) {
    fail(sprintf("Gene '%s' not found in DE results.", gene_name))
  }
  lfc <- de[[logfc_col]][de[[id_col]] == gene_name]
  sign_ok <- if (expected_sign == "positive") lfc > 0 else lfc < 0
  if (!sign_ok) {
    fail(sprintf("Gene '%s' logFC=%.3f, expected %s.", gene_name, lfc, expected_sign))
  }
  cat(sprintf("PASS logFC direction: %s=%.3f (%s)\n", gene_name, lfc, expected_sign))
}

# ── Check critical output files exist ────────────────────────────────────────

expected_outputs <- get_opt_val("--expect-outputs")
if (!is.na(expected_outputs)) {
  for (fpattern in strsplit(expected_outputs, ",")[[1L]]) {
    found <- list.files(run_dir, pattern = fpattern, recursive = TRUE, full.names = TRUE)
    if (length(found) == 0L) {
      fail(sprintf("Expected output matching '%s' not found.", fpattern))
    }
    for (fp in found) {
      if (file.info(fp)$size == 0) fail(sprintf("Output file '%s' is empty.", fp))
    }
  }
  cat(sprintf("PASS expected outputs found: %s\n", expected_outputs))
}

# ── Output integrity check ───────────────────────────────────────────────────

integrity_path <- file.path(run_dir, "tables", "output_integrity.tsv")
if (file.exists(integrity_path)) {
  integ <- utils::read.delim(integrity_path, stringsAsFactors = FALSE)
  missing <- integ[integ$status != "OK", ]
  if (nrow(missing) > 0L) {
    cat("WARNING: Some outputs are missing or empty:\n")
    print(missing)
  } else {
    cat("PASS output_integrity: all outputs present and non-empty\n")
  }
}

cat(sprintf("ALL_ASSERTIONS_PASS for execution_state=%s\n", expected_execution))
