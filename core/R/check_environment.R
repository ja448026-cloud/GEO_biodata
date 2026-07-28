#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
output_file <- if (length(args) >= 1L) args[[1L]] else "environment.tsv"

packages <- data.frame(
  package = c(
    "GEOquery", "Biobase", "httr2", "jsonlite", "digest", "yaml",
    "limma", "edgeR", "DESeq2", "ggplot2", "pheatmap",
    "fgsea", "clusterProfiler", "msigdbr",
    "Seurat", "Matrix", "patchwork"
  ),
  capability = c(
    "core", "core", "core", "core", "core", "core",
    "bulk", "bulk", "bulk", "figures", "figures",
    "gsea", "gsea", "gsea",
    "scrna", "scrna", "figures"
  ),
  required = c(
    TRUE, TRUE, TRUE, TRUE, TRUE, TRUE,
    FALSE, FALSE, FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE,
    FALSE, FALSE, FALSE
  ),
  stringsAsFactors = FALSE
)

packages$installed <- vapply(
  packages$package,
  requireNamespace,
  logical(1),
  quietly = TRUE
)
packages$version <- vapply(
  packages$package,
  function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
    as.character(utils::packageVersion(pkg))
  },
  character(1)
)

output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
utils::write.table(
  packages,
  file = output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

cat("R version:", R.version.string, "\n")
cat("Environment table:", normalizePath(output_file, winslash = "/", mustWork = TRUE), "\n")

entrez_email <- Sys.getenv("ENTREZ_EMAIL", unset = "")
if (!nzchar(entrez_email)) {
  entrez_email <- getOption("Entrez.email", default = "")
}
if (!nzchar(entrez_email)) {
  cat("Warning: ENTREZ_EMAIL or options('Entrez.email') is not set. Set one before repeated GEO/NCBI requests.\n")
}

missing_core <- packages$package[packages$required & !packages$installed]
if (length(missing_core) > 0L) {
  cat("Missing core packages:", paste(missing_core, collapse = ", "), "\n")
  cat("Install Bioconductor packages with BiocManager and CRAN packages with install.packages().\n")
  quit(status = 2L)
}

cat("CORE_ENVIRONMENT_READY\n")
