#!/usr/bin/env Rscript

# Deprecated scRNA analysis template.
#
# This file is intentionally retained as a traceable compatibility marker, but
# it must not be used as an automatic analysis entry point. Earlier versions
# mixed raw-count processing with author-object review and included fallback
# paths that could densify large H5AD matrices. Those behaviors are unsafe for
# an agent-driven one-accession workflow.
#
# Use manifest-driven routes instead:
#   - scrna_raw_counts for count-bearing MTX/H5/H5AD/Seurat inputs;
#   - scrna_author_object for author-processed objects, preserving existing
#     metadata, labels, embeddings, and clusters by default.

stop(
  paste(
    "analyze_scrna_template.R is deprecated and fail-closed.",
    "Use run_manifest.yaml plus scrna_raw_counts or scrna_author_object drivers.",
    "Author objects must not be reclustered or re-embedded unless explicitly requested.",
    sep = "\n"
  ),
  call. = FALSE
)
