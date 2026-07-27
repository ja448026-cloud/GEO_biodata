#!/usr/bin/env Rscript

# Deprecated bulk analysis template.
#
# This file is intentionally retained as a traceable compatibility marker, but
# it must not be used as an automatic analysis entry point. Earlier versions
# mixed raw-count, normalized, and microarray routing in one script. That made
# it too easy for an agent to choose a statistical model from heuristics instead
# of an explicit, reviewed run manifest.
#
# Use a manifest-driven route instead:
#   1. create run_manifest.yaml from templates/run_manifest.example.yaml;
#   2. validate it with scripts/validate_manifest.R;
#   3. run a route-specific driver for bulk_raw_counts, bulk_normalized, or
#      microarray_series_matrix.

stop(
  paste(
    "analyze_bulk_template.R is deprecated and fail-closed.",
    "Use run_manifest.yaml plus a route-specific bulk driver.",
    "Do not choose DESeq2, edgeR, or limma from value-range heuristics.",
    sep = "\n"
  ),
  call. = FALSE
)
