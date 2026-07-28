#!/usr/bin/env Rscript
file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
this_file <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
source(file.path(dirname(dirname(this_file)), "shared", "run_legacy.R"))
run_legacy_script(file.path("drivers", "run_bulk_counts.R"))
