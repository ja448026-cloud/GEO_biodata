#!/usr/bin/env Rscript
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
src <- if (is.na(arg)) file.path("skills", "geo-biodata-workflow", "scripts", "drivers", "run_bulk_counts.R") else sub("^--file=", "", arg)
root <- normalizePath(file.path(dirname(normalizePath(src, mustWork = TRUE)), "..", "..", "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "core", "R", "drivers", "run_bulk_counts.R"))
