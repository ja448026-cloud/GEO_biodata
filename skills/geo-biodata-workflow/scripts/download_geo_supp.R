#!/usr/bin/env Rscript

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg) > 0L) {
  normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)
} else {
  normalizePath(file.path("skills", "geo-biodata-workflow", "scripts", "download_geo_supp.R"), mustWork = TRUE)
}

repo_root <- normalizePath(file.path(dirname(script_path), "..", "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(repo_root, "core", "R", "download_reviewed_files.R"))
