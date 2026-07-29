#!/usr/bin/env Rscript
arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)[1]
src <- if (is.na(arg)) file.path("skills", "geo-biodata-workflow", "scripts", "generate_download_plan.R") else sub("^--file=", "", arg)
root <- normalizePath(file.path(dirname(normalizePath(src, mustWork = TRUE)), "..", "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(root, "core", "R", "generate_download_plan.R"))
