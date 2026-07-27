#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
usage <- paste(
  "Usage:",
  "  bootstrap_environment.R --profile core --check",
  "  bootstrap_environment.R --profile bulk --plan",
  "  bootstrap_environment.R --profile scrna --install",
  sep = "\n"
)
if (length(args) < 2L) stop(usage, call. = FALSE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
get_opt <- function(flag, default = NA_character_) {
  hit <- which(args == flag)
  if (length(hit) == 0L || hit[[1L]] == length(args)) return(default)
  args[[hit[[1L]] + 1L]]
}
mode <- intersect(args, c("--check", "--plan", "--install"))
if (length(mode) != 1L) stop("Specify exactly one of --check, --plan, or --install.", call. = FALSE)
mode <- mode[[1L]]
profile_name <- get_opt("--profile", "core")

find_repo_file <- function(starts, relative_path) {
  for (start in starts) {
    current <- normalizePath(start, winslash = "/", mustWork = FALSE)
    for (i in seq_len(8L)) {
      candidate <- file.path(current, relative_path)
      if (file.exists(candidate)) return(candidate)
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
  }
  NA_character_
}
cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_dir <- if (length(file_arg) > 0L) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), mustWork = TRUE)) else getwd()
repo_profile <- find_repo_file(c(getwd(), script_dir), "dependency_profiles.yaml")
if (is.na(repo_profile)) stop("Could not find dependency_profiles.yaml.", call. = FALSE)
if (!requireNamespace("yaml", quietly = TRUE)) stop("Missing required package: yaml", call. = FALSE)

profiles <- yaml::read_yaml(repo_profile)$profiles
if (is.null(profiles[[profile_name]])) {
  stop("Unknown profile: ", profile_name, ". Available: ", paste(names(profiles), collapse = ", "), call. = FALSE)
}

collect_profile <- function(name, seen = character()) {
  if (name %in% seen) stop("Dependency profile cycle detected: ", paste(c(seen, name), collapse = " -> "), call. = FALSE)
  profile <- profiles[[name]]
  inherited <- profile$extends %||% character()
  collected <- lapply(inherited, collect_profile, seen = c(seen, name))
  cran <- unique(unlist(c(lapply(collected, `[[`, "cran"), profile$packages$cran %||% character()), use.names = FALSE))
  bioc <- unique(unlist(c(lapply(collected, `[[`, "bioc"), profile$packages$bioc %||% character()), use.names = FALSE))
  list(cran = cran, bioc = bioc)
}
deps <- collect_profile(profile_name)
dep_table <- rbind(
  data.frame(profile = profile_name, source = "cran", package = deps$cran, stringsAsFactors = FALSE),
  data.frame(profile = profile_name, source = "bioc", package = deps$bioc, stringsAsFactors = FALSE)
)
if (nrow(dep_table) == 0L) dep_table <- data.frame(profile = profile_name, source = character(), package = character())
dep_table$installed <- vapply(dep_table$package, requireNamespace, logical(1), quietly = TRUE)
dep_table$version <- vapply(dep_table$package, function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return("")
  as.character(utils::packageVersion(pkg))
}, character(1))

out_path <- file.path(getwd(), paste0("environment_profile_", profile_name, ".tsv"))
utils::write.table(dep_table, out_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

if (mode == "--plan") {
  cat("DEPENDENCY_PLAN\n")
  cat("Profile: ", profile_name, "\n", sep = "")
  if (length(deps$cran) > 0L) cat("CRAN: ", paste(deps$cran, collapse = ", "), "\n", sep = "")
  if (length(deps$bioc) > 0L) cat("Bioconductor: ", paste(deps$bioc, collapse = ", "), "\n", sep = "")
  cat(normalizePath(out_path, winslash = "/", mustWork = TRUE), "\n", sep = "")
  quit(status = 0L)
}

missing <- dep_table$package[!dep_table$installed]
if (mode == "--check") {
  cat(if (length(missing) == 0L) "DEPENDENCIES_READY" else "DEPENDENCIES_MISSING", "\n", sep = "")
  if (length(missing) > 0L) cat("Missing: ", paste(missing, collapse = ", "), "\n", sep = "")
  cat(normalizePath(out_path, winslash = "/", mustWork = TRUE), "\n", sep = "")
  quit(status = if (length(missing) == 0L) 0L else 1L)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
if (length(deps$cran) > 0L) {
  install.packages(setdiff(deps$cran, rownames(utils::installed.packages())), repos = "https://cloud.r-project.org")
}
if (length(deps$bioc) > 0L) {
  BiocManager::install(deps$bioc, ask = FALSE, update = FALSE)
}
cat("DEPENDENCY_INSTALL_ATTEMPTED\n")
cat(normalizePath(out_path, winslash = "/", mustWork = TRUE), "\n", sep = "")
