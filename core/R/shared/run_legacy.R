run_legacy_script <- function(relative_script) {
  args <- commandArgs(trailingOnly = TRUE)
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  this_file <- if (length(file_arg) > 0L) {
    normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
  } else {
    normalizePath(sys.frames()[[1L]]$ofile %||% "", winslash = "/", mustWork = FALSE)
  }

  find_repo_root <- function(start_dir) {
    current <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
    for (i in seq_len(12L)) {
      if (file.exists(file.path(current, "dependency_profiles.yaml")) &&
          dir.exists(file.path(current, "skills", "geo-biodata-workflow", "scripts"))) {
        return(current)
      }
      parent <- dirname(current)
      if (identical(parent, current)) break
      current <- parent
    }
    stop("Could not find GEO_biodata repository root.", call. = FALSE)
  }

  repo_root <- find_repo_root(dirname(this_file))
  target <- file.path(repo_root, "skills", "geo-biodata-workflow", "scripts", relative_script)
  if (!file.exists(target)) {
    stop("Legacy implementation script does not exist: ", target, call. = FALSE)
  }
  message("NOTE: core/R wrapper forwarding to compatibility implementation: ", relative_script)
  status <- system2("Rscript", c(shQuote(target), shQuote(args)))
  quit(status = as.integer(status %||% 0L))
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
