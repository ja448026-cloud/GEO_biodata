#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: validate_manifest.R /path/to/run_manifest.yaml", call. = FALSE)
}

manifest_path <- args[[1L]]
if (!file.exists(manifest_path)) stop("Manifest file does not exist: ", manifest_path, call. = FALSE)
if (!requireNamespace("yaml", quietly = TRUE)) stop("Missing required package: yaml", call. = FALSE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

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
manifest_dir <- dirname(normalizePath(manifest_path, winslash = "/", mustWork = TRUE))
schema_path <- find_repo_file(c(getwd(), manifest_dir, script_dir), file.path("schemas", "run_manifest.schema.yaml"))
if (is.na(schema_path)) stop("Could not find schemas/run_manifest.schema.yaml.", call. = FALSE)

manifest <- yaml::read_yaml(manifest_path)
schema <- yaml::read_yaml(schema_path)

errors <- character()
warnings <- character()
add_error <- function(...) errors <<- c(errors, paste0(...))
add_warning <- function(...) warnings <<- c(warnings, paste0(...))

is_true <- function(x) isTRUE(x) || identical(tolower(as.character(x)), "true")
has_field <- function(x, field) !is.null(x[[field]]) && length(x[[field]]) > 0L && !identical(x[[field]], "")
resolve_manifest_path <- function(path) {
  if (!nzchar(path %||% "")) return("")
  if (grepl("^[A-Za-z]:[\\\\/]|^/", path)) return(path)
  file.path(manifest_dir, path)
}

for (field in schema$required_top_level) {
  if (!has_field(manifest, field)) add_error("Missing top-level field: ", field)
}

accession <- toupper(manifest$accession %||% "")
if (!grepl("^GSE[0-9]+$", accession)) add_error("accession must be a GEO Series accession like GSE123456.")

route <- manifest$route %||% ""
route_rules <- schema$route_ontology[[route]]
allowed_routes <- names(schema$route_ontology)
if (is.null(route_rules)) {
  add_error("route must be one of: ", paste(allowed_routes, collapse = ", "))
  route_rules <- list(
    input_types = character(),
    requires_input_file = TRUE,
    requires_sample_mapping = TRUE,
    requires_design = TRUE,
    requires_contrast = TRUE
  )
}

for (field in schema$required_input_fields) {
  if (!has_field(manifest$input, field)) add_error("Missing input field: input.", field)
}
input_type <- manifest$input$input_type %||% ""
if (!input_type %in% (route_rules$input_types %||% character())) {
  add_error(
    "input.input_type '", input_type, "' is not compatible with route '", route,
    "'. Allowed: ", paste(route_rules$input_types %||% character(), collapse = ", ")
  )
}

location_type <- manifest$input$location_type %||% "local"
if (!location_type %in% schema$allowed_location_types) {
  add_error("input.location_type must be one of: ", paste(schema$allowed_location_types, collapse = ", "))
}

input_file <- manifest$input$file %||% ""
input_path <- resolve_manifest_path(input_file)
if (isTRUE(route_rules$requires_input_file)) {
  if (!has_field(manifest$input, "file")) {
    add_error("Missing input field: input.file")
  } else if (identical(location_type, "local") && !file.exists(input_path)) {
    add_error("Local input.file does not exist: ", input_file)
  } else if (identical(location_type, "remote") && !grepl("^(https?|ftp|s3|gs)://", input_file, ignore.case = TRUE)) {
    add_error("Remote input.file must be a URI: ", input_file)
  }
}

review_mode <- manifest$review$mode %||% "manual"
if (!review_mode %in% names(schema$review_modes)) {
  add_error("review.mode must be one of: ", paste(names(schema$review_modes), collapse = ", "))
}
if (identical(review_mode, "manual")) {
  for (field in schema$review_modes$manual$required_flags) {
    if (!has_field(manifest$review, field) || !is_true(manifest$review[[field]])) {
      add_error("Manual review gate is not confirmed: review.", field)
    }
  }
} else if (identical(review_mode, "agent_adjudicated")) {
  for (field in schema$review_modes$agent_adjudicated$required_fields) {
    if (!has_field(manifest$review, field)) add_error("Missing agent adjudication field: review.", field)
  }
  if (has_field(manifest$review, "unresolved_conflicts")) {
    unresolved <- suppressWarnings(as.integer(manifest$review$unresolved_conflicts))
    if (!is.finite(unresolved) || unresolved != 0L) add_error("review.unresolved_conflicts must be 0 for automatic analysis.")
  }
  if (has_field(manifest$review, "evidence_threshold_passed") && !is_true(manifest$review$evidence_threshold_passed)) {
    add_error("review.evidence_threshold_passed must be true for automatic analysis.")
  }
  decision_record <- manifest$review$decision_record %||% ""
  if (nzchar(decision_record) && !file.exists(resolve_manifest_path(decision_record))) {
    add_error("review.decision_record does not exist: ", decision_record)
  }
}

autonomy_mode <- manifest$autonomy$mode %||% "balanced"
if (!autonomy_mode %in% schema$autonomy_modes) {
  add_error("autonomy.mode must be one of: ", paste(schema$autonomy_modes, collapse = ", "))
}

sample_map <- NULL
sample_path <- ""
sample_id_col <- manifest$sample_mapping$sample_id_column %||% ""
group_col <- manifest$sample_mapping$group_column %||% ""
if (isTRUE(route_rules$requires_sample_mapping)) {
  for (field in schema$required_sample_mapping_fields) {
    if (!has_field(manifest$sample_mapping, field)) add_error("Missing sample_mapping field: sample_mapping.", field)
  }
  sample_file <- manifest$sample_mapping$file %||% ""
  sample_path <- resolve_manifest_path(sample_file)
  if (nzchar(sample_file) && !file.exists(sample_path)) {
    add_error("sample_mapping.file does not exist: ", sample_file)
  }
}

if (nzchar(sample_path) && file.exists(sample_path)) {
  sample_map <- utils::read.delim(sample_path, stringsAsFactors = FALSE, check.names = FALSE)
  for (col in c(sample_id_col, group_col)) {
    if (!col %in% names(sample_map)) add_error("Sample mapping column not found: ", col)
  }
  if (sample_id_col %in% names(sample_map)) {
    ids <- sample_map[[sample_id_col]]
    if (any(is.na(ids) | !nzchar(as.character(ids)))) add_error("Sample IDs contain missing or empty values.")
    if (any(duplicated(ids))) add_error("Sample IDs are duplicated.")
  }
}

if (isTRUE(route_rules$requires_design)) {
  for (field in schema$required_design_fields) {
    if (!has_field(manifest$design, field)) add_error("Missing design field: design.", field)
  }
}
if (isTRUE(route_rules$requires_contrast)) {
  for (field in schema$required_contrast_fields) {
    if (!has_field(manifest$design$contrast, field)) add_error("Missing contrast field: design.contrast.", field)
  }
}

design_formula <- manifest$design$formula %||% ""
contrast <- manifest$design$contrast
if (nzchar(design_formula) && !is.null(sample_map)) {
  formula_vars <- all.vars(stats::as.formula(design_formula))
  missing_formula_vars <- setdiff(formula_vars, names(sample_map))
  if (length(missing_formula_vars) > 0L) {
    add_error("Design formula fields missing from sample mapping: ", paste(missing_formula_vars, collapse = ", "))
  }
  if (length(missing_formula_vars) == 0L && length(formula_vars) > 0L) {
    design_matrix <- tryCatch(
      stats::model.matrix(stats::as.formula(design_formula), data = sample_map),
      error = function(e) {
        add_error("Design matrix could not be built: ", conditionMessage(e))
        NULL
      }
    )
    if (!is.null(design_matrix) && qr(design_matrix)$rank < ncol(design_matrix)) add_error("Design matrix is not full rank.")
  }
  if (isTRUE(route_rules$requires_contrast) && !is.null(contrast)) {
    contrast_factor <- contrast$factor %||% group_col
    if (!contrast_factor %in% names(sample_map)) {
      add_error("Contrast factor missing from sample mapping: ", contrast_factor)
    } else {
      levels_present <- unique(as.character(sample_map[[contrast_factor]]))
      if (!contrast$numerator %in% levels_present) add_error("Contrast numerator not present: ", contrast$numerator)
      if (!contrast$denominator %in% levels_present) add_error("Contrast denominator not present: ", contrast$denominator)
      if (!contrast_factor %in% formula_vars) add_error("Contrast factor is not included in design formula: ", contrast_factor)
      unit_col <- manifest$sample_mapping$biological_unit_column %||% sample_id_col
      if (unit_col %in% names(sample_map)) {
        units_by_group <- aggregate(sample_map[[unit_col]], list(group = sample_map[[contrast_factor]]), function(x) length(unique(x)))
        names(units_by_group)[2L] <- "n_units"
        low <- units_by_group$group %in% c(contrast$numerator, contrast$denominator) & units_by_group$n_units < 2L
        if (any(low)) add_error("Each contrast level needs at least two biological units for basic DE: ", paste(units_by_group$group[low], collapse = ", "))
      }
    }
  }
}

read_matrix_input <- function(path) {
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
}

if (identical(location_type, "local") && nzchar(input_path) && file.exists(input_path) &&
    route %in% c("bulk_raw_counts", "bulk_normalized", "microarray_series_matrix")) {
  mat_df <- tryCatch(read_matrix_input(input_path), error = function(e) {
    add_error("Could not read input matrix as tabular text: ", conditionMessage(e))
    NULL
  })
  if (!is.null(mat_df)) {
    if (ncol(mat_df) < 2L) {
      add_error("Input matrix must have a feature column and at least one sample column.")
    } else {
      sample_cols <- colnames(mat_df)[-1L]
      value_mat <- as.matrix(mat_df[, -1L, drop = FALSE])
      suppressWarnings(storage.mode(value_mat) <- "numeric")
      if (any(!is.finite(value_mat))) add_error("Input matrix contains NA, NaN, Inf, or non-numeric values.")
      if (route == "bulk_raw_counts") {
        if (any(value_mat < 0, na.rm = TRUE)) add_error("Raw counts contain negative values.")
        if (any(abs(value_mat - round(value_mat)) > 1e-8, na.rm = TRUE)) add_error("Raw counts must be integer-like.")
      }
      if (!is.null(sample_map) && sample_id_col %in% names(sample_map)) {
        mapped <- as.character(sample_map[[sample_id_col]])
        missing_in_mapping <- setdiff(sample_cols, mapped)
        extra_in_mapping <- setdiff(mapped, sample_cols)
        if (length(missing_in_mapping) > 0L) add_error("Matrix columns missing from sample mapping: ", paste(missing_in_mapping, collapse = ", "))
        if (length(extra_in_mapping) > 0L) add_error("Sample mapping IDs missing from matrix columns: ", paste(extra_in_mapping, collapse = ", "))
      }
    }
  }
}

# ── Conditional requirements (route + analysis.intent) ───────────────────────

if (!is.null(route_rules$conditional_requirements)) {
  analysis_intent <- manifest$analysis$intent %||% "differential_expression"
  for (cr in route_rules$conditional_requirements) {
    cond <- cr$condition
    # Check if this condition's analysis.intent matches
    if (!is.null(cond$analysis.intent) && !identical(cond$analysis.intent, analysis_intent)) next

    # Apply overrides
    if (isTRUE(cr$requires_design) && !isTRUE(route_rules$requires_design)) {
      route_rules$requires_design <- TRUE
      if (!has_field(manifest$design, "formula")) {
        add_error("route + analysis.intent=", analysis_intent, " requires design.formula.")
      }
      if (!has_field(manifest$design, "contrast")) {
        add_error("route + analysis.intent=", analysis_intent, " requires design.contrast.")
      }
    }
    if (isTRUE(cr$requires_contrast) && !isTRUE(route_rules$requires_contrast)) {
      route_rules$requires_contrast <- TRUE
    }

    # Validate required scale fields
    if (!is.null(cr$required_scale_fields) && length(cr$required_scale_fields) > 0L) {
      scale_config <- manifest$input$scale %||% list()
      for (sf in cr$required_scale_fields) {
        if (!has_field(scale_config, sf)) {
          add_error("route + analysis.intent=", analysis_intent, " requires input.scale.", sf)
        }
      }
      # Scale value rules
      if (!is.null(cr$scale_value_rules)) {
        for (svr in cr$scale_value_rules) {
          sv_field <- svr$field
          applies <- svr$applies_when
          if (!is.null(applies$input_type)) {
            if (!input_type %in% unlist(applies$input_type)) next
          }
          val <- scale_config[[sv_field]]
          if (!is.null(svr$must_equal) && !identical(val, svr$must_equal)) {
            add_error("input.scale.", sv_field, " must be ", svr$must_equal, " for ", input_type, ", got: ", val %||% "NULL")
          }
          if (!is.null(svr$must_start_with) && !grepl(paste0("^", svr$must_start_with), val %||% "")) {
            add_error("input.scale.", sv_field, " must start with '", svr$must_start_with, "' for ", input_type, ", got: ", val %||% "NULL")
          }
        }
      }
    }
  }
}

# ── Scale contract validation ────────────────────────────────────────────────
if (route %in% c("bulk_normalized", "microarray_series_matrix")) {
  scale_config <- manifest$input$scale %||% list()
  analysis_intent <- manifest$analysis$intent %||% "differential_expression"

  if (!analysis_intent %in% c("eda_only", "differential_expression")) {
    add_warning("analysis.intent should be 'eda_only' or 'differential_expression', got: ", analysis_intent)
  }

  if (route == "bulk_normalized") {
    input_type <- manifest$input$input_type %||% ""
    if (identical(input_type, "normalized_expression") &&
        !isTRUE(scale_config$evidence_source) && !nzchar(scale_config$evidence_source %||% "")) {
      add_warning("normalized_expression without scale evidence: only EDA is safe. Set analysis.intent=eda_only or provide input.scale.evidence_source.")
    }
    if (identical(input_type, "tpm") &&
        (!isTRUE(scale_config$transformed) || !nzchar(scale_config$transform %||% ""))) {
      add_warning("TPM without confirmed log-transform: limma DE may be inappropriate. Set input.scale.transformed=true and input.scale.transform.")
    }
  }
}

# ── Agent adjudication gate ──────────────────────────────────────────────────
decision_record <- manifest$review$decision_record %||% ""
decision_path <- if (nzchar(decision_record)) resolve_manifest_path(decision_record) else ""
# Only run agent gate for agent_adjudicated review mode
if (identical(manifest$review$mode %||% "", "agent_adjudicated")) {
  if (!nzchar(decision_path) || !file.exists(decision_path)) {
    add_error("review.decision_record '", decision_record, "' does not exist or is empty.")
  } else {
    decision <- tryCatch(
      utils::read.delim(decision_path, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL
    )
    if (is.null(decision) || nrow(decision) == 0L) {
      add_error("analysis_decisions.tsv is empty or unreadable.")
    } else {
      # Required columns (with aliases)
      required_map <- list(
        selected_value = c("selected_value"),
        confidence = c("confidence"),
        source_tier = c("source_tier"),
        decision_rule = c("decision_rule"),
        requires_user_input = c("requires_user_input"),
        conflict_status = c("conflict_status"),
        evidence_sources = c("evidence_sources"),
        retrieved_at = c("retrieved_at", "retrieved_at_utc"),
        decided_at = c("decided_at")
      )
      for (col_name in names(required_map)) {
        found <- any(required_map[[col_name]] %in% names(decision))
        if (!found) {
          add_error("analysis_decisions missing required column: ", col_name,
            " (accepted: ", paste(required_map[[col_name]], collapse = ", "), ")")
        }
      }

      # Multi-row: select the correct row by decision_type or decision_id
      if (nrow(decision) > 1L) {
        target_id <- manifest$review$decision_id %||% ""
        if (nzchar(target_id) && "decision_id" %in% names(decision)) {
          decision <- decision[decision$decision_id == target_id, , drop = FALSE]
          if (nrow(decision) == 0L) {
            add_error("analysis_decisions has no row with decision_id '", target_id, "'.")
          }
        } else if ("decision_type" %in% names(decision)) {
          sel <- decision[grepl("route", tolower(decision$decision_type)), , drop = FALSE]
          if (nrow(sel) == 0L) sel <- decision  # fallback: use first row
          decision <- sel[1L, , drop = FALSE]
        }
      }
      if (nrow(decision) == 0L) {
        add_error("analysis_decisions could not resolve a valid decision row.")
      }

      # selected_value must match manifest route
      if ("selected_value" %in% names(decision)) {
        if (decision$selected_value[[1L]] != route) {
          add_error("analysis_decisions selected_value '", decision$selected_value[[1L]],
            "' does not match manifest route '", route, "'.")
        }
      }

      # requires_user_input must be FALSE
      if ("requires_user_input" %in% names(decision)) {
        if (isTRUE(decision$requires_user_input[[1L]]) ||
            identical(tolower(as.character(decision$requires_user_input[[1L]])), "true")) {
          add_error("analysis_decisions requires_user_input is TRUE.")
        }
      }

      # Confidence threshold (hard error)
      autonomy <- manifest$autonomy$mode %||% "balanced"
      min_confidence <- switch(autonomy, conservative = 0.90, balanced = 0.75, autonomous = 0.60)
      if ("confidence" %in% names(decision)) {
        dec_conf <- as.numeric(decision$confidence[[1L]])
        if (is.na(dec_conf) || dec_conf < min_confidence) {
          add_error(sprintf("Decision confidence %.2f below autonomy '%s' threshold %.2f.",
            dec_conf %||% 0, autonomy, min_confidence))
        }
      }

      # Conflict status must be structured (none/resolved/unresolved)
      if ("conflict_status" %in% names(decision)) {
        cs <- tolower(as.character(decision$conflict_status[[1L]] %||% ""))
        if (!cs %in% c("none", "resolved")) {
          add_error("analysis_decisions conflict_status must be 'none' or 'resolved', got: ", cs)
        }
      } else if ("conflicting_evidence" %in% names(decision)) {
        # Legacy: structured conflict_status column is required
        add_error("analysis_decisions missing required column: conflict_status. Add conflict_status: none|resolved|unresolved.")
      }

      # Source tier must be in registry (hard error)
      if ("source_tier" %in% names(decision)) {
        source_registry <- find_repo_file(c(getwd(), manifest_dir, script_dir), file.path("knowledge", "source_registry.yaml"))
        valid_tiers <- if (!is.na(source_registry)) {
          src <- yaml::read_yaml(source_registry)
          names(src$source_tiers) %||% character()
        } else {
          character()
        }
        if (length(valid_tiers) > 0L && !decision$source_tier[[1L]] %in% valid_tiers) {
          add_error("analysis_decisions source_tier '", decision$source_tier[[1L]], "' not in registered tiers.")
        }
      }

      # Decision rule must exist in knowledge/decision_rules/
      if ("decision_rule" %in% names(decision)) {
        dr <- decision$decision_rule[[1L]] %||% ""
        if (!nzchar(dr) || grepl("^rules_only", dr)) {
          add_error("decision_rule is empty or 'rules_only'.")
        }
        # Verify rule exists in decision_rules directory
        rules_dir <- find_repo_file(c(getwd(), manifest_dir, script_dir), file.path("knowledge", "decision_rules"))
        if (!is.na(rules_dir)) {
          rule_files <- list.files(rules_dir, pattern = "\\.yaml$", full.names = TRUE, recursive = TRUE)
          rule_ids <- character()
          for (rf in rule_files) {
            ry <- tryCatch(yaml::read_yaml(rf), error = function(e) NULL)
            if (!is.null(ry$rule_id)) rule_ids <- c(rule_ids, ry$rule_id)
          }
          if (length(rule_ids) > 0L && !dr %in% rule_ids) {
            add_error("decision_rule '", dr, "' not found in knowledge/decision_rules/.")
          }
        }
      }

      # Timestamps: must be parseable and decided_at >= retrieved_at
      parse_ts <- function(val) {
        if (is.null(val) || is.na(val) || !nzchar(as.character(val))) return(NA_real_)
        as.numeric(tryCatch(as.POSIXct(val, tz = "UTC"), error = function(e) NA_real_))
      }
      ret_ts <- parse_ts(decision$retrieved_at[[1L]] %||% decision$retrieved_at_utc[[1L]])
      dec_ts <- parse_ts(decision$decided_at[[1L]])
      if (is.na(ret_ts) && is.na(dec_ts)) {
        add_error("analysis_decisions has unparseable or missing timestamps.")
      } else if (!is.na(ret_ts) && !is.na(dec_ts) && dec_ts < ret_ts) {
        add_error("analysis_decisions decided_at is before retrieved_at.")
      }
    }
  }
}

out_dir <- dirname(manifest_path)
status <- if (length(errors) == 0L) "MANIFEST_VALIDATED" else "MANIFEST_INVALID"
validation <- data.frame(
  status = c(rep("ERROR", length(errors)), rep("WARNING", length(warnings))),
  message = c(errors, warnings),
  stringsAsFactors = FALSE
)
if (nrow(validation) == 0L) validation <- data.frame(status = "OK", message = "Manifest validation passed.", stringsAsFactors = FALSE)
utils::write.table(validation, file.path(out_dir, "manifest_validation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

workflow_status <- data.frame(
  state = status,
  updated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  note = if (length(errors) == 0L) "Manifest structure, route requirements, input checks, and review gates passed." else paste(errors, collapse = " | "),
  stringsAsFactors = FALSE
)
utils::write.table(workflow_status, file.path(out_dir, "workflow_status.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "")

cat(status, "\n", sep = "")
cat(normalizePath(file.path(out_dir, "manifest_validation.tsv"), winslash = "/", mustWork = TRUE), "\n")
if (length(errors) > 0L) quit(status = 1L)
