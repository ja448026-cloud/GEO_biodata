#!/usr/bin/env Rscript
#
# gene_id_mapping.R
# Gene-ID mapping contract for GEO_biodata enrichment.
# Source this file; do not run standalone.
#
# Supports:
#   - GEO platform (GPL) annotation files
#   - org.*.eg.db Bioconductor annotation packages
#   - Custom two-column mapping tables
#
# The mapping contract ensures:
#   1. Probe/gene ID → gene symbol mapping with dedup strategy documented
#   2. Universe (background) definition from the platform or reference
#   3. Coverage and dedup statistics recorded in an audit table

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

clean_gene_symbols <- function(x, min_symbol_nchar = 2L) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  split_symbols <- unlist(strsplit(x, "\\s*(///|;|,)\\s*"), use.names = FALSE)
  split_symbols <- trimws(split_symbols)
  split_symbols <- split_symbols[nzchar(split_symbols) & nchar(split_symbols) >= min_symbol_nchar]
  unique(split_symbols)
}

# ── GPL annotation reader ─────────────────────────────────────────────────────

read_gpl_annotation <- function(gpl_path) {
  is_gz <- grepl("\\.gz$", gpl_path, ignore.case = TRUE)

  # On Windows ucrt R builds, gzfile() and file() connections with readLines(n=1)
  # can segfault on large files. Use shell decompression when available.
  read_path <- gpl_path
  tmp_file <- NULL
  if (is_gz) {
    gzip_bin <- Sys.which("gzip")
    if (nzchar(gzip_bin)) {
      tmp_file <- tempfile(fileext = ".annot")
      ret <- system2(gzip_bin, c("-dc", shQuote(gpl_path)), stdout = tmp_file)
      if (ret == 0L && file.exists(tmp_file) && file.info(tmp_file)$size > 0) {
        read_path <- tmp_file
        is_gz <- FALSE
      }
    }
    if (is_gz) stop("Cannot decompress .gz file and gzfile() is unreliable on this R build. Install gzip.")
  }
  on.exit(if (!is.null(tmp_file) && file.exists(tmp_file)) try(unlink(tmp_file), silent = TRUE))

  # Read all lines into memory (avoids connection-based read that segfaults on ucrt)
  all_lines <- readLines(read_path, warn = FALSE)
  table_begin <- which(all_lines == "!platform_table_begin")
  if (length(table_begin) == 0L) stop("GPL annotation missing !platform_table_begin marker.")
  table_begin <- table_begin[[1L]]
  table_end <- which(seq_along(all_lines) > table_begin & all_lines == "!platform_table_end")
  if (length(table_end) == 0L) stop("GPL annotation missing !platform_table_end marker.")
  table_end <- table_end[[1L]]
  if (table_end <= table_begin + 1L) stop("GPL annotation has no table header or data rows.")

  header_lines <- if (table_begin > 1L) all_lines[seq_len(table_begin - 1L)] else character()
  data_lines <- all_lines[(table_begin + 1L):(table_end - 1L)]

  if (length(data_lines) < 2L) stop("GPL annotation has no data rows.")

  col_header <- strsplit(data_lines[1L], "\t")[[1L]]
  data_lines <- data_lines[-1L]

  # Parse via tempfile to avoid textConnection memory issues
  data_tmp <- tempfile(fileext = ".tsv")
  on.exit(try(unlink(data_tmp), silent = TRUE), add = TRUE)
  writeLines(data_lines, data_tmp)
  annot <- utils::read.delim(data_tmp, header = FALSE, stringsAsFactors = FALSE,
                              col.names = col_header, check.names = FALSE,
                              quote = "", comment.char = "", na.strings = "")

  # Extract platform metadata
  meta <- list()
  for (line in header_lines) {
    if (grepl("^!.*=.*", line)) {
      parts <- strsplit(sub("^!", "", line), " = ")[[1L]]
      key <- trimws(parts[1L])
      val <- if (length(parts) > 1L) trimws(paste(parts[-1L], collapse = " = ")) else ""
      meta[[key]] <- val
    }
  }

  id_col <- "ID"
  symbol_col <- "Gene symbol"
  entrez_col <- "Gene ID"

  if (!id_col %in% names(annot)) {
    stop("GPL annotation missing ID column.")
  }

  list(
    annotation = annot,
    meta = meta,
    platform_id = meta[["Annotation_platform"]] %||% "unknown",
    id_col = id_col,
    symbol_col = if (symbol_col %in% names(annot)) symbol_col else NULL,
    entrez_col = if (entrez_col %in% names(annot)) entrez_col else NULL
  )
}

# ── Probe-to-symbol mapping ───────────────────────────────────────────────────

build_probe_symbol_map <- function(gpl_annotation, min_symbol_nchar = 2L) {
  annot <- gpl_annotation$annotation
  id_col <- gpl_annotation$id_col
  symbol_col <- gpl_annotation$symbol_col

  if (is.null(symbol_col)) stop("GPL annotation has no Gene symbol column.")

  split_rows <- lapply(seq_len(nrow(annot)), function(i) {
    symbols <- clean_gene_symbols(annot[[symbol_col]][[i]], min_symbol_nchar = min_symbol_nchar)
    if (!length(symbols)) return(NULL)
    data.frame(
      probe_id = as.character(annot[[id_col]][[i]]),
      gene_symbol = symbols,
      stringsAsFactors = FALSE
    )
  })
  df <- do.call(rbind, split_rows)
  if (is.null(df) || nrow(df) == 0L) stop("GPL annotation has no usable gene symbols.")
  df <- df[!duplicated(df), ]

  n_probes_mapped <- length(unique(df$probe_id))
  n_single <- sum(!duplicated(df$probe_id, fromLast = FALSE) & !duplicated(df$probe_id, fromLast = TRUE))
  n_multi <- n_probes_mapped - n_single
  n_probes_total <- length(unique(as.character(annot[[id_col]])))

  map <- df
  attr(map, "stats") <- list(
    total_probes = n_probes_total,
    mapped_probes = n_probes_mapped,
    unique_genes = length(unique(map$gene_symbol)),
    single_map_probes = n_single,
    multi_map_probes = n_multi,
    coverage = round(n_probes_mapped / n_probes_total * 100, 1),
    platform_id = gpl_annotation$platform_id %||% "unknown"
  )

  map
}

# ── Universe from GPL ─────────────────────────────────────────────────────────

build_universe_from_gpl <- function(gpl_annotation) {
  annot <- gpl_annotation$annotation
  symbol_col <- gpl_annotation$symbol_col

  if (is.null(symbol_col)) {
    # Fallback: use all probe IDs as universe
    return(unique(as.character(annot[[gpl_annotation$id_col]])))
  }

  clean_gene_symbols(annot[[symbol_col]])
}

# ── Dedup strategy for multi-mapping probes ───────────────────────────────────

map_probes_to_symbols <- function(de_result, probe_symbol_map, method = c("max_abs_logfc", "min_pvalue", "mean")) {
  method <- match.arg(method)
  map <- probe_symbol_map

  if (!"feature_id" %in% names(de_result)) {
    de_result$feature_id <- rownames(de_result)
  }

  merged <- merge(de_result, map, by.x = "feature_id", by.y = "probe_id", all.x = TRUE)
  merged <- merged[!is.na(merged$gene_symbol) & nzchar(merged$gene_symbol), ]

  logfc_col <- if ("logFC" %in% names(merged)) "logFC" else names(merged)[grep("logFC|log2FoldChange", names(merged))][1L]
  pval_col <- if ("P.Value" %in% names(merged)) "P.Value" else names(merged)[grep("P\\.Value|pvalue", names(merged))][1L]

  dup_genes <- merged$gene_symbol[duplicated(merged$gene_symbol)]
  if (length(dup_genes) > 0L) {
    dedup_groups <- split(merged, merged$gene_symbol)
    dedup_records <- list()
    for (gene in names(dedup_groups)) {
      group <- dedup_groups[[gene]]
      if (nrow(group) == 1L) {
        dedup_records[[gene]] <- group
      } else {
        if (method == "max_abs_logfc" && !is.na(logfc_col)) {
          best <- group[which.max(abs(group[[logfc_col]])), ][1L, ]
        } else if (method == "min_pvalue" && !is.na(pval_col)) {
          best <- group[which.min(group[[pval_col]]), ][1L, ]
        } else {
          best <- group[1L, ]
        }
        dedup_records[[gene]] <- best
      }
    }
    merged <- do.call(rbind, dedup_records)
    rownames(merged) <- NULL
  }

  merged
}

# ── Mapping audit ─────────────────────────────────────────────────────────────

write_mapping_audit <- function(gpl_stats, probe_symbol_map, de_result, out_path) {
  if (!"feature_id" %in% names(de_result)) {
    de_result$feature_id <- rownames(de_result)
  }
  n_de_probes <- nrow(de_result)
  mapped_de <- sum(de_result$feature_id %in% probe_symbol_map$probe_id)
  n_unique_after_dedup <- length(unique(de_result$feature_id))

  audit <- data.frame(
    field = c(
      "platform_id",
      "total_probes_on_platform",
      "probes_mapped_to_symbol",
      "platform_mapping_coverage_pct",
      "unique_gene_symbols",
      "de_probes_input",
      "de_probes_mapped",
      "de_mapping_coverage_pct",
      "unique_de_features_input",
      "dedup_method",
      "timestamp_utc"
    ),
    value = c(
      gpl_stats$platform_id %||% "unknown",
      as.character(gpl_stats$total_probes %||% NA_integer_),
      as.character(gpl_stats$mapped_probes %||% NA_integer_),
      as.character(gpl_stats$coverage %||% NA_real_),
      as.character(gpl_stats$unique_genes %||% NA_integer_),
      as.character(n_de_probes),
      as.character(mapped_de),
      as.character(round(mapped_de / max(1, n_de_probes) * 100, 1)),
      as.character(n_unique_after_dedup),
      "max_abs_logfc",
      format(Sys.time(), tz = "UTC", usetz = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  utils::write.table(audit, out_path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  audit
}
