---
name: geo-biodata
description: Unified GEO_biodata entry point. Use when a user gives a GEO accession or project-level GEO_biodata request and needs routing across intake, reviewed downloads, bulk analysis, preranked GSEA, or read-only scRNA intake.
---

# GEO_biodata

Use this as the default entry point. This project is a skill library plus small R helpers, not a full software framework. Route to the smallest module that helps the user finish the GEO task.

## Module Selection

| User intent | Module | Executor |
|---|---|---|
| GEO metadata, DOI/PMID/PMCID, open supplements, reviewed download, file audit, manifest draft | `geo-biodata-intake` | `core/R/*intake*`, `core/R/download_reviewed_files.R` |
| Bulk raw-count, normalized, or microarray DE/QC | `geo-biodata-bulk` | `core/R/drivers/run_bulk_*.R`, `core/R/drivers/run_microarray.R` |
| Preranked GSEA from an existing DE/rank table and GMT | `geo-biodata-enrichment` | `core/R/enrichment/run_preranked_gsea.R` |
| Seurat/H5AD/RDS object inventory without reclustering | `geo-biodata-scrna-intake` | `core/R/scrna/inspect_object.R` |
| Single-purpose plot planning, source linkage, figure generation guidance, and QA | `geo-biodata-figure` | independent playbook; agent-authored R/Python plotting from source-linked outputs |

## Defaults

Start with intake unless the user already provides a validated manifest or a DE/rank table.

Do not run statistical analysis from filenames alone. Preserve raw downloads unchanged. Keep handoffs file-based and short.

For scRNA clustering, marker annotation, cell communication, trajectory, or condition DE, use mature single-cell tools/skills after `geo-biodata-scrna-intake` confirms object content and metadata.
