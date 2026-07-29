---
name: geo-biodata-workflow
description: >
  Legacy-compatible GEO_biodata entry point. Use for existing prompts that name
  geo-biodata-workflow; otherwise prefer the slimmer current modules:
  geo-biodata, geo-biodata-intake, geo-biodata-bulk, geo-biodata-enrichment,
  geo-biodata-scrna, or geo-biodata-figure.
---

# GEO_biodata Workflow

Deprecated compatibility entry point. Keep it available for old prompts; do not add new workflow rules here.

This compatibility skill keeps the original script bundle available. For new work, route to `geo-biodata` first and then use the smallest current module.

## Choose A Narrow Skill

| Task | Prefer |
|---|---|
| GEO metadata, DOI/PMID/PMCID, publication supplements, download plan, reviewed download, file-level QC, manifest draft | `geo-biodata-intake` |
| bulk raw-count, bulk normalized, or microarray DE/QC from a validated manifest | `geo-biodata-bulk` |
| preranked GSEA from an existing rank table and local GMT file | `geo-biodata-enrichment` |
| Seurat/H5AD/RDS/MTX object inventory without reclustering | `geo-biodata-scrna` |
| figure planning, plot source linkage, and QA guidance | `geo-biodata-figure` |

## Shared Script Bundle

Stable public executors now live under `core/R`. These legacy paths remain available for wrappers and old prompts:

- `scripts/check_environment.R`
- `scripts/bootstrap_environment.R`
- `scripts/discover_geo.R`
- `scripts/generate_download_plan.R`
- `scripts/download_geo_supp.R`
- `scripts/validate_manifest.R`
- `scripts/drivers/run_bulk_counts.R`
- `scripts/drivers/run_bulk_normalized.R`
- `scripts/drivers/run_microarray.R`
- `scripts/drivers/inspect_scrna_object.R`
- `scripts/run_gsea_template.R`
- `scripts/marker_utilities.R`

Do not read method-heavy references unless the selected current module requires them.

## Stop States

Keep blocked states explicit:

- `RESOURCE_INVENTORY_COMPLETE`
- `DISCOVERY_PARTIAL`
- `MANIFEST_VALIDATED`
- `MANIFEST_INVALID`
- `EXECUTION_COMPLETE`
- `EXECUTION_FAILED`
- `BLOCKED`
- `RAW_COMPUTE_REQUIRED`

Do not silently convert an inventory/download result into an analysis result.
