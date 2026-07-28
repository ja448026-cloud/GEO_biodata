---
name: geo-biodata-workflow
description: >
  Legacy-compatible GEO_biodata entry point. Use for existing prompts that name
  geo-biodata-workflow; otherwise prefer narrower skills including geo-biodata-router,
  geo-biodata-intake-download, geo-biodata-bulk-de, geo-biodata-bulk-enrichment,
  geo-biodata-scrna-intake, geo-biodata-scrna-cluster-marker, geo-biodata-figure,
  or geo-biodata-maintainer.
---

# GEO_biodata Workflow

This compatibility skill keeps the original script bundle available. For new work, route to a narrower skill first.

## Choose A Narrow Skill

| Task | Prefer |
|---|---|
| GEO metadata, DOI/PMID/PMCID, publication supplements, download plan, reviewed download, file-level QC | `geo-biodata-intake-download` |
| bulk raw-count, bulk normalized, or microarray DE/QC from a validated manifest | `geo-biodata-bulk-de` |
| GO/KEGG/Reactome/MSigDB/ORA/GSEA from an existing DE table or rank vector | `geo-biodata-bulk-enrichment` |
| Seurat/H5AD/RDS/MTX object inventory without reclustering | `geo-biodata-scrna-intake` |
| scRNA reclustering, marker ranking, marker review, or label sanity checks | `geo-biodata-scrna-cluster-marker` |
| sample/donor-level pseudobulk condition DE within scRNA cell types | `geo-biodata-scrna-pseudobulk` |
| one simple source-table-linked diagnostic plot | `geo-biodata-figure` |
| repository structure, smoke tests, CI, fixtures, or release hygiene | `geo-biodata-maintainer` |

## Shared Script Bundle

The current stable scripts remain here:

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

Do not read method-heavy references unless the selected narrow skill requires them.

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
