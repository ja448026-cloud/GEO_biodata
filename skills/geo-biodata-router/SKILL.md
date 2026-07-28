---
name: geo-biodata-router
description: Route GEO_biodata tasks to the narrow skill needed for the user's current goal. Use when a user gives a GEO accession or GEO_biodata project request and may need download-only intake, bulk DE, enrichment, scRNA object intake, scRNA clustering/marker review, plotting, or maintainer validation.
---

# GEO_biodata Router

Use this as the lightweight entry point. Do not load analysis-specific references until the task requires them.

## Route by Intent

| User intent | Use skill | Heavy dependencies |
|---|---|---|
| Find GEO files, metadata, DOI, supplements, or download reviewed files | `geo-biodata-intake-download` | no |
| Validate a manifest and run bulk/microarray DE/QC | `geo-biodata-bulk-de` | yes, bulk only |
| Run GO/KEGG/Reactome/GSEA from a DE/rank table | `geo-biodata-bulk-enrichment` | yes, enrichment only |
| Inspect Seurat/H5AD/RDS/MTX object contents | `geo-biodata-scrna-intake` | moderate |
| Recluster cells or review cluster markers | `geo-biodata-scrna-cluster-marker` | yes, scRNA only |
| Test condition DE within scRNA cell types | `geo-biodata-scrna-pseudobulk` | yes, scRNA plus bulk counts |
| Make one diagnostic plot from an existing source table | `geo-biodata-figure` | depends on backend |
| Edit project structure, CI, fixtures, or release state | `geo-biodata-maintainer` | no by default |

## Default Behavior

Start with intake/download unless the user explicitly asks for analysis or provides a validated manifest and matrix.

Keep raw downloads unchanged. Do not infer statistical models from file names alone. Use route-specific manifests for analysis.

If the user only asks to download or inventory a GEO accession, stop after the resource/download contract and do not read bulk, scRNA, marker, enrichment, or figure references.
