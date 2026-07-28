---
name: geo-biodata-scrna-intake
description: GEO_biodata single-cell object intake workflow. Use when the user has Seurat, H5AD, RDS, SingleCellExperiment, MTX/H5, or author-provided scRNA objects and needs a read-only inventory of assays, layers, reductions, metadata fields, author labels, sample/donor fields, and raw-count availability without clustering or marker annotation.
---

# GEO_biodata scRNA Intake

This skill inventories single-cell objects. It does not recluster or assign cell types by default.

## Required Context

Read only scRNA intake rules:

- `../geo-biodata-workflow/references/basic-scrna.md`
- `../../knowledge/decision_rules/scrna_author_object.yaml`
- `../../knowledge/decision_rules/scrna_object_intake.yaml`
- `../../knowledge/decision_rules/pseudobulk_requirement.yaml`

## Driver

Validate the manifest, then run:

```powershell
Rscript skills\geo-biodata-workflow\scripts\drivers\inspect_scrna_object.R runs\GSE000000\run_manifest.yaml
```

## Outputs

Expect object inventory, metadata fields, assay/layer availability, author-label fields, sample/donor candidates, and a clear status.

If raw counts and sample-level biological replication are present and the user asks for condition DE, route to pseudobulk rather than cell-level DE.

