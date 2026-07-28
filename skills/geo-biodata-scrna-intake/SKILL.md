---
name: geo-biodata-scrna-intake
description: GEO_biodata single-cell object intake workflow. Use when the user has Seurat, H5AD, RDS, SingleCellExperiment, MTX/H5, or author-provided scRNA objects and needs a read-only inventory of assays, layers, reductions, metadata fields, author labels, sample/donor fields, and raw-count availability without clustering or marker annotation.
---

# GEO_biodata scRNA Intake

This skill inventories single-cell objects. It does not recluster, annotate new cell types, compute markers, or run condition DE.

## Required Context

Read only the intake rules needed to inspect object content and decide whether downstream single-cell skills are appropriate:

- `../geo-biodata-workflow/references/basic-scrna.md`
- `../../knowledge/decision_rules/scrna_author_object.yaml`
- `../../knowledge/decision_rules/scrna_object_intake.yaml`
- `../../knowledge/decision_rules/pseudobulk_requirement.yaml`
- `../../docs/handoffs/scrna-clustering.md`
- `../../docs/handoffs/scrna-pseudobulk.md`

## Driver

Validate the manifest, then run:

```powershell
Rscript core\R\validate_manifest.R runs\GSE000000\run_manifest.yaml
Rscript core\R\scrna\inspect_object.R runs\GSE000000\run_manifest.yaml
```

## Outputs

Expect object inventory, metadata fields, assay/layer availability, author-label fields, sample/donor candidates, and a clear status.

If the user asks for clustering, marker review, annotation, trajectory, communication, or condition DE, stop after intake and hand off to dedicated local single-cell skills. Use the handoff docs to preserve:

- object path and format.
- raw-count assay/layer availability.
- sample/donor/case fields.
- author labels and reductions.
- suggested biological replication unit.

Condition DE must use sample/donor-level pseudobulk when raw counts and biological replicates are available. Cell-level tests are not valid biological replication.
