---
name: geo-biodata-scrna
description: GEO_biodata single-cell workflow for read-only object intake, QC planning, label provenance, and pseudobulk handoff from Seurat, H5AD, RDS, SingleCellExperiment, MTX, or H5 inputs.
---

# GEO_biodata scRNA

Use this skill for scRNA object inventory and conservative handoff planning. It may inspect Seurat/H5AD/RDS/MTX/H5 objects, but it does not recluster, relabel, run marker claims, or perform condition DE by itself.

## Driver

Validate the manifest, then run:

```powershell
Rscript core\R\validate_manifest.R runs\GSE000000\run_manifest.yaml
Rscript core\R\scrna\inspect_object.R runs\GSE000000\run_manifest.yaml
```

Supported routes are `scrna_author_object` and `scrna_raw_counts`.

## Intake Inventory

Before QC or annotation, record object format, matrix orientation, sparse/dense storage, raw-count layer or assay, normalized and scaled layers, cell and feature metadata fields, embeddings, graphs, author labels, cluster labels, sample/donor/capture fields, filtering notes, and conversion history.

Keep the original object or matrix unchanged. Avoid dense conversion of large sparse matrices. Do not assume AnnData and Seurat use the same matrix orientation. Do not use a converted object until slot/layer inventories have been compared. Preserve author labels in their original field and write any harmonized labels to a separate field.

If raw counts are absent, reproduce metadata, labels, and embeddings only. Do not describe the result as raw-count QC.

## QC And Analysis Boundary

Post-count work must be sample-aware and provenance-preserving:

1. Preserve raw counts; normalize only into a separate assay or layer.
2. Review detected genes, total UMI, mitochondrial fraction, and batch/capture partitions per sample before selecting thresholds.
3. Record threshold justification plots and retained or removed cell counts.
4. Assess ambient RNA and call doublets per sample or capture when supported; if unavailable, report the blocker.
5. Do not normalize already-normalized data.
6. Do not reflexively regress total counts, mitochondrial fraction, or cell cycle.
7. Use UMAP for visualization, not distance-based inference.

If clustering or annotation is requested after intake, use dedicated local single-cell tools. Cluster in PCA/neighbor graph space, sweep resolution, rank markers descriptively, compare inferred labels to author labels and supplements, and allow unknown or ambiguous labels. Treat marker tests as descriptive ranking, not proof of a new cell type or treatment effect.

## Pseudobulk Boundary

Condition-level scRNA DE must aggregate by biological unit and cell state when possible; cells are not independent replicates. Formal DE requires donor or sample ID, condition, cell type or reviewed cluster labels, and sufficient cells per unit-cell type. If donor/sample fields are missing, stop before condition-level inference.

Safe outputs before formal statistics include `cell_qc_metrics.tsv`, `qc_thresholds.tsv`, `sample_composition`, `cluster_resolution_summary.tsv`, `cluster_markers_descriptive.tsv`, and `annotation_confidence.tsv`.

Until GEO_biodata has dedicated pseudobulk executors, aggregate raw counts by sample x cell type with local single-cell tools, then pass each cell-type matrix to the bulk count DE route.
