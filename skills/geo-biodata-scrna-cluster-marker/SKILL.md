---
name: geo-biodata-scrna-cluster-marker
description: GEO_biodata single-cell clustering and marker-review workflow. Use when the user explicitly asks to recluster cells, choose resolution, build UMAP/tSNE/PCA views, find cluster markers, score marker panels, or review cell-type labels. Do not use for bulk transcriptome DE.
---

# GEO_biodata scRNA Cluster And Marker Review

Use only after scRNA intake confirms a usable count-bearing object or reviewed processed object.

## Required External Skill Guidance

Use the local single-cell skills as method references when available:

- `bio-single-cell-clustering`
- `bio-single-cell-markers-annotation`
- `bio-single-cell-preprocessing`
- `bio-single-cell-doublet-detection` when doublet risk is relevant

## Rules

- Cluster in PCA/neighbor-graph space, not UMAP coordinates.
- Sweep resolution when reclustering; do not tune one value until labels look convenient.
- Treat marker tests as descriptive ranking for annotation.
- Do not report cluster marker p-values as proof that a cluster is a real biological population.
- Do not test treatment-vs-control effects cell-by-cell. Use pseudobulk with sample/donor as the unit.

## Marker Content Boundary

Generic marker panels and marker utilities are optional. Load marker rules only when the user asks for labels, markers, or signature scores.

