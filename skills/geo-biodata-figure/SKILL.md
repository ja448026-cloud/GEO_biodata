---
name: geo-biodata-figure
description: GEO_biodata figure reference and plotting router. Use when the user wants simple source-linked visualizations from GEO_biodata outputs, figure planning, figure QA, or one-step diagnostic plot polishing. This skill delegates plotting style and QA to existing local figure skills instead of reimplementing a full plotting framework.
---

# GEO_biodata Figure

Use this as a light visualization layer. It is not a statistical module.

## Execution Boundary

Prefer plots already written by the active route driver:

- bulk library size, PCA, p-value histogram, mean-variance, sample correlation.
- GSEA NES dot/bar plot from `core/R/enrichment/run_preranked_gsea.R`.
- scRNA inventory plots only after a dedicated scRNA executor supports them.

For custom plotting, use local figure skills as the executor:

- `figure-planner`: decide the single purpose of the plot.
- `nature-figure`: backend-specific R/Python scientific plotting and export.
- `omics-figure-qa`: source-table linkage, overcrowding, pathway-label readability, and panel QA.

## Contract

Every GEO_biodata figure should record:

- source table or object path.
- script or driver that generated it.
- one-sentence purpose.
- filtering rule.
- output path.

Do not combine unrelated panels by default. One simple figure per step is preferred.

