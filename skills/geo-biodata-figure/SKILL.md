---
name: geo-biodata-figure
description: GEO_biodata lightweight diagnostic plotting workflow. Use when the user wants one simple source-table-linked visualization from GEO_biodata outputs, such as PCA, sample correlation heatmap, library size, volcano, MA, p-value histogram, GO dotplot, GSEA NES barplot, UMAP, or marker dotplot.
---

# GEO_biodata Figure

Use this to make one diagnostic figure from an existing table or object. Do not build large composite figures by default.

## Figure Contract

Before plotting, record:

- source table or object path.
- one-sentence purpose of the plot.
- plot type.
- filtering rule.
- output path.

## Style Defaults

- One plot, one message.
- Prefer PDF/SVG for editable scientific plots and PNG only for quick preview.
- Keep labels readable; shorten long pathway names in the plot and keep full names in the source table.
- Use restrained, consistent colors.
- Do not change analysis values while polishing plots.

## Reference Skills

Use local figure skills as needed:

- `figure-planner` for claim/panel logic.
- `nature-figure` for backend-specific publication plotting and export standards.
- `omics-figure-qa` for source-table linkage, overcrowding, and panel usability checks.

