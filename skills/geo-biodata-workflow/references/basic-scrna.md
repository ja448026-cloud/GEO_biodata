# Basic single-cell route

Read `scrna-reusable-rules-and-markers.md` before assigning or comparing cell labels. Use the marker panels there as first-pass diagnostic checks only.

## Intake checks

1. Identify the raw-count layer and confirm non-negative integer-like values.
2. Map each cell to a sample or capture channel; retain donor/subject when available.
3. Import author-provided cell labels without overwriting them.
4. Check species-specific mitochondrial gene naming.
5. Record whether ambient-RNA and doublet assessment can be performed.

If the processed object lacks raw counts, reproduce its metadata and embeddings only; do not describe the result as raw-count QC.

## Minimum analysis

1. Calculate per-cell detected genes, total counts, and mitochondrial fraction.
2. Plot distributions per sample before selecting thresholds.
3. Record proposed thresholds and retained-cell counts; use `REVIEW_REQUIRED` when judgment is needed.
4. Preserve counts, then normalize in a separate layer/assay.
5. Select variable features, run PCA, choose a documented PC range, build a neighbor graph, cluster, and calculate UMAP.
6. Find cluster markers for descriptive annotation.
7. Compare inferred labels with author labels and supplements; allow unknown or ambiguous labels.

Do not use cells as independent replicates for donor-level differential claims.

## Minimum figures

- per-sample QC distributions and before/after retention;
- PCA variance or elbow diagnostic;
- UMAP by cluster and by sample/donor;
- cluster-size table or bar plot;
- marker dot plot or small marker heatmap;
- annotation agreement/confidence summary when labels are available.

## Method references

- [Seurat standard workflow](https://satijalab.org/seurat/articles/pbmc3k_tutorial)
- [Seurat visualization](https://satijalab.org/seurat/articles/visualization_vignette)
- [SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment/)
