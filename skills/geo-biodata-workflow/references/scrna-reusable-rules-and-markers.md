# Reusable scRNA Rules And Markers

Use this reference for generic scRNA-seq GEO analyses. Keep it lightweight: the goal is to help an agent avoid common mistakes and produce first-pass diagnostic figures, not to make final biological labels.

## Rules to reuse

1. Inspect input type before analysis: 10X MTX/H5, h5ad, Seurat RDS, SingleCellExperiment, or author tables.
2. Preserve author-provided labels in separate metadata fields. Do not overwrite them with inferred labels.
3. Calculate QC metrics before filtering: detected genes, counts/UMIs, mitochondrial percentage, ribosomal percentage when possible.
4. Plot QC distributions by sample or donor before setting thresholds. Fixed cutoffs are defaults for demonstration only, not final filtering decisions.
5. Record whether raw counts are available. If only normalized data or embeddings are available, do not claim full QC rerun.
6. For multiple samples, keep `sample`, `donor`, `condition`, and `batch` distinct. Never treat cells as independent biological replicates for donor-level inference.
7. If suspected doublets can be evaluated, record scores or a recommended tool route. Do not silently remove doublets without saving the score table and rule.
8. Explore clustering resolution lightly. One default resolution is acceptable for a first pass, but final labels require marker review and, where possible, author/reference agreement. Prefer Leiden clustering when the local Seurat installation supports it; record fallback to Louvain/default.
9. Use marker panels as sanity checks. Marker presence or dot plots can support annotation, but they do not replace manual review or reference-based annotation.
10. Stop at `REVIEW_REQUIRED` when sample mapping, donor structure, author labels, or QC decisions are ambiguous.
11. Batch integration is not automatic. Use it only after checking whether batch structure is technical, biological, or confounded with condition.
12. If author cell annotations exist, compare cluster-by-label tables before proposing harmonized labels.

## Generic human marker panels

These markers are intentionally broad and useful across many tumor or immune microenvironment datasets.

| Cell class | Markers |
|---|---|
| Epithelial | EPCAM, KRT8, KRT18, KRT19, KRT7 |
| T cell | CD3D, CD3E, TRAC, IL7R |
| CD8 T cell | CD8A, CD8B, GZMK |
| Treg | FOXP3, IL2RA, CTLA4 |
| NK or cytotoxic | NKG7, GNLY, PRF1, GZMB |
| B cell | MS4A1, CD79A, CD79B, CD19 |
| Plasma cell | MZB1, JCHAIN, XBP1, IGKC |
| Myeloid or monocyte | LST1, LYZ, S100A8, S100A9, FCN1 |
| Macrophage | C1QA, C1QB, APOE, CD68 |
| Dendritic cell | FCER1A, CLEC10A, CST3, LILRA4 |
| Fibroblast | COL1A1, COL1A2, DCN, LUM |
| Endothelial | PECAM1, VWF, KDR, PLVAP |
| Platelet | PPBP, PF4 |
| Mast cell | TPSAB1, TPSB2, KIT, CPA3 |
| Cycling | MKI67, TOP2A, STMN1, HMGB2 |
| Stress or dissociation | FOS, JUN, JUNB, DUSP1, ATF3, HSPA1A |

## Generic mouse marker notes

Use mouse-case equivalents where possible: epithelial `Epcam/Krt8/Krt18/Krt19`, T cell `Cd3d/Cd3e/Trac`, B cell `Ms4a1/Cd79a`, myeloid `Lyz2/S100a8/S100a9`, macrophage `C1qa/C1qb/Apoe`, fibroblast `Col1a1/Dcn/Lum`, endothelial `Pecam1/Vwf/Kdr`, cycling `Mki67/Top2a`.

## Minimal figures

Produce only diagnostic figures in the public workflow:

- QC violin plots by sample.
- PCA elbow plot.
- UMAP by cluster and sample.
- Clustering resolution summary.
- Dot plot for generic marker panels.
- Cluster marker table and top-marker dot plot.
- Author-label comparison table when labels exist.
