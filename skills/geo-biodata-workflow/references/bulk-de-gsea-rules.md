# Bulk DE And GSEA Rules

Use this reference for ordinary bulk RNA-seq, expression matrix, and microarray GEO analyses. Keep the implementation thin and rely on established R packages.

## Differential expression rules

1. Start from input readiness: rows are features, columns are samples, sample IDs match metadata, and the biological unit is clear.
2. Separate raw counts from normalized expression. Do not pass normalized or microarray values to a raw-count model.
3. Set factor levels and contrast direction explicitly before testing. Record the numerator and denominator group.
4. Include known batch, pair, donor, or subject structure in the design when the metadata supports it. Do not subtract batch and then run DE on the corrected matrix.
5. Use DESeq2 or edgeR for raw integer counts; use limma or limma-voom for microarray, normalized expression, or voom-weighted count data.
6. Use design-aware low-expression filtering for counts. For edgeR, prefer `filterByExpr`; for DESeq2, pre-filter only for speed and rely on independent filtering at results time.
7. Treat very small sample designs as exploratory. Do not overstate FDR-stable biology when replication is weak.
8. For visualization, use VST, rlog, logCPM, or an appropriate log-scale matrix. Never use raw counts for PCA or heatmaps.
9. Save the design table, contrast, package versions, DE table, and source data for every diagnostic plot.
10. For DESeq2 effect-size reporting, consider `lfcShrink()` after the primary Wald or LRT result when stable ranking or volcano labels matter.
11. Add a p-value histogram and MA plot as sanity checks. Strong asymmetry, spikes, or all-flat p-values usually mean the design or input scale needs review.
12. For STAR/RSEM-style files, verify strandedness and whether TPM-derived quantities such as `lengthScaledTPM` were converted back to count-like values. Do not treat scaled TPM as raw counts.
13. Harmonize pseudoautosomal annotations explicitly; for example, avoid letting `PAR_Y` suffix variants create duplicate male chromosome genes without review.

## GSEA rules

1. Run GSEA only when a signed ranked vector is available for most genes. Use ORA instead for a small unranked gene list.
2. Prefer a signed, variance-calibrated statistic: DESeq2 Wald `stat`, limma moderated `t`, or a signed p-value transform when no statistic exists.
3. Do not rank by raw p-value alone; it erases direction. Avoid bare log2FC ranking when low-count noise dominates.
4. Deduplicate gene IDs before running GSEA. One gene should contribute one statistic.
5. Sort the named numeric vector decreasing before testing.
6. Record the gene-set source, ID type, date or package version, min/max gene-set size, exponent, seed, and multiple-testing method.
7. For preranked fgsea/clusterProfiler results, state that the null is gene permutation. If a matrix and design are available and correlation-aware testing matters, use CAMERA/fry/roast from limma.
8. Interpret NES direction relative to the contrast: positive means enriched at the numerator/up end of the ranking.
9. Inspect leading-edge genes before treating a pathway as meaningful.
10. For ORA, define the universe as tested genes after filtering, not all genes in the genome.
11. For GO over-representation, reduce redundant terms with a documented semantic-similarity or parent-child rule when many related terms dominate.
12. For KEGG or online annotation sources, cache the exact gene-set table or package version used for the run.

## Minimal outputs

- `sample_mapping_used.tsv`
- `de_results_<contrast>.tsv`
- `pca_coordinates.tsv`
- `bulk_pca.pdf`
- `bulk_ma_<contrast>.pdf`
- `bulk_pvalue_histogram_<contrast>.pdf`
- `bulk_volcano_<contrast>.pdf`
- `gsea_results_<collection>.tsv` when GSEA is run
- `gsea_top_<collection>.pdf` or another compact diagnostic plot

## Package preference

- Raw counts: DESeq2 or edgeR.
- Microarray or normalized expression: limma.
- Pre-ranked GSEA with a local gene-set list or GMT: fgsea.
- GO/KEGG/Reactome/MSigDB with annotation support: clusterProfiler and related annotation packages.
- Correlation-aware gene-set testing with a matrix and design: limma CAMERA, fry, or roast.
