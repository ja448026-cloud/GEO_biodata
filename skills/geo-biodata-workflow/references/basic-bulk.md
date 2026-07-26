# Basic bulk or microarray route

Read `bulk-de-gsea-rules.md` before running differential expression or GSEA. Use it as the compact policy layer for model choice, contrast direction, batch handling, and ranked-list enrichment.

## Intake checks

1. Determine whether values are raw counts, normalized expression, or microarray intensities.
2. Confirm that rows are genes/probes and columns are samples.
3. Preserve original feature identifiers and record any annotation mapping.
4. Build `sample_mapping.tsv` from GEO characteristics and author supplements.
5. Confirm the experimental unit, groups, pairs, batches, and requested contrast.

If group or biological-unit mapping is unresolved, produce EDA only and set `BLOCKED_METADATA`.

## Minimum analysis

### Raw counts

- Filter very lowly expressed features with a recorded rule.
- Plot library sizes and detected-feature counts.
- Apply a documented variance-stabilizing or log-CPM transform for visualization.
- Plot sample correlation and PCA/MDS, colored by group and labeled by sample.
- Use DESeq2 or edgeR for differential analysis only when the design matrix is valid.

### Normalized or microarray expression

- Inspect value distributions and missingness.
- Plot sample distributions, correlation heatmap, and PCA.
- Use limma only after confirming the scale and design.
- Do not pass normalized values to a raw-count model.

## Minimum figures

- sample/library QC;
- sample correlation heatmap;
- PCA or MDS;
- if differential analysis is valid: MA or mean-difference plot, volcano plot, and a small top-feature heatmap.

Save the numerical source of every plot. Avoid decorative composites and biological claims not supported by the study design.

## Method references

- [DESeq2 documentation](https://bioconductor.org/packages/DESeq2/)
- [edgeR documentation](https://bioconductor.org/packages/edgeR/)
- [limma documentation](https://bioconductor.org/packages/limma/)
