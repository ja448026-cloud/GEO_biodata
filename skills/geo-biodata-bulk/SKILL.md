---
name: geo-biodata-bulk
description: GEO_biodata bulk and microarray analysis workflow. Use when a validated manifest needs bulk raw-count DESeq2 analysis, normalized-expression limma analysis, microarray limma analysis, sample QC, design/contrast validation, diagnostic plots, output integrity, and ranked DE tables. Do not use for scRNA clustering, marker annotation, or cell-level condition DE.
---

# GEO_biodata Bulk

Use only after intake identifies a matrix and `run_manifest.yaml` passes validation.

## Required References

Read only the bulk-relevant references:

- `skills/geo-biodata-workflow/references/basic-bulk.md`
- `skills/geo-biodata-workflow/references/bulk-de-gsea-rules.md`
- `knowledge/route_maturity.yaml`

## Executors

```powershell
Rscript core\R\validate_manifest.R runs\GSE000000\run_manifest.yaml
Rscript core\R\drivers\run_bulk_counts.R runs\GSE000000\run_manifest.yaml
Rscript core\R\drivers\run_bulk_normalized.R runs\GSE000000\run_manifest.yaml
Rscript core\R\drivers\run_microarray.R runs\GSE000000\run_manifest.yaml
```

Choose exactly one driver by manifest route:

- `bulk_raw_counts`: DESeq2.
- `bulk_normalized`: limma for confirmed transformed normalized matrices; EDA-only or blocked when scale is not DE-safe.
- `microarray_series_matrix`: limma with platform/mapping audit.

Do not pass normalized values to raw-count models. For paired designs, model the patient/donor term explicitly.

## Rules

Biological inference uses the donor, patient, animal, or independently prepared sample as the unit, not cells, spots, reads, or technical sections. Record `sample_mapping.biological_unit`, count independent biological units per contrast level, and stop with `BLOCKED_METADATA` or `REVIEW_REQUIRED` when the unit cannot be reconstructed or a contrast level has fewer than two biological units.

Choose the statistical model from verified input scale, not filenames. `bulk_raw_counts` requires a nonnegative integer-like gene-by-sample matrix and uses DESeq2 or edgeR-style count models; `bulk_normalized` requires finite normalized/log expression or abundance and uses limma or exploratory QC; `microarray_series_matrix` requires a platform-reviewed expression matrix or ExpressionSet and uses limma.

Verify matrix columns exactly match reviewed sample mapping. Reject negative or fractional values for raw-count routes. Do not infer raw-count status from names such as count, matrix, expression, or data.

Paired or blocked designs must be explicit before analysis: add the blocking factor to sample mapping, include it in `design.formula`, and confirm the contrast factor is estimable after blocking. Stop when pair identifiers are implied but unavailable or the model matrix is not full rank.

Enrichment is downstream. Use `geo-biodata-enrichment` only after a ranked table exists.
