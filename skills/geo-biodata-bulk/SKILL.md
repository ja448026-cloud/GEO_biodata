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
- `knowledge/decision_rules/bulk_input_scale.yaml`
- `knowledge/decision_rules/paired_design.yaml`
- `knowledge/decision_rules/biological_replication.yaml`
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

Enrichment is downstream. Use `geo-biodata-enrichment` only after a ranked table exists.
