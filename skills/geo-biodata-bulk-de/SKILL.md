---
name: geo-biodata-bulk-de
description: GEO_biodata bulk and microarray differential-expression workflow. Use when a validated GEO matrix needs bulk raw-count DESeq2/edgeR-style analysis, normalized-expression limma analysis, microarray limma analysis, sample QC, design/contrast validation, and ranked DE tables. Do not use for scRNA clustering or marker annotation.
---

# GEO_biodata Bulk DE

Use only after intake/download has identified a matrix and sample mapping.

## Required Context

Read only the bulk-relevant references:

- `../geo-biodata-workflow/references/basic-bulk.md`
- `../geo-biodata-workflow/references/bulk-de-gsea-rules.md`
- `../../knowledge/decision_rules/bulk_input_scale.yaml`
- `../../knowledge/decision_rules/paired_design.yaml`
- `../../knowledge/decision_rules/biological_replication.yaml`

## Drivers

Validate the manifest before running any driver:

```powershell
Rscript skills\geo-biodata-workflow\scripts\validate_manifest.R runs\GSE000000\run_manifest.yaml
```

Then choose exactly one route:

- `bulk_raw_counts`: `drivers/run_bulk_counts.R`
- `bulk_normalized`: `drivers/run_bulk_normalized.R`
- `microarray_series_matrix`: `drivers/run_microarray.R`

Do not pass normalized values to raw-count models. For paired designs, model the patient/donor term explicitly.

## Outputs

Expect sample mapping, design matrix, contrast matrix, DE table, QC checks, fallback events, output integrity, library/PCA/sample-correlation diagnostics, and session information.

Enrichment is downstream. Use `geo-biodata-bulk-enrichment` only after a DE or ranked table exists.

