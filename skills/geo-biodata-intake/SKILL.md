---
name: geo-biodata-intake
description: GEO_biodata intake workflow. Use when the user wants GEO metadata, sample tables, DOI/PMID/PMCID, open publication supplements, supplementary-file inventory, reviewed downloads, file-level QC, intake_handoff.yaml, or a run_manifest draft without running DE, enrichment, clustering, or marker analysis.
---

# GEO_biodata Intake

This module turns a GEO accession into a reviewed resource inventory and a file-based handoff.

## Executors

- `core/R/check_environment.R`
- `core/R/bootstrap_environment.R --profile intake --plan|--check|--install`
- `core/R/discover_geo.R <GSE> <run-dir>`
- `core/R/generate_dataset_report.R <run-dir>`
- `core/R/generate_download_plan.R <run-dir> <regex> [reason]`
- `core/R/download_reviewed_files.R <download_plan.tsv> <raw-dir>`
- `core/R/inspect_downloaded_input.R <run-dir> [input-file]`
- `core/R/generate_manifest_draft.R <intake_handoff.yaml> [run_manifest.draft.yaml]`

## Rules

- Use `discover_geo.R` first; it records series-level and sample-level supplements when GEO exposes both.
- Use `generate_dataset_report.R` for DOI/PMCID/abstract enrichment and a human-readable dataset summary. DOI lookup is helpful context, not a blocker.
- Download only rows that were selected and reviewed in `plans/download_plan.tsv`.
- Use the default download order: GEOquery, then aria2c, then curl, then `utils::download.file`.
- Preserve downloaded files unchanged; inspect and analyze derived copies.
- Do not confirm scale automatically. The handoff may suggest a route, but analysis still needs manual review flags in `run_manifest.yaml`.

## Outputs

Write:

- `resources/publication_links.tsv`
- `resources/publication_supplements.tsv`
- `dataset_report.md`
- `plans/download_plan.tsv`
- `raw/download_manifest.tsv`
- `tables/file_inventory.tsv`
- `tables/matrix_intake_audit.tsv`
- `tables/sample_overlap.tsv`
- `intake_handoff.yaml`
