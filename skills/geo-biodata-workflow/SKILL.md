---
name: geo-biodata-workflow
description: >
  Inspect a GEO Series accession (GSE), collect its public metadata, supplementary files,
  publication links, and author-provided annotations, then run a conservative basic bulk
  transcriptome or single-cell RNA-seq workflow. Use when a user gives a GSE accession and
  wants reproducible download, resource inventory, input routing, basic QC, clustering or
  differential analysis, and minimal diagnostic figures without project-specific assumptions.
---

# geo_biodata_workflow

Treat this skill as a reusable execution path — a coordination layer over existing public tools. It is NOT a universal analysis engine. Prefer author-provided processed matrices, objects, metadata, and supplementary tables before considering raw FASTQ processing.

## Required order

1. Create a new accession-specific run directory under `runs/<GSE>/`. Never mix two GSE studies. Reserve `validation/runs/` for project-maintainer smoke tests only.
2. Run `scripts/check_environment.R` and record missing optional capabilities.
3. Run `scripts/discover_geo.R <GSE> <run-dir>` before downloading expression files.
4. Review the resulting GEO metadata, sample index, supplement index, publication links, species, assay type, sample count, and likely biological unit.
5. Read `references/routing-and-resources.md` and select exactly one input route.
6. Generate a route-specific download plan with `scripts/generate_download_plan.R <run-dir> <regex> [reason]`.
7. Review `plans/download_plan.tsv`, set `reviewed=TRUE` only for files confirmed to match the selected route, then download with `scripts/download_geo_supp.R <download-plan.tsv> <raw-dir>`.
8. Preserve downloaded files unchanged. Record hashes and analyze copies or derived objects.
9. Create `run_manifest.yaml` from `templates/run_manifest.example.yaml`; the manifest is the analysis fact source for route, input type, biological unit, design, and contrast.
10. Validate the manifest with `scripts/validate_manifest.R <run_manifest.yaml>`. Do not run statistical analysis unless the state is `MANIFEST_VALIDATED`.
11. For bulk or microarray data, read `references/basic-bulk.md` plus `references/bulk-de-gsea-rules.md`; use route-specific manifest-driven drivers when available.
12. For single-cell data, read `references/basic-scrna.md` plus `references/scrna-reusable-rules-and-markers.md`; keep `scrna_raw_counts` and `scrna_author_object` separate.
13. Save the input inventory, sample mapping, code, tables, figures, package versions, warnings, and final status.

## Allowed outcomes

- `RESOURCE_INVENTORY_COMPLETE`: resources were collected but analysis has not started.
- `DISCOVERY_PARTIAL`: public-resource discovery ran, but one or more discovery steps failed or returned incomplete metadata.
- `READY_FOR_BASIC_ANALYSIS`: input and metadata support a basic workflow.
- `MANIFEST_VALIDATED`: input, route, sample mapping, design, contrast, and review gates passed structural validation.
- `MANIFEST_INVALID`: manifest validation failed; analysis must not run.
- `BASIC_ANALYSIS_COMPLETE`: required outputs and checks exist.
- `REVIEW_REQUIRED`: a human must confirm sample mapping, QC thresholds, or labels.
- `BLOCKED_INPUT`: no supported processed input is available.
- `BLOCKED_METADATA`: biological units or groups cannot be reconstructed safely.
- `RAW_COMPUTE_REQUIRED`: only FASTQ or another raw sequencing route is available.

Do not silently convert a blocked state into an analysis result.

## Minimum output contract

```text
<GSE>/
  resources/
    series_metadata.tsv
    sample_index.tsv
    sample_characteristics.tsv   (extracted GSM characteristics)
    supplement_index.tsv
    publication_links.tsv
    publication_supplements.tsv  (supplementary links from publisher)
    routing_evidence.tsv
  plans/
    download_plan.tsv
  raw/
    download_manifest.tsv
  derived/
  tables/
  figures/
  logs/
  scripts/
  environment.tsv
  manifest_validation.tsv
  workflow_status.tsv
  workflow_events.tsv
  summary.md
```

Every figure must have a source table or documented object and the script that generated it. Produce diagnostic figures only; they are not publication-ready.

## Guardrails

- Verify species, tissue, assay, sample identifiers, donor or subject identifiers, condition, and replicate structure.
- Review SuperSeries/SubSeries relationships before deciding which GSE is the analysis unit.
- Prefer GEO/sample metadata and author supplements over filename inference.
- Record whether matrices are raw counts, normalized expression, TPM/FPKM/CPM, microarray intensities, or author-processed object layers in `run_manifest.yaml` before selecting DE tools.
- Preserve author labels in a separate field; never overwrite them with inferred labels.
- Do not apply fixed single-cell QC cutoffs without inspecting per-sample distributions.
- Do not run donor-level inference by treating cells as independent replicates.
- Do not scrape paywalled papers. Record DOI/PMID/PMCID and use lawful open resources or user-supplied files.
- Do not download every supplementary file when a reviewed filename filter is sufficient. Generate a download plan first and require `reviewed=TRUE` before transfer.
- Stop at `RAW_COMPUTE_REQUIRED` rather than implementing alignment or quantification ad hoc.

## Analysis template scripts

Deprecated template scripts in `scripts/` are retained as traceable compatibility markers:

- `analyze_bulk_template.R` — deprecated and fail-closed; do not use as an automatic analysis entry point.
- `run_gsea_template.R` — lightweight fgsea preranked GSEA template for a named rank vector and local gene sets.
- `analyze_scrna_template.R` — deprecated and fail-closed; keep raw-count and author-object routes separate.
- `marker_utilities.R` — Shared gene signature scoring, generic marker panels, marker-presence checks, and marker table helpers.
- `generate_download_plan.R` — creates a reviewable `plans/download_plan.tsv` from the supplement index before any file transfer.
- `validate_manifest.R` — validates `run_manifest.yaml` before analysis and writes `manifest_validation.tsv`.

The deprecated analysis templates intentionally stop immediately. Route-specific drivers should read a validated manifest and must not infer statistical models from filename patterns or expression-value ranges.

## Stable interface

- Entry point: one GEO Series accession number.
- Core language: R.
- Required packages: GEOquery, Biobase, httr2, jsonlite, digest, yaml.
- Optional packages: limma, edgeR, DESeq2, ggplot2, pheatmap, fgsea, clusterProfiler, msigdbr, Seurat, Matrix, patchwork.
- Required outputs: series metadata, sample index, supplement index, publication links, download manifest when files are downloaded, workflow status, summary, and diagnostic figures when analysis runs.

## Upstream tools

This workflow coordinates established tools; it does not vendor them:

- [GEOquery](https://bioconductor.org/packages/GEOquery/) for GEO access;
- [nf-core/fetchngs](https://nf-co.re/fetchngs/latest/) for optional raw-data handoff;
- [Seurat](https://satijalab.org/seurat/) for common scRNA route;
- [DESeq2](https://bioconductor.org/packages/DESeq2/), [edgeR](https://bioconductor.org/packages/edgeR/), and [limma](https://bioconductor.org/packages/limma/) for input-appropriate bulk routes;
- [Europe PMC](https://europepmc.org/RestfulWebService) for publication metadata and open resources.
