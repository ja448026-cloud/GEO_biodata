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
6. Download only the files needed for that route with `scripts/download_geo_supp.R <GSE> <raw-dir> <regex>`.
7. Preserve downloaded files unchanged. Record hashes and analyze copies or derived objects.
8. For bulk or microarray data, read `references/basic-bulk.md` plus `references/bulk-de-gsea-rules.md`, then optionally use `scripts/analyze_bulk_template.R` as a starting point.
9. For single-cell data, read `references/basic-scrna.md` plus `references/scrna-reusable-rules-and-markers.md`, then optionally use `scripts/analyze_scrna_template.R` as a starting point.
10. Save the input inventory, sample mapping, code, tables, figures, package versions, warnings, and final status.

## Allowed outcomes

- `RESOURCE_INVENTORY_COMPLETE`: resources were collected but analysis has not started.
- `READY_FOR_BASIC_ANALYSIS`: input and metadata support a basic workflow.
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
  raw/
    download_manifest.tsv
  derived/
  tables/
  figures/
  logs/
  scripts/
  environment.tsv
  workflow_status.tsv
  summary.md
```

Every figure must have a source table or documented object and the script that generated it. Produce diagnostic figures only; they are not publication-ready.

## Guardrails

- Verify species, tissue, assay, sample identifiers, donor or subject identifiers, condition, and replicate structure.
- Review SuperSeries/SubSeries relationships before deciding which GSE is the analysis unit.
- Prefer GEO/sample metadata and author supplements over filename inference.
- Record whether matrices are raw counts, normalized expression, TPM/FPKM/CPM, microarray intensities, or author-processed object layers before selecting DE tools.
- Preserve author labels in a separate field; never overwrite them with inferred labels.
- Do not apply fixed single-cell QC cutoffs without inspecting per-sample distributions.
- Do not run donor-level inference by treating cells as independent replicates.
- Do not scrape paywalled papers. Record DOI/PMID/PMCID and use lawful open resources or user-supplied files.
- Do not download every supplementary file when a reviewed filename filter is sufficient.
- Stop at `RAW_COMPUTE_REQUIRED` rather than implementing alignment or quantification ad hoc.

## Analysis template scripts

Template scripts in `scripts/` provide a starting point for the agent:

- `analyze_bulk_template.R` — DESeq2/limma workflow skeleton for bulk data; use the rules reference when writing an edgeR-specific driver.
- `run_gsea_template.R` — lightweight fgsea preranked GSEA template for a named rank vector and local gene sets.
- `analyze_scrna_template.R` — Seurat workflow skeleton for single-cell data.
- `marker_utilities.R` — Shared gene signature scoring, generic marker panels, marker-presence checks, and marker table helpers.

These templates are intentionally skeletal. The agent must customize thresholds, group definitions, contrast specifications, and figure parameters per dataset. They exist to reduce repetitive boilerplate and to enforce consistent output conventions, not to replace dataset-specific judgment.

Set `template_dir` before sourcing templates from outside the `scripts/` directory; otherwise the templates look for `marker_utilities.R` in the current working directory.

## Stable interface

- Entry point: one GEO Series accession number.
- Core language: R.
- Required packages: GEOquery, Biobase, httr2, jsonlite, digest.
- Optional packages: limma, edgeR, DESeq2, ggplot2, pheatmap, fgsea, clusterProfiler, msigdbr, Seurat, Matrix, patchwork.
- Required outputs: series metadata, sample index, supplement index, publication links, download manifest when files are downloaded, workflow status, summary, and diagnostic figures when analysis runs.

## Upstream tools

This workflow coordinates established tools; it does not vendor them:

- [GEOquery](https://bioconductor.org/packages/GEOquery/) for GEO access;
- [nf-core/fetchngs](https://nf-co.re/fetchngs/latest/) for optional raw-data handoff;
- [Seurat](https://satijalab.org/seurat/) for common scRNA route;
- [DESeq2](https://bioconductor.org/packages/DESeq2/), [edgeR](https://bioconductor.org/packages/edgeR/), and [limma](https://bioconductor.org/packages/limma/) for input-appropriate bulk routes;
- [Europe PMC](https://europepmc.org/RestfulWebService) for publication metadata and open resources.
