# Public MVP Plan

Date: 2026-07-27

## Goal

Build `geo_biodata_workflow` as a lightweight public workflow for agents. A user should be able to provide a GEO Series accession and have an agent follow a repeatable path for public resource discovery, environment checks, selective downloads, basic QC, simple transcriptome or scRNA analysis, and diagnostic figures.

## Design Boundary

This project should be a workflow and skill, not a heavy R package.

Keep deterministic operations as small scripts:

- environment check;
- GEO metadata and supplement discovery;
- publication supplement lookup;
- reviewable download plans, selective download, and hashes;
- manifest validation before analysis;
- a manifest-driven `bulk_raw_counts` driver plus optional templates for GSEA and future scRNA first-pass analysis.

Keep method judgment as short references:

- routing and resource triage;
- bulk DE and GSEA principles;
- scRNA QC, clustering, annotation, and marker-panel principles;
- stop states and reporting expectations.

Do not copy whole external skill libraries into this repo. Use their reusable patterns to write concise local guidance.

## Release Boundary For v0.1

Include:

- `skills/geo-biodata-workflow/SKILL.md`;
- `agents/openai.yaml`;
- lightweight R scripts;
- concise reference files for routing, bulk/DE/GSEA, and scRNA;
- validation notes for syntax, skill validation, resource discovery, and download guard behavior.

Do not include:

- local absolute paths;
- downloaded GEO data;
- generated run directories;
- project-specific biology;
- FASTQ alignment/quantification;
- heavy single-cell integration, trajectory, or communication workflows;
- claims that every GEO record can be analyzed automatically.

## Method Modules

| Module | Public route | Main dependency |
|---|---|---|
| GEO discovery | metadata, GSM index, supplements, publication links | GEOquery, httr2 |
| Resource collection | reviewed download plan, manifest, hashes | GEOquery, digest |
| Analysis readiness | run manifest validation and review gates | yaml |
| `bulk_raw_counts` | QC, PCA, DE, p-value histogram | DESeq2 |
| `bulk_normalized` / `microarray_series_matrix` | QC, PCA, limma-compatible DE when implemented | limma |
| Enrichment | ranked-list GSEA and guidance for GO/KEGG/Reactome/MSigDB | fgsea, clusterProfiler |
| scRNA | QC metrics, clustering, UMAP, markers, author-label comparison | Seurat |

## Promotion Criteria

- A clean checkout contains no private path or generated data.
- The skill folder passes the host agent's skill validator.
- R scripts parse without executing downloads.
- The README explains the one-accession workflow without overstating automation.
- The workflow stops explicitly when metadata or inputs are insufficient.
- Download transfer requires a route-specific plan and reviewed file rows.
- Statistical analysis requires a validated `run_manifest.yaml`.
- Deprecated all-in-one analysis templates fail closed and the `bulk_raw_counts` route has a manifest-driven driver.
- CI runs deterministic smoke checks for parseability, guards, manifest validation, and public-tree hygiene.
- The skill can be installed by another agent without hidden local context.

## Later Work

Add real examples only after the structure is stable. Prefer small, public, processed examples and keep downloaded files out of Git. Extract an R package only if repeated use reveals stable functions that are awkward as scripts.
