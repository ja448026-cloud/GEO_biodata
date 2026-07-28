# GEO_biodata

`GEO_biodata` is a lightweight, agent-oriented workflow for public GEO biodata.

Give an agent one GEO Series accession such as `GSE000000`; the project guides it through:

1. environment checks;
2. GEO metadata and resource discovery;
3. publication and supplementary-material lookup;
4. selective download of reviewed files;
5. manifest validation before analysis;
6. basic analysis only after route-specific drivers and review gates support it;
7. diagnostic tables, logs, and a clear stop state.

This is not a new bioinformatics framework and not a full R package. It is a reusable execution path over established tools such as GEOquery, limma, DESeq2, edgeR, fgsea, clusterProfiler, and Seurat.

## Why this exists

Many GEO studies already provide processed matrices, Seurat/H5AD objects, sample metadata, cell annotations, or supplementary tables. A useful agent workflow should not immediately start from FASTQ or invent labels from filenames. It should first collect public resources, decide the safest analysis route, and only run the basic analysis that the available inputs justify.

The core idea is:

```text
GSE accession -> resource inventory -> input routing -> selective download -> manifest validation -> basic analysis or explicit stop state
```

## What this project provides

- A Codex/Claude-compatible skill: `skills/geo-biodata-workflow/`
- R scripts for deterministic steps:
  - environment check;
  - dependency profile planning/checking;
  - GEO metadata/sample/supplement discovery;
  - publication supplement lookup;
  - reviewable download-plan generation;
  - guarded GEO supplementary-file download;
  - manifest validation before analysis;
  - manifest-driven bulk raw-count, bulk-normalized, microarray, and scRNA object-intake drivers;
  - a lightweight GSEA template plus deprecated fail-closed bulk/scRNA compatibility markers.
- Route ontology, dependency profiles, and generic decision rules for repeatable agent adjudication.
- `knowledge/skill_integration_map.yaml`, which records which reusable bio/ngs method patterns were folded into this workflow and which heavier analyses remain explicit handoffs.
- Short reference notes for common decisions:
  - input routing;
  - bulk expression and DE/GSEA rules;
  - single-cell QC, clustering, marker, and annotation principles.
- A consistent output contract for raw inputs, hashes, logs, tables, figures, and workflow status.

## Scope

Included:

- GEO Series (`GSE`) resource inventory.
- GEO/GSM metadata and sample-characteristics extraction.
- GEO supplementary-file index and selective download.
- Publication identifiers and open-resource links where available.
- Manifest-driven readiness checks for bulk RNA-seq, expression matrix, microarray, and scRNA routes.
- Route-specific drivers for validated raw-count matrices, normalized matrices, microarray matrices, and read-only scRNA object intake.
- Deprecated all-in-one bulk/scRNA scripts that fail closed and are not automatic entry points.
- Basic preranked GSEA with `fgsea`.
- GO/KEGG/Reactome/MSigDB enrichment guidance.
- Single-cell QC, clustering, author-object review, and marker diagnostics as route-specific design guidance.

Not included:

- Project-specific biology or private analysis assumptions.
- Automatic final cell-type labels or phenotype claims.
- FASTQ alignment, quantification, or full nf-core execution.
- Heavy batch integration, trajectory inference, cell communication, spatial transcriptomics, or multi-omics modeling.
- Paywalled article scraping.
- A guarantee that every GEO record is machine-readable or analysis-ready.

## Requirements

Minimum R packages for resource discovery:

```r
install.packages(c("httr2", "jsonlite", "digest", "yaml"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("GEOquery", "Biobase"))
```

Profile-based installation (recommended):

```r
# Core: resource discovery and manifest validation
Rscript skills/geo-biodata-workflow/scripts/bootstrap_environment.R --profile core --install

# Bulk: DESeq2/edgeR/limma for all bulk routes
Rscript skills/geo-biodata-workflow/scripts/bootstrap_environment.R --profile bulk --install

# Enrichment: GSEA and pathway analysis
Rscript skills/geo-biodata-workflow/scripts/bootstrap_environment.R --profile enrichment --install

# scRNA: Seurat, SingleCellExperiment, and H5AD readers
Rscript skills/geo-biodata-workflow/scripts/bootstrap_environment.R --profile scrna --install
```

Manual per-package installation:

```r
BiocManager::install(c("limma", "edgeR", "DESeq2", "fgsea", "clusterProfiler"))
install.packages(c("ggplot2", "pheatmap", "Seurat", "Matrix", "patchwork"))
```

H5AD support (optional — needed only for `inspect_scrna_object.R` on .h5ad files):

```r
BiocManager::install("anndataR")    # pure R, recommended
# or
BiocManager::install("zellkonverter")  # requires Python anndata
```

For repeated NCBI/GEO requests, set an email address before running discovery:

```r
Sys.setenv(ENTREZ_EMAIL = "your_email@example.org")
```

or in your shell:

```bash
export ENTREZ_EMAIL="your_email@example.org"
```

On Windows PowerShell:

```powershell
$env:ENTREZ_EMAIL = "your_email@example.org"
```

Dependency profiles are declared in `dependency_profiles.yaml`. To inspect a profile without installing heavy packages:

```powershell
Rscript skills\geo-biodata-workflow\scripts\bootstrap_environment.R --profile core --plan
Rscript skills\geo-biodata-workflow\scripts\bootstrap_environment.R --profile bulk --check
```

## Quick start without installing the skill

From the repository root:

```bash
mkdir -p runs/GSE000000
Rscript skills/geo-biodata-workflow/scripts/check_environment.R runs/GSE000000/environment.tsv
Rscript skills/geo-biodata-workflow/scripts/discover_geo.R GSE000000 runs/GSE000000
```

Then inspect:

```text
runs/GSE000000/summary.md
runs/GSE000000/workflow_status.tsv
runs/GSE000000/resources/series_metadata.tsv
runs/GSE000000/resources/sample_index.tsv
runs/GSE000000/resources/sample_characteristics.tsv
runs/GSE000000/resources/supplement_index.tsv
runs/GSE000000/resources/routing_hint.tsv
runs/GSE000000/resources/routing_evidence.tsv
runs/GSE000000/resources/route_candidates.tsv
runs/GSE000000/resources/analysis_decisions.tsv
runs/GSE000000/resources/publication_links.tsv
runs/GSE000000/resources/sra_links.tsv
runs/GSE000000/workflow_events.tsv
```

Generate a download plan for route-relevant files, review it, then download only reviewed rows:

```powershell
Rscript skills/geo-biodata-workflow/scripts/generate_download_plan.R `
  runs/GSE000000 `
  "counts|matrix|metadata|h5ad|rds|mtx" `
  "candidate processed expression or metadata file"

# Edit runs/GSE000000/plans/download_plan.tsv and set reviewed=TRUE
# only for files confirmed to match the selected route.
Rscript skills/geo-biodata-workflow/scripts/download_geo_supp.R `
  runs/GSE000000/plans/download_plan.tsv `
  runs/GSE000000/raw
```

The plan generator and downloader intentionally refuse unrestricted patterns such as `.*`. The downloader also refuses selected rows that have not been reviewed. This prevents accidental downloads of every supplementary file in a large GEO record.

For slow GEO transfers, the downloader retries failed plan-mode downloads and uses a longer timeout than base R. Override with:

```powershell
$env:GEO_BIODATA_DOWNLOAD_TIMEOUT_SEC = "600"
$env:GEO_BIODATA_DOWNLOAD_RETRIES = "5"
```

Before analysis, create and validate a manifest:

```powershell
Copy-Item templates\run_manifest.example.yaml runs\GSE000000\run_manifest.yaml
# Edit run_manifest.yaml so input files, sample mapping, route, design, and contrast are explicit.
Rscript skills\geo-biodata-workflow/scripts/validate_manifest.R runs/GSE000000/run_manifest.yaml
```

## Typical agent prompt

After installing or exposing the skill to an agent:

```text
Use GEO_biodata to analyze GSE000000.
First inventory the GEO resources, then choose the safest route.
Do not download all supplements blindly.
Stop with a clear status if metadata or input files are insufficient.
```

Expected behavior:

1. Create one accession-specific run directory.
2. Check the R environment.
3. Run GEO/resource discovery before downloading expression data.
4. Review metadata, sample characteristics, supplements, route candidates, analysis decisions, article links, species, assay type, and biological unit.
5. Select one route:
   - `bulk_raw_counts`;
   - `bulk_normalized`;
   - `microarray_series_matrix`;
   - `scrna_raw_counts`;
   - `scrna_author_object`;
   - `metadata_only`;
   - `raw_fastq_handoff`.
6. Generate a download plan, review it, and download only `reviewed=TRUE` route-relevant files.
7. Preserve raw downloads unchanged and record hashes.
8. Create `run_manifest.yaml` and validate it before any statistical analysis.
9. Save status, summary, scripts, tables, figures, logs, and session information.

## Skill installation

### Claude Code

Copy or symlink:

```text
skills/geo-biodata-workflow/
```

into your Claude skill directory, then invoke the skill by name:

```text
/geo-biodata-workflow GSE000000
```

### Codex/OpenAI

Use the included agent metadata:

```text
skills/geo-biodata-workflow/agents/openai.yaml
```

The public project name is `GEO_biodata`. The installable skill id is `geo-biodata-workflow` because current skill systems use lowercase hyphen-case names.

## Repository layout

The public tree is split by role:

- root files are project-level documentation, license, dependency profiles, and CI configuration;
- `skills/` is the installable agent workflow;
- `schemas/`, `templates/`, and `knowledge/` are reusable contracts and decision rules;
- `validation/` contains deterministic fixtures plus smoke and real-data validation notes;
- `containers/` provides optional reproducible runtime definitions.

```text
GEO_biodata/
  .github/
    workflows/
      smoke.yml
  AGENTS.md
  README.md
  dependency_profiles.yaml
  containers/
  schemas/
    route_ontology.yaml
    run_manifest.schema.yaml
  knowledge/
    decision_rules/
      bulk_input_scale.yaml
      scrna_post_count_qc.yaml
    source_registry.yaml
    skill_integration_map.yaml
  templates/
    run_manifest.example.yaml
  skills/
    geo-biodata-workflow/
      SKILL.md
      agents/
        openai.yaml
      references/
        routing-and-resources.md
        basic-bulk.md
        bulk-de-gsea-rules.md
        basic-scrna.md
        scrna-reusable-rules-and-markers.md
      scripts/
        bootstrap_environment.R
        check_environment.R
        discover_geo.R
        collect_publication_supplements.R
        generate_download_plan.R
        download_geo_supp.R
        validate_manifest.R
        drivers/
          run_bulk_counts.R
          run_bulk_normalized.R
          run_microarray.R
          inspect_scrna_object.R
        bulk_limma_common.R
        analyze_bulk_template.R
        run_gsea_template.R
        analyze_scrna_template.R
        marker_utilities.R
  validation/
    fixtures/
    REAL_DATA_VALIDATION_20260728.md
    SMOKE_TEST_20260727.md
    expected_download_refusal.log
    run_smoke_checks.R
```

Generated data and analysis outputs should be written under `runs/` and kept out of Git. `validation/runs/` is reserved for maintainer smoke tests and should not contain committed data.

## Output contract

Each GEO accession should produce a self-contained run directory:

```text
<GSE>/
  resources/
    series_metadata.tsv
    sample_index.tsv
    sample_characteristics.tsv
    supplement_index.tsv
    publication_links.tsv
    publication_supplements.tsv
    routing_hint.tsv
    routing_evidence.tsv
    route_candidates.tsv
    analysis_decisions.tsv
    sra_links.tsv
  plans/
    download_plan.tsv
  raw/
    download_manifest.tsv
  derived/
  tables/
    fallback_events.tsv
  figures/
  logs/
  scripts/
  environment.tsv
  manifest_validation.tsv
  workflow_status.tsv
  workflow_events.tsv
  summary.md
```

Every figure should have a source script and either a source table or a documented analysis object. Figures from this workflow are diagnostic by default, not publication-ready.

## Workflow states

| State | Meaning |
|---|---|
| `DISCOVERY_PARTIAL` | Public-resource discovery ran, but one or more discovery steps failed or returned incomplete metadata |
| `RESOURCE_INVENTORY_COMPLETE` | GEO/public resources were collected; analysis has not started |
| `REVIEW_REQUIRED` | Human review is needed for sample characteristics, SuperSeries/SubSeries relationships, metadata, QC thresholds, labels, or conflicting resources |
| `MANIFEST_VALIDATED` | Input, route, sample mapping, design, contrast, and review gates passed structural validation |
| `MANIFEST_INVALID` | Manifest validation failed; analysis must not run |
| `READY_FOR_BASIC_ANALYSIS` | Inputs and metadata support one basic route |
| `BASIC_ANALYSIS_COMPLETE` | Required tables, figures, logs, and status exist |
| `BLOCKED_INPUT` | No supported processed input is available |
| `BLOCKED_METADATA` | Biological groups or units cannot be reconstructed safely |
| `RAW_COMPUTE_REQUIRED` | Only FASTQ/raw sequencing data are available |

Do not silently turn a blocked state into an analysis result.

## Route selection guide

| Evidence found | Recommended route |
|---|---|
| Integer-like gene-by-sample matrix with clear groups | `bulk_raw_counts` |
| Log-normalized expression, TPM, FPKM, CPM, or similar matrix | `bulk_normalized` |
| ExpressionSet or microarray series matrix | `microarray_series_matrix` |
| MTX/H5 count matrix, count-bearing H5AD, or count-bearing Seurat object | `scrna_raw_counts` |
| Author Seurat/H5AD/RDS object with labels but uncertain raw counts | `scrna_author_object` |
| SRA/FASTQ only | `raw_fastq_handoff` |
| Missing groups, donor identity, or sample mapping | `metadata_only`, `BLOCKED_METADATA`, or EDA-only review |

When a GEO record is a SuperSeries or SubSeries, review related accessions before choosing the analysis unit.

## Basic analysis outputs

The previous all-in-one bulk and scRNA templates are deprecated and fail closed. They remain in the tree as compatibility markers. Use route-specific manifest-driven drivers instead.

### Bulk raw-count route (`scripts/drivers/run_bulk_counts.R`)

- `sample_mapping_used.tsv`
- `design_matrix.tsv`
- `library_sizes.tsv`
- `pca_coordinates.tsv`
- `de_results_<contrast>.tsv`
- `qc_checks.tsv`
- `bulk_library_sizes.pdf`
- `bulk_pca.pdf`
- `bulk_pvalue_histogram_<contrast>.pdf`

### Bulk normalized route (`scripts/drivers/run_bulk_normalized.R`)

Uses limma via shared `bulk_limma_common.R`. For TPM, FPKM, CPM, or log-normalized matrices.

- `sample_mapping_used.tsv`
- `design_matrix.tsv`
- `library_sizes.tsv`
- `pca_coordinates.tsv`
- `sample_correlation.tsv`
- `de_results_<contrast>.tsv`
- `qc_checks.tsv`
- `bulk_library_sizes.pdf`
- `bulk_pca.pdf`
- `sample_correlation_heatmap.pdf`
- `bulk_pvalue_histogram_<contrast>.pdf`
- `bulk_meanvar_<contrast>.pdf`

### Microarray route (`scripts/drivers/run_microarray.R`)

Uses limma via shared `bulk_limma_common.R`. Handles Series Matrix parsing, platform annotation, and probe-to-gene mapping.

- Same outputs as bulk normalized, plus:
- `probe_mapping.tsv`
- `processing_notes.txt`

### scRNA object intake (`scripts/drivers/inspect_scrna_object.R`)

Read-only inventory for Seurat RDS and H5AD objects. Never modifies the original object.

- `inventory.tsv`
- `cell_metadata_fields.tsv`
- `feature_metadata_fields.tsv`
- `inventory_summary.md`

### GSEA route

- `gsea_results_<collection>.tsv`
- `gsea_top_<collection>.pdf`

## Guardrails for agents

- Start with resource discovery, not downloading.
- Prefer author-provided processed matrices, objects, metadata, and labels when they are adequate.
- Do not infer sample groups only from filenames when GEO/GSM metadata or supplements disagree.
- Record whether expression values are raw counts, normalized expression, TPM/FPKM/CPM, microarray intensities, or author-processed object layers.
- Do not pass normalized values into raw-count DE models.
- Define contrast direction explicitly.
- For ORA, use the tested/filtered gene universe, not the whole genome.
- Preserve author cell labels and write harmonized labels separately.
- Do not apply fixed scRNA QC thresholds without reviewing per-sample distributions.
- Do not treat cells as independent biological replicates for donor-level inference.
- Do not remove doublets silently; save scores and rules first.
- Do not scrape paywalled papers.

## Validation

The current public workflow has smoke and real-data validation notes at:

```text
validation/SMOKE_TEST_20260727.md
validation/REAL_DATA_VALIDATION_20260728.md
```

The smoke test verifies:

- skill folder structure;
- R script parseability;
- GEO discovery without expression-data download;
- guarded refusal of unrestricted supplementary download and unreviewed download plans;
- manifest validation success and failure gates;
- dependency profile planning;
- optional bulk raw-count driver execution when DESeq2 is available;
- absence of deprecated dangerous analysis patterns;
- absence of local/private path strings in the public project tree.

It does not validate biological conclusions for any specific GEO dataset.

## Design principle

Keep the workflow simple:

- deterministic operations live in scripts;
- method judgment lives in short references;
- dataset-specific decisions stay in the run directory;
- blocked inputs produce explicit stop states;
- established R/Bioconductor packages do the analysis.
