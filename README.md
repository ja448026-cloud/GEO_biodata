# geo_biodata_workflow

`geo_biodata_workflow` is a lightweight, agent-oriented workflow for public GEO biodata.

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
  - GEO metadata/sample/supplement discovery;
  - publication supplement lookup;
  - reviewable download-plan generation;
  - guarded GEO supplementary-file download;
  - manifest validation before analysis;
  - a lightweight GSEA template plus deprecated fail-closed bulk/scRNA compatibility markers.
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
- Deprecated all-in-one bulk/scRNA scripts that fail closed until route-specific drivers are implemented.
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

Optional analysis packages:

```r
BiocManager::install(c("limma", "edgeR", "DESeq2", "fgsea", "clusterProfiler"))
install.packages(c("ggplot2", "pheatmap", "Seurat", "Matrix", "patchwork"))
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

Before analysis, create and validate a manifest:

```powershell
Copy-Item templates\run_manifest.example.yaml runs\GSE000000\run_manifest.yaml
# Edit run_manifest.yaml so input files, sample mapping, route, design, and contrast are explicit.
Rscript skills\geo-biodata-workflow/scripts/validate_manifest.R runs/GSE000000/run_manifest.yaml
```

## Typical agent prompt

After installing or exposing the skill to an agent:

```text
Use geo_biodata_workflow to analyze GSE000000.
First inventory the GEO resources, then choose the safest route.
Do not download all supplements blindly.
Stop with a clear status if metadata or input files are insufficient.
```

Expected behavior:

1. Create one accession-specific run directory.
2. Check the R environment.
3. Run GEO/resource discovery before downloading expression data.
4. Review metadata, sample characteristics, supplements, article links, species, assay type, and biological unit.
5. Select one route:
   - bulk counts;
   - normalized/microarray bulk;
   - scRNA count/object;
   - author-processed object;
   - metadata-only;
   - raw-only handoff.
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

The public project name is `geo_biodata_workflow`. The installable skill id is `geo-biodata-workflow` because current skill systems use lowercase hyphen-case names.

## Repository layout

```text
geo_biodata_workflow/
  .github/
    workflows/
      smoke.yml
  AGENTS.md
  README.md
  PUBLIC_MVP_PLAN.md
  schemas/
    run_manifest.schema.yaml
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
        check_environment.R
        discover_geo.R
        collect_publication_supplements.R
        generate_download_plan.R
        download_geo_supp.R
        validate_manifest.R
        analyze_bulk_template.R
        run_gsea_template.R
        analyze_scrna_template.R
        marker_utilities.R
  validation/
    fixtures/
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
    sra_links.tsv
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
| Integer-like gene-by-sample matrix with clear groups | Bulk counts |
| ExpressionSet, microarray series matrix, or log-normalized expression | Normalized bulk or microarray |
| MTX/H5 count matrix, count-bearing H5AD, or count-bearing Seurat object | scRNA count workflow |
| Author Seurat/H5AD object with labels but uncertain raw counts | Author-processed object review |
| SRA/FASTQ only | `RAW_COMPUTE_REQUIRED` handoff |
| Missing groups, donor identity, or sample mapping | `BLOCKED_METADATA` or EDA only |

When a GEO record is a SuperSeries or SubSeries, review related accessions before choosing the analysis unit.

## Basic analysis outputs

The previous all-in-one bulk and scRNA templates are deprecated and fail closed. They remain in the tree as compatibility markers while route-specific manifest-driven drivers are developed. Do not use them as automatic entry points.

Bulk route:

- `sample_mapping_used.tsv`
- `pca_coordinates.tsv`
- `de_results_<contrast>.tsv`
- `bulk_library_sizes.pdf`
- `bulk_pca.pdf`
- `bulk_sample_correlation.pdf`
- `bulk_pvalue_histogram_<contrast>.pdf`
- `bulk_ma_<contrast>.pdf`
- `bulk_volcano_<contrast>.pdf`
- `bulk_top_de_heatmap_<contrast>.pdf`

GSEA route:

- `gsea_results_<collection>.tsv`
- `gsea_top_<collection>.pdf`

scRNA route:

- `cell_qc_metrics.tsv`
- `cluster_resolution_summary.tsv`
- `generic_marker_presence.tsv`
- `cluster_markers.tsv`
- `cluster_composition.tsv`
- `scrna_qc_prefilter_<metric>.pdf`
- `scrna_pca_elbow.pdf`
- `scrna_umap_cluster.pdf`
- `scrna_umap_sample.pdf`
- `scrna_generic_marker_dotplot.pdf`
- `scrna_marker_dotplot.pdf`
- processed Seurat object under `derived/`

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

The current public prototype has a smoke-test note at:

```text
validation/SMOKE_TEST_20260727.md
```

The smoke test verifies:

- skill folder structure;
- R script parseability;
- GEO discovery without expression-data download;
- guarded refusal of unrestricted supplementary download and unreviewed download plans;
- manifest validation success and failure gates;
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
