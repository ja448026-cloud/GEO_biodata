# GEO_biodata

`GEO_biodata` is an agent-oriented workflow for public GEO biodata.

It starts from a GEO Series accession, inventories public resources, builds reviewed download plans, validates an explicit manifest, and runs only the analysis route supported by the available data.

It is not a full R package and not a FASTQ processing framework. Raw data, downloaded supplements, sample mapping, design choices, logs, and status files are kept traceable.

## Choose The Narrow Skill

Use the smallest skill that matches the current task:

| Need | Skill |
|---|---|
| Decide which module to use | `geo-biodata` |
| GEO metadata, DOI/PMID/PMCID, publication supplements, file inventory, reviewed download, file-level QC, manifest draft | `geo-biodata-intake` |
| Bulk raw-count, bulk-normalized, or microarray DE/QC from a validated manifest | `geo-biodata-bulk` |
| Preranked GSEA from an existing DE/rank table and local GMT file | `geo-biodata-enrichment` |
| Read-only scRNA object inventory for Seurat/H5AD/RDS/MTX | `geo-biodata-scrna-intake` |
| Independent source-table-linked figure planning, generation guidance, and QA | `geo-biodata-figure` |
| Legacy prompts naming the old entry point | `geo-biodata-workflow` |

Stable executor paths live under `core/R/`. The original script bundle remains under `skills/geo-biodata-workflow/scripts/` as a compatibility implementation while the public entry points move to the slimmer module layout.

## Fast Download-Only Start

From the repository root:

```powershell
New-Item -ItemType Directory -Force runs\GSE000000 | Out-Null
Rscript core\R\check_environment.R runs\GSE000000\environment.tsv
Rscript core\R\bootstrap_environment.R --profile intake --check
Rscript core\R\discover_geo.R GSE000000 runs\GSE000000
```

Review:

```text
runs/GSE000000/resources/series_metadata.tsv
runs/GSE000000/resources/sample_index.tsv
runs/GSE000000/resources/sample_characteristics.tsv
runs/GSE000000/resources/supplement_index.tsv
runs/GSE000000/resources/publication_links.tsv
runs/GSE000000/resources/publication_supplements.tsv
runs/GSE000000/resources/route_candidates.tsv
runs/GSE000000/resources/routing_evidence.tsv
```

Generate and review a selective download plan:

```powershell
Rscript core\R\generate_download_plan.R `
  runs\GSE000000 `
  "counts|matrix|expression|fpkm|tpm|h5ad|rds|mtx" `
  "candidate processed expression or metadata file"

# Edit runs\GSE000000\plans\download_plan.tsv.
# Set reviewed=TRUE only for selected rows that match the intended route.
Rscript core\R\download_reviewed_files.R `
  runs\GSE000000\plans\download_plan.tsv `
  runs\GSE000000\raw

Rscript core\R\inspect_downloaded_input.R runs\GSE000000
Rscript core\R\generate_manifest_draft.R runs\GSE000000\intake_handoff.yaml
```

Download-only QC includes transfer status, hashes, file size, readable compression/table checks when practical, sample-ID overlap, and route hints. It does not run differential expression, enrichment, clustering, or marker analysis.

For slow GEO transfers:

```powershell
$env:GEO_BIODATA_DOWNLOAD_TIMEOUT_SEC = "600"
$env:GEO_BIODATA_DOWNLOAD_RETRIES = "5"
```

### Download fallback strategy

When `utils::download.file` exhausts all retries, the downloader falls back to a backup method:

1. **Primary**: R's built-in `utils::download.file` (libcurl backend) with configurable retries and exponential backoff
2. **Backup (curl)**: If `curl` is on PATH and `GEO_BIODATA_DOWNLOAD_BACKUP=curl` (default), the downloader retries via the `curl` CLI with its own retry and timeout
3. **Last resort (GEOquery)**: For NCBI FTP/GEO URLs specifically, `GEOquery::getGEOSuppFiles` is tried as a final alternative

Configure backup behavior:

```powershell
# Disable backup (primary-only)
$env:GEO_BIODATA_DOWNLOAD_BACKUP = "none"

# Increase backup timeout (default 1200s = 20 min)
$env:GEO_BIODATA_BACKUP_TIMEOUT_SEC = "3600"
```

Linux/macOS shells use the same `Rscript` commands with `/` paths, for example:

```bash
mkdir -p runs/GSE000000
Rscript core/R/check_environment.R runs/GSE000000/environment.tsv
Rscript core/R/bootstrap_environment.R --profile intake --check
Rscript core/R/discover_geo.R GSE000000 runs/GSE000000
```

## Analysis Start

Before analysis, create and validate a manifest:

```powershell
Copy-Item templates\run_manifest.example.yaml runs\GSE000000\run_manifest.yaml
Rscript core\R\validate_manifest.R runs\GSE000000\run_manifest.yaml
```

Then select exactly one route:

| Manifest route | Driver |
|---|---|
| `bulk_raw_counts` | `core\R\drivers\run_bulk_counts.R` |
| `bulk_normalized` | `core\R\drivers\run_bulk_normalized.R` |
| `microarray_series_matrix` | `core\R\drivers\run_microarray.R` |
| `scrna_author_object` | `core\R\scrna\inspect_object.R` |

Use `geo-biodata-enrichment` only after a DE or ranked table exists. The current executable enrichment driver is preranked GSEA with a local GMT file; GO/KEGG/Reactome/ORA remain guidance handoffs until gene-ID mapping and background-universe contracts are supplied.

Use `geo-biodata-scrna-intake` for object inventory only. If reclustering, marker review, annotation, or condition DE is needed, use the handoff notes in `docs/handoffs/` with dedicated local single-cell skills. Cluster markers are descriptive; condition DE in scRNA should use sample/donor-level pseudobulk.

Use `geo-biodata-figure` after source tables, DE/rank tables, route-driver diagnostics, or inspected objects exist. It is an independent omics visualization playbook, not an external skill router and not a statistical executor. Agents may use currently available R, Python, or user-specified plotting tools, but figures must stay source-linked and must not change upstream statistics.

## Dependency Profiles

Profiles are declared in `dependency_profiles.yaml`.

```powershell
Rscript core\R\bootstrap_environment.R --profile manifest --plan
Rscript core\R\bootstrap_environment.R --profile intake --plan
Rscript core\R\bootstrap_environment.R --profile bulk_limma --plan
Rscript core\R\bootstrap_environment.R --profile bulk_counts --plan
Rscript core\R\bootstrap_environment.R --profile enrichment_gsea --plan
Rscript core\R\bootstrap_environment.R --profile scrna_intake --plan
```

Default use should be `--plan` or `--check`. Install heavy profiles only when that route is actually needed.

## Output Contract

Typical run directory:

```text
runs/GSE000000/
  resources/
  plans/
  raw/
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

Preserve downloaded files unchanged. Analyze derived copies or explicit local paths recorded in `run_manifest.yaml`.

Formal figures may add `tables/figure_sources/`, `scripts/figures/`, `figure_plans/`, `figure_qa/`, and `captions/`. Diagnostic and exploratory figures require at minimum a source table, plotting script, and figure file; manuscript-facing figures require the full plan/source/script/PDF/PNG/QA/caption contract.

## Validation

Quick local check:

```powershell
Rscript validation\run_smoke_checks.R --quick
```

Full local smoke, including optional driver checks when dependencies are installed:

```powershell
Rscript validation\run_smoke_checks.R
```

GitHub Actions default smoke is intentionally light. Heavy bulk, scRNA, enrichment, and real-data validation should be run manually or by route-specific workflows.

## Current Real-Data Validation

`validation/REAL_DATA_VALIDATION_20260728.md` records a real bulk-normalized gastric cancer validation on `GSE214293`, including selective download, manifest validation, paired limma DE, QC status, and reverse-contrast consistency.
