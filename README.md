# GEO_biodata

`GEO_biodata` is an agent-oriented workflow for public GEO biodata.

It starts from a GEO Series accession, inventories public resources, builds reviewed download plans, validates an explicit manifest, and runs only the analysis route supported by the available data.

It is not a full R package and not a FASTQ processing framework. Raw data, downloaded supplements, sample mapping, design choices, logs, and status files are kept traceable.

## Choose The Narrow Skill

Use the smallest skill that matches the current task:

| Need | Skill |
|---|---|
| Decide which module to use | `geo-biodata-router` |
| GEO metadata, DOI/PMID/PMCID, publication supplements, file inventory, reviewed download, file-level QC | `geo-biodata-intake-download` |
| Bulk raw-count, bulk-normalized, or microarray DE/QC from a validated manifest | `geo-biodata-bulk-de` |
| GO/KEGG/Reactome/MSigDB/ORA/GSEA from an existing DE table or rank vector | `geo-biodata-bulk-enrichment` |
| Read-only scRNA object inventory for Seurat/H5AD/RDS/MTX | `geo-biodata-scrna-intake` |
| scRNA reclustering, resolution review, marker ranking, or label sanity checks | `geo-biodata-scrna-cluster-marker` |
| scRNA condition DE within reviewed cell types using sample/donor replicates | `geo-biodata-scrna-pseudobulk` |
| One simple source-table-linked diagnostic plot | `geo-biodata-figure` |
| Repo structure, CI, fixtures, release checks | `geo-biodata-maintainer` |
| Legacy prompts naming the old entry point | `geo-biodata-workflow` |

The original script bundle remains under `skills/geo-biodata-workflow/scripts/` for compatibility. New skill entries reference that shared implementation instead of duplicating code.

## Fast Download-Only Start

From the repository root:

```powershell
New-Item -ItemType Directory -Force runs\GSE000000 | Out-Null
Rscript skills\geo-biodata-workflow\scripts\check_environment.R runs\GSE000000\environment.tsv
Rscript skills\geo-biodata-workflow\scripts\bootstrap_environment.R --profile core --check
Rscript skills\geo-biodata-workflow\scripts\discover_geo.R GSE000000 runs\GSE000000
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
Rscript skills\geo-biodata-workflow\scripts\generate_download_plan.R `
  runs\GSE000000 `
  "counts|matrix|expression|fpkm|tpm|h5ad|rds|mtx" `
  "candidate processed expression or metadata file"

# Edit runs\GSE000000\plans\download_plan.tsv.
# Set reviewed=TRUE only for selected rows that match the intended route.
Rscript skills\geo-biodata-workflow\scripts\download_geo_supp.R `
  runs\GSE000000\plans\download_plan.tsv `
  runs\GSE000000\raw
```

Download-only QC includes transfer status, hashes, file size, readable compression/table checks when practical, sample-ID overlap, and route hints. It does not run differential expression, enrichment, clustering, or marker analysis.

For slow GEO transfers:

```powershell
$env:GEO_BIODATA_DOWNLOAD_TIMEOUT_SEC = "600"
$env:GEO_BIODATA_DOWNLOAD_RETRIES = "5"
```

## Analysis Start

Before analysis, create and validate a manifest:

```powershell
Copy-Item templates\run_manifest.example.yaml runs\GSE000000\run_manifest.yaml
Rscript skills\geo-biodata-workflow\scripts\validate_manifest.R runs\GSE000000\run_manifest.yaml
```

Then select exactly one route:

| Manifest route | Driver |
|---|---|
| `bulk_raw_counts` | `skills\geo-biodata-workflow\scripts\drivers\run_bulk_counts.R` |
| `bulk_normalized` | `skills\geo-biodata-workflow\scripts\drivers\run_bulk_normalized.R` |
| `microarray_series_matrix` | `skills\geo-biodata-workflow\scripts\drivers\run_microarray.R` |
| `scrna_author_object` | `skills\geo-biodata-workflow\scripts\drivers\inspect_scrna_object.R` |

Use `geo-biodata-bulk-enrichment` only after a DE or ranked table exists.

Use `geo-biodata-scrna-cluster-marker` only after scRNA intake confirms that reclustering or marker review is appropriate. Cluster markers are descriptive; condition DE in scRNA should use sample/donor-level pseudobulk.

## Dependency Profiles

Profiles are declared in `dependency_profiles.yaml`.

```powershell
Rscript skills\geo-biodata-workflow\scripts\bootstrap_environment.R --profile core --plan
Rscript skills\geo-biodata-workflow\scripts\bootstrap_environment.R --profile bulk_limma --plan
Rscript skills\geo-biodata-workflow\scripts\bootstrap_environment.R --profile bulk_counts --plan
Rscript skills\geo-biodata-workflow\scripts\bootstrap_environment.R --profile enrichment --plan
Rscript skills\geo-biodata-workflow\scripts\bootstrap_environment.R --profile scrna_intake --plan
Rscript skills\geo-biodata-workflow\scripts\bootstrap_environment.R --profile scrna --plan
Rscript skills\geo-biodata-workflow\scripts\bootstrap_environment.R --profile scrna_pseudobulk --plan
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
