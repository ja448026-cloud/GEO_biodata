# Smoke Test

Date: 2026-07-27

## Current Verified Boundary

- Skill structure validation passed before the rename.
- Core environment check passed under the local R environment.
- Live discovery test with GSE2553 completed without expression-data download.
- NCBI summary and GEO supplement inventory completed.
- Sample index contained 181 records.
- No series-level supplementary files were reported for that test accession, and the workflow recorded `NO_SUPPLEMENT_FILES`.
- Unfiltered supplementary download request (`.*`) was correctly refused with a non-zero exit status.

## After Rename

The project was renamed to `geo_biodata_workflow`, with installable skill id `geo-biodata-workflow`.

- Skill validation after rename: PASS.
- R syntax parse after rename: PASS for all 8 scripts.
- Local absolute path scan after rename: PASS, no Windows drive-style absolute paths, private project strings, or old skill-id strings found in the public project tree.

## After Review-Driven Simplification

The 2026-07-27 review was applied selectively. The project remains a lightweight workflow/skill, not a full R package.

Accepted changes:

- Added safe `on.exit()` handling for analysis-template log sinks.
- Added GEO Entrez-email warning/header support and SuperSeries/SubSeries routing flags.
- Added bulk p-value histogram output, Okabe-Ito-compatible DE plot colors, robust heatmap clipping, `ward.D2` heatmap clustering, session info, and workflow status output.
- Clarified DE/GSEA rules: input scale, edgeR `filterByExpr`, DESeq2 shrinkage note, ORA universe, GO simplification, KEGG/gene-set provenance, leading-edge counts.
- Updated scRNA template to prefer Leiden clustering when available, fall back loudly when not, record a cluster-resolution summary, document doublet-scoring expectations, and tolerate Seurat v4/v5 assay-access differences.
- Expanded generic marker and QC helper sets while keeping them broad and non-project-specific.

Validation after simplification:

- Skill validation: PASS.
- R syntax parse: PASS for all 8 scripts.
- Public-tree string scan: PASS, no Windows drive-style absolute paths, private project strings, or old skill-id strings found.

## Pre-GitHub Cleanup

- Validation directory was reduced to `SMOKE_TEST_20260727.md` and `expected_download_refusal.log`.
- Developer-only validation drivers and local environment TSV files were removed from the public tree.
- Public user run directory was standardized as `runs/<GSE>/`; `validation/runs/` is reserved for maintainer smoke tests.
- MVP promotion criteria no longer imply that this repository vendors the host skill validator.

## Plan-First Hardening

The 2026-07-27 follow-up hardening addressed the first operational gap found during local use: discovery can be partially successful while sample-characteristic extraction or supplement lookup fails.

Accepted changes:

- `discover_geo.R` now writes `workflow_events.tsv` and can emit `DISCOVERY_PARTIAL`, `REVIEW_REQUIRED`, or `BLOCKED_METADATA` instead of always declaring `RESOURCE_INVENTORY_COMPLETE`.
- `supplement_index.tsv` now includes explicit `supplement_url` and `file_name` fields when GEO reports downloadable supplements.
- Added `generate_download_plan.R`, which creates `plans/download_plan.tsv` and refuses broad patterns such as `.*`.
- `download_geo_supp.R` now supports plan-first downloads and only transfers rows with `selected=TRUE` and `reviewed=TRUE`.
- Added `validation/run_smoke_checks.R` for deterministic local checks: R parseability, broad-regex refusal, unreviewed-plan refusal, invalid-URL refusal, and public-tree path scanning.

Validation after plan-first hardening:

- `validation/run_smoke_checks.R`: PASS, including parse checks for all 9 R scripts.
- Live lightweight discovery smoke with `GSE2553 --no-characteristics`: PASS.
- The lightweight smoke wrote `REVIEW_REQUIRED`, recorded `sample_characteristics=SKIPPED`, and preserved the route hint `bulk_microarray` without expression-data download.

## Manifest Alignment Hardening

The next-round Codex alignment added an analysis-readiness gate without claiming that route-specific bulk or scRNA drivers are complete.

Accepted changes:

- Added `schemas/run_manifest.schema.yaml` and `templates/run_manifest.example.yaml`.
- Added `validate_manifest.R`, which validates route/input compatibility, required review gates, sample mapping fields, design formula fields, contrast levels, duplicate sample IDs, and design-matrix rank when sample metadata are available.
- `discover_geo.R` now writes `resources/routing_evidence.tsv` so route hints are auditable evidence rather than implicit permission to analyze.
- Deprecated `analyze_bulk_template.R` and `analyze_scrna_template.R` now fail closed instead of running risky all-in-one analysis logic.
- Added a minimal GitHub Actions smoke workflow.
- `validation/run_smoke_checks.R` now also validates a fixed manifest fixture, rejects an unconfirmed review gate, and scans for legacy dangerous patterns.

Validation after manifest alignment:

- `validation/run_smoke_checks.R`: PASS, including parse checks for all 10 R scripts and manifest validation gates.

## Boundary

This smoke test validates the environment, skill structure, public-resource discovery, output creation, manifest validation, and download guards. It does not claim biological validity for bulk DE, GSEA, or scRNA clustering. The old all-in-one bulk and scRNA templates are deprecated and fail closed until route-specific manifest-driven drivers are implemented.
