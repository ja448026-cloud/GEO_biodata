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

## Boundary

This smoke test validates the environment, skill structure, public-resource discovery, output creation, and download guard. It does not claim biological validity for bulk DE, GSEA, or scRNA clustering. Those templates remain lightweight workflow aids that require accession-specific metadata review.
