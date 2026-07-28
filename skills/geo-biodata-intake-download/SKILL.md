---
name: geo-biodata-intake-download
description: GEO_biodata download-only and resource-intake workflow. Use when the user wants GEO metadata, sample table, DOI/PMID/PMCID, publication supplementary links, supplementary-file inventory, download planning, reviewed downloads, hashes, and file-level QC without running differential expression, enrichment, clustering, or marker analysis.
---

# GEO_biodata Intake And Download

This skill is intentionally light. It inventories public resources and downloads only reviewed files.

## Steps

1. Create `runs/<GSE>/` or another accession-specific run directory.
2. Run `../geo-biodata-workflow/scripts/check_environment.R`.
3. Use `../geo-biodata-workflow/scripts/bootstrap_environment.R --profile core --plan` or `--check`; install only if needed.
4. Run `../geo-biodata-workflow/scripts/discover_geo.R <GSE> <run-dir>`.
5. Review `resources/series_metadata.tsv`, `sample_index.tsv`, `sample_characteristics.tsv`, `supplement_index.tsv`, `publication_links.tsv`, `publication_supplements.tsv`, `route_candidates.tsv`, and `routing_evidence.tsv`.
6. Generate a selective plan with `../geo-biodata-workflow/scripts/generate_download_plan.R <run-dir> <regex> [reason]`.
7. Set `reviewed=TRUE` only for route-relevant files, then run `../geo-biodata-workflow/scripts/download_geo_supp.R <plan> <raw-dir>`.

## QC Boundary

Do include file/resource QC:

- URL and transfer status.
- downloaded size and SHA-256 hash.
- compression/readability checks when practical.
- matrix dimensions, row/column names, and sample-ID overlap with GEO metadata.
- whether the file looks like raw counts, normalized expression, microarray intensity, scRNA object, or raw sequencing handoff.

Do not run biological/statistical QC, DE, enrichment, clustering, or marker analysis in this skill.

## Publication Resources

Record DOI, PMID, PMCID, article title, open article links, and open supplementary-material links when discoverable. Do not scrape paywalled content.

