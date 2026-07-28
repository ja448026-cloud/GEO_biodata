---
name: geo-biodata-scrna-pseudobulk
description: GEO_biodata single-cell pseudobulk differential-expression bridge. Use when scRNA data have raw counts, sample or donor metadata, cell-type labels or reviewed clusters, and the user asks for treatment-vs-control or condition DE within cell types. Do not use cell-level marker tests for population DE.
---

# GEO_biodata scRNA Pseudobulk

Use this for condition DE in scRNA data.

## Preconditions

- Raw counts are available.
- Cells have sample/donor identifiers.
- Cell types or clusters are reviewed.
- Each contrast has enough biological replication at the sample/donor level.

## Route

Aggregate raw counts by `sample_id x cell_type` or `sample_id x reviewed_cluster`, then hand each cell-type count matrix to the bulk DE engine.

Do not run `FindMarkers` or `rank_genes_groups` for treatment-vs-control population claims. Those are marker-ranking tools, not replicate-aware condition DE.

## Dependencies

Use `scrna_pseudobulk` only when this route is needed. It extends scRNA intake plus bulk count dependencies.

