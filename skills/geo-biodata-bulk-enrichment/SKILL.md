---
name: geo-biodata-bulk-enrichment
description: GEO_biodata downstream enrichment workflow. Use when bulk or microarray DE results already exist and the user asks for GO, KEGG, Reactome, MSigDB, ORA, preranked GSEA, pathway tables, enrichment plots, or ranked-list interpretation.
---

# GEO_biodata Bulk Enrichment

Use this only after a DE table or signed ranked vector exists.

## Preconditions

- Species is known.
- Gene ID type is known or mapped.
- The DE table has an effect-size/statistic column and adjusted p-values.
- The background universe is recorded for ORA.
- The rank metric is recorded for GSEA.

## Methods

- Use ORA for a filtered gene list with a clear background universe.
- Use preranked GSEA when most genes have a signed statistic.
- Prefer `fgsea` for local GMT/ranked-vector workflows.
- Use `clusterProfiler`/Reactome/MSigDB only when the required annotation packages are installed.

## Boundary

Do not rerun DE in this skill. Do not use enrichment to rescue an underpowered or invalid DE design.

Write source tables and one simple plot per analysis type, such as a GO dotplot or GSEA NES barplot. Use `geo-biodata-figure` for plot polishing.

