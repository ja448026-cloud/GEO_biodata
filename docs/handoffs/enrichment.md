# Enrichment Handoff

`geo-biodata-enrichment` currently has one executable path: preranked GSEA with `fgsea` and a local GMT file.

GO, KEGG, Reactome, MSigDB ORA, and online enrichment require explicit gene ID type, mapping provenance, and background universe. GEO_biodata includes local GMT-based GSEA and ORA helpers; use broader pathway skills only when a method outside those helpers is needed.

Minimum future ORA driver contract:

- input gene list.
- background universe.
- species.
- gene ID type.
- mapping table and mapping loss.
- database/source version.
- result table and source-linked figure.
