---
name: geo-biodata-enrichment
description: GEO_biodata executable preranked GSEA workflow. Use when a bulk or microarray DE table or rank table already exists and the user asks for preranked GSEA with a local GMT gene-set file. GO/KEGG/Reactome/ORA are planned handoffs unless their gene-ID mapping and universe are explicitly provided.
---

# GEO_biodata Enrichment

This module currently provides one deterministic executor: preranked GSEA with `fgsea` and a local GMT file.

## Executor

```powershell
Rscript core\R\enrichment\run_preranked_gsea.R `
  --rank-table runs\GSE000000\tables\de_results_tumor_vs_normal.tsv `
  --gene-column feature_id `
  --rank-column t `
  --gmt gene_sets\hallmark.gmt `
  --out-dir runs\GSE000000 `
  --collection hallmark
```

## Rules

- Use the full ranked gene list when available.
- Prefer limma `t`, DESeq2 `stat`, or another signed statistic.
- Positive ranks mean numerator/up from the originating contrast.
- Deduplicate gene IDs before enrichment; keep the largest absolute rank.
- Record GMT path, hash, size filters, rank column, and session info.

GO/KEGG/Reactome/ORA require explicit ID type, mapping provenance, and background universe. Until a dedicated driver exists, use `docs/handoffs/enrichment.md`.
