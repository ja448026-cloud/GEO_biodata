---
name: geo-biodata-enrichment
description: GEO_biodata enrichment analysis. Provides smoke-validated preranked GSEA and alpha ORA with a required gene-ID mapping and universe contract.
---

# GEO_biodata Enrichment

Two deterministic executors are available: preranked GSEA and over-representation analysis (ORA). Treat ORA and gene-ID mapping as alpha routes until validated on the target study's real annotation and gene-set files.

## Preranked GSEA

```powershell
Rscript core\R\enrichment\run_preranked_gsea.R `
  --rank-table runs\GSE000000\tables\de_results_tumor_vs_normal.tsv `
  --gene-column feature_id `
  --rank-column t `
  --gmt gene_sets\hallmark.gmt `
  --out-dir runs\GSE000000 `
  --collection hallmark
```

Rules:
- Use the full ranked gene list when available.
- Prefer limma `t`, DESeq2 `stat`, or another signed statistic.
- Positive ranks mean numerator/up from the originating contrast.
- Deduplicate gene IDs; keep the largest absolute rank.
- Confirm gene-ID space before interpretation; the executor writes `gsea_id_overlap_gate_*.tsv` and stops when rank IDs barely overlap the GMT.
- Record GMT path, hash, size filters, rank column, and session info.
- Pathway redundancy is automatically reduced via leading-edge Jaccard clustering; results in `gsea_pathway_redundancy_*.tsv`.

## Over-Representation Analysis (ORA)

```powershell
Rscript core\R\enrichment\run_ora_enrichment.R `
  --gene-list runs\GSE000000\tables\de_genes_mapped.tsv `
  --gene-column gene_symbol `
  --universe runs\GSE000000\tables\universe.txt `
  --gmt-dir gene_sets\msigdb `
  --out-dir runs\GSE000000 `
  --collections h.all,c2.cp.kegg_medicus,c5.go `
  --fdr-threshold 0.05
```

### ORA prerequisites

Before running ORA, you must:
1. Map probe/feature IDs to gene symbols using `core/R/mapping/gene_id_mapping.R`
2. Define the universe (background): platform-detected, expression-filtered, reference, or custom
3. Save a mapping audit table (`mapping_audit.tsv`)

Example mapping workflow for microarray (GPL annotation):
```r
source("core/R/mapping/gene_id_mapping.R")
gpl <- read_gpl_annotation("GPL6884.annot.gz")
pmap <- build_probe_symbol_map(gpl)
universe <- build_universe_from_gpl(gpl)
writeLines(universe, "universe.txt")
de_mapped <- map_probes_to_symbols(de_results, pmap, method = "max_abs_logfc")
write.table(de_mapped, "de_genes_mapped.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
write_mapping_audit(attr(pmap, "stats"), pmap, de_results, "mapping_audit.tsv")
```

### Rules

- Use Fisher's exact test for ORA (one-tailed, greater).
- Define the universe explicitly; do not use the entire genome unless the platform supports it.
- Fail closed when the universe has fewer than 100 genes or more than 10% of query genes are outside the universe.
- Deduplicate gene symbols with a documented strategy (default: max_abs_logfc).
- Pathway redundancy is automatically reduced via overlap Jaccard clustering.
- Record GMT sources, universe size, FDR threshold, and session info.

## Output Contract

GSEA produces `tables/gsea_results_*.tsv`, `tables/gsea_id_overlap_gate_*.tsv`, optional `tables/gsea_pathway_redundancy_*.tsv`, `figures/gsea_top_*.pdf`, `tables/gsea_manifest_*.tsv`, `gsea_status.tsv`, and `logs/sessionInfo_gsea_*.txt`.

ORA produces `tables/ora_results.tsv`, optional `tables/ora_pathway_redundancy.tsv`, `figures/ora_top_dot.pdf`, `tables/ora_manifest.tsv`, `ora_status.tsv`, and `logs/sessionInfo_ora.txt`.
