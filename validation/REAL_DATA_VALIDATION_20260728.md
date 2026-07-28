# Real Data Validation

Date: 2026-07-28

## GSE214293 Bulk Normalized Validation

Purpose: validate the `bulk_normalized` route on a real gastric cancer GEO Series.

Dataset:

- Accession: `GSE214293`
- GEO title: Construction and validation of a scoring system for perineural invasion-mediated inflammation in gastric cancer
- Assay: expression profiling by high throughput sequencing
- Species: Homo sapiens
- GEO samples: 50
- Selected input: `GSE214293_FPKM_processed_file.csv.gz`
- Selected route: `bulk_normalized`

Route decision:

- The author-provided supplementary file is an FPKM processed matrix.
- Raw downloaded CSV was preserved unchanged.
- Annotation columns were removed in a derived matrix.
- Untransformed FPKM values were transformed to `log2(FPKM + 1)` before limma.
- Only complete tumor/adjacent-normal patient pairs were used for paired DE.

Analysis manifest:

- Route: `bulk_normalized`
- Input type: `fpkm`
- Analysis intent: `differential_expression`
- Biological unit: patient
- Samples used: 46
- Patients used: 23 complete pairs
- Design formula: `~ patient_id + group`
- Contrast: `tumor` vs `adjacent_normal`

Validation result:

- Manifest validator: `MANIFEST_VALIDATED`
- Driver execution: `EXECUTION_COMPLETE`
- Contract state: `VALID`
- Technical QC: `REVIEW_REQUIRED`
- Result signal: `STRONG_SIGNAL`
- Genes in input: 18,605
- Genes tested after filtering: 15,154
- Genes with adjusted P < 0.05: 7,812
- Fallback events: none; empty `fallback_events.tsv` contract table was written.
- QC review reason: PCA outlier `N25`.

Contrast reversal:

- Reverse contrast: `adjacent_normal` vs `tumor`
- Genes compared: 15,154
- Maximum absolute `logFC_main + logFC_reverse`: `1.065814e-14`
- Maximum absolute P-value delta: `2.109424e-15`
- Maximum absolute adjusted P-value delta: `2.997602e-15`
- Result: PASS

Workflow fixes prompted by this validation:

- Plan-mode downloads now use configurable timeout/retry settings and remove partial files before retry.
- GEO supplement discovery now backfills `supplement_url` from GEOquery's `url` column.
- Download-plan generation can recover `supplement_url` from legacy `url` columns.
- limma-based normalized/microarray DE outputs now include `sample_correlation.tsv` and `sample_correlation_heatmap.pdf`.
- limma QC output now records outlier sample IDs, not only outlier counts.
- Drivers now always write `fallback_events.tsv`, including no-fallback runs.

Boundary:

This validation supports beta-level confidence for the manifest-driven `bulk_normalized` route on one real gastric cancer FPKM matrix with reviewed paired sample mapping. It does not validate every normalized matrix format, every cancer type, or the microarray route.

