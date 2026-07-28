# scRNA Pseudobulk Handoff

Use pseudobulk for condition DE in scRNA data.

Preconditions:

- raw counts are available.
- sample or donor field is present.
- cell type or reviewed cluster labels are present.
- each contrast has biological replication at sample/donor level.

Future GEO_biodata executors should implement:

- `aggregate_pseudobulk.R`
- `validate_pseudobulk_design.R`
- `write_celltype_manifests.R`

Until then, aggregate raw counts by sample x cell type with local single-cell skills, then pass each cell-type matrix to the bulk count DE route.
