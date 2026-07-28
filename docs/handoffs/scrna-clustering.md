# scRNA Clustering And Marker Handoff

GEO_biodata scRNA intake is read-only. Recluster or annotate only after object inventory confirms usable counts, sample metadata, and author-label provenance.

Use local skills:

- `bio-single-cell-preprocessing`
- `bio-single-cell-clustering`
- `bio-single-cell-markers-annotation`
- `single-cell-atlas` for broader workflow planning

Boundary:

- Cluster in PCA/neighbor graph space, not UMAP coordinates.
- Sweep resolution; choose defensible granularity.
- Treat marker tests as descriptive ranking, not proof of a real cell type.
- Do not use cluster-marker tests for treatment-vs-control claims.
