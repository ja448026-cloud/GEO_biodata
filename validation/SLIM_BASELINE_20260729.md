# Slim Baseline 2026-07-29

Baseline: `e3054ab499625e2de6c8d5f80d94f8edf73fc740`

PowerShell equivalent of `wc -c validation/run_smoke_checks.R`:

```text
34459
```

PowerShell equivalent of `find knowledge schemas docs/handoffs -type f | sort`:

```text
docs/handoffs/enrichment.md
docs/handoffs/scrna-clustering.md
docs/handoffs/scrna-pseudobulk.md
knowledge/decision_rules/biological_replication.yaml
knowledge/decision_rules/bulk_input_scale.yaml
knowledge/decision_rules/enrichment_universe.yaml
knowledge/decision_rules/gene_id_mapping.yaml
knowledge/decision_rules/paired_design.yaml
knowledge/decision_rules/pseudobulk_requirement.yaml
knowledge/decision_rules/raw_fastq_handoff.yaml
knowledge/decision_rules/scrna_author_object.yaml
knowledge/decision_rules/scrna_object_intake.yaml
knowledge/decision_rules/scrna_post_count_qc.yaml
knowledge/decision_rules/superseries_selection.yaml
knowledge/platform_registry.tsv
knowledge/route_maturity.yaml
knowledge/skill_integration_map.yaml
knowledge/source_registry.yaml
knowledge/visualization_source_registry.yaml
schemas/intake_handoff.schema.yaml
schemas/route_ontology.yaml
schemas/run_manifest.schema.yaml
```

PowerShell equivalent of `ls skills/`:

```text
geo-biodata
geo-biodata-bulk
geo-biodata-enrichment
geo-biodata-figure
geo-biodata-intake
geo-biodata-scrna-intake
geo-biodata-workflow
```

`git log --oneline -3`:

```text
e3054ab Harden GEO route execution gates
7b37266 Refocus GEO intake download helpers
ccaaaef Stabilize enrichment and mapping contracts
```
