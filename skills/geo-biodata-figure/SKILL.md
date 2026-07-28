---
name: geo-biodata-figure
description: Independent GEO_biodata omics visualization playbook for choosing, designing, generating, and auditing source-linked diagnostic, exploratory, claim-supporting, or manuscript figures from existing GEO_biodata outputs. Use when an agent needs figure guidance without relying on another local figure skill.
---

# GEO_biodata Figure

Mode:

- `skill_mode`: playbook
- `execution_model`: agent-authored
- `deterministic_driver`: false
- `external_skill_dependency`: none
- `maturity`: methodological-guidance

This skill turns an existing analysis result into a defensible figure. It is not a statistical analysis module and does not rerun differential expression, enrichment, clustering, or count generation.

## Required Inputs

Before plotting, identify:

- `figure question`: one sentence stating what the reader must learn.
- `stage`: `diagnostic`, `exploratory`, `claim-supporting`, `manuscript`, or `quick`.
- `source`: table path or object path plus the exact layer/columns used.
- `unit`: biological sample, patient, donor, animal, cell, gene, pathway, or other explicit unit.
- `statistics source`: upstream result table or `none`.
- `output purpose`: audit, internal exploration, claim support, manuscript, or user preview.

If there is no inspectable source table/object, return `FIGURE_BLOCKED` or first create the figure-specific source table.

## Standard Workflow

1. Write `Figure question: ...`.
2. Confirm the figure stage and output purpose.
3. Confirm source path, data scale, biological unit, grouping variables, and any paired/blocking structure.
4. Read `references/figure-decision-matrix.md` when choosing among plot types.
5. Read only the recipe file needed for the chosen plot family.
6. Export a figure-specific source table before plotting. Do not plot directly from a complex object without saving the data actually used.
7. Use currently available R, Python, or user-specified tools to create the simplest valid plot.
8. Run the relevant scientific and visual QA checks.
9. Save the required output contract for the figure stage.
10. Return one completion state: `FIGURE_COMPLETE`, `FIGURE_REVIEW_REQUIRED`, or `FIGURE_BLOCKED`.

## Quick Figure Mode

Use Quick Figure Mode when the user asks to quickly inspect a fresh result or says to draw one plot first.

Quick Mode rules:

- Answer one figure question.
- Read one source table or one object layer.
- Generate one plot by default.
- Save the source table, script, and figure.
- Check unit, direction, sample count, filtering rule, and label readability.
- Do not call the result a manuscript figure.

## Lazy References

Load only what is needed:

- `references/principles-and-qa.md`: core visualization principles, scientific QA, visual QA, and integrity QA.
- `references/figure-decision-matrix.md`: choose a plot type from the data structure and figure question.
- `references/omics-plot-recipes.md`: PCA/MDS, MA, volcano, heatmap, enrichment, scRNA, composition, network, and overview recipes.
- `references/manuscript-figure-planning.md`: multi-panel planning, caption requirements, export expectations, and claim boundaries.

Use `../../knowledge/visualization_source_registry.yaml` to inspect source provenance and license/use notes when adding or revising principles. Summarize sources in your own words.

## Output Contract

Diagnostic or exploratory figures require:

- `tables/figure_sources/<figure_id>_source.tsv`
- `scripts/figures/<figure_id>.R` or `scripts/figures/<figure_id>.py`
- `figures/<figure_id>.pdf` or `figures/<figure_id>.png`

Claim-supporting or manuscript figures require:

- `figure_plans/<figure_id>_plan.yaml`
- `tables/figure_sources/<figure_id>_source.tsv`
- `scripts/figures/<figure_id>.R` or `scripts/figures/<figure_id>.py`
- `figures/<figure_id>.pdf`
- `figures/<figure_id>.png`
- `figure_qa/<figure_id>_qa.tsv`
- `captions/<figure_id>_caption.md`

Quick Mode requires the diagnostic/exploratory three-file contract plus the five quick QA checks recorded in the response or a QA file.

## Templates

Use these optional templates for formal figures:

- `templates/figure_plan.yaml`
- `templates/figure_spec.yaml`
- `templates/figure_qa.tsv`
- `templates/figure_caption.md`

## Hard Boundaries

- Do not change upstream statistics to make a figure look stronger.
- Do not hide excluded samples, failed groups, missing values, or inconvenient directions.
- Do not treat cells, fields of view, spots, or technical replicates as independent patients/donors.
- Do not label a marker plot as condition differential expression.
- Do not make a volcano plot the only evidence for a DE result.
- Do not copy published figures, proprietary templates, icons, or long code/text blocks.
- Do not require another local skill or a fixed plotting executor.
- Do not build unrelated multi-panel figures by default.

## Completion States

Return `FIGURE_COMPLETE` when the source, script, figure, and required stage contract are present and QA has no blocking issue.

Return `FIGURE_REVIEW_REQUIRED` when the figure was created but has a scientific or visual issue that needs user review, such as weak labels, crowded pathways, small n, ambiguous unit, or unresolved batch/paired structure.

Return `FIGURE_BLOCKED` when the source data, biological unit, statistics source, or requested claim is not available or is scientifically inappropriate for plotting.
