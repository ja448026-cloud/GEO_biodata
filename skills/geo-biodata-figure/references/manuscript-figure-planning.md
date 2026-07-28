# Manuscript Figure Planning

Use this file only when the figure is claim-supporting or manuscript-facing.

## Planning Rules

- Write the claim before selecting panels.
- Include only panels that support the same claim or a necessary diagnostic bridge.
- Keep diagnostic panels separate unless the manuscript needs them to justify interpretation.
- Prefer audited individual panels before building a composite.
- Use panel letters only after the panel set is stable.

## Minimum Manuscript Contract

Each manuscript-facing figure needs:

- `figure_plans/<figure_id>_plan.yaml`
- `tables/figure_sources/<figure_id>_source.tsv`
- `scripts/figures/<figure_id>.R` or `.py`
- `figures/<figure_id>.pdf`
- `figures/<figure_id>.png`
- `figure_qa/<figure_id>_qa.tsv`
- `captions/<figure_id>_caption.md`

## Caption Requirements

Caption drafts must state:

- data source and GEO accession when relevant.
- displayed unit and sample count.
- transformation or scale.
- contrast direction and statistics source.
- filtering or label rule.
- whether the figure is diagnostic, exploratory, or claim-supporting.

## Multi-Panel QA

Before completion, verify:

- all panels answer one figure-level question.
- panel order follows the analysis logic.
- shared legends and color meanings are consistent.
- labels remain readable at final export size.
- no panel implies stronger evidence than the upstream analysis supports.
