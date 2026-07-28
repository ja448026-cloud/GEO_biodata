# Omics Plot Recipes

Use these recipes after defining the figure question, source, unit, and stage.

## PCA Or MDS

- Use independent samples as points.
- Show condition plus one major technical or biological covariate when available.
- Label explained variance for PCA axes when computed.
- Use ellipses only when n and interpretation are defensible.
- Record any outlier labeling rule.

## MA Plot

- Use average abundance on x and effect size on y.
- Color by adjusted threshold when the upstream result provides one.
- Keep the full point cloud visible.
- Use this as a core DE diagnostic before relying on a volcano summary.

## Volcano Plot

- Use explicit log fold-change direction on x.
- Label whether y is raw p value or adjusted p value.
- Predeclare the label rule: largest effect, strongest statistic, known candidate list, or top rank.
- State any p-value truncation in the caption.

## Heatmap And Color

- Record why rows were selected: DE threshold, variance, curated signature, pathway membership, or predefined gene list.
- Label the value scale: raw, log-normalized, centered, row z-score, rank, percentage, or score.
- Do not compare absolute expression across genes when the displayed value is row z-score.
- Include group, batch, donor, tissue, or paired annotations when they affect interpretation.
- Use ordered, perceptual color maps for continuous values and limited categorical palettes for groups.

## Dot, Box, Violin, And Paired Plots

- Show individual biological observations when feasible.
- Use paired lines for paired tumor-normal, before-after, or repeated-measure designs.
- Label n per group or make it recoverable from the source table.
- State the error-bar type if summaries are displayed.

## Enrichment And Gene-Set Plots

- Prefer directional effect/NES plus FDR and gene-set size.
- Record pathway filtering and redundancy handling.
- Collapse or group highly similar pathways when many terms are significant.
- Use leading-edge displays only when they explain the mechanism.

## scRNA Figures

- Treat UMAP/t-SNE as visual embedding, not quantitative distance evidence.
- Inspect donor/sample composition for clusters before making biological claims.
- Compare cell proportions at sample/donor level.
- Distinguish expression fraction from average expression in marker dot plots.
- Use sample/donor-level pseudobulk for condition DE claims.

## Composition Plots

- Keep denominator explicit: sample total, donor total, tissue region, spot count, or cell count.
- Use sample-level values for group comparisons.
- Avoid stacked bars when too many categories obscure the main message; group rare classes with a stated rule.

## Network Figures

- Define node type, edge meaning, edge source, and filtering rule before layout.
- Reduce nodes to the message-bearing subset.
- Group modules visually when modules are part of the interpretation.
- Do not use dense networks as decorative evidence.

## Overview Schematics

- Make the schematic original and source-linked.
- Show analysis order, input types, output contracts, and decision gates.
- Do not copy a paper diagram, icon set, or proprietary visual system.
