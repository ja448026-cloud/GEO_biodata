# Principles And QA

Use these rules for GEO_biodata diagnostic, exploratory, claim-supporting, and manuscript figures.

## Core Principles

- Start from the figure question, not a preferred plot type.
- Keep one main message per figure unless a multi-panel layout is needed to support one coherent claim.
- Preserve the correct biological unit. Patients, donors, animals, and independent samples are not interchangeable with cells, spots, fields, or technical replicates.
- Prefer effect size, direction, confidence, and uncertainty over decorative significance markers.
- Make filtering rules explicit for selected genes, pathways, labels, heatmap rows, cells, or samples.
- Treat color as data encoding. Use distinguishable palettes and redundant encodings when group identity matters.
- Separate diagnostic figures from claim-supporting or manuscript figures.
- Keep figures traceable to source data, code, parameters, and output files.

## Scientific QA

Check before completion:

- Does every point, bar, line, cell, or tile represent the intended unit?
- Is the contrast direction identical to the upstream result table?
- Are p values, FDR values, ranks, or effect sizes read from the correct source?
- Are paired, blocked, repeated-measure, donor, or batch structures visible or documented?
- Are excluded samples, missing values, and failed groups recorded?
- Are normalization, log scale, z-score, percentage, or rank transformations labeled?
- Does the figure avoid causal language unless the upstream design supports it?

## Visual QA

Check before completion:

- Can the main message be understood quickly at the intended size?
- Are axes, legends, group labels, units, and titles clear?
- Are labels readable without overlap?
- Are continuous color scales perceptually ordered and not rainbow-like?
- Are categorical colors distinguishable in color-blind and grayscale contexts when possible?
- Are axis ranges honest and not cropped to exaggerate effects?
- Are multi-panel figures aligned around one claim rather than a collection of unrelated outputs?

## Integrity QA

For diagnostic/exploratory figures, confirm the source table, script, and figure exist.

For claim-supporting/manuscript figures, confirm the plan, source table, script, PDF, PNG, QA file, and caption draft exist.

Record `FIGURE_REVIEW_REQUIRED` if any scientific or visual QA item cannot be resolved without user judgment.
