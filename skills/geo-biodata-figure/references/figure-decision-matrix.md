# Figure Decision Matrix

Choose the simplest plot that answers the figure question while preserving the correct unit.

| Figure question | Preferred plot | Avoid |
|---|---|---|
| Do samples separate by condition, batch, tissue, or donor? | PCA or MDS scatter with shape/color for key covariates and explained variance when available | Treating separation as proof of causality |
| Is there a global DE signal and stable model behavior? | MA plot, p-value histogram, mean-variance plot, sample-correlation heatmap | Volcano as the only diagnostic |
| Which genes have large effects and strong statistical support? | Volcano or ranked effect plot with predeclared label rule | Hand-picked labels that only support the story |
| How does one gene or score vary across independent samples? | Dot, box, violin, or paired plot; show individual observations | Bar plot hiding n and paired structure |
| How do multiple genes, samples, or pathways form a matrix pattern? | Heatmap with explicit row selection and scale | Unexplained top rows or confusing row z-score with absolute expression |
| Which enriched pathways are strongest and directional? | NES/effect dot or bar plot with FDR and gene-set size | Long unfiltered pathway list |
| How do cell-type proportions differ? | Sample-level composition plot, dot/box summaries by group | Pooling all cells as one group-level value |
| How does a marker pattern support annotation? | Dot plot or heatmap showing expression fraction and average expression | Treating marker plot as condition DE |
| What relationship exists between two continuous variables? | Scatter with method-labeled correlation or model fit when justified | Dual y-axis or unexplained smoothing |
| What network or mechanism should be communicated? | Reduced network with defined node/edge meaning and modules | Dense all-edge network with no message |
| What is the overall workflow or design? | Original overview schematic with explicit inputs, outputs, and analysis order | Copying a published figure style or diagram |

If no row fits, state the figure question and propose a new plot type with source, unit, and QA implications.
