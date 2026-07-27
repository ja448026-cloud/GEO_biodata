# Input routing and public resources

## 1. Inventory first

Collect and keep separate:

- GEO series metadata and sample characteristics;
- GEO supplementary-file index;
- author-provided expression matrices or objects;
- author-provided sample/cell annotations;
- DOI, PMID, PMCID, and lawful article or supplement links;
- linked SRA/ENA accessions when processed files are insufficient.

Use GEO metadata plus file signatures. A filename alone is not evidence of assay type, count scale, or sample identity.

Reject non-Series accessions early. This workflow starts from `GSE` records; `GDS`, `GPL`, `GSM`, BioProject, SRA, or FASTQ-only accessions are secondary resources that need a routed handoff.

Check whether the record is a SuperSeries or SubSeries. If related GSE records exist, decide the analysis unit before download: parent SuperSeries for complete inventory, or the specific SubSeries that matches the assay and biological question.

## 2. Choose one route

| Evidence | Route | Outcome |
|---|---|---|
| Integer-like gene-by-sample matrix plus reviewed sample groups | Bulk counts | Basic bulk QC; differential analysis only if design is clear |
| ExpressionSet or normalized matrix from a microarray/series matrix | Normalized bulk or microarray | EDA and limma-compatible analysis |
| MTX/H5 count matrix, count-bearing RDS, or count-bearing H5AD | scRNA counts | Per-sample QC, normalization, PCA, clustering, UMAP, markers |
| Author Seurat/H5AD object with metadata but uncertain raw layer | Author processed object | Reproduce labels/embeddings and inspect resources; do not claim raw-count QC |
| FASTQ/SRA only | Raw sequencing | `RAW_COMPUTE_REQUIRED`; optionally prepare an nf-core/fetchngs handoff |
| Ambiguous sample mapping or missing biological groups | Any | `BLOCKED_METADATA` or EDA only |

Prefer the smallest author-provided processed resource that preserves the quantities needed for the requested analysis.

## 3. Download plan before transfer

Do not pass a broad regular expression directly to a downloader. After routing, generate `plans/download_plan.tsv`, then review each selected row.

Minimum fields to check before setting `reviewed=TRUE`:

- `file_name` and `supplement_url` match the selected route;
- selected files are not redundant raw archives when a smaller processed matrix/object is available;
- expected size is acceptable for the local device and network;
- the selection reason is recorded in route-specific language.

The downloader should transfer only rows where both `selected=TRUE` and `reviewed=TRUE`.

## 4. Minimum provenance checks

- Record whether expression values are raw counts, log-normalized expression, TPM/FPKM/CPM, microarray intensities, or author-processed object layers.
- Record evidence for species, tissue, donor/sample identity, disease/control labels, and batch fields.
- Prefer author metadata and supplementary tables over inferred file-name parsing.
- If multiple processed resources disagree, stop at `REVIEW_REQUIRED` and document the conflict.

## 5. Manifest before analysis

`routing_hint.tsv` and `routing_evidence.tsv` are advisory. They must not directly trigger statistical analysis.

Before analysis, create `run_manifest.yaml` and validate it with `scripts/validate_manifest.R`. The manifest must explicitly record:

- route and input type;
- species and gene ID type when known;
- reviewed sample mapping file;
- biological unit;
- design formula;
- structured contrast;
- confirmed review gates.

If the manifest is not `MANIFEST_VALIDATED`, stop before differential analysis, clustering, or enrichment.

## 6. Resource priority

1. GEO series and GSM characteristics.
2. Supplementary sample or cell metadata tables.
3. Metadata embedded in the author-provided object.
4. Main article Methods and figure legends.
5. Open supplementary article files.
6. Inferred mapping, recorded explicitly with confidence and evidence.

Never replace original labels. Store `author_label`, `harmonized_label`, `label_source`, `confidence`, and `review_note`.

## 7. Recommended public components

- GEO access and supplements: [GEOquery](https://bioconductor.org/packages/GEOquery/)
- Raw accession/FASTQ handoff: [nf-core/fetchngs](https://nf-co.re/fetchngs/latest/)
- Publication identifiers and open resources: [Europe PMC](https://europepmc.org/RestfulWebService)
- Optional accession graph: [ffq](https://github.com/pachterlab/ffq)

These are upstream tools. This workflow coordinates them and records decisions; it does not vendor their code.
