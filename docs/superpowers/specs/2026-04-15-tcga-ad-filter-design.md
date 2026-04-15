# TCGA AD Filter Stage Design

Date: 2026-04-15

## Goal

Insert a new workflow stage after `workflow/02-Select_Preclinical.R` that filters `selected_genes` by cross-system distribution concordance against matched TCGA or TARGET expression data.

The new stage follows first principles:

- Only test genes that already passed the existing preclinical significance and intersection logic.
- Compare underlying expression distributions directly instead of relying only on meta-analysis summary statistics.
- Standardize each cohort independently with `z-score` before testing so the test reflects shape and relative distribution, not absolute scale.
- Keep genes whose preclinical and TCGA distributions are not significantly different by Anderson-Darling.

## Stage Layout

The workflow becomes:

- `02-Select_Preclinical.R`
  - Keep current behavior.
  - Output `selected_genes.csv`.
- `03-TCGA_AD_Filter.R`
  - New stage.
  - Read `selected_genes.csv`.
  - Map `tumor_type` to TCGA or TARGET cohort.
  - Pull raw expression vectors for each selected gene from `cell`, `pdcpdx`, and matched TCGA/TARGET data.
  - `z-score` each cohort-specific vector before comparison.
  - Run Anderson-Darling tests for `cell vs TCGA` and `pdcpdx vs TCGA`.
  - Output full statistics and filtered genes.
- `04-Clinical_Validation.R`
  - Renamed from current `03-Clinical_Validation.R`.
  - Read AD-filtered genes from stage 03 instead of reading `selected_genes.csv` directly.

This keeps each stage focused on one decision:

- Stage 02: significance and overlap
- Stage 03: preclinical-to-TCGA distribution concordance
- Stage 04: clinical validation

## Input Data

### Preclinical genes

The new stage consumes:

- `workflow/Output/<drug>/<tumor_slug>/selected_genes.csv`

Expected key field:

- `name`: gene symbol used throughout downstream filtering

### TCGA/TARGET cohort mapping

Use the following mapping from DROMA tumor type to cohort identifier:

```r
droma_to_tcga_tumor_type <- c(
  "haematopoietic/lymphoid cancer" = "TCGA-LAML",
  "nervous system cancer" = "TCGA-GBM",
  "sarcoma" = "TCGA-SARC",
  "lung cancer" = "TCGA-LUAD",
  "breast cancer" = "TCGA-BRCA",
  "prostate cancer" = "TCGA-PRAD",
  "stomach cancer" = "TCGA-STAD",
  "bladder cancer" = "TCGA-BLCA",
  "skin cancer" = "TCGA-SKCM",
  "ovarian cancer" = "TCGA-OV",
  "kidney cancer" = "TCGA-KIRC",
  "thyroid cancer" = "TCGA-THCA",
  "aerodigestive tract cancer" = "TCGA-HNSC",
  "vulvar cancer" = NA,
  "endometrial cancer" = "TCGA-UCEC",
  "non-cancer" = NA,
  "colon cancer" = "TCGA-COAD",
  "pancreatic cancer" = "TCGA-PAAD",
  "cervical cancer" = "TCGA-CESC",
  "liver cancer" = "TCGA-LIHC",
  "choriocarcinoma" = NA,
  "uterine cancer" = "TCGA-UCS",
  "testicular cancer" = "TCGA-TGCT",
  "nasopharyngeal cancer" = NA,
  "retinoblastoma" = "TARGET-RT"
)
```

If a tumor type maps to `NA`, the stage should fail fast with a clear message that no supported TCGA/TARGET cohort exists for AD filtering.

### TCGA/TARGET expression resources

Use:

- `tcga_rna_counts_dir <- "/Users/peng/Library/CloudStorage/OneDrive-Personal/28PHD_peng/250301-DROMA_project/archive260314/251112-DROMA_align/benchmark_mini/Input/TCGA/"`
- `gene_probe_map_path <- "/Users/peng/Desktop/Project/DROMA/Data/gencode.human.v49.annotation.gene.probeMap"`

The implementation must resolve a selected gene symbol to the corresponding TCGA expression row through the probe map before extracting the cohort expression vector.

## Statistical Rule

For each selected gene:

1. Extract raw expression values from:
   - matched `cell` samples
   - matched `pdcpdx` samples
   - matched TCGA/TARGET cohort
2. Remove missing and non-finite values.
3. Apply `z-score` independently within each vector.
4. Run Anderson-Darling k-sample test for:
   - `cell_z` vs `tcga_z`
   - `pdcpdx_z` vs `tcga_z`
5. Define concordance as:
   - `p >= ad_p_t`

This means the retained null is "the two samples are not detectably different in distribution" at the configured threshold.

## Retention Rule

Use a strict retention rule:

- Keep a gene only if both comparisons pass:
  - `cell_vs_tcga_concordant == TRUE`
  - `pdcpdx_vs_tcga_concordant == TRUE`

Genes with insufficient data for either comparison are not retained.

This is stricter than allowing one-sided support, but it matches the intent of requiring translation consistency across both preclinical systems before clinical validation.

## Script And Function Boundaries

### New workflow script

Add:

- `workflow/03-TCGA_AD_Filter.R`

Responsibilities:

- load `selected_genes.csv`
- resolve matched TCGA/TARGET cohort from `tumor_type`
- load the required preclinical and TCGA/TARGET expression data
- call reusable helper functions
- write stage outputs
- print concise stage summary

Rename:

- current `workflow/03-Clinical_Validation.R` -> `workflow/04-Clinical_Validation.R`

Update the renamed script to read:

- `selected_genes_ad_filtered.csv`

instead of:

- `selected_genes.csv`

### New helper file

Add a new helper under `R/FuncXX`, with the exact subdirectory chosen to match existing project organization.

Recommended contents:

- tumor-type to cohort mapping helper
- TCGA expression loading helper
- per-gene vector extraction helper
- safe `z-score` helper
- safe Anderson-Darling wrapper
- batch runner over `selected_genes`
- output shaping helper

The workflow script should remain orchestration-only. Statistical logic and data extraction logic should live in the helper file.

## Proposed Function Responsibilities

Recommended function set:

- `.get_tcga_cohort_for_tumor_type(tumor_type)`
  - returns a single cohort identifier or throws a clear error
- `.safe_zscore(x)`
  - returns standardized numeric vector
  - returns empty vector when input has fewer than 2 usable values or zero variance
- `.run_ad_test(x, y)`
  - performs Anderson-Darling comparison
  - returns statistic, p-value, pass flag, sample counts, and status message
- `.extract_tcga_gene_vector(gene, cohort_id, tcga_rna_counts_dir, gene_probe_map_path)`
  - returns the matched TCGA/TARGET expression vector for one gene
- `.extract_preclinical_gene_vectors(gene, drug, tumor_type, ...)`
  - returns `cell` and `pdcpdx` raw expression vectors for one gene from the corresponding preclinical data source
- `run_tcga_ad_filter(selected_genes, drug, tumor_type, ...)`
  - loops over genes
  - computes pairwise AD results
  - returns both a stats table and a filtered table

Names can change to fit repository style, but the separation of responsibilities should remain.

## Output Files

Stage 03 should write at least:

- `selected_genes_ad_stats.csv`
  - one row per selected gene
  - includes test statistics, p-values, sample sizes, pass flags, and status columns
- `selected_genes_ad_filtered.csv`
  - subset of genes that pass the strict AD concordance rule

Recommended columns in `selected_genes_ad_stats.csv`:

- inherited gene identity columns from `selected_genes`
- `tcga_cohort`
- `cell_n`
- `pdcpdx_n`
- `tcga_n`
- `cell_vs_tcga_ad_stat`
- `cell_vs_tcga_ad_p`
- `cell_vs_tcga_concordant`
- `pdcpdx_vs_tcga_ad_stat`
- `pdcpdx_vs_tcga_ad_p`
- `pdcpdx_vs_tcga_concordant`
- `ad_concordant`
- `ad_status`

`ad_status` should explain failures such as:

- no cohort mapping
- gene not found in probe map
- gene not found in TCGA matrix
- insufficient cell data
- insufficient pdcpdx data
- insufficient TCGA data
- zero variance after cleaning

## Error Handling

Handle expected failure modes explicitly:

- Unsupported tumor type mapping: stop with a clear error before running the stage.
- Missing or unreadable TCGA directory or probe map: stop with a clear error.
- Missing gene mapping or missing expression row: record per-gene failure in stats output, do not crash the whole stage.
- Too few usable values after filtering: mark the gene as failed for that comparison.
- Constant vectors after filtering: mark as failed because `z-score` is undefined or uninformative.

The stage should fail only for configuration or input resource problems. Per-gene data issues should be recorded and skipped.

## Testing Strategy

Validation should cover:

- tumor type correctly maps to the expected TCGA/TARGET cohort
- unsupported tumor types fail fast
- `z-score` helper rejects degenerate vectors correctly
- AD wrapper returns stable columns for success and failure cases
- genes with one failed comparison do not pass final filtering
- renamed stage 04 consumes the AD-filtered output without changing its own statistical thresholds

At minimum, a smoke test should run the new stage for:

- `drug <- "Paclitaxel"`
- `tumor_type <- "breast cancer"`

and verify:

- stage 03 writes both output files
- stage 04 can read stage 03 filtered output

## Non-Goals

This change does not:

- modify the significance thresholds in stage 02
- alter the clinical validation method in stage 04
- relax the requirement to use matched TCGA/TARGET cohorts
- introduce one-sided AD retention logic

## Recommended Implementation Order

1. Add reusable helper file under `R/FuncXX`.
2. Add `workflow/03-TCGA_AD_Filter.R`.
3. Rename `workflow/03-Clinical_Validation.R` to `workflow/04-Clinical_Validation.R`.
4. Update stage 04 input to `selected_genes_ad_filtered.csv`.
5. Run a smoke test on the Paclitaxel breast cancer path.

