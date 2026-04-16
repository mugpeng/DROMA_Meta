# ============================================================================
# 00-Eligible_Drug_Tumor.R
# Export valid drugs and tumor types that satisfy project-count thresholds
# ============================================================================

library(data.table)
suppressPackageStartupMessages(library(DROMA.Meta))
library(DROMA.Set)
library(DROMA.R)

db_path <- "/Users/peng/Desktop/Project/DROMA/Data/droma.sqlite"
project_root <- normalizePath(file.path(getwd(), "..", "Meta_Example"), mustWork = TRUE)
# project_root <- file.path(normalizePath(getwd(), mustWork = TRUE), "Meta_Example")
defaults <- getMetaWorkflowDefaults(project_root = project_root)
output_base <- defaults$output_base
dir.create(output_base, showWarnings = FALSE, recursive = TRUE)

cell_n_datasets_t <- 3
pdcpdx_n_datasets_t <- 1

cat("\n=== 00: Eligible Drug And Tumor Type ===\n")

con <- connectDROMADatabase(db_path)

project_anno <- DROMA.Set::listDROMAProjects()
drug_anno <- getDROMAAnnotation("drug")
sample_anno <- getDROMAAnnotation("sample")

valid_inputs <- getValidDrugsAndTumorTypes(
  project_anno = project_anno,
  drug_anno = drug_anno,
  sample_anno = sample_anno,
  cell_n_datasets_t = cell_n_datasets_t,
  pdcpdx_n_datasets_t = pdcpdx_n_datasets_t
)

valid_drugs <- valid_inputs$valid_drugs
valid_tumor_types <- valid_inputs$valid_tumor_types

candidate_pairs <- data.table::CJ(
  drug = valid_drugs,
  tumor_type = valid_tumor_types,
  unique = TRUE
)

cell_eligible_pairs <- candidate_pairs[
  vapply(
    seq_len(.N),
    function(i) {
      drug_i <- candidate_pairs$drug[[i]]
      tumor_type_i <- candidate_pairs$tumor_type[[i]]

      tryCatch(
        {
          filterProjectsForDrugTumor(
            project_anno = project_anno,
            drug_anno = drug_anno,
            sample_anno = sample_anno,
            dataset_types = c("CellLine", "PDC"),
            drug = drug_i,
            tumor_type = tumor_type_i,
            min_project_count = cell_n_datasets_t
          )
          TRUE
        },
        error = function(e) FALSE
      )
    },
    logical(1)
  )
]

eligible_pairs <- cell_eligible_pairs[
  vapply(
    seq_len(.N),
    function(i) {
      drug_i <- cell_eligible_pairs$drug[[i]]
      tumor_type_i <- cell_eligible_pairs$tumor_type[[i]]

      tryCatch(
        {
          filterProjectsForDrugTumor(
            project_anno = project_anno,
            drug_anno = drug_anno,
            sample_anno = sample_anno,
            dataset_types = c("PDO", "PDX"),
            drug = drug_i,
            tumor_type = tumor_type_i,
            min_project_count = pdcpdx_n_datasets_t
          )
          TRUE
        },
        error = function(e) FALSE
      )
    },
    logical(1)
  )
]

eligible_pairs_path <- file.path(output_base, "eligible_drug_tumor_pairs.csv")
fwrite(eligible_pairs, eligible_pairs_path)

cat(sprintf("  OK cell_eligible_pairs: %d\n", nrow(cell_eligible_pairs)))
cat(sprintf("  OK eligible_pairs (PDO/PDX min_project_count = %d): %d\n", pdcpdx_n_datasets_t, nrow(eligible_pairs)))
cat(sprintf("  Wrote: %s\n", eligible_pairs_path))
