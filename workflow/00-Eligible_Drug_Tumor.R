# ============================================================================
# 00-Eligible_Drug_Tumor.R
# Export valid drugs and tumor types that satisfy project-count thresholds
# ============================================================================

library(data.table)
suppressPackageStartupMessages(library(DROMA.Meta))
library(DROMA.Set)
library(DROMA.R)

db_path <- "/Users/peng/Desktop/Project/DROMA/Data/droma.sqlite"
output_base <- "/Users/peng/Desktop/Project/DROMA/Meta_project3/workflow/Output"
dir.create(output_base, showWarnings = FALSE, recursive = TRUE)

cell_n_datasets_t <- 3
pdcpdx_n_datasets_t <- 2

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

eligible_pairs <- data.table::CJ(
  drug = valid_drugs,
  tumor_type = valid_tumor_types,
  unique = TRUE
)

eligible_pairs <- eligible_pairs[
  vapply(
    seq_len(.N),
    function(i) {
      drug_i <- eligible_pairs$drug[[i]]
      tumor_type_i <- eligible_pairs$tumor_type[[i]]

      cell_ok <- tryCatch(
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

      if (!cell_ok) {
        return(FALSE)
      }

      pdcpdx_ok <- tryCatch(
        {
          filterProjectsForDrugTumor(
            project_anno = project_anno,
            drug_anno = drug_anno,
            sample_anno = sample_anno,
            dataset_types = c("PDO", "PDX"),
            drug = drug_i,
            tumor_type = tumor_type_i,
            min_project_count = 1
          )
          TRUE
        },
        error = function(e) FALSE
      )

      pdcpdx_ok
    },
    logical(1)
  )
]

valid_drugs <- unique(eligible_pairs$drug)
valid_tumor_types <- unique(eligible_pairs$tumor_type)

fwrite(
  data.frame(drug = valid_drugs, stringsAsFactors = FALSE),
  file.path(output_base, "valid_drugs.csv")
)
fwrite(
  data.frame(tumor_type = valid_tumor_types, stringsAsFactors = FALSE),
  file.path(output_base, "valid_tumor_types.csv")
)
fwrite(
  eligible_pairs,
  file.path(output_base, "eligible_drug_tumor_pairs.csv")
)

cat(sprintf("  OK valid_drugs: %d\n", length(valid_drugs)))
cat(sprintf("  OK valid_tumor_types: %d\n", length(valid_tumor_types)))
cat(sprintf("  OK eligible_pairs: %d\n", nrow(eligible_pairs)))
