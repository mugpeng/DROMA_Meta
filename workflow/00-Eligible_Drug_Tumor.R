# ============================================================================
# 00-Eligible_Drug_Tumor.R
# Export valid drugs and tumor types that satisfy project-count thresholds
# ============================================================================

library(data.table)
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

valid_inputs <- .get_valid_drugs_and_tumor_types(
  project_anno = project_anno,
  drug_anno = drug_anno,
  sample_anno = sample_anno,
  cell_n_datasets_t = cell_n_datasets_t,
  pdcpdx_n_datasets_t = pdcpdx_n_datasets_t
)

valid_drugs <- valid_inputs$valid_drugs
valid_tumor_types <- valid_inputs$valid_tumor_types

writeLines(valid_drugs, file.path(output_base, "valid_drugs.txt"))
writeLines(valid_tumor_types, file.path(output_base, "valid_tumor_types.txt"))

fwrite(
  data.frame(drug = valid_drugs, stringsAsFactors = FALSE),
  file.path(output_base, "valid_drugs.csv")
)
fwrite(
  data.frame(tumor_type = valid_tumor_types, stringsAsFactors = FALSE),
  file.path(output_base, "valid_tumor_types.csv")
)

cat(sprintf("  OK valid_drugs: %d\n", length(valid_drugs)))
cat(sprintf("  OK valid_tumor_types: %d\n", length(valid_tumor_types)))
