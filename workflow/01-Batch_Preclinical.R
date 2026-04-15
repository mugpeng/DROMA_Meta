# ============================================================================
# 01-Batch_Preclinical.R
# BatchFindSignificantFeatures on cell_sets and pdcpdx_sets for one drug/tumor
# ============================================================================

library(data.table)
library(DROMA.Set)
library(DROMA.R)
source("/Users/peng/Desktop/Project/DROMA/Meta_project3/R/FuncHelper.R", local = FALSE)
source("/Users/peng/Desktop/Project/DROMA/Meta_project3/R/FuncValidCheck.R", local = FALSE)

db_path <- "/Users/peng/Desktop/Project/DROMA/Data/droma.sqlite"
drug <- "Paclitaxel"
tumor_type <- "breast cancer"
feature2_type <- "mRNA"
data_type <- "all"
cores <- 3

cell_min_intersected_cells <- 20
pdcpdx_min_intersected_cells <- 8

output_base <- "/Users/peng/Desktop/Project/DROMA/Meta_project3/workflow/Output"
tumor_type_slug <- sanitizeName(tumor_type)
output_dir <- file.path(output_base, drug, tumor_type_slug)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 01: Batch Preclinical ===\n")
cat("drug:", drug, "\n")
cat("tumor_type:", tumor_type, "\n")

con <- connectDROMADatabase(db_path)

project_anno <- DROMA.Set::listDROMAProjects()
drug_anno <- getDROMAAnnotation("drug")
sample_anno <- getDROMAAnnotation("sample")

valid_inputs <- getValidDrugsAndTumorTypes(
  project_anno = project_anno,
  drug_anno = drug_anno,
  sample_anno = sample_anno,
  cell_n_datasets_t = 3,
  pdcpdx_n_datasets_t = 2
)
valid_drugs <- valid_inputs$valid_drugs
valid_tumor_types <- valid_inputs$valid_tumor_types

if (!drug %in% valid_drugs) {
  stop("drug does not satisfy both cell_sets and pdcpdx_sets project-count requirements: ", drug, call. = FALSE)
}
if (!tumor_type %in% valid_tumor_types) {
  stop("tumor_type does not satisfy both cell_sets and pdcpdx_sets project-count requirements: ", tumor_type, call. = FALSE)
}

cell_names_all <- project_anno[project_anno$dataset_type %in% c("CellLine", "PDC"), ]$project_name
cell_drug_projects <- unique(as.character(
  drug_anno$ProjectID[!is.na(drug_anno$DrugName) & drug_anno$DrugName == drug]
))
cell_tumor_projects <- unique(as.character(
  sample_anno$ProjectID[!is.na(sample_anno$TumorType) & sample_anno$TumorType == tumor_type]
))
cell_names <- intersect(cell_names_all, intersect(cell_drug_projects, cell_tumor_projects))
if (length(cell_names) < 3) {
  stop(
    sprintf(
      "cell_sets eligible projects for drug '%s' and tumor_type '%s' < 3; cannot satisfy n_datasets_t",
      drug, tumor_type
    ),
    call. = FALSE
  )
}
cell_sets <- createMultiDromaSetFromAllProjects(
  db_path = db_path,
  include_projects = cell_names,
  con = con
)

pdcpdx_names_all <- project_anno[project_anno$dataset_type %in% c("PDO", "PDX"), ]$project_name
pdcpdx_drug_projects <- unique(as.character(
  drug_anno$ProjectID[!is.na(drug_anno$DrugName) & drug_anno$DrugName == drug]
))
pdcpdx_tumor_projects <- unique(as.character(
  sample_anno$ProjectID[!is.na(sample_anno$TumorType) & sample_anno$TumorType == tumor_type]
))
pdcpdx_names <- intersect(pdcpdx_names_all, intersect(pdcpdx_drug_projects, pdcpdx_tumor_projects))
if (length(pdcpdx_names) < 2) {
  stop(
    sprintf(
      "pdcpdx_sets eligible projects for drug '%s' and tumor_type '%s' < 2; cannot satisfy n_datasets_t",
      drug, tumor_type
    ),
    call. = FALSE
  )
}
pdcpdx_sets <- createMultiDromaSetFromAllProjects(
  db_path = db_path,
  include_projects = pdcpdx_names,
  con = con
)

writeLines(cell_names, file.path(output_dir, "cell_sets_projects.txt"))
writeLines(pdcpdx_names, file.path(output_dir, "pdcpdx_sets_projects.txt"))

cat("  ", drug, " on cell_sets (", feature2_type, ")...\n", sep = "")
batch_cell <- batchFindSignificantFeatures(
  cell_sets,
  feature1_type = "drug",
  feature1_name = drug,
  feature2_type = feature2_type,
  data_type = data_type,
  tumor_type = tumor_type,
  overlap_only = FALSE,
  cores = cores,
  min_intersected_cells = cell_min_intersected_cells
)
fwrite(batch_cell, file.path(output_dir, "batch_cell_sets_mRNA.csv"))
cat(sprintf("  OK cell_sets: %d rows\n", nrow(batch_cell)))

cat("  ", drug, " on pdcpdx_sets (", feature2_type, ")...\n", sep = "")
batch_pdcpdx <- batchFindSignificantFeatures(
  pdcpdx_sets,
  feature1_type = "drug",
  feature1_name = drug,
  feature2_type = feature2_type,
  data_type = data_type,
  tumor_type = tumor_type,
  overlap_only = FALSE,
  cores = cores,
  min_intersected_cells = pdcpdx_min_intersected_cells
)
fwrite(batch_pdcpdx, file.path(output_dir, "batch_pdcpdx_sets_mRNA.csv"))
cat(sprintf("  OK pdcpdx_sets: %d rows\n", nrow(batch_pdcpdx)))
