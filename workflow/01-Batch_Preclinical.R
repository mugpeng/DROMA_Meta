# ============================================================================
# 01-Batch_Preclinical.R
# BatchFindSignificantFeatures on cell_sets and pdcpdx_sets for one drug/tumor
# ============================================================================

library(data.table)
suppressPackageStartupMessages(library(DROMA.Meta))
library(DROMA.Set)
library(DROMA.R)

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

cell_names_all <- project_anno[project_anno$dataset_type %in% c("CellLine", "PDC"), ]$project_name
cell_drug_projects <- unique(as.character(
  drug_anno$ProjectID[!is.na(drug_anno$DrugName) & drug_anno$DrugName == drug]
))

cell_sets <- createMultiDromaSetFromAllProjects(
  db_path = db_path,
  include_projects = cell_names_all,
  con = con
)

pdcpdx_names_all <- project_anno[project_anno$dataset_type %in% c("PDO", "PDX"), ]$project_name
pdcpdx_sets <- createMultiDromaSetFromAllProjects(
  db_path = db_path,
  include_projects = pdcpdx_names_all,
  con = con
)

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
