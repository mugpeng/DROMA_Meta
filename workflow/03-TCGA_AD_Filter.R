# ============================================================================
# 03-TCGA_AD_Filter.R
# Calculate TCGA Anderson-Darling statistics for selected_genes
# ============================================================================

library(data.table)
library(DROMA.Set)
library(DROMA.R)

source("/Users/peng/Desktop/Project/DROMA/Meta_project3/R/FuncHelper.R", local = FALSE)
source("/Users/peng/Desktop/Project/DROMA/Meta_project3/R/FuncValidCheck.R", local = FALSE)
source("/Users/peng/Desktop/Project/DROMA/Meta_project3/R/FuncTcgaAD.R", local = FALSE)

db_path <- "/Users/peng/Desktop/Project/DROMA/Data/droma.sqlite"
drug <- "Paclitaxel"
tumor_type <- "breast cancer"
feature2_type <- "mRNA"
data_type <- "all"

tcga_ad_p_t <- 0.01
tcga_rna_counts_dir <- "/Users/peng/Library/CloudStorage/OneDrive-Personal/28PHD_peng/250301-DROMA_project/archive260314/251112-DROMA_align/benchmark_mini/Input/TCGA/rna_counts"
gene_probe_map_path <- "/Users/peng/Desktop/Project/DROMA/Data/gencode.human.v49.annotation.gene.probeMap"

output_base <- "/Users/peng/Desktop/Project/DROMA/Meta_project3/workflow/Output"
tumor_type_slug <- sanitizeName(tumor_type)
output_dir <- file.path(output_base, drug, tumor_type_slug)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 03: TCGA AD Filter ===\n")

selected_genes <- fread(file.path(output_dir, "selected_genes.csv"))

cat(sprintf("  loaded selected_genes: %d biomarkers\n", nrow(selected_genes)))

con <- connectDROMADatabase(db_path)
on.exit(try(close(con), silent = TRUE), add = TRUE)

project_anno <- listDROMAProjects()
drug_anno <- getDROMAAnnotation("drug")
sample_anno <- getDROMAAnnotation("sample")

cell_names <- filterProjectsForDrugTumor(
  project_anno = project_anno,
  drug_anno = drug_anno,
  sample_anno = sample_anno,
  dataset_types = c("CellLine", "PDC"),
  drug = drug,
  tumor_type = tumor_type,
  min_project_count = 3
)
pdcpdx_names <- filterProjectsForDrugTumor(
  project_anno = project_anno,
  drug_anno = drug_anno,
  sample_anno = sample_anno,
  dataset_types = c("PDO", "PDX"),
  drug = drug,
  tumor_type = tumor_type,
  min_project_count = 2
)

cell_sets <- createMultiDromaSetFromAllProjects(
  db_path = db_path,
  include_projects = cell_names,
  con = con
)
pdcpdx_sets <- createMultiDromaSetFromAllProjects(
  db_path = db_path,
  include_projects = pdcpdx_names,
  con = con
)

tcga_tumor_type <- getMatchedTcgaTumorType(tumor_type)
cat("  matched TCGA/TARGET cohort:", tcga_tumor_type, "\n")

ad_stats <- batchFindTcgaADConcordantFeatures(
  selected_features = selected_genes,
  cell_set = cell_sets,
  pdcpdx_set = pdcpdx_sets,
  tumor_type = tumor_type,
  tcga_rna_counts_dir = tcga_rna_counts_dir,
  gene_probe_map_path = gene_probe_map_path,
  feature_type = feature2_type,
  data_type = data_type,
  p_t = tcga_ad_p_t
)

fwrite(ad_stats, file.path(output_dir, "selected_genes_ad_stats.csv"))
saveRDS(ad_stats, file.path(output_dir, "selected_genes_ad_stats.rds"))

cat(sprintf("  OK selected_genes_ad_stats: %d biomarkers\n", nrow(ad_stats)))
