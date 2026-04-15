# ============================================================================
# 03-TCGA_AD_Filter.R
# Calculate TCGA Anderson-Darling statistics for selected_genes
# ============================================================================

library(data.table)
library(DROMA.Set)
library(DROMA.R)

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

ccle_set <- createMultiDromaSetFromAllProjects(
  db_path = db_path,
  include_projects = "CCLE",
  con = con
)

tcga_tumor_type <- getMatchedTcgaTumorType(tumor_type)
cat("  matched TCGA/TARGET cohort:", tcga_tumor_type, "\n")
cat("  preclinical cohort: CCLE\n")

ad_stats <- batchFindTcgaADConcordantFeatures(
  selected_features = selected_genes,
  preclinical_set = ccle_set,
  tumor_type = tumor_type,
  tcga_rna_counts_dir = tcga_rna_counts_dir,
  gene_probe_map_path = gene_probe_map_path,
  feature_type = feature2_type,
  data_type = data_type,
  p_t = tcga_ad_p_t,
  preclinical_label = "ccle"
)

fwrite(ad_stats, file.path(output_dir, "selected_genes_ad_stats.csv"))
# saveRDS(ad_stats, file.path(output_dir, "selected_genes_ad_stats.rds"))

cat(sprintf("  OK selected_genes_ad_stats: %d biomarkers\n", nrow(ad_stats)))

selected_genes_ad_filtered <- fread(file.path(output_dir, "selected_genes_ad_filtered.csv"))
