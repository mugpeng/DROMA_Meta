# ============================================================================
# 03-TCGA_AD_Filter.R
# Filter selected_genes by TCGA concordance using Anderson-Darling tests
# ============================================================================

library(data.table)
library(DROMA.Set)
library(DROMA.R)

source("/Users/peng/Desktop/Project/DROMA/Meta_project3/R/run_drug_tumor_biomarker_workflow.R", local = FALSE)

db_path <- "/Users/peng/Desktop/Project/DROMA/Data/droma.sqlite"
drug <- "Paclitaxel"
tumor_type <- "breast cancer"
feature2_type <- "mRNA"
data_type <- "all"

tcga_ad_p_t <- 0.05
tcga_rna_counts_dir <- "/Users/peng/Library/CloudStorage/OneDrive-Personal/28PHD_peng/250301-DROMA_project/archive260314/251112-DROMA_align/benchmark_mini/Input/TCGA/rna_counts"
gene_probe_map_path <- "/Users/peng/Desktop/Project/DROMA/Data/gencode.human.v49.annotation.gene.probeMap"

output_base <- "/Users/peng/Desktop/Project/DROMA/Meta_project3/workflow/Output"
tumor_type_slug <- .sanitize_name(tumor_type)
output_dir <- file.path(output_base, drug, tumor_type_slug)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 03: TCGA AD Filter ===\n")

config <- .resolve_biomarker_workflow_config(
  drug = drug,
  tumor_type = tumor_type,
  db_path = db_path,
  feature2_type = feature2_type,
  data_type = data_type,
  output_base = output_base,
  tcga_ad_p_t = tcga_ad_p_t,
  tcga_rna_counts_dir = tcga_rna_counts_dir,
  gene_probe_map_path = gene_probe_map_path,
  workflow_dir = file.path("/Users/peng/Desktop/Project/DROMA/Meta_project3", "workflow")
)

selected_genes <- fread(file.path(output_dir, "selected_genes.csv"))
ad_results <- .run_tcga_ad(
  config = config,
  selected_genes = selected_genes
)

cat(sprintf("  OK selected_genes_ad_stats: %d biomarkers\n", nrow(ad_results$ad_stats)))
cat(sprintf("  OK selected_genes_ad_filtered: %d biomarkers\n", nrow(ad_results$ad_filtered)))
