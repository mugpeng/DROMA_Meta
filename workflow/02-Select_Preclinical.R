# ============================================================================
# 02-Select_Preclinical.R
# Select significant mRNA biomarkers from cell_sets and pdcpdx_sets
# ============================================================================

library(data.table)
library(DROMA.Set)
library(DROMA.R)
source("/Users/peng/Desktop/Project/DROMA/Meta_project3/R/FuncHelper.R", local = FALSE)

drug <- "Paclitaxel"
tumor_type <- "breast cancer"
output_base <- "/Users/peng/Desktop/Project/DROMA/Meta_project3/workflow/Output"
tumor_type_slug <- sanitizeName(tumor_type)
output_dir <- file.path(output_base, drug, tumor_type_slug)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 02: Select Preclinical ===\n")

batch_cell <- fread(file.path(output_dir, "batch_cell_sets_mRNA.csv"))
batch_pdcpdx <- fread(file.path(output_dir, "batch_pdcpdx_sets_mRNA.csv"))

mRNA_cell_sig <- getSignificantFeatures(
  batch_cell,
  es_t = 0.1,
  P_t = 0.05,
  use_p_value = TRUE,
  n_datasets_t = 3
)

mRNA_pdcpdx_sig <- getSignificantFeatures(
  batch_pdcpdx,
  es_t = 0.1,
  P_t = 0.05,
  use_p_value = TRUE,
  n_datasets_t = 2
)

selected_genes <- getIntersectSignificantFeatures(
  cell = mRNA_cell_sig,
  pdcpdx = mRNA_pdcpdx_sig
)

fwrite(mRNA_cell_sig, file.path(output_dir, "mRNA_cell_sig.csv"))
fwrite(mRNA_pdcpdx_sig, file.path(output_dir, "mRNA_pdcpdx_sig.csv"))
fwrite(selected_genes, file.path(output_dir, "selected_genes.csv"))

cat(sprintf("  OK mRNA_cell_sig: %d biomarkers\n", nrow(mRNA_cell_sig)))
cat(sprintf("  OK mRNA_pdcpdx_sig: %d biomarkers\n", nrow(mRNA_pdcpdx_sig)))
cat(sprintf("  OK selected_genes (intersection): %d biomarkers\n", nrow(selected_genes)))
