# ============================================================================
# 21-Breast_Shared_AD_Signature.R
# Generate a breast-cancer shared response signature from selected_genes_ad_filtered.csv
# ============================================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(DROMA.Meta.visualization))

project_root <- getVisProjectRoot()
defaults <- getVisWorkflowDefaults(project_root = project_root)
output_base <- defaults$output_base
vis_output <- defaults$vis_output
dir.create(vis_output, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 21: Breast Shared AD-Filtered Signature ===\n")

selected_ad <- collectAllSelectedAdFiltered(output_base)
fwrite(selected_ad, file.path(vis_output, "collected_selected_genes_ad_filtered.csv"))
cat(sprintf("  OK collected_selected_genes_ad_filtered: %d rows\n", nrow(selected_ad)))

signature <- prepareSharedSignatureMatrix(
  selected_ad,
  tumor_type = "breast_cancer",
  effect_col = "effect_size_cell",
  min_drugs_per_gene = 4L,
  top_n_genes = 200L,
  min_genes_per_drug = 10L
)

matrix_dt <- as.data.table(signature$matrix, keep.rownames = "name")
fwrite(matrix_dt, file.path(vis_output, "breast_shared_ad_signature_matrix.csv"))
fwrite(signature$genes, file.path(vis_output, "breast_shared_ad_signature_genes.csv"))
fwrite(signature$drugs, file.path(vis_output, "breast_shared_ad_signature_drugs.csv"))

ht <- plotSharedSignatureHeatmap(
  signature,
  title = "Breast cancer shared AD-filtered response signature",
  row_font_size = if (nrow(signature$matrix) > 120L) 3.5 else 5
)

pdf_path <- file.path(vis_output, "breast_shared_ad_signature_heatmap.pdf")
saveComplexHeatmapPdf(
  ht,
  pdf_path,
  width = max(7.2, ncol(signature$matrix) * 0.55 + 3),
  height = max(7.2, nrow(signature$matrix) * 0.045 + 3)
)

cat(sprintf(
  "  OK breast_shared_ad_signature_heatmap.pdf: %d genes x %d drugs\n",
  nrow(signature$matrix),
  ncol(signature$matrix)
))
cat("\nDone.\n")
