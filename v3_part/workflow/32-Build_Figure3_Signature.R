# ============================================================================
# 32-Build_Figure3_Signature.R
# Build breast shared-state signature matrix and cluster assignments
# ============================================================================

fallback_if_missing <- function(x, y) {
  if (is.null(x) || length(x) == 0 || !nzchar(x)) y else x
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- sub("^--file=", "", fallback_if_missing(file_arg[1], ""))
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path, mustWork = TRUE)) else getwd()
root_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(grid))

source(file.path(root_dir, "visualization", "R", "FuncVisHelper.R"))
source(file.path(root_dir, "visualization", "R", "FuncVisBreastSignature.R"))
source(file.path(root_dir, "v3_part", "R", "FuncFigure3Helper.R"))

paths <- getFigure3Paths()
dir.create(paths$output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 32: Build Figure 3 Signature ===\n")

selected_ad <- fread(paths$input_selected_ad)

signature <- prepareSharedSignatureMatrix(
  selected_ad,
  tumor_type = "breast_cancer",
  effect_col = "effect_size_cell",
  min_drugs_per_gene = 4L,
  top_n_genes = 240L,
  min_genes_per_drug = 10L
)
gene_clusters <- assignSharedSignatureClusters(signature$matrix, k = 3L)
drug_clusters <- assignSharedSignatureDrugClusters(signature$matrix, k = 3L)

fwrite(as.data.table(signature$matrix, keep.rownames = "name"), paths$signature_matrix_csv)
fwrite(signature$genes, paths$signature_genes_csv)
fwrite(signature$drugs, paths$signature_drugs_csv)
fwrite(gene_clusters, paths$gene_clusters_csv)
fwrite(drug_clusters, paths$drug_clusters_csv)

ht <- plotSharedSignatureHeatmap(
  signature,
  title = "A. Breast shared-state response signature",
  row_font_size = 3.5
)
saveComplexHeatmapPdf(
  ht,
  paths$fig3a_pdf,
  width = max(7.2, ncol(signature$matrix) * 0.4 + 3),
  height = max(7.2, nrow(signature$matrix) * 0.045 + 3)
)

cat(sprintf("  OK Figure3A heatmap: %d genes x %d drugs\n", nrow(signature$matrix), ncol(signature$matrix)))
cat(sprintf("  OK gene clusters: %d\n", uniqueN(gene_clusters$gene_cluster)))
cat(sprintf("  OK drug clusters: %d\n", uniqueN(drug_clusters$drug_cluster)))

