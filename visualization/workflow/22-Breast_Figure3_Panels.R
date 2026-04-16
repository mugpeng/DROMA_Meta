# ============================================================================
# 22-Breast_Figure3_Panels.R
# Figure 3-like breast cancer shared AD-filtered response program panels.
# ============================================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(DROMA.Meta.visualization))

project_root <- getVisProjectRoot()
defaults <- getVisWorkflowDefaults(project_root = project_root)
output_base <- defaults$output_base
vis_output <- defaults$vis_output
figure_dir <- file.path(vis_output, "breast_figure3")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 22: Breast Figure 3-like Panels ===\n")

selected_ad <- collectAllSelectedAdFiltered(output_base)
final_biomarkers <- collectAllFinalBiomarkers(output_base)

breast_universe <- unique(selected_ad[tumor_type %in% c("breast_cancer", "breast cancer"), name])
fwrite(data.table(name = breast_universe), file.path(figure_dir, "breast_selected_ad_universe.csv"))

signature <- prepareSharedSignatureMatrix(
  selected_ad,
  tumor_type = "breast_cancer",
  effect_col = "effect_size_cell",
  min_drugs_per_gene = 4L,
  top_n_genes = 200L,
  min_genes_per_drug = 10L
)
matrix_dt <- as.data.table(signature$matrix, keep.rownames = "name")
fwrite(matrix_dt, file.path(figure_dir, "Figure3A_breast_shared_signature_matrix.csv"))
fwrite(signature$genes, file.path(figure_dir, "Figure3A_breast_shared_signature_genes.csv"))
fwrite(signature$drugs, file.path(figure_dir, "Figure3A_breast_shared_signature_drugs.csv"))

gene_clusters <- assignSharedSignatureClusters(signature$matrix, k = 3L)
drug_clusters <- assignSharedSignatureDrugClusters(signature$matrix, k = 3L)
fwrite(gene_clusters, file.path(figure_dir, "Figure3A_gene_clusters.csv"))
fwrite(drug_clusters, file.path(figure_dir, "Figure3A_drug_clusters.csv"))

ht <- plotSharedSignatureHeatmap(
  signature,
  title = "A. Breast shared AD-filtered response signature",
  row_font_size = 3.5
)
saveComplexHeatmapPdf(
  ht,
  file.path(figure_dir, "Figure3A_breast_shared_signature_heatmap.pdf"),
  width = max(7.2, ncol(signature$matrix) * 0.55 + 3),
  height = max(7.2, nrow(signature$matrix) * 0.045 + 3)
)
cat(sprintf("  OK Figure3A heatmap: %d genes x %d drugs\n", nrow(signature$matrix), ncol(signature$matrix)))

go_dt <- runClusterGoEnrichment(
  gene_clusters,
  universe = breast_universe,
  pvalue_cutoff = 1,
  qvalue_cutoff = 1
)
fwrite(go_dt, file.path(figure_dir, "Figure3B_GO_enrichment_by_gene_cluster.csv"))
if (nrow(go_dt)) {
  p_go <- plotClusterGoEnrichment(go_dt, top_n = 8L) +
    ggplot2::labs(title = "B. GO Biological Process enrichment")
  saveMetaVisPdf(
    p_go,
    file.path(figure_dir, "Figure3B_GO_enrichment_by_gene_cluster.pdf"),
    width = 8.4,
    height = 7.2
  )
  cat(sprintf("  OK Figure3B GO enrichment: %d terms\n", nrow(go_dt)))
} else {
  cat("  SKIP Figure3B GO plot: no significant terms\n")
}

marker_dt <- buildBreastMarkerOverlay(
  selected_ad,
  tumor_type = "breast_cancer",
  effect_col = "effect_size_cell"
)
fwrite(marker_dt, file.path(figure_dir, "Figure3C_breast_marker_overlay.csv"))
if (nrow(marker_dt)) {
  p_marker <- plotBreastMarkerOverlay(marker_dt) +
    ggplot2::labs(title = "C. Breast subtype and response marker overlay")
  saveMetaVisPdf(
    p_marker,
    file.path(figure_dir, "Figure3C_breast_marker_overlay.pdf"),
    width = 10.5,
    height = 4.8
  )
  cat(sprintf("  OK Figure3C marker overlay: %d marker-drug entries\n", nrow(marker_dt)))
} else {
  cat("  SKIP Figure3C marker overlay: no marker entries in selected AD-filtered genes\n")
}

final_overlay <- buildFinalBiomarkerClusterOverlay(final_biomarkers, gene_clusters)
fwrite(final_overlay, file.path(figure_dir, "Figure3D_final_biomarker_cluster_overlay.csv"))
if (nrow(final_overlay)) {
  p_final <- plotFinalBiomarkerClusterOverlay(final_overlay) +
    ggplot2::labs(title = "D. Final biomarkers mapped to shared signature clusters")
  saveMetaVisPdf(
    p_final,
    file.path(figure_dir, "Figure3D_final_biomarker_cluster_overlay.pdf"),
    width = 7.2,
    height = 4.8
  )
  cat(sprintf("  OK Figure3D final overlay: %d cluster-drug entries\n", nrow(final_overlay)))
} else {
  cat("  SKIP Figure3D final overlay: no final biomarkers overlap shared signature genes\n")
}

writeLines(c(
  "# Breast Figure 3-like panels",
  "",
  "These panels adapt the lung-cancer Figure 3 logic to breast cancer using `selected_genes_ad_filtered.csv` rather than final biomarkers.",
  "",
  "- `Figure3A_breast_shared_signature_heatmap.pdf`: shared AD-filtered gene-by-drug effect-size signature.",
  "- `Figure3B_GO_enrichment_by_gene_cluster.pdf`: clusterProfiler GO Biological Process enrichment by shared gene cluster.",
  "- `Figure3C_breast_marker_overlay.pdf`: breast subtype/response marker effect-size overlay.",
  "- `Figure3D_final_biomarker_cluster_overlay.pdf`: sparse final biomarkers mapped back to shared AD-filtered clusters.",
  "",
  "Interpretation note: the shared program is defined at the AD-filtered selected-gene layer; final biomarkers remain drug-specific and are shown as an overlay rather than the clustering input."
), file.path(figure_dir, "legend.md"), useBytes = TRUE)

cat("\nDone.\n")
