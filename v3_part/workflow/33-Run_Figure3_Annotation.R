# ============================================================================
# 33-Run_Figure3_Annotation.R
# Run GO enrichment and breast-state / pseudo-atlas annotations
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
suppressPackageStartupMessages(library(ggplot2))

source(file.path(root_dir, "visualization", "R", "FuncVisHelper.R"))
source(file.path(root_dir, "visualization", "R", "FuncVisBreastSignature.R"))
source(file.path(root_dir, "v3_part", "R", "FuncFigure3Helper.R"))
source(file.path(root_dir, "v3_part", "R", "FuncFigure3Plot.R"))

paths <- getFigure3Paths()
dir.create(paths$output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 33: Run Figure 3 Annotation ===\n")

selected_ad <- fread(paths$input_selected_ad)
final_biomarkers <- fread(paths$input_final_biomarkers)
breast_universe <- fread(paths$universe_csv)$name
gene_clusters <- fread(paths$gene_clusters_csv)

go_dt <- tryCatch(
  runClusterGoEnrichment(
    gene_clusters,
    universe = breast_universe,
    pvalue_cutoff = 1,
    qvalue_cutoff = 1
  ),
  error = function(e) {
    warning("GO enrichment skipped: ", conditionMessage(e), call. = FALSE)
    data.table()
  }
)
fwrite(go_dt, paths$go_csv)
if (nrow(go_dt)) {
  p_go <- plotClusterGoEnrichment(go_dt, top_n = 8L) +
    ggplot2::labs(title = "B. GO Biological Process enrichment")
  saveMetaVisPdf(p_go, paths$fig3b_pdf, width = 8.8, height = 7.2)
  cat(sprintf("  OK Figure3B GO enrichment: %d rows\n", nrow(go_dt)))
} else {
  cat("  SKIP Figure3B GO plot: no enrichment result\n")
}

marker_gene_overlay <- buildFigure3ProgramGeneOverlay(selected_ad)
marker_program_overlay <- buildFigure3ProgramOverlay(marker_gene_overlay)
fwrite(marker_gene_overlay, paths$marker_gene_overlay_csv)
fwrite(marker_program_overlay, paths$marker_program_overlay_csv)
if (nrow(marker_gene_overlay)) {
  p_marker <- plotFigure3ProgramGeneOverlay(marker_gene_overlay) +
    ggplot2::labs(title = "C. PAM50 and breast-state marker overlay")
  saveMetaVisPdf(p_marker, paths$fig3c_pdf, width = 12.5, height = 8.6)
  cat(sprintf("  OK Figure3C marker overlay: %d entries\n", nrow(marker_gene_overlay)))
} else {
  cat("  SKIP Figure3C marker overlay: no marker overlap\n")
}

cluster_program_annotation <- annotateFigure3GeneClusters(gene_clusters)
pseudo_atlas_annotation <- annotateFigure3PseudoAtlas(gene_clusters)
fwrite(cluster_program_annotation, paths$cluster_program_csv)
fwrite(pseudo_atlas_annotation, paths$pseudo_atlas_csv)

if (nrow(cluster_program_annotation)) {
  p_program <- plotFigure3ClusterProgramAnnotation(cluster_program_annotation, top_n = 3L) +
    ggplot2::labs(title = "D. Shared clusters map to breast transcriptional programs")
  saveMetaVisPdf(p_program, paths$fig3d_pdf, width = 8.8, height = 5.8)
  cat(sprintf("  OK Figure3D cluster-program annotation: %d rows\n", nrow(cluster_program_annotation)))
}

if (nrow(pseudo_atlas_annotation)) {
  p_atlas <- plotFigure3PseudoAtlasAnnotation(pseudo_atlas_annotation, top_n = 3L) +
    ggplot2::labs(title = "E. Marker-guided pseudo-atlas annotation")
  saveMetaVisPdf(p_atlas, paths$fig3e_pdf, width = 8.8, height = 5.8)
  cat(sprintf("  OK Figure3E pseudo-atlas annotation: %d rows\n", nrow(pseudo_atlas_annotation)))
}

final_overlay <- buildFinalBiomarkerClusterOverlay(final_biomarkers, gene_clusters)
fwrite(final_overlay, paths$final_overlay_csv)

