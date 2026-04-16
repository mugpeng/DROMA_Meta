# ============================================================================
# 13-Heatmap_Biomarkers.R
# Generate biomarker effect-size heatmap
# ============================================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(DROMA.Meta.visualization))

project_root <- getVisProjectRoot()
defaults   <- getVisWorkflowDefaults(project_root = project_root)
vis_output <- defaults$vis_output
dir.create(vis_output, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 13: Heatmap Biomarkers ===\n")

final_path <- file.path(vis_output, "collected_final_biomarkers.csv")
if (!file.exists(final_path)) {
  stop("Run 10-Collect_Results.R first", call. = FALSE)
}

final_biomarkers <- fread(final_path)

if (nrow(final_biomarkers) == 0) {
  cat("  SKIP heatmap: no final biomarkers\n")
} else {
  ht <- plotBiomarkerHeatmap(final_biomarkers,
    es_col = "effect_size_ctrdb",
    title = "Final biomarker effect sizes",
    show_cell_text = nrow(final_biomarkers) <= 40
  )

  pdf_path <- file.path(vis_output, "biomarker_heatmap.pdf")
  saveComplexHeatmapPdf(ht, pdf_path,
    width = max(7.2, length(unique(paste(final_biomarkers$drug, final_biomarkers$tumor_type))) * 2 + 3),
    height = max(5, nrow(final_biomarkers) * 0.3 + 3)
  )
  cat(sprintf("  OK biomarker_heatmap.pdf\n"))
}

cat("\nDone.\n")
