# ============================================================================
# 18-Forest_Plots.R
# Generate stage-wise forest plots for top biomarkers
# ============================================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DROMA.Meta))

defaults   <- getVisWorkflowDefaults()
vis_output <- defaults$vis_output
dir.create(vis_output, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 18: Forest Plots ===\n")

final_path <- file.path(vis_output, "collected_final_biomarkers.csv")
if (!file.exists(final_path)) {
  stop("Run 10-Collect_Results.R first", call. = FALSE)
}

final_biomarkers <- fread(final_path)

if (nrow(final_biomarkers) == 0) {
  cat("  SKIP forest: no final biomarkers\n")
} else {
  top_n <- min(15, nrow(final_biomarkers))
  p_forest <- plotStageForest(final_biomarkers, top_n = top_n)
  saveMetaVisPdf(p_forest, file.path(vis_output, "biomarker_forest.pdf"),
    width = 7.2, height = max(5, top_n * 0.4 + 2)
  )
  cat(sprintf("  OK biomarker_forest.pdf (top %d)\n", top_n))
}

cat("\nDone.\n")
