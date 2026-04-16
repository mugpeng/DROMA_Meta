# ============================================================================
# 16-Summary_Plots.R
# Generate direction distribution and stage comparison figures
# ============================================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DROMA.Meta))

defaults   <- getVisWorkflowDefaults()
vis_output <- defaults$vis_output
dir.create(vis_output, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 16: Summary Plots ===\n")

final_path   <- file.path(vis_output, "collected_final_biomarkers.csv")
summary_path <- file.path(vis_output, "collected_summary.csv")

if (!file.exists(summary_path)) {
  stop("Run 10-Collect_Results.R first", call. = FALSE)
}

final_biomarkers <- fread(final_path)
summary_dt       <- fread(summary_path)

# 1) Direction distribution
if (nrow(final_biomarkers) > 0 && "direction" %in% names(final_biomarkers)) {
  p_dir <- plotDirectionSummary(final_biomarkers)
  saveMetaVisPdf(p_dir, file.path(vis_output, "direction_summary.pdf"),
    width = 7.2, height = 4.5
  )
  cat("  OK direction_summary.pdf\n")
} else {
  cat("  SKIP direction_summary: no direction data\n")
}

# 2) Dumbbell attrition chart
if (nrow(summary_dt) > 0) {
  p_stage <- plotStageCountComparison(summary_dt)
  saveMetaVisPdf(p_stage, file.path(vis_output, "stage_attrition_dumbbell.pdf"),
    width = 7.2, height = 4.5
  )
  cat("  OK stage_attrition_dumbbell.pdf\n")
} else {
  cat("  SKIP stage_attrition: no summary data\n")
}

cat("\nDone.\n")
