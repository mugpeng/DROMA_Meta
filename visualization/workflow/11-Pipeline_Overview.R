# ============================================================================
# 11-Pipeline_Overview.R
# Generate pipeline funnel and pair-level summary figures
# ============================================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DROMA.Meta))

defaults   <- getVisWorkflowDefaults()
vis_output <- defaults$vis_output
dir.create(vis_output, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 11: Pipeline Overview ===\n")

summary_path <- file.path(vis_output, "collected_summary.csv")
final_path   <- file.path(vis_output, "collected_final_biomarkers.csv")

if (!file.exists(summary_path)) {
  stop("Run 10-Collect_Results.R first", call. = FALSE)
}

summary_dt      <- fread(summary_path)
final_biomarkers <- fread(final_path)

# 1) Pipeline funnel (aggregate across all pairs)
p_funnel <- plotPipelineFunnel(summary_dt)
saveMetaVisPdf(p_funnel, file.path(vis_output, "pipeline_funnel.pdf"),
  width = 7.2, height = 4.5
)
cat("  OK pipeline_funnel.pdf\n")

# 2) Biomarker count per pair (direction stacked)
if (nrow(final_biomarkers) > 0) {
  p_count <- plotBiomarkerCountByPair(final_biomarkers)
  saveMetaVisPdf(p_count, file.path(vis_output, "biomarker_count_by_pair.pdf"),
    width = 7.2, height = 4.5
  )
  cat("  OK biomarker_count_by_pair.pdf\n")
} else {
  cat("  SKIP biomarker_count_by_pair: no final biomarkers\n")
}

cat("\nDone.\n")
