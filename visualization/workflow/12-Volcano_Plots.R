# ============================================================================
# 12-Volcano_Plots.R
# Generate volcano plots for each drug-tumor pair
# ============================================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DROMA.Meta))

defaults   <- getVisWorkflowDefaults()
vis_output <- defaults$vis_output
volcano_dir <- file.path(vis_output, "volcano")
dir.create(volcano_dir, showWarnings = FALSE, recursive = TRUE)

es_t <- 0.1
P_t <- 0.05
use_p_value <- FALSE

cat("\n=== 12: Volcano Plots ===\n")

batch_cell_path <- file.path(vis_output, "collected_batch_cell.csv")
if (!file.exists(batch_cell_path)) {
  stop("Run 10-Collect_Results.R first", call. = FALSE)
}

batch_cell <- fread(batch_cell_path)

batchPlotVolcano(
  batch_all  = batch_cell,
  output_dir = volcano_dir,
  es_t       = es_t,
  P_t        = P_t,
  use_p_value = use_p_value
)

cat("\nDone.\n")
