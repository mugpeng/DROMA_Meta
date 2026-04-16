# ============================================================================
# 10-Collect_Results.R
# Scan all output directories and build master summary tables for visualization
# ============================================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DROMA.Meta))

# --- Configuration ---
project_root <- file.path(normalizePath(getwd(), mustWork = TRUE), "Meta_Example")
defaults    <- getVisWorkflowDefaults(project_root = project_root)
output_base <- defaults$output_base
vis_output  <- defaults$vis_output
dir.create(vis_output, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 10: Collect Results ===\n")

# 1) Pair-level summary (from clinical_validation_info.csv)
summary_dt <- collectWorkflowResults(output_base)
fwrite(summary_dt, file.path(vis_output, "collected_summary.csv"))
cat(sprintf("  OK collected_summary: %d pairs\n", nrow(summary_dt)))

# 2) All final biomarkers merged
final_biomarkers <- collectAllFinalBiomarkers(output_base)
fwrite(final_biomarkers, file.path(vis_output, "collected_final_biomarkers.csv"))
cat(sprintf("  OK collected_final_biomarkers: %d biomarkers\n", nrow(final_biomarkers)))

# 3) All batch results (cell and PDC/PDX separately)
batch_cell <- collectAllBatchResults(output_base, type = "cell")
fwrite(batch_cell, file.path(vis_output, "collected_batch_cell.csv"))
cat(sprintf("  OK collected_batch_cell: %d rows\n", nrow(batch_cell)))

batch_pdcpdx <- collectAllBatchResults(output_base, type = "pdcpdx")
fwrite(batch_pdcpdx, file.path(vis_output, "collected_batch_pdcpdx.csv"))
cat(sprintf("  OK collected_batch_pdcpdx: %d rows\n", nrow(batch_pdcpdx)))

# 4) AD statistics
ad_stats <- collectAllAdStats(output_base)
fwrite(ad_stats, file.path(vis_output, "collected_ad_stats.csv"))
cat(sprintf("  OK collected_ad_stats: %d rows\n", nrow(ad_stats)))

# 5) Significant features (cell and PDC/PDX)
cell_sig <- collectAllSigFeatures(output_base, type = "cell")
fwrite(cell_sig, file.path(vis_output, "collected_cell_sig.csv"))
cat(sprintf("  OK collected_cell_sig: %d rows\n", nrow(cell_sig)))

pdcpdx_sig <- collectAllSigFeatures(output_base, type = "pdcpdx")
fwrite(pdcpdx_sig, file.path(vis_output, "collected_pdcpdx_sig.csv"))
cat(sprintf("  OK collected_pdcpdx_sig: %d rows\n", nrow(pdcpdx_sig)))

cat("\nDone.\n")
