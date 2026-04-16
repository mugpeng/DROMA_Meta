# ============================================================================
# 23-Compute_Target_Distance.R
# Compute biomarker-to-target and background distance distributions
# ============================================================================

suppressPackageStartupMessages(library(data.table))

root_dir <- normalizePath(file.path(getwd(), "v2_part"), mustWork = TRUE)
source(file.path(root_dir, "R", "FuncFigure2Helper.R"))

paths <- getFigure2Paths()

cat("\n=== 23: Compute Target Distance ===\n")

required <- c(paths$input_features_csv, paths$drug_target_csv, paths$reactome_edges_csv)
missing_files <- required[!file.exists(required)]
if (length(missing_files)) {
  stop(sprintf("Missing required inputs:\n%s", paste(missing_files, collapse = "\n")), call. = FALSE)
}

feature_dt <- fread(paths$input_features_csv)
target_dt <- fread(paths$drug_target_csv)
edge_dt <- fread(paths$reactome_edges_csv)

cat("  computing biomarker-to-target distances...\n")
distance_dt <- computeFigure2Distances(feature_dt, target_dt, edge_dt)
cat("  summarizing distance distributions...\n")
summary_dt <- summarizeFigure2Distances(distance_dt)

fwrite(distance_dt, paths$distance_results_csv)
fwrite(summary_dt, paths$distance_summary_csv)

cat(sprintf("  OK distance_results: %d rows\n", nrow(distance_dt)))
cat(sprintf("  OK distance_summary_by_stage: %d rows\n", nrow(summary_dt)))
cat("\nDone.\n")
