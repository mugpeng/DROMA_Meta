# ============================================================================
# 21-Collect_Figure2_Inputs.R
# Collect and normalize Figure 2 biomarker inputs from meta_batch outputs
# ============================================================================

suppressPackageStartupMessages(library(data.table))

root_dir <- normalizePath(file.path(getwd(), "v2_part"), mustWork = TRUE)
source(file.path(root_dir, "R", "FuncFigure2Helper.R"))

paths <- getFigure2Paths()

cat("\n=== 21: Collect Figure 2 Inputs ===\n")

feature_dt <- collectFigure2Inputs()
fwrite(feature_dt, paths$input_features_csv)

stage_summary <- feature_dt[
  ,
  .(
    n_rows = .N,
    n_unique_pairs = uniqueN(paste(drug, tumor_type)),
    n_unique_genes = uniqueN(name),
    n_non_missing_effect = sum(!is.na(effect_size))
  ),
  by = .(stage)
]
fwrite(stage_summary, paths$stage_summary_csv)

cat(sprintf("  OK figure2_input_features: %d rows\n", nrow(feature_dt)))
cat(sprintf("  OK figure2_stage_summary: %d stages\n", nrow(stage_summary)))
cat("\nDone.\n")
