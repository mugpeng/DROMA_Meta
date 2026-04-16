# ============================================================================
# 31-Collect_Figure3_Inputs.R
# Collect breast cancer inputs for Figure 3 shared-state analyses
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

source(file.path(root_dir, "visualization", "R", "FuncVisHelper.R"))
source(file.path(root_dir, "visualization", "R", "FuncVisBreastSignature.R"))
source(file.path(root_dir, "v3_part", "R", "FuncFigure3Helper.R"))

paths <- getFigure3Paths()
dir.create(paths$output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 31: Collect Figure 3 Inputs ===\n")

selected_ad <- collectFigure3SelectedAdFiltered(paths$output_base)
final_biomarkers <- collectFigure3FinalBiomarkers(paths$output_base)

fwrite(selected_ad, paths$input_selected_ad)
fwrite(final_biomarkers, paths$input_final_biomarkers)

breast_universe <- unique(selected_ad[tumor_type %in% c("breast_cancer", "breast cancer"), name])
fwrite(data.table(name = sort(breast_universe)), paths$universe_csv)

cat(sprintf("  OK selected AD-filtered rows: %d\n", nrow(selected_ad)))
cat(sprintf("  OK final biomarker rows: %d\n", nrow(final_biomarkers)))
cat(sprintf("  OK breast universe genes: %d\n", length(breast_universe)))

