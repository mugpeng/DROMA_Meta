# ============================================================================
# run_all.R
# Run the full visualization workflow from collected meta-workflow outputs.
#
# Usage:
#   source("visualization/workflow/run_all.R")
#
# Required inputs:
#   OUTPUT_BASE  - path to the meta workflow batch output directory
#   VIS_OUTPUT   - path to the final output directory. Intermediate visualization
#                  files are written to a temporary staging directory and the
#                  final reorganized article package is written here.
# ============================================================================

suppressPackageStartupMessages(library(data.table))

resolve_env_path <- function(...) {
  candidates <- c(...)
  for (name in candidates) {
    value <- Sys.getenv(name, unset = NA_character_)
    if (!is.na(value) && nzchar(value)) {
      return(normalizePath(value, mustWork = FALSE))
    }
  }
  NA_character_
}

script_path <- tryCatch(sys.frames()[[1]]$ofile, error = function(e) NULL)
if (is.null(script_path) || !nzchar(script_path)) {
  script_path <- normalizePath(file.path(getwd(), "visualization", "workflow", "run_all.R"),
    mustWork = FALSE
  )
}
repo_root <- normalizePath(file.path(dirname(script_path), "../.."), mustWork = TRUE)

output_base <- resolve_env_path("OUTPUT_BASE", "DROMA_META_OUTPUT_BASE")
final_output <- resolve_env_path("VIS_OUTPUT", "DROMA_META_VIS_ARTICLE_OUTPUT", "DROMA_META_VIS_OUTPUT")

if (is.na(output_base)) {
  stop("Set OUTPUT_BASE (or DROMA_META_OUTPUT_BASE) before running run_all.R", call. = FALSE)
}
if (is.na(final_output)) {
  stop("Set VIS_OUTPUT before running run_all.R", call. = FALSE)
}

staging_root <- tempfile(pattern = "droma_meta_vis_")
dir.create(staging_root, recursive = TRUE, showWarnings = FALSE)
vis_output <- file.path(staging_root, "visualization")
article_output <- final_output

Sys.setenv(
  DROMA_META_OUTPUT_BASE = output_base,
  DROMA_META_VIS_OUTPUT = vis_output,
  DROMA_META_VIS_INPUT = vis_output,
  DROMA_META_VIS_ARTICLE_OUTPUT = article_output
)

cat("run_all.R\n")
cat("  repo_root: ", repo_root, "\n", sep = "")
cat("  OUTPUT_BASE: ", output_base, "\n", sep = "")
cat("  STAGING_OUTPUT: ", vis_output, "\n", sep = "")
cat("  FINAL_OUTPUT: ", article_output, "\n", sep = "")

scripts <- c(
  "visualization/workflow/10-Collect_Results.R",
  "visualization/workflow/11-Pipeline_Overview.R",
  "visualization/workflow/12-Volcano_Plots.R",
  "visualization/workflow/13-Heatmap_Biomarkers.R",
  "visualization/workflow/14-Concordance_Scatter.R",
  "visualization/workflow/15-Upset_Intersection.R",
  "visualization/workflow/16-Summary_Plots.R",
  "visualization/workflow/17-TCGA_AD_Plots.R",
  "visualization/workflow/18-Forest_Plots.R",
  "visualization/workflow/19-Figure1_Panel.R",
  "visualization/workflow/20-Organize_Article_Figures.R"
)

for (script in scripts) {
  script_file <- file.path(repo_root, script)
  cat("\n>>> ", script, "\n", sep = "")
  source(script_file, local = FALSE)
}

if (dir.exists(staging_root)) {
  unlink(staging_root, recursive = TRUE)
}

cat("\nVisualization workflow complete.\n")
