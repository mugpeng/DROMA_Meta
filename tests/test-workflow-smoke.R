`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || !nzchar(x)) y else x
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- sub("^--file=", "", file_arg[1] %||% "")
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path, mustWork = TRUE)) else getwd()
root_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

required_files <- c(
  file.path(root_dir, "R", "run_drug_tumor_biomarker_workflow.R"),
  file.path(root_dir, "workflow", "00-Eligible_Drug_Tumor.R"),
  file.path(root_dir, "workflow", "01-Batch_Preclinical.R"),
  file.path(root_dir, "workflow", "02-Select_Preclinical.R"),
  file.path(root_dir, "workflow", "03-Clinical_Validation.R")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(sprintf("Missing required files:\n%s", paste(missing_files, collapse = "\n")))
}

source(required_files[[1]], local = FALSE)

if (!exists("run_drug_tumor_biomarker_workflow", mode = "function")) {
  stop("run_drug_tumor_biomarker_workflow() is not defined")
}

if (!exists(".resolve_biomarker_workflow_config", mode = "function")) {
  stop(".resolve_biomarker_workflow_config() is not defined")
}

cfg <- .resolve_biomarker_workflow_config(
  drug = "Paclitaxel",
  tumor_type = "breast cancer",
  output_base = tempfile("biomarker-output-"),
  workflow_dir = file.path(root_dir, "workflow")
)

stopifnot(cfg$feature2_type == "mRNA")
stopifnot(cfg$cell_min_intersected_cells == 20L)
stopifnot(cfg$pdcpdx_min_intersected_cells == 8L)
stopifnot(identical(cfg$output_dir, file.path(cfg$output_base, "Paclitaxel", "breast_cancer")))

for (path in required_files[-1]) {
  parse(file = path)
}

cat("workflow smoke test passed\n")
