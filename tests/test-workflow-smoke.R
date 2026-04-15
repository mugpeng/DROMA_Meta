`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || !nzchar(x)) y else x
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- sub("^--file=", "", file_arg[1] %||% "")
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path, mustWork = TRUE)) else getwd()
root_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

workflow_files <- c(
  file.path(root_dir, "workflow", "00-Eligible_Drug_Tumor.R"),
  file.path(root_dir, "workflow", "01-Batch_Preclinical.R"),
  file.path(root_dir, "workflow", "02-Select_Preclinical.R"),
  file.path(root_dir, "workflow", "03-TCGA_AD_Filter.R"),
  file.path(root_dir, "workflow", "04-Clinical_Validation.R")
)

missing_files <- workflow_files[!file.exists(workflow_files)]
if (length(missing_files) > 0) {
  stop(sprintf("Missing workflow files:\n%s", paste(missing_files, collapse = "\n")))
}

for (path in workflow_files) {
  parse(file = path)
}

source(file.path(root_dir, "R", "FuncHelper.R"), local = FALSE)
source(file.path(root_dir, "R", "FuncValidCheck.R"), local = FALSE)
source(file.path(root_dir, "R", "FuncTcgaAD.R"), local = FALSE)
source(file.path(root_dir, "R", "FuncMetaWorkflow.R"), local = FALSE)

stopifnot(exists("runMetaWorkflowBatch", mode = "function"))
stopifnot(exists("runMetaWorkflowOne", mode = "function"))

cat("workflow smoke test passed\n")
