fallback_if_missing <- function(x, y) {
  if (is.null(x) || length(x) == 0 || !nzchar(x)) y else x
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- sub("^--file=", "", fallback_if_missing(file_arg[1], ""))
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path, mustWork = TRUE)) else getwd()
v2_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

workflow_files <- c(
  file.path(v2_root, "workflow", "21-Collect_Figure2_Inputs.R"),
  file.path(v2_root, "workflow", "22-Build_Target_Network.R"),
  file.path(v2_root, "workflow", "23-Compute_Target_Distance.R"),
  file.path(v2_root, "workflow", "24-Plot_Figure2A_Distance_Distribution.R"),
  file.path(v2_root, "workflow", "25-Plot_Figure2B_Distance_vs_Effect.R")
)

missing_files <- workflow_files[!file.exists(workflow_files)]
if (length(missing_files) > 0) {
  stop(sprintf("Missing figure2 workflow files:\n%s", paste(missing_files, collapse = "\n")))
}

for (path in workflow_files) {
  parse(file = path)
}

cat("v2 figure2 workflow smoke test passed\n")
