`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || !nzchar(x)) y else x
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- sub("^--file=", "", file_arg[1] %||% "")
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path, mustWork = TRUE)) else getwd()
root_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

required_files <- c(
  file.path(root_dir, "DESCRIPTION"),
  file.path(root_dir, "NAMESPACE"),
  file.path(root_dir, "one.R"),
  file.path(root_dir, "R", "FuncHelper.R"),
  file.path(root_dir, "R", "FuncValidCheck.R"),
  file.path(root_dir, "R", "FuncTcgaAD.R"),
  file.path(root_dir, "R", "FuncMetaWorkflow.R")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(sprintf("Missing required files:\n%s", paste(missing_files, collapse = "\n")))
}

desc <- read.dcf(file.path(root_dir, "DESCRIPTION"))
stopifnot(desc[1, "Package"] == "DROMA.Meta")

ns_lines <- readLines(file.path(root_dir, "NAMESPACE"), warn = FALSE)
stopifnot(any(grepl("^export\\(runMetaWorkflowOne\\)$", ns_lines)))
stopifnot(any(grepl("^export\\(runMetaWorkflowBatch\\)$", ns_lines)))

one_lines <- readLines(file.path(root_dir, "one.R"), warn = FALSE)
stopifnot(any(grepl("valid_drugs_csv", one_lines, fixed = TRUE)))
stopifnot(any(grepl("valid_tumor_types_csv", one_lines, fixed = TRUE)))
stopifnot(any(grepl("summary_csv", one_lines, fixed = TRUE)))

source(file.path(root_dir, "R", "FuncHelper.R"), local = FALSE)
source(file.path(root_dir, "R", "FuncValidCheck.R"), local = FALSE)
source(file.path(root_dir, "R", "FuncTcgaAD.R"), local = FALSE)
source(file.path(root_dir, "R", "FuncMetaWorkflow.R"), local = FALSE)

stopifnot(exists("runMetaWorkflowOne", mode = "function"))
stopifnot(exists("runMetaWorkflowBatch", mode = "function"))
stopifnot(exists("buildDrugTumorGrid", mode = "function"))

defaults <- getMetaWorkflowDefaults(root_dir)
stopifnot(is.null(defaults$valid_drugs_csv))
stopifnot(is.null(defaults$valid_tumor_types_csv))
stopifnot(is.null(defaults$batch_summary_csv))

grid <- buildDrugTumorGrid(
  file.path(root_dir, "workflow", "Output", "valid_drugs.csv"),
  file.path(root_dir, "workflow", "Output", "valid_tumor_types.csv")
)
stopifnot(nrow(grid) > 0)
stopifnot(all(c("drug", "tumor_type") %in% colnames(grid)))

cat("package smoke test passed\n")
