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
stopifnot(any(grepl("^export\\(runMetaWorkflow\\)$", ns_lines)))

one_lines <- readLines(file.path(root_dir, "one.R"), warn = FALSE)
stopifnot(any(grepl("library(DROMA.Meta)", one_lines, fixed = TRUE)))
stopifnot(!any(grepl("source(\"R/", one_lines, fixed = TRUE)))
stopifnot(any(grepl("drugs <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("tumor_types <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("for (drug in drugs)", one_lines, fixed = TRUE)))
stopifnot(any(grepl("for (tumor_type in tumor_types)", one_lines, fixed = TRUE)))
stopifnot(any(grepl("db_path <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("tcga_rna_counts_dir <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("gene_probe_map_path <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("data_type <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("cores <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("cell_min_intersected_cells <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("pdcpdx_min_intersected_cells <-", one_lines, fixed = TRUE)))

suppressPackageStartupMessages(library(DROMA.Meta))

stopifnot(exists("runMetaWorkflow", mode = "function"))
stopifnot(exists("buildDrugTumorGrid", mode = "function"))

defaults <- getMetaWorkflowDefaults(root_dir)
stopifnot(is.null(defaults$valid_drugs_csv))
stopifnot(is.null(defaults$valid_tumor_types_csv))
stopifnot(is.null(defaults$batch_summary_csv))
stopifnot(identical(defaults$droma_db_path, "/Users/peng/Desktop/Project/DROMA/Data/droma.sqlite"))

grid <- buildDrugTumorGrid(
  file.path(root_dir, "workflow", "Output", "valid_drugs.csv"),
  file.path(root_dir, "workflow", "Output", "valid_tumor_types.csv")
)
stopifnot(nrow(grid) > 0)
stopifnot(all(c("drug", "tumor_type") %in% colnames(grid)))

cat("package smoke test passed\n")
