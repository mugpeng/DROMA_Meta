get_workflow_root <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  candidates <- c(
    sub("^--file=", "", file_arg),
    if (!is.null(sys.frames()[[1]]$ofile)) sys.frames()[[1]]$ofile else character(0),
    file.path(getwd(), "workflow", "00-Workflow_Common.R"),
    file.path(getwd(), "00-Workflow_Common.R")
  )
  candidates <- unique(candidates[nzchar(candidates)])
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop("Could not locate workflow root")
  }
  dirname(normalizePath(existing[[1]], mustWork = TRUE))
}

workflow_root <- get_workflow_root()
meta_project_root <- normalizePath(dirname(workflow_root), mustWork = TRUE)
repo_root <- normalizePath(dirname(meta_project_root), mustWork = TRUE)
droma_r2_root <- file.path(meta_project_root, "DROMA_R2")

library(data.table)
library(DROMA.Set)
library(DROMA.R)

source_dir <- function(path) {
  files <- sort(list.files(path, pattern = "\\.[Rr]$", full.names = TRUE))
  invisible(lapply(files, source, local = FALSE))
}

source(file.path(repo_root, "DB_project", "DROMA_R", "R", "FuncPairDataLoading.R"), local = FALSE)
source(file.path(repo_root, "DB_project", "DROMA_R", "R", "FuncPairDataPairing.R"), local = FALSE)
source(file.path(repo_root, "DB_project", "DROMA_R", "R", "FuncPairMetaAnalysis.R"), local = FALSE)
source(file.path(repo_root, "DB_project", "DROMA_R", "R", "FuncPairBatchFeature.R"), local = FALSE)
source(file.path(repo_root, "DB_project", "DROMA_R", "R", "FuncClinical.R"), local = FALSE)
source_dir(file.path(droma_r2_root, "R"))

workflow_config <- buildWorkflowConfig(
  repo_root = repo_root,
  output_base = file.path(workflow_root, "Output"),
  drug_names = "Paclitaxel",
  tumor_types = "breast cancer",
  feature_type = "mRNA",
  es_t = 0.1,
  fdr_t = 0.1,
  pair_cor_t = 0.2,
  pair_p_t = 0.05,
  tcga_fdr_t = 0.01,
  requested_cores = 3L
)

dir.create(workflow_config$output_base, recursive = TRUE, showWarnings = FALSE)

save_stage <- function(object, subdir, filename) {
  out_dir <- file.path(workflow_config$output_base, subdir)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(object, file = file.path(out_dir, filename))
}

read_stage <- function(subdir, filename) {
  readRDS(file.path(workflow_config$output_base, subdir, filename))
}
