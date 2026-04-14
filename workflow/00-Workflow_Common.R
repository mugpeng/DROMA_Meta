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

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || !nzchar(x)) y else x
}

get_current_script_dir <- function(default_name) {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  candidate <- sub("^--file=", "", file_arg[1] %||% "")
  if (nzchar(candidate) && file.exists(candidate)) {
    return(dirname(normalizePath(candidate, mustWork = TRUE)))
  }
  if (!is.null(sys.frames()[[1]]$ofile) && file.exists(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile, mustWork = TRUE)))
  }
  normalizePath(file.path(getwd(), "workflow"), mustWork = FALSE)
}

workflow_root <- get_workflow_root()
find_repo_root <- function(start_dir) {
  current <- normalizePath(start_dir, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(current, "Data", "droma.sqlite"))) {
      return(current)
    }
    if (dir.exists(file.path(current, "DB_project")) && dir.exists(file.path(current, "Meta_project"))) {
      return(current)
    }
    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate repository root from workflow path (expect ancestor with Data/droma.sqlite or DB_project+Meta_project)")
    }
    current <- parent
  }
}

repo_root <- find_repo_root(workflow_root)
meta_project_root <- normalizePath(dirname(workflow_root), mustWork = TRUE)
droma_r2_root <- normalizePath(file.path(meta_project_root, "DROMA_R2"), mustWork = FALSE)
if (!dir.exists(file.path(droma_r2_root, "R")) ||
    length(list.files(file.path(droma_r2_root, "R"), pattern = "\\.[Rr]$", full.names = TRUE)) == 0) {
  droma_r2_root <- normalizePath(file.path(repo_root, "Meta_project", "cellline_meta", "_staging", "DROMA_R2"), mustWork = FALSE)
}

library(data.table)
library(DROMA.Set)
library(DROMA.R)

source_dir <- function(path) {
  files <- sort(list.files(path, pattern = "\\.[Rr]$", full.names = TRUE))
  invisible(lapply(files, source, local = FALSE))
}

source_dir(file.path(droma_r2_root, "R"))

workflow_config <- buildWorkflowConfig(
  repo_root = repo_root,
  tcga_dir = "/Users/peng/Library/CloudStorage/OneDrive-Personal/28PHD_peng/250301-DROMA_project/archive260314/251112-DROMA_align/benchmark_mini/Input/TCGA/rna_counts",
  output_base = file.path(workflow_root, "Output"),
  drug_names = c("Paclitaxel"),
  tumor_types = c("breast cancer"),
  es_t = 0.1,
  fdr_t = 0.1,
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
