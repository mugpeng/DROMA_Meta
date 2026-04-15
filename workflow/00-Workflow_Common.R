droma_sqlite_path <- ""

get_project_root <- function() {
  wd <- normalizePath(getwd(), mustWork = TRUE)
  if (basename(wd) == "workflow") {
    project_root <- dirname(wd)
  } else {
    project_root <- wd
  }
  if (!file.exists(file.path(project_root, "workflow", "00-Workflow_Common.R"))) {
    stop("Run from Meta_project2 root or Meta_project2/workflow")
  }
  project_root
}

meta_project_root <- get_project_root()
workflow_root <- file.path(meta_project_root, "workflow")
repo_root <- normalizePath(dirname(meta_project_root), mustWork = TRUE)
droma_r2_root <- file.path(meta_project_root, "DROMA_R2")

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
  db_path = if (nzchar(droma_sqlite_path)) {
    droma_sqlite_path
  } else {
    file.path(repo_root, "Data", "droma.sqlite")
  },
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
