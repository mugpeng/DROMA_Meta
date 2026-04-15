# Run with getwd() == project root (parent of workflow/ and DROMA_R2/).
# Override main DROMA SQLite: set absolute path, or "" for ../Data/droma.sqlite (relative to project root).
droma_sqlite_path <- ""
# droma_sqlite_path <- "/home/data/denglab/bigData/DROMA/droma.sqlite"

# Override CTRDB SQLite: set absolute path, or "" for ../Data/ctrdb.sqlite.
ctrdb_sqlite_path <- ""
# ctrdb_sqlite_path <- "/home/data/denglab/bigData/DROMA/ctrdb.sqlite"

# TCGA RNA counts directory: absolute path, or "" for ../Data/TCGA/rna_counts.
tcga_rna_counts_dir <- "/Users/peng/Library/CloudStorage/OneDrive-Personal/28PHD_peng/250301-DROMA_project/archive260314/251112-DROMA_align/benchmark_mini/Input/TCGA/"
# tcga_rna_counts_dir <- "/home/data/denglab/bigData/DROMA/rna_counts"

workflow_root <- file.path(".", "workflow")
meta_project_root <- normalizePath(".", mustWork = TRUE)
repo_root <- meta_project_root
droma_r2_root <- file.path(".", "DROMA_R2")

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
    file.path(repo_root, "..", "Data", "droma.sqlite")
  },
  ctrdb_path = if (nzchar(ctrdb_sqlite_path)) {
    ctrdb_sqlite_path
  } else {
    file.path(repo_root, "..", "Data", "ctrdb.sqlite")
  },
  tcga_dir = if (nzchar(tcga_rna_counts_dir)) {
    tcga_rna_counts_dir
  } else {
    file.path(repo_root, "..", "Data", "TCGA", "rna_counts")
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

# One `{ ... }` block so `on.exit` is tied to a single eval frame under `source()`:
# otherwise each top-level line is a separate `eval()`, and `on.exit` runs as soon
# as the `on.exit(...)` line finishes, closing `con` before `listDROMAProjects`.
{
  cat("\n=== 00: Project grouping ===\n")
  con <- DROMA.Set::connectDROMADatabase(workflow_config$db_path)
  on.exit(DROMA.Set::closeDROMADatabase(con), add = TRUE)
  dir.create(file.path(workflow_config$output_base, "01-projects"), recursive = TRUE, showWarnings = FALSE)

  project_anno <- listDROMAProjects(connection = con)
  project_groups <- createWorkflowProjectGroups(project_anno)
  group_sets <- lapply(project_groups, function(project_names) {
    DROMA.Set::createMultiDromaSetFromAllProjects(
      db_path = workflow_config$db_path,
      include_projects = sort(unique(project_names)),
      con = con
    )
  })
  shared_features <- getSharedGroupFeatures(group_sets, workflow_config$feature_type)

  save_stage(project_anno, "01-projects", "project_anno.rds")
  save_stage(project_groups, "01-projects", "project_groups.rds")
  save_stage(group_sets, "01-projects", "group_sets.rds")
  save_stage(shared_features, "01-projects", "shared_features.rds")

  fwrite(project_anno, file.path(workflow_config$output_base, "01-projects", "project_anno.csv"))
  for (group_name in names(project_groups)) {
    fwrite(
      data.table(group = group_name, project_name = project_groups[[group_name]]),
      file.path(workflow_config$output_base, "01-projects", paste0(group_name, "_projects.csv"))
    )
  }
  fwrite(
    data.table(feature = shared_features),
    file.path(workflow_config$output_base, "01-projects", "shared_features.csv")
  )
}
