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

# Gene symbol -> Ensembl id map for TCGA row lookup (TSV: id, gene, ...). "" -> ../Data/gencode.v22.annotation.gene.probeMap
gene_probe_map_path <- "/Users/peng/Desktop/Project/DROMA/Data/gencode.human.v49.annotation.gene.probeMap"
# gene_probe_map_path <- "/Users/peng/Desktop/Project/DROMA/Data/gencode.v22.annotation.gene.probeMap"

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
  gene_probe_map = if (nzchar(gene_probe_map_path)) {
    gene_probe_map_path
  } else {
    file.path(repo_root, "..", "Data", "gencode.v22.annotation.gene.probeMap")
  },
  output_base = file.path(workflow_root, "Output"),
  drug_names = "Paclitaxel",
  tumor_types = "breast cancer",
  feature_type = "mRNA",
  # Meta: |effect_size| >= es_t, and either p_value < meta_p_t (meta_use_p_value TRUE) or q_value < fdr_t (FALSE)
  es_t = 0.1,
  fdr_t = 0.1,
  meta_p_t = 0.05,
  meta_use_p_value = TRUE,
  # TCGA translation: raw AD/KS p > tcga_p_t => supported
  tcga_p_t = 0.1,
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

# Filled by workflow/00b-Project_Grouping.R in-memory (no required 01-projects files).
.droma_wf_01_projects_rds <- c(
  project_anno = "project_anno.rds",
  project_groups = "project_groups.rds",
  group_sets = "group_sets.rds",
  shared_features = "shared_features.rds"
)

read_project_grouping <- function(name) {
  if (!name %in% names(.droma_wf_01_projects_rds)) {
    stop("read_project_grouping: unknown name: ", name, call. = FALSE)
  }
  cache <- get0(".droma_wf_01_projects", envir = .GlobalEnv, ifnotfound = NULL)
  if (!is.null(cache) && !is.null(cache[[name]])) {
    return(cache[[name]])
  }
  read_stage("01-projects", .droma_wf_01_projects_rds[[name]])
}
