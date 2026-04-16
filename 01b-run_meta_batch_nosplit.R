# Batch driver for DROMA.Meta ----

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DROMA.Meta))
suppressPackageStartupMessages(library(DROMA.Set))
suppressPackageStartupMessages(library(DROMA.R))

# project_root <- normalizePath(getwd(), mustWork = TRUE)
project_root <- normalizePath(file.path(getwd(), "Meta_Example"), mustWork = TRUE)
defaults <- getMetaWorkflowDefaults(project_root = project_root)

# Driver-level input and output locations. Keep these runtime choices here so
# the files under R/ remain reusable pure function definitions.
eligible_pairs_csv <- file.path(defaults$output_base, "eligible_drug_pairs_filter.csv")

output_base_batch <- file.path(defaults$output_base, "meta_batch")
dir.create(output_base_batch, recursive = TRUE, showWarnings = FALSE)

summary_csv <- file.path(output_base_batch, "meta_workflow_batch_summary.csv")

# Batch-loop inputs. Edit these values in the script when you want to change
# how the workflow is run.
eligible_pairs <- data.table::fread(eligible_pairs_csv, data.table = TRUE)

# suppressPackageStartupMessages(library(DROMA.Meta.visualization))
# project_root <- "/Users/peng/Desktop/Project/DROMA/Meta_project/Meta_Example"
# defaults   <- getVisWorkflowDefaults(project_root = project_root)
# vis_output <- defaults$vis_output
# summary_path <- file.path(vis_output, "collected_summary.csv")
# summary_dt      <- fread(summary_path)
# eligible_pairs_bk <- eligible_pairs
# summary_dt$pair <- paste0(summary_dt$drug, "_", summary_dt$tumor_type)
# eligible_pairs$pair <- paste0(eligible_pairs$drug, "_", eligible_pairs$tumor_type)
# eligible_pairs[, pair := gsub("([a-zA-Z0-9'-]+_)([a-z]+) ([a-z]+)$", "\\1\\2_\\3", pair)]
# eligible_pairs <- eligible_pairs[eligible_pairs$pair %in% summary_dt$pair,]
# defaults <- getMetaWorkflowDefaults(project_root = project_root)

eligible_pairs
eligible_pairs <- unique(eligible_pairs)
if ("n_datasets_pdcpdx_invitro" %in% names(eligible_pairs)) {
  eligible_pairs[
    ,
    pdcpdx_ge_2_pair := !is.na(n_datasets_pdcpdx_invitro) & n_datasets_pdcpdx_invitro >= 2
  ]
} else {
  eligible_pairs[, pdcpdx_ge_2_pair := NA]
}
eligible_pair_flags <- eligible_pairs[
  ,
  .(
    pdcpdx_ge_2_pair = {
      values <- pdcpdx_ge_2_pair[!is.na(pdcpdx_ge_2_pair)]
      if (length(values)) any(values) else NA
    }
  ),
  by = .(drug, tumor_type)
]
eligible_pairs_run <- unique(eligible_pairs[, .(drug, tumor_type)])

# Run all eligible pairs directly without splitting ----
eligible_pairs_run <- eligible_pairs_run[order(drug, tumor_type)]

# db_path <- defaults$droma_db_path
# ctrdb_path <- defaults$ctrdb_path
# tcga_rna_counts_dir <- defaults$tcga_rna_counts_dir
# gene_probe_map_path <- defaults$gene_probe_map_path

db_path <- "/home/data/denglab/bigData/DROMA/droma.sqlite"
ctrdb_path <- "/home/data/denglab/bigData/DROMA/ctrdb.sqlite"
tcga_rna_counts_dir <- "/home/data/denglab/bigData/DROMA/rna_counts"
gene_probe_map_path <- "/home/data/denglab/bigData/DROMA/gencode.human.v49.annotation.gene.probeMap"

output_base <- defaults$output_base

feature2_type <- "mRNA"
data_type <- "all"
cores <- 3

cell_min_intersected_cells <- 10
pdcpdx_min_intersected_cells <- 5

cell_es_t <- 0.1
cell_P_t <- 0.05
cell_n_datasets_t <- 3

pdcpdx_es_t <- 0.1
pdcpdx_P_t <- 0.1
pdcpdx_n_datasets_t <- 1

tcga_ad_p_t <- 0.01

clinical_es_t <- 0.1
clinical_P_t <- 0.1
clinical_n_datasets_t <- NULL

batch_results <- list()
for (i in seq_len(nrow(eligible_pairs_run))) {
  pair <- eligible_pairs_run[i]
  batch_results[[i]] <- runMetaWorkflow(
    drug = pair$drug[[1]],
    tumor_type = pair$tumor_type[[1]],
    feature2_type = feature2_type,
    data_type = data_type,
    cores = cores,
    cell_min_intersected_cells = cell_min_intersected_cells,
    pdcpdx_min_intersected_cells = pdcpdx_min_intersected_cells,
    cell_es_t = cell_es_t,
    cell_P_t = cell_P_t,
    cell_n_datasets_t = cell_n_datasets_t,
    pdcpdx_es_t = pdcpdx_es_t,
    pdcpdx_P_t = pdcpdx_P_t,
    pdcpdx_n_datasets_t = pdcpdx_n_datasets_t,
    tcga_ad_p_t = tcga_ad_p_t,
    clinical_es_t = clinical_es_t,
    clinical_P_t = clinical_P_t,
    clinical_n_datasets_t = clinical_n_datasets_t,
    db_path = db_path,
    ctrdb_path = ctrdb_path,
    tcga_rna_counts_dir = tcga_rna_counts_dir,
    gene_probe_map_path = gene_probe_map_path,
    output_base = output_base_batch,
    override = F,
    verbose = TRUE
  )
}

summary_dt <- if (length(batch_results)) {
  data.table::rbindlist(batch_results, fill = TRUE)
} else {
  data.table::data.table()
}
summary_dt <- merge(
  summary_dt,
  eligible_pair_flags,
  by = c("drug", "tumor_type"),
  all.x = TRUE
)
data.table::setcolorder(
  summary_dt,
  c(
    "drug",
    "tumor_type",
    "pdcpdx_ge_2_pair",
    "ctrdb_status",
    setdiff(names(summary_dt), c("drug", "tumor_type", "pdcpdx_ge_2_pair", "ctrdb_status"))
  )
)
data.table::fwrite(summary_dt, summary_csv)
print(summary_dt)
