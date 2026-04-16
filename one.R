# Batch driver for DROMA.Meta ----

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DROMA.Meta))
suppressPackageStartupMessages(library(DROMA.Set))
suppressPackageStartupMessages(library(DROMA.R))

project_root <- file.path(normalizePath(getwd(), mustWork = TRUE), "Meta_Example")
# project_root <- normalizePath(getwd(), mustWork = TRUE)
defaults <- getMetaWorkflowDefaults(project_root = project_root)

# Driver-level input and output locations. Keep these runtime choices here so
# the files under R/ remain reusable pure function definitions.
eligible_pairs_csv <- file.path(defaults$output_base, "eligible_drug_tumor_pairs.csv")

output_base_batch <- file.path(defaults$output_base, "meta_batch")
dir.create(output_base_batch, recursive = TRUE, showWarnings = FALSE)

summary_csv <- file.path(output_base_batch, "meta_workflow_batch_summary.csv")

# Batch-loop inputs. Edit these values in the script when you want to change
# how the workflow is run.
eligible_pairs <- data.table::fread(eligible_pairs_csv, data.table = TRUE)
eligible_pairs <- unique(eligible_pairs[, .(drug, tumor_type)])

db_path <- defaults$droma_db_path
# db_path <- "/home/data/denglab/bigData/DROMA/droma.sqlite"
ctrdb_path <- defaults$ctrdb_path
# ctrdb_path <- "/home/data/denglab/bigData/DROMA/ctrdb.sqlite"

tcga_rna_counts_dir <- defaults$tcga_rna_counts_dir
# tcga_rna_counts_dir <- "/home/data/denglab/bigData/DROMA/rna_counts"

gene_probe_map_path <- defaults$gene_probe_map_path
# gene_probe_map_path <- "/Users/peng/Desktop/Project/DROMA/Data/gencode.human.v49.annotation.gene.probeMap"
output_base <- defaults$output_base

feature2_type <- "mRNA"
data_type <- "all"
cores <- 4

cell_min_intersected_cells <- 10
pdcpdx_min_intersected_cells <- 5

cell_es_t <- 0.1
cell_P_t <- 0.05
cell_n_datasets_t <- 3

pdcpdx_es_t <- 0.1
pdcpdx_P_t <- 0.05
pdcpdx_n_datasets_t <- 1

tcga_ad_p_t <- 0.01

clinical_es_t <- 0.1
clinical_P_t <- 0.05
clinical_n_datasets_t <- NULL

batch_results <- list()
for (i in seq_len(nrow(eligible_pairs))) {
  pair <- eligible_pairs[i]

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

summary_dt <- data.table::rbindlist(batch_results, fill = TRUE)
data.table::fwrite(summary_dt, summary_csv)
print(summary_dt)
