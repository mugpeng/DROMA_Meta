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
eligible_pairs_pdcpdx_ge_2_csv <- file.path(defaults$output_base, "eligible_drug_tumor_pairs_pdcpdx_ge_2.csv")
eligible_pairs_pdcpdx_eq_1_csv <- file.path(defaults$output_base, "eligible_drug_tumor_pairs_pdcpdx_eq_1.csv")

output_base_1 <- file.path(defaults$output_base, "meta_batch_pdcpdx_ge_2")
output_base_2 <- file.path(defaults$output_base, "meta_batch_pdcpdx_eq_1")
dir.create(output_base_1, recursive = TRUE, showWarnings = FALSE)
dir.create(output_base_2, recursive = TRUE, showWarnings = FALSE)

summary_csv_1 <- file.path(output_base_1, "meta_workflow_batch_summary.csv")
summary_csv_2 <- file.path(output_base_2, "meta_workflow_batch_summary.csv")

# Batch-loop inputs. Edit these values in the script when you want to change
# how the workflow is run.
eligible_pairs_1 <- data.table::fread(eligible_pairs_pdcpdx_ge_2_csv, data.table = TRUE)
eligible_pairs_1 <- unique(eligible_pairs_1[, .(drug, tumor_type)])

eligible_pairs_2 <- data.table::fread(eligible_pairs_pdcpdx_eq_1_csv, data.table = TRUE)
eligible_pairs_2 <- unique(eligible_pairs_2[, .(drug, tumor_type)])

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
pdcpdx_n_datasets_t <- 2

tcga_ad_p_t <- 0.01

clinical_es_t <- 0.1
clinical_P_t <- 0.05
clinical_n_datasets_t <- NULL

for (batch in list(
  list(pairs = eligible_pairs_1, output_base = output_base_1, summary_csv = summary_csv_1)
  # list(pairs = eligible_pairs_2, output_base = output_base_2, summary_csv = summary_csv_2)
)) {
  batch_results <- list()
  for (i in seq_len(nrow(batch$pairs))) {
    pair <- batch$pairs[i]

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
      output_base = batch$output_base,
      override = F,
      verbose = TRUE
    )
  }

  summary_dt <- data.table::rbindlist(batch_results, fill = TRUE)
  data.table::fwrite(summary_dt, batch$summary_csv)
  print(summary_dt)
}
