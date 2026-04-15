# Batch driver for DROMA.Meta ----

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DROMA.Meta))
suppressPackageStartupMessages(library(DROMA.Set))
suppressPackageStartupMessages(library(DROMA.R))

project_root <- normalizePath(getwd(), mustWork = TRUE)
defaults <- getMetaWorkflowDefaults(project_root = project_root)

# Driver-level input and output locations. Keep these runtime choices here so
# the files under R/ remain reusable pure function definitions.
valid_drugs_csv <- file.path(project_root, "workflow", "Output", "valid_drugs.csv")
valid_tumor_types_csv <- file.path(project_root, "workflow", "Output", "valid_tumor_types.csv")
summary_csv <- file.path(project_root, "workflow", "Output", "meta_workflow_batch_summary.csv")

# Batch-loop inputs. Edit these values in the script when you want to change
# how the workflow is run.
drugs <- readSingleColumnCsv(valid_drugs_csv, "drug")
tumor_types <- readSingleColumnCsv(valid_tumor_types_csv, "tumor_type")

db_path <- defaults$droma_db_path
ctrdb_path <- defaults$ctrdb_path
tcga_rna_counts_dir <- defaults$tcga_rna_counts_dir
gene_probe_map_path <- defaults$gene_probe_map_path
output_base <- defaults$output_base

feature2_type <- "mRNA"
data_type <- "all"
cores <- 3

cell_min_intersected_cells <- 20
pdcpdx_min_intersected_cells <- 8

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

results <- list()
result_id <- 1L

# drug <- drugs[1]
# tumor_type <- tumor_types[1]

for (drug in drugs) {
  for (tumor_type in tumor_types) {
    results[[result_id]] <- runMetaWorkflow(
      drug = drug,
      tumor_type = tumor_type,
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
      output_base = output_base,
      verbose = TRUE
    )
    result_id <- result_id + 1L
  }
}

summary_dt <- data.table::rbindlist(results, fill = TRUE)
data.table::fwrite(summary_dt, summary_csv)

print(summary_dt)
