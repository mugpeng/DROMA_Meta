# Batch driver for DROMA.Meta ----

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DROMA.Meta))
suppressPackageStartupMessages(library(DROMA.Set))
suppressPackageStartupMessages(library(DROMA.R))

# project_root <- normalizePath(getwd(), mustWork = TRUE)
project_root <- normalizePath(file.path(getwd(), "..", "Meta_Example"), mustWork = TRUE)
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

# Split into four nearly equal batches ----
n_batch_parts <- 4L
eligible_pairs_run <- eligible_pairs_run[order(drug, tumor_type)]
eligible_pairs_run[
  ,
  batch_part := rep(
    seq_len(n_batch_parts),
    each = ceiling(.N / n_batch_parts)
  )[seq_len(.N)]
]
eligible_pairs_run_part_1 <- eligible_pairs_run[batch_part == 1L][, batch_part := NULL]
eligible_pairs_run_part_2 <- eligible_pairs_run[batch_part == 2L][, batch_part := NULL]
eligible_pairs_run_part_3 <- eligible_pairs_run[batch_part == 3L][, batch_part := NULL]
eligible_pairs_run_part_4 <- eligible_pairs_run[batch_part == 4L][, batch_part := NULL]

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
cores <- 8

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

batch_part_1_rds <- file.path(output_base_batch, "meta_workflow_batch_part_1.rds")
batch_part_2_rds <- file.path(output_base_batch, "meta_workflow_batch_part_2.rds")
batch_part_3_rds <- file.path(output_base_batch, "meta_workflow_batch_part_3.rds")
batch_part_4_rds <- file.path(output_base_batch, "meta_workflow_batch_part_4.rds")

# Part 1 — copy this block alone to run batch 1 ----
{
  eligible_pairs_run <- eligible_pairs_run_part_1
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
  part_dt <- if (length(batch_results)) {
    data.table::rbindlist(batch_results, fill = TRUE)
  } else {
    data.table::data.table()
  }
  saveRDS(part_dt, batch_part_1_rds)
}

# Part 2 — copy this block alone to run batch 2 ----
{
  eligible_pairs_run <- eligible_pairs_run_part_2
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
  part_dt <- if (length(batch_results)) {
    data.table::rbindlist(batch_results, fill = TRUE)
  } else {
    data.table::data.table()
  }
  saveRDS(part_dt, batch_part_2_rds)
}

# Part 3 — copy this block alone to run batch 3 ----
{
  eligible_pairs_run <- eligible_pairs_run_part_3
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
  part_dt <- if (length(batch_results)) {
    data.table::rbindlist(batch_results, fill = TRUE)
  } else {
    data.table::data.table()
  }
  saveRDS(part_dt, batch_part_3_rds)
}

# Part 4 — copy this block alone to run batch 4 ----
{
  eligible_pairs_run <- eligible_pairs_run_part_4
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
  part_dt <- if (length(batch_results)) {
    data.table::rbindlist(batch_results, fill = TRUE)
  } else {
    data.table::data.table()
  }
  saveRDS(part_dt, batch_part_4_rds)
}

# Merge part RDS files into summary (run after all four parts finished) ----
summary_dt <- data.table::rbindlist(
  list(
    readRDS(batch_part_1_rds),
    readRDS(batch_part_2_rds),
    readRDS(batch_part_3_rds),
    readRDS(batch_part_4_rds)
  ),
  fill = TRUE
)
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
