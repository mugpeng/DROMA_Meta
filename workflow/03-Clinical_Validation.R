# ============================================================================
# 03-Clinical_Validation.R
# Clinical validation on CTRDB using preclinical-selected mRNA biomarkers
# ============================================================================

library(data.table)
library(DROMA.Set)
library(DROMA.R)

ctrdb_path <- "/Users/peng/Desktop/Project/DROMA/Data/ctrdb.sqlite"
drug <- "Paclitaxel"
tumor_type <- "breast cancer"
data_type <- "all"
cores <- 3

clinical_es_t <- 0.05
clinical_P_t <- 0.1
clinical_n_datasets_t <- NULL

output_base <- "/Users/peng/Desktop/Project/DROMA/Meta_project3/workflow/Output"
tumor_type_slug <- .sanitize_name(tumor_type)
output_dir <- file.path(output_base, drug, tumor_type_slug)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 03: Clinical Validation ===\n")

mRNA_cell_sig <- fread(file.path(output_dir, "mRNA_cell_sig.csv"))
selected_genes <- fread(file.path(output_dir, "selected_genes.csv"))

connectCTRDatabase(ctrdb_path)

clinical_batch <- tryCatch(
  batchFindClinicalSigResponse(
    select_omics = unique(selected_genes$name),
    select_drugs = drug,
    data_type = data_type,
    tumor_type = tumor_type,
    cores = cores
  ),
  error = function(e) {
    warning("Clinical validation returned no usable result: ", e$message, call. = FALSE)
    .empty_meta_df()
  }
)

if (nrow(clinical_batch) > 0) {
  clinical_sig <- getSignificantFeatures(
    clinical_batch,
    es_t = clinical_es_t,
    P_t = clinical_P_t,
    use_p_value = TRUE,
    n_datasets_t = clinical_n_datasets_t
  )
} else {
  clinical_sig <- data.frame(name = character(0), stringsAsFactors = FALSE)
}

if (nrow(selected_genes) > 0 && nrow(clinical_sig) > 0) {
  final_biomarkers <- getIntersectSignificantFeatures(
    preclinical = as.data.frame(selected_genes),
    ctrdb = as.data.frame(clinical_sig),
    direction_cols = c(preclinical = "direction_pdcpdx", ctrdb = "direction")
  )

  if (nrow(final_biomarkers) > 0) {
    final_biomarkers$drug <- drug
    final_biomarkers$tumor_type <- tumor_type
    final_biomarkers$direction <- coalesce(final_biomarkers$direction_ctrdb, final_biomarkers$direction_pdcpdx_preclinical, final_biomarkers$direction_cell_preclinical)
    final_biomarkers$direction_ctrdb <- NULL
    final_biomarkers$direction_pdcpdx_preclinical <- NULL
    final_biomarkers$direction_cell_preclinical <- NULL
  }
} else {
  final_biomarkers <- data.frame(name = character(0), stringsAsFactors = FALSE)
}

fwrite(clinical_batch, file.path(output_dir, "clinical_batch_mRNA.csv"))
fwrite(clinical_sig, file.path(output_dir, "clinical_sig_mRNA.csv"))
fwrite(final_biomarkers, file.path(output_dir, "final_biomarkers.csv"))

cat(sprintf("  OK clinical_sig: %d biomarkers\n", nrow(clinical_sig)))
cat(sprintf("  OK final_biomarkers: %d biomarkers\n", nrow(final_biomarkers)))
