# ============================================================================
# 03-Clinical_Validation.R
# Clinical validation on CTRDB using pdcpdx significant mRNA biomarkers
# ============================================================================

library(data.table)
library(DROMA.Set)
library(DROMA.R)

source("/Users/peng/Desktop/Project/DROMA/Meta_project3/R/run_drug_tumor_biomarker_workflow.R", local = FALSE)

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
mRNA_pdcpdx_sig <- fread(file.path(output_dir, "mRNA_pdcpdx_sig.csv"))

connectCTRDatabase(ctrdb_path)

clinical_batch <- tryCatch(
  batchFindClinicalSigResponse(
    select_omics = unique(mRNA_pdcpdx_sig$name),
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

if (nrow(mRNA_pdcpdx_sig) > 0 && nrow(clinical_sig) > 0) {
  mRNA_pdcpdx_sig2 <- as.data.frame(mRNA_pdcpdx_sig)
  names(mRNA_pdcpdx_sig2)[names(mRNA_pdcpdx_sig2) == "direction"] <- "direction_pdcpdx"

  final_biomarkers <- getIntersectSignificantFeatures(
    pdcpdx = mRNA_pdcpdx_sig2,
    ctrdb = as.data.frame(clinical_sig),
    direction_cols = c(pdcpdx = "direction_pdcpdx", ctrdb = "direction")
  )

  if (nrow(final_biomarkers) > 0) {
    final_biomarkers$drug <- drug
    final_biomarkers$tumor_type <- tumor_type
    final_biomarkers$cell_supported <- final_biomarkers$name %in% mRNA_cell_sig$name
  }
} else {
  final_biomarkers <- data.frame(name = character(0), stringsAsFactors = FALSE)
}

fwrite(clinical_batch, file.path(output_dir, "clinical_batch_mRNA.csv"))
fwrite(clinical_sig, file.path(output_dir, "clinical_sig_mRNA.csv"))
fwrite(final_biomarkers, file.path(output_dir, "final_biomarkers.csv"))
saveRDS(clinical_batch, file.path(output_dir, "clinical_batch_mRNA.rds"))
saveRDS(clinical_sig, file.path(output_dir, "clinical_sig_mRNA.rds"))
saveRDS(final_biomarkers, file.path(output_dir, "final_biomarkers.rds"))

cat(sprintf("  OK clinical_sig: %d biomarkers\n", nrow(clinical_sig)))
cat(sprintf("  OK final_biomarkers: %d biomarkers\n", nrow(final_biomarkers)))
