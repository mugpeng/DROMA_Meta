# ============================================================================
# 04-Clinical_Validation.R
# Clinical validation on CTRDB using TCGA-AD-filtered mRNA biomarkers
# ============================================================================

library(data.table)
library(DROMA.Set)
library(DROMA.R)
source("/Users/peng/Desktop/Project/DROMA/Meta_project3/R/FuncHelper.R", local = FALSE)

ctrdb_path <- "/Users/peng/Desktop/Project/DROMA/Data/ctrdb.sqlite"
drug <- "Paclitaxel"
tumor_type <- "breast cancer"
data_type <- "all"
cores <- 3

clinical_es_t <- 0.05
clinical_P_t <- 0.1
clinical_n_datasets_t <- NULL

output_base <- "/Users/peng/Desktop/Project/DROMA/Meta_project3/workflow/Output"
tumor_type_slug <- sanitizeName(tumor_type)
output_dir <- file.path(output_base, drug, tumor_type_slug)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 04: Clinical Validation ===\n")

selected_genes_ad_filtered <- fread(file.path(output_dir, "selected_genes_ad_filtered.csv"))

connectCTRDatabase(ctrdb_path)

clinical_batch <- tryCatch(
  batchFindClinicalSigResponse(
    select_omics = unique(selected_genes_ad_filtered$name),
    select_drugs = drug,
    data_type = data_type,
    tumor_type = tumor_type,
    cores = cores
  ),
  error = function(e) {
    warning("Clinical validation returned no usable result: ", e$message, call. = FALSE)
    createEmptyMetaDf()
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

if (nrow(selected_genes_ad_filtered) > 0 && nrow(clinical_sig) > 0) {
  final_biomarkers <- getIntersectSignificantFeatures(
    pdcpdx = as.data.frame(selected_genes_ad_filtered),
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

cat(sprintf("  OK clinical_sig: %d biomarkers\n", nrow(clinical_sig)))
cat(sprintf("  OK final_biomarkers: %d biomarkers\n", nrow(final_biomarkers)))
