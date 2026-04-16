# ============================================================================
# 04-Clinical_Validation.R
# Clinical validation on CTRDB using TCGA-AD-filtered mRNA biomarkers
# ============================================================================

library(data.table)
suppressPackageStartupMessages(library(DROMA.Meta))
library(DROMA.Set)
library(DROMA.R)

ctrdb_path <- "/Users/peng/Desktop/Project/DROMA/Data/ctrdb.sqlite"
drug <- "Paclitaxel"
tumor_type <- "breast cancer"
data_type <- "all"
cores <- 3

clinical_es_t <- 0.1
clinical_P_t <- 0.05
clinical_n_datasets_t <- NULL

output_base <- "/Users/peng/Desktop/Project/DROMA/Meta_project3/workflow/Output"
tumor_type_slug <- sanitizeName(tumor_type)
output_dir <- file.path(output_base, drug, tumor_type_slug)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 04: Clinical Validation ===\n")

# Need DROMA database connection for annotation fetching used internally
droma_path <- "/Users/peng/Desktop/Project/DROMA/Data/DROMA.sqlite"
connectDROMADatabase(droma_path)

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
    invitro = as.data.frame(selected_genes_ad_filtered),
    ctrdb = as.data.frame(clinical_sig),
    direction_cols = c(invitro = "direction_pdcpdx", ctrdb = "direction")
  )

  if (nrow(final_biomarkers) > 0) {
    final_biomarkers$drug <- drug
    final_biomarkers$tumor_type <- tumor_type

    # Consolidate direction columns
    dir_cols <- grep("direction", names(final_biomarkers), value = TRUE)
    if (length(dir_cols) > 0) {
      final_biomarkers$direction <- final_biomarkers[[tail(dir_cols, 1)]]
      # Remove other direction columns
      cols_to_remove <- setdiff(dir_cols, "direction")
      if (length(cols_to_remove) > 0) {
        final_biomarkers[, cols_to_remove] <- NULL
      }
    }
  }
} else {
  final_biomarkers <- data.frame(name = character(0), stringsAsFactors = FALSE)
}

# fwrite(clinical_batch, file.path(output_dir, "clinical_batch_mRNA.csv"))
fwrite(clinical_sig, file.path(output_dir, "clinical_sig_mRNA.csv"))
fwrite(final_biomarkers, file.path(output_dir, "final_biomarkers.csv"))

cat(sprintf("  OK clinical_sig: %d biomarkers\n", nrow(clinical_sig)))
cat(sprintf("  OK final_biomarkers: %d biomarkers\n", nrow(final_biomarkers)))
