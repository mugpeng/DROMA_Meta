# ============================================================================
# 04-Clinical_Validation.R
# Clinical validation on CTRDB using TCGA-AD-filtered mRNA biomarkers
# ============================================================================

library(data.table)
suppressPackageStartupMessages(library(DROMA.Meta))
library(DROMA.Set)
library(DROMA.R)

ctrdb_path <- "/Users/peng/Desktop/Project/DROMA/Data/ctrdb.sqlite"
project_root <- normalizePath(file.path(getwd(), "..", "Meta_Example"), mustWork = TRUE)
defaults <- getMetaWorkflowDefaults(project_root = project_root)
drug <- "Paclitaxel"
tumor_type <- "breast cancer"
data_type <- "all"
cores <- 3

clinical_es_t <- 0.1
clinical_P_t <- 0.05
clinical_n_datasets_t <- NULL

output_base <- defaults$output_base
tumor_type_slug <- sanitizeName(tumor_type)
output_dir <- file.path(output_base, drug, tumor_type_slug)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
clinical_info_path <- file.path(output_dir, "clinical_validation_info.csv")

cat("\n=== 04: Clinical Validation ===\n")

# Need DROMA database connection for annotation fetching used internally
droma_path <- "/Users/peng/Desktop/Project/DROMA/Data/DROMA.sqlite"
connectDROMADatabase(droma_path)

selected_genes_ad_filtered <- fread(file.path(output_dir, "selected_genes_ad_filtered.csv"))

connectCTRDatabase(ctrdb_path)

fetchClinicalBatch <- function(query_tumor_type) {
  tryCatch(
    batchFindClinicalSigResponse(
      select_omics = unique(selected_genes_ad_filtered$name),
      select_drugs = drug,
      data_type = data_type,
      tumor_type = query_tumor_type,
      cores = cores
    ),
    error = function(e) {
      createEmptyMetaDf()
    }
  )
}

clinical_batch <- fetchClinicalBatch(tumor_type)
ctrdb_status <- "tumor_type_matched"
ctrdb_fallback <- FALSE
clinical_query_tumor_type <- tumor_type

if (nrow(clinical_batch) == 0L) {
  clinical_batch <- fetchClinicalBatch("all")
  clinical_query_tumor_type <- "all"

  if (nrow(clinical_batch) > 0L) {
    ctrdb_status <- "fallback_all"
    ctrdb_fallback <- TRUE
  } else {
    ctrdb_status <- "no_ctrdb_data"
    ctrdb_fallback <- FALSE
  }
}

if (ncol(clinical_batch) == 0L) {
  clinical_batch <- createEmptyMetaDf()
}

if (nrow(clinical_batch) > 0) {
  clinical_sig <- getSignificantFeatures(
    clinical_batch,
    es_t = clinical_es_t,
    P_t = clinical_P_t,
    use_p_value = TRUE,
    n_datasets_t = clinical_n_datasets_t
  )
} else {
  clinical_sig <- createEmptyMetaDf()
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

    final_biomarkers$ctrdb_fallback <- ctrdb_fallback
    final_biomarkers$ctrdb_status <- ctrdb_status
  } else {
    final_biomarkers <- createEmptyFinalBiomarkersDf()
  }
} else {
  final_biomarkers <- createEmptyFinalBiomarkersDf()
}

clinical_info <- data.table::data.table(
  drug = drug,
  tumor_type = tumor_type,
  ctrdb_status = ctrdb_status,
  ctrdb_fallback = ctrdb_fallback,
  clinical_query_tumor_type = clinical_query_tumor_type,
  n_clinical_batch = nrow(clinical_batch),
  n_clinical_sig = nrow(clinical_sig),
  n_final_biomarkers = nrow(final_biomarkers)
)

fwrite(clinical_batch, file.path(output_dir, "clinical_batch.csv"))
fwrite(clinical_sig, file.path(output_dir, "clinical_sig_mRNA.csv"))
fwrite(final_biomarkers, file.path(output_dir, "final_biomarkers.csv"))
fwrite(clinical_info, clinical_info_path)

cat(sprintf("  OK clinical_batch: %d biomarkers\n", nrow(clinical_batch)))
cat(sprintf("  OK clinical_sig: %d biomarkers\n", nrow(clinical_sig)))
cat(sprintf("  OK final_biomarkers: %d biomarkers\n", nrow(final_biomarkers)))
