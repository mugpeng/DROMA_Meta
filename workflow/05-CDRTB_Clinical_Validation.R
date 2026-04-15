source(file.path(".", "workflow", "00-Setup.R"), local = FALSE)

cat("\n=== 05: CDRTB clinical validation ===\n")
dir.create(file.path(workflow_config$output_base, "06-clinical-validation"), recursive = TRUE, showWarnings = FALSE)

merged_candidates <- read_stage("05-preclinical-merge", "merged_candidates.rds")
clinical_rows <- list()

if (nrow(merged_candidates) > 0) {
  split_candidates <- split(merged_candidates, interaction(merged_candidates$drug, merged_candidates$tumor_type, drop = TRUE))
  for (candidate_df in split_candidates) {
    re <- runClinicalValidationForCandidates(
      candidate_df = candidate_df,
      drug_name = candidate_df$drug[1],
      tumor_type = candidate_df$tumor_type[1],
      ctrdb_path = workflow_config$ctrdb_path,
      cores = .get_safe_cores(workflow_config$requested_cores),
      es_t = workflow_config$es_t,
      fdr_t = workflow_config$fdr_t
    )
    clinical_rows[[length(clinical_rows) + 1L]] <- re
  }
}

clinical_results <- if (length(clinical_rows) > 0) rbindlist(clinical_rows, fill = TRUE) else merged_candidates

save_stage(clinical_results, "06-clinical-validation", "clinical_results.rds")
fwrite(clinical_results, file.path(workflow_config$output_base, "06-clinical-validation", "clinical_results.csv"))
