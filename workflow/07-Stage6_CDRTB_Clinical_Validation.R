workflow_bootstrap_ok <- FALSE
for (d in c(file.path(getwd(), "workflow"), getwd())) {
  boot <- file.path(d, "00-resolve_workflow_dir.R")
  if (file.exists(boot)) {
    source(boot, local = FALSE)
    workflow_bootstrap_ok <- TRUE
    break
  }
}
if (!workflow_bootstrap_ok) {
  stop("Could not locate workflow/00-resolve_workflow_dir.R; setwd() to Meta_project or workflow.", call. = FALSE)
}
script_dir <- resolve_workflow_script_dir()
source(file.path(script_dir, "00-Workflow_Common.R"), local = FALSE)

cat("\n=== Step 7: Stage6 CDRTB Clinical Validation ===\n")
merged <- read_stage("06-preclinical-merge", "merged_candidates.rds")
clinical_results <- list()
dir.create(file.path(workflow_config$output_base, "07-clinical-validation"), recursive = TRUE, showWarnings = FALSE)

for (drug_name in workflow_config$drug_names) {
  candidate_df <- merged
  clinical_results[[drug_name]] <- runClinicalValidation(
    candidate_df = candidate_df,
    drug_name = drug_name,
    tumor_type = "all",
    db_path = workflow_config$db_path,
    ctrdb_path = workflow_config$ctrdb_path,
    cores = .get_safe_cores(workflow_config$requested_cores),
    es_t = workflow_config$es_t,
    fdr_t = workflow_config$fdr_t
  )

  fwrite(
    clinical_results[[drug_name]],
    file.path(workflow_config$output_base, "07-clinical-validation", paste0(.sanitize_name(drug_name), "_clinical.csv"))
  )
}

save_stage(clinical_results, "07-clinical-validation", "clinical_results.rds")
