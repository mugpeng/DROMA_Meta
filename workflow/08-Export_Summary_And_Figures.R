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

cat("\n=== Step 8: Export Summary ===\n")
meta_results <- read_stage("04-meta", "meta_results.rds")
tcga_results <- read_stage("05-tcga-filter", "tcga_results.rds")
merged <- read_stage("06-preclinical-merge", "merged_candidates.rds")
clinical_results <- read_stage("07-clinical-validation", "clinical_results.rds")
dir.create(file.path(workflow_config$output_base, "08-summary"), recursive = TRUE, showWarnings = FALSE)

summary_dt <- rbindlist(list(
  data.table(stage = "meta_invitro", n = nrow(meta_results$invitro)),
  data.table(stage = "meta_invivo", n = nrow(meta_results$invivo)),
  data.table(stage = "tcga_invitro", n = nrow(tcga_results$invitro)),
  data.table(stage = "tcga_invivo", n = nrow(tcga_results$invivo)),
  data.table(stage = "merged_preclinical", n = nrow(merged)),
  data.table(stage = "clinical_validation", n = sum(vapply(clinical_results, nrow, integer(1))))
), fill = TRUE)

fwrite(summary_dt, file.path(workflow_config$output_base, "08-summary", "run_summary.csv"))
fwrite(merged, file.path(workflow_config$output_base, "08-summary", "final_preclinical_candidates.csv"))

if (length(clinical_results) > 0) {
  clinical_merged <- rbindlist(clinical_results, fill = TRUE, idcol = "drug")
  fwrite(clinical_merged, file.path(workflow_config$output_base, "08-summary", "final_clinical_validation.csv"))
  fwrite(
    clinical_merged[clinical_supported == TRUE],
    file.path(workflow_config$output_base, "08-summary", "final_clinically_supported_candidates.csv")
  )
  fwrite(
    clinical_merged[retained == TRUE],
    file.path(workflow_config$output_base, "08-summary", "final_retained_candidates.csv")
  )
}
