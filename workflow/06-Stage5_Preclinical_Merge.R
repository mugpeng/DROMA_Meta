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

cat("\n=== Step 6: Stage5 Preclinical Merge ===\n")
tcga_results <- read_stage("05-tcga-filter", "tcga_results.rds")
dir.create(file.path(workflow_config$output_base, "06-preclinical-merge"), recursive = TRUE, showWarnings = FALSE)

invitro_post <- tcga_results$invitro[is.na(tcga_supported) | tcga_supported == TRUE]
invivo_post <- tcga_results$invivo[is.na(tcga_supported) | tcga_supported == TRUE]

merged <- mergePreclinicalCandidates(
  invitro_meta = invitro_post,
  invivo_meta = invivo_post,
  fdr_t = workflow_config$fdr_t,
  es_t = workflow_config$es_t
)

save_stage(merged, "06-preclinical-merge", "merged_candidates.rds")
fwrite(merged, file.path(workflow_config$output_base, "06-preclinical-merge", "merged_candidates.csv"))
