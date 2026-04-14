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

cat("\n=== Step 4: Stage3 Group Meta ===\n")
within_results <- read_stage("03-within-study", "within_results.rds")
meta_results <- list()
dir.create(file.path(workflow_config$output_base, "04-meta"), recursive = TRUE, showWarnings = FALSE)

for (group_name in names(within_results)) {
  meta_results[[group_name]] <- runGroupedMetaAnalysis(
    within_study_dt = within_results[[group_name]],
    model_group = group_name,
    fdr_t = workflow_config$fdr_t,
    es_t = workflow_config$es_t,
    min_studies = workflow_config$groups[[group_name]]$min_studies
  )

  fwrite(
    meta_results[[group_name]],
    file.path(workflow_config$output_base, "04-meta", paste0(group_name, "_meta.csv"))
  )
}

save_stage(meta_results, "04-meta", "meta_results.rds")
