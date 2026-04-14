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

cat("\n=== Step 5: Stage4 TCGA Translation Filter ===\n")
group_sets <- read_stage("01-projects", "group_sets.rds")
meta_results <- read_stage("04-meta", "meta_results.rds")
tcga_results <- list()
dir.create(file.path(workflow_config$output_base, "05-tcga-filter"), recursive = TRUE, showWarnings = FALSE)

for (group_name in names(meta_results)) {
  tcga_results[[group_name]] <- runTcgaTranslationFilter(
    meta_candidates = meta_results[[group_name]],
    group_set = group_sets[[group_name]],
    tcga_dir = workflow_config$tcga_dir,
    feature_type = workflow_config$feature_type,
    fdr_t = workflow_config$fdr_t,
    es_t = workflow_config$es_t
  )

  fwrite(
    tcga_results[[group_name]],
    file.path(workflow_config$output_base, "05-tcga-filter", paste0(group_name, "_tcga.csv"))
  )
}

save_stage(tcga_results, "05-tcga-filter", "tcga_results.rds")
