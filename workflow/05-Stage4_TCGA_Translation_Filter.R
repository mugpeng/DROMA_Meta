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
  stop("Could not locate workflow/00-resolve_workflow_dir.R", call. = FALSE)
}
script_dir <- resolve_workflow_script_dir()
source(file.path(script_dir, "00-Workflow_Common.R"), local = FALSE)

cat("\n=== Step 5: Stage4 TCGA Translation Filter ===\n")
dir.create(file.path(workflow_config$output_base, "05-tcga-filter"), recursive = TRUE, showWarnings = FALSE)

group_sets <- read_stage("01-projects", "group_sets.rds")
preclinical_candidates <- read_stage("04-meta", "preclinical_candidates.rds")

tcga_results <- runTcgaTranslationFilter(
  preclinical_candidates = preclinical_candidates,
  cellline_set = group_sets$cellline,
  pdcpdx_set = group_sets$pdcpdx,
  tcga_dir = workflow_config$tcga_dir,
  feature_type = workflow_config$feature_type,
  fdr_t = workflow_config$tcga_fdr_t,
  cores = .get_safe_cores(workflow_config$requested_cores)
)

save_stage(tcga_results, "05-tcga-filter", "tcga_results.rds")
fwrite(tcga_results, file.path(workflow_config$output_base, "05-tcga-filter", "tcga_results.csv"))
