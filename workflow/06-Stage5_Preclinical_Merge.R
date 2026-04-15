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

cat("\n=== Step 6: Stage5 Preclinical Merge ===\n")
dir.create(file.path(workflow_config$output_base, "06-preclinical-merge"), recursive = TRUE, showWarnings = FALSE)

meta_results <- read_stage("04-meta", "meta_results.rds")
tcga_results <- read_stage("05-tcga-filter", "tcga_results.rds")

merged_candidates <- mergePreclinicalCandidates(
  cellline_meta = meta_results$cellline,
  pdcpdx_meta = meta_results$pdcpdx,
  tcga_results = tcga_results,
  fdr_t = workflow_config$fdr_t,
  es_t = workflow_config$es_t
)

merged_candidates <- merged_candidates[direction_concordant == TRUE]

save_stage(merged_candidates, "06-preclinical-merge", "merged_candidates.rds")
fwrite(merged_candidates, file.path(workflow_config$output_base, "06-preclinical-merge", "merged_candidates.csv"))
