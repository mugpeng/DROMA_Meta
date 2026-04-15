source(file.path(if (basename(getwd()) == "workflow") "." else "workflow", "00-Workflow_Common.R"), local = FALSE)

cat("\n=== Step 5: Stage4 Preclinical Merge ===\n")
dir.create(file.path(workflow_config$output_base, "05-preclinical-merge"), recursive = TRUE, showWarnings = FALSE)

meta_results <- read_stage("03-meta", "meta_results.rds")
tcga_results <- read_stage("04-tcga-filter", "tcga_results.rds")

merged_candidates <- mergePreclinicalCandidates(
  cellline_meta = meta_results$cellline,
  pdcpdx_meta = meta_results$pdcpdx,
  tcga_results = tcga_results,
  fdr_t = workflow_config$fdr_t,
  es_t = workflow_config$es_t
)

merged_candidates <- merged_candidates[direction_concordant == TRUE]

save_stage(merged_candidates, "05-preclinical-merge", "merged_candidates.rds")
fwrite(merged_candidates, file.path(workflow_config$output_base, "05-preclinical-merge", "merged_candidates.csv"))
