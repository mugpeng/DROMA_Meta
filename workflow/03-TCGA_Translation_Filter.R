source(file.path(".", "workflow", "00-Setup.R"), local = FALSE)

cat("\n=== 03: TCGA translation filter ===\n")
dir.create(file.path(workflow_config$output_base, "04-tcga-filter"), recursive = TRUE, showWarnings = FALSE)

group_sets <- read_stage("01-projects", "group_sets.rds")
preclinical_candidates <- read_stage("03-meta", "preclinical_candidates.rds")
preclinical_candidates2 <- preclinical_candidates[preclinical_candidates$direction_concordant %in% TRUE,]

tcga_results <- runTcgaTranslationFilter(
  preclinical_candidates = preclinical_candidates2,
  cellline_set = group_sets$cellline,
  pdcpdx_set = group_sets$pdcpdx,
  tcga_dir = workflow_config$tcga_dir,
  feature_type = workflow_config$feature_type,
  fdr_t = workflow_config$tcga_fdr_t,
  cores = .get_safe_cores(workflow_config$requested_cores),
  gene_probe_map = workflow_config$gene_probe_map
)

save_stage(tcga_results, "04-tcga-filter", "tcga_results.rds")
fwrite(tcga_results, file.path(workflow_config$output_base, "04-tcga-filter", "tcga_results.csv"))
