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

cat("\n=== Step 4: Stage3 Group Meta ===\n")
dir.create(file.path(workflow_config$output_base, "04-meta"), recursive = TRUE, showWarnings = FALSE)

group_sets <- read_stage("01-projects", "group_sets.rds")
shared_features <- read_stage("01-projects", "shared_features.rds")
coverage_results <- read_stage("02-coverage", "coverage_results.rds")
meta_results <- list()
sig_results <- list()

for (group_name in names(group_sets)) {
  group_meta_rows <- list()
  runtime_dt <- coverage_results[[group_name]]$candidates_runtime

  if (nrow(runtime_dt) > 0) {
    for (i in seq_len(nrow(runtime_dt))) {
      drug_name <- runtime_dt$drug[[i]]
      tumor_type <- runtime_dt$tumor_type[[i]]
      meta_dt <- runGroupedMetaAnalysis(
        group_set = group_sets[[group_name]],
        group_name = group_name,
        drug_name = drug_name,
        tumor_type = tumor_type,
        feature_names = shared_features,
        feature_type = workflow_config$feature_type,
        cores = .get_safe_cores(workflow_config$requested_cores),
        preloaded = TRUE
      )
      if (nrow(meta_dt) > 0) {
        group_meta_rows[[length(group_meta_rows) + 1L]] <- meta_dt
      }
    }
  }

  meta_results[[group_name]] <- if (length(group_meta_rows) > 0) {
    rbindlist(group_meta_rows, fill = TRUE)
  } else {
    .empty_meta()
  }
  sig_results[[group_name]] <- meta_results[[group_name]][
    q_value < workflow_config$fdr_t & abs(effect_size) >= workflow_config$es_t
  ]

  fwrite(meta_results[[group_name]], file.path(workflow_config$output_base, "04-meta", paste0(group_name, "_meta.csv")))
  fwrite(sig_results[[group_name]], file.path(workflow_config$output_base, "04-meta", paste0(group_name, "_meta_sig.csv")))
}

preclinical_candidates <- mergePreclinicalCandidates(
  cellline_meta = meta_results$cellline,
  pdcpdx_meta = meta_results$pdcpdx,
  fdr_t = workflow_config$fdr_t,
  es_t = workflow_config$es_t
)

save_stage(meta_results, "04-meta", "meta_results.rds")
save_stage(sig_results, "04-meta", "sig_results.rds")
save_stage(preclinical_candidates, "04-meta", "preclinical_candidates.rds")
fwrite(preclinical_candidates, file.path(workflow_config$output_base, "04-meta", "preclinical_candidates.csv"))
