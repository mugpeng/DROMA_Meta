source(file.path(".", "workflow", "00-Setup.R"), local = FALSE)

cat("\n=== 02: Group meta ===\n")
dir.create(file.path(workflow_config$output_base, "03-meta"), recursive = TRUE, showWarnings = FALSE)

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
      filtered_projects <- runtime_dt$filtered_projects[[i]]
      keep_projects <- intersect(unique(as.character(filtered_projects)), names(group_sets[[group_name]]@DromaSets))
      filtered_group_set <- DROMA.Set::createMultiDromaSetFromObjects(
        group_sets[[group_name]]@DromaSets[keep_projects],
        project_names = keep_projects
      )
      meta_dt <- batchFindSignificantFeatures(
        dromaset_object = filtered_group_set,
        feature1_type = "drug",
        feature1_name = drug_name,
        feature2_type = workflow_config$feature_type,
        feature2_name = shared_features,
        data_type = "all",
        tumor_type = tumor_type,
        overlap_only = FALSE,
        cores = .get_safe_cores(workflow_config$requested_cores),
        show_progress = TRUE,
        preloaded = TRUE,
        verbose = TRUE
      )
      if (nrow(meta_dt) > 0) {
        meta_dt <- data.table::as.data.table(meta_dt)
        if (!"q_value" %in% colnames(meta_dt)) {
          meta_dt[, q_value := stats::p.adjust(p_value, method = "BH")]
        }
        meta_dt[, `:=`(
          drug = drug_name,
          tumor_type = tumor_type,
          model_group = group_name,
          heterogeneity_p = if ("pval.Q" %in% colnames(meta_dt)) pval.Q else NA_real_,
          i2 = if ("I2" %in% colnames(meta_dt)) I2 else NA_real_,
          direction = ifelse(effect_size >= 0, "Up", "Down")
        )]
        keep_cols <- c(
          "drug", "tumor_type", "name", "p_value", "q_value", "effect_size",
          "n_datasets", "model_group", "heterogeneity_p", "i2", "direction"
        )
        for (col in keep_cols) {
          if (!col %in% colnames(meta_dt)) {
            meta_dt[, (col) := NA]
          }
        }
        meta_dt <- unique(meta_dt[, ..keep_cols])
        group_meta_rows[[length(group_meta_rows) + 1L]] <- meta_dt
      }
    }
  }

  meta_results[[group_name]] <- if (length(group_meta_rows) > 0) {
    rbindlist(group_meta_rows, fill = TRUE)
  } else {
    .empty_meta()
  }
  sig_results[[group_name]] <- if (isTRUE(workflow_config$meta_use_p_value)) {
    meta_results[[group_name]][
      !is.na(p_value) & p_value < workflow_config$meta_p_t & abs(effect_size) >= workflow_config$es_t
    ]
  } else {
    meta_results[[group_name]][
      !is.na(q_value) & q_value < workflow_config$fdr_t & abs(effect_size) >= workflow_config$es_t
    ]
  }

  fwrite(meta_results[[group_name]], file.path(workflow_config$output_base, "03-meta", paste0(group_name, "_meta.csv")))
  fwrite(sig_results[[group_name]], file.path(workflow_config$output_base, "03-meta", paste0(group_name, "_meta_sig.csv")))
}

preclinical_candidates <- mergePreclinicalCandidates(
  cellline_meta = meta_results$cellline,
  pdcpdx_meta = meta_results$pdcpdx,
  fdr_t = workflow_config$fdr_t,
  es_t = workflow_config$es_t,
  use_p_value = workflow_config$meta_use_p_value,
  p_t = workflow_config$meta_p_t
)

save_stage(meta_results, "03-meta", "meta_results.rds")
save_stage(sig_results, "03-meta", "sig_results.rds")
save_stage(preclinical_candidates, "03-meta", "preclinical_candidates.rds")
fwrite(preclinical_candidates, file.path(workflow_config$output_base, "03-meta", "preclinical_candidates.csv"))
