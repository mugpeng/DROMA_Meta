source(file.path(".", "workflow", "00a-Setup.R"), local = FALSE)

cat("\n=== 05: CDRTB clinical validation ===\n")
dir.create(file.path(workflow_config$output_base, "05-clinical-validation"), recursive = TRUE, showWarnings = FALSE)

merged_candidates <- read_stage("04-preclinical-merge", "merged_candidates.rds")
merged_candidates <- merged_candidates[merged_candidates$tcga_supported %in% TRUE, ]

clinical_results <- data.table::copy(merged_candidates)

if (nrow(merged_candidates) > 0) {
  DROMA.Set::connectCTRDatabase(workflow_config$ctrdb_path)
  on.exit({
    if (exists("ctrdb_connection", envir = .GlobalEnv)) {
      rm("ctrdb_connection", envir = .GlobalEnv)
    }
  }, add = TRUE)

  clinical_rows <- list()
  split_candidates <- split(merged_candidates, interaction(merged_candidates$drug, merged_candidates$tumor_type, drop = TRUE))

  for (candidate_df in split_candidates) {
    run_one_scope <- function(scope) {
      tryCatch(
        DROMA.R::batchFindClinicalSigResponse(
          select_omics = candidate_df$name,
          select_drugs = candidate_df$drug[[1]],
          data_type = "all",
          tumor_type = scope,
          cores = .get_safe_cores(workflow_config$requested_cores)
        ),
        error = function(e) NULL
      )
    }

    clinical_scope <- candidate_df$tumor_type[[1]]
    fallback_to_all <- FALSE
    cli_df <- run_one_scope(clinical_scope)

    if (is.null(cli_df) || !is.data.frame(cli_df) || nrow(cli_df) == 0) {
      cli_df <- run_one_scope("all")
      clinical_scope <- "all"
      fallback_to_all <- TRUE
    }

    if (is.null(cli_df) || !is.data.frame(cli_df) || nrow(cli_df) == 0) {
      out <- data.table::copy(candidate_df)
      out$clinical_scope <- clinical_scope
      out$fallback_to_all <- fallback_to_all
      out$clinical_supported <- FALSE
      out$direction_clinical <- NA_character_
      out$direction_concordant <- FALSE
      out$retained <- FALSE
      clinical_rows[[length(clinical_rows) + 1L]] <- out
      next
    }

    clinical_candidates <- DROMA.R::getSignificantFeatures(
      meta_df = cli_df,
      es_t = workflow_config$es_t,
      P_t = workflow_config$fdr_t,
      use_p_value = TRUE
    )

    if (!"direction" %in% colnames(cli_df) && "effect_size" %in% colnames(cli_df)) {
      cli_df$direction <- ifelse(cli_df$effect_size >= 0, "Up", "Down")
    }

    merged <- merge(
      candidate_df,
      cli_df[, c("name", "p_value", "q_value", "effect_size", "n_datasets", "direction"), with = FALSE],
      by = "name",
      all.x = TRUE,
      suffixes = c("_preclinical", "_clinical")
    )

    merged$clinical_scope <- clinical_scope
    merged$fallback_to_all <- fallback_to_all
    merged$clinical_supported <- merged$name %in% clinical_candidates$name
    merged$direction_clinical <- merged$direction
    merged$direction_concordant <- merged$clinical_supported &
      !is.na(merged$direction_clinical) &
      merged$preclinical_direction == merged$direction_clinical
    merged$retained <- merged$clinical_supported & merged$direction_concordant

    clinical_rows[[length(clinical_rows) + 1L]] <- merged
  }

  clinical_results <- data.table::rbindlist(clinical_rows, fill = TRUE)
}

save_stage(clinical_results, "05-clinical-validation", "clinical_results.rds")
fwrite(clinical_results, file.path(workflow_config$output_base, "05-clinical-validation", "clinical_results.csv"))
