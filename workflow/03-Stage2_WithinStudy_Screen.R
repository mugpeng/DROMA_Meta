# ============================================================================
# 03-Stage2_WithinStudy_Screen.R
# Grouped biomarker screening using DROMA.R batch meta-analysis
# Reuses: DROMA.R::batchFindSignificantFeatures()
# Logic: drug-response projects are limited by Stage 1 eligibility, while
#        expression can come from any matched sample/model within the same group
#        (same behavior as batchFindSignificantFeatures on MultiDromaSet)
# Statistic: Spearman correlation
# ============================================================================

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

cat("\n=== Step 3: Stage2 Within-Study Screen ===\n")

# Load grouped DromaSet objects and Stage 1 eligible drug/tumor combinations.
group_sets <- read_stage("01-projects", "group_sets.rds")
coverage_results <- read_stage("02-coverage", "coverage_results.rds")
within_results <- list()
dir.create(file.path(workflow_config$output_base, "03-within-study"), recursive = TRUE, showWarnings = FALSE)

for (group_name in names(group_sets)) {
  group_dt <- coverage_results[[group_name]]$candidates_runtime
  within_rows <- list()

  # Skip empty groups early to keep stage output shape stable.
  if (nrow(group_dt) == 0) {
    within_results[[group_name]] <- .empty_within_study()
    next
  }

  # Each row is one drug/tumor candidate that passed Stage 1 coverage rules.
  for (i in seq_len(nrow(group_dt))) {
    row <- group_dt[i]
    dt <- runWithinStudyScreen(
      group_set = group_sets[[group_name]],
      drug_name = row$drug[[1]],
      tumor_type = row$tumor_type[[1]],
      feature_type = workflow_config$feature_type,
      projects = row$eligible_projects[[1]],
      cores = .get_safe_cores(workflow_config$requested_cores),
      preloaded = TRUE,
      verbose = FALSE
    )
    if (nrow(dt) > 0) {
      within_rows[[i]] <- dt
    }
  }

  within_results[[group_name]] <- if (length(within_rows) > 0) {
    rbindlist(within_rows, fill = TRUE)
  } else {
    .empty_within_study()
  }

  # Export one CSV per workflow group for quick inspection.
  fwrite(
    within_results[[group_name]],
    file.path(workflow_config$output_base, "03-within-study", paste0(group_name, "_within_study.csv"))
  )
}

# Save the full object for the stage-3 formatter.
save_stage(within_results, "03-within-study", "within_results.rds")
