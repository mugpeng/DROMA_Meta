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
group_sets <- read_stage("01-projects", "group_sets.rds")
coverage_results <- read_stage("02-coverage", "coverage_results.rds")
within_results <- list()
dir.create(file.path(workflow_config$output_base, "03-within-study"), recursive = TRUE, showWarnings = FALSE)

for (group_name in names(group_sets)) {
  group_dt <- coverage_results[[group_name]]$candidates_runtime
  within_rows <- list()
  if (nrow(group_dt) == 0) {
    within_results[[group_name]] <- .empty_within_study()
    next
  }

  for (i in seq_len(nrow(group_dt))) {
    row <- group_dt[i]
    dt <- runWithinStudyScreen(
      group_set = group_sets[[group_name]],
      drug_name = row$drug[[1]],
      tumor_type = row$tumor_type[[1]],
      feature_type = workflow_config$feature_type,
      projects = row$eligible_projects[[1]]
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

  fwrite(
    within_results[[group_name]],
    file.path(workflow_config$output_base, "03-within-study", paste0(group_name, "_within_study.csv"))
  )
}

save_stage(within_results, "03-within-study", "within_results.rds")
