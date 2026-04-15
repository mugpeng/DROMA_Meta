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

cat("\n=== Step 2: Stage1 Coverage Filter ===\n")
con <- DROMA.Set::connectDROMADatabase(workflow_config$db_path)
on.exit(DROMA.Set::closeDROMADatabase(con), add = TRUE)
dir.create(file.path(workflow_config$output_base, "02-coverage"), recursive = TRUE, showWarnings = FALSE)

project_groups <- read_stage("01-projects", "project_groups.rds")
coverage_results <- list()

for (group_name in names(project_groups)) {
  group_cfg <- workflow_config$groups[[group_name]]
  coverage_results[[group_name]] <- computeCoverageForGroup(
    project_names = project_groups[[group_name]],
    drug_names = workflow_config$drug_names,
    tumor_types = workflow_config$tumor_types,
    db_path = workflow_config$db_path,
    con = con,
    feature_type = workflow_config$feature_type,
    dataset_label = group_name,
    min_overlap_per_pair = group_cfg$min_overlap_per_pair,
    min_pairs = group_cfg$min_pairs
  )
  fwrite(
    coverage_results[[group_name]]$coverage,
    file.path(workflow_config$output_base, "02-coverage", paste0(group_name, "_coverage.csv"))
  )
  fwrite(
    coverage_results[[group_name]]$candidates_save,
    file.path(workflow_config$output_base, "02-coverage", paste0(group_name, "_candidates.csv"))
  )
}

save_stage(coverage_results, "02-coverage", "coverage_results.rds")
