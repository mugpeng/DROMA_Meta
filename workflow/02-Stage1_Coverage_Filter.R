# ============================================================================
# 02-Stage1_Coverage_Filter.R
# Coverage prefilter for grouped biomarker screening
# Drug response is anchored to each focal project, while expression/model
# availability is traversed across all projects in the same workflow group.
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

cat("\n=== Step 2: Stage1 Coverage Filter ===\n")
con <- connectDROMADatabase(workflow_config$db_path)
on.exit(closeDROMADatabase(con), add = TRUE)

# Load project grouping prepared in Step 1.
project_groups <- read_stage("01-projects", "project_groups.rds")
coverage_results <- list()
dir.create(file.path(workflow_config$output_base, "02-coverage"), recursive = TRUE, showWarnings = FALSE)

for (group_name in names(project_groups)) {
  group_cfg <- workflow_config$groups[[group_name]]

  # For each workflow group, require enough focal drug-response projects whose
  # samples can be matched to expression data found anywhere within the group.
  coverage_results[[group_name]] <- computeCoverageForGroup(
    project_names = project_groups[[group_name]],
    drug_names = workflow_config$drug_names,
    tumor_types = workflow_config$tumor_types,
    db_path = workflow_config$db_path,
    con = con,
    feature_type = workflow_config$feature_type,
    dataset_label = group_name,
    min_overlap_per_study = group_cfg$min_overlap_per_study,
    min_studies = group_cfg$min_studies
  )

  # Export detailed coverage and candidate tables for manual inspection.
  fwrite(
    coverage_results[[group_name]]$coverage,
    file.path(workflow_config$output_base, "02-coverage", paste0(group_name, "_coverage.csv"))
  )
  fwrite(
    coverage_results[[group_name]]$candidates_save,
    file.path(workflow_config$output_base, "02-coverage", paste0(group_name, "_candidates.csv"))
  )
}

# Save stage output for downstream screening.
save_stage(coverage_results, "02-coverage", "coverage_results.rds")
