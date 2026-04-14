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

cat("\n=== Step 1: Setup And Project Grouping ===\n")
con <- connectDROMADatabase(workflow_config$db_path)
on.exit(closeDROMADatabase(con), add = TRUE)

project_anno <- as.data.table(DROMA.Set::listDROMAProjects(connection = con))
project_groups <- createWorkflowProjectGroups(project_anno)

group_sets <- lapply(project_groups, function(project_names) {
  createWorkflowGroupSet(workflow_config$db_path, con, project_names)
})

save_stage(project_anno, "01-projects", "project_anno.rds")
save_stage(project_groups, "01-projects", "project_groups.rds")
save_stage(group_sets, "01-projects", "group_sets.rds")
dir.create(file.path(workflow_config$output_base, "01-projects"), recursive = TRUE, showWarnings = FALSE)

fwrite(project_anno, file.path(workflow_config$output_base, "01-projects", "project_anno.csv"))
for (group_name in names(project_groups)) {
  fwrite(
    data.table(group = group_name, project_name = project_groups[[group_name]]),
    file.path(workflow_config$output_base, "01-projects", paste0(group_name, "_projects.csv"))
  )
}
