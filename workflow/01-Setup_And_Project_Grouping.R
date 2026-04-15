source(file.path(if (basename(getwd()) == "workflow") "." else "workflow", "00-Workflow_Common.R"), local = FALSE)

cat("\n=== Step 1: Setup And Project Grouping ===\n")
con <- DROMA.Set::connectDROMADatabase(workflow_config$db_path)
on.exit(DROMA.Set::closeDROMADatabase(con), add = TRUE)
dir.create(file.path(workflow_config$output_base, "01-projects"), recursive = TRUE, showWarnings = FALSE)

project_anno <- data.table::as.data.table(DROMA.Set::listDROMAProjects(connection = con))
project_groups <- createWorkflowProjectGroups(project_anno)
group_sets <- lapply(project_groups, function(project_names) {
  createWorkflowGroupSet(workflow_config$db_path, con, project_names)
})
shared_features <- getSharedGroupFeatures(group_sets, workflow_config$feature_type)

save_stage(project_anno, "01-projects", "project_anno.rds")
save_stage(project_groups, "01-projects", "project_groups.rds")
save_stage(group_sets, "01-projects", "group_sets.rds")
save_stage(shared_features, "01-projects", "shared_features.rds")

fwrite(project_anno, file.path(workflow_config$output_base, "01-projects", "project_anno.csv"))
for (group_name in names(project_groups)) {
  fwrite(
    data.table(group = group_name, project_name = project_groups[[group_name]]),
    file.path(workflow_config$output_base, "01-projects", paste0(group_name, "_projects.csv"))
  )
}
fwrite(
  data.table(feature = shared_features),
  file.path(workflow_config$output_base, "01-projects", "shared_features.csv")
)
