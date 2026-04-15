source(file.path(".", "workflow", "00a-Setup.R"), local = FALSE)

# One `{ ... }` block so `on.exit` is tied to a single eval frame under `source()`:
# otherwise each top-level line is a separate `eval()`, and `on.exit` runs as soon
# as the `on.exit(...)` line finishes, closing `con` before `listDROMAProjects`.
{
  cat("\n=== 00b: Project grouping ===\n")
  con <- DROMA.Set::connectDROMADatabase(workflow_config$db_path)
  on.exit(DROMA.Set::closeDROMADatabase(con), add = TRUE)

  project_anno <- listDROMAProjects(connection = con)
  project_groups <- createWorkflowProjectGroups(project_anno)
  group_sets <- lapply(project_groups, function(project_names) {
    DROMA.Set::createMultiDromaSetFromAllProjects(
      db_path = workflow_config$db_path,
      include_projects = sort(unique(project_names)),
      con = con
    )
  })
  shared_features <- getSharedGroupFeatures(group_sets, workflow_config$feature_type)

  assign(
    ".droma_wf_01_projects",
    list(
      project_anno = project_anno,
      project_groups = project_groups,
      group_sets = group_sets,
      shared_features = shared_features
    ),
    envir = .GlobalEnv
  )
}
