# Back-compat entry: loads libraries, workflow_config, save_stage/read_stage, read_project_grouping.
# Run workflow/00b-Project_Grouping.R before 01+ (same R session); it caches grouping in .GlobalEnv.
source(file.path(".", "workflow", "00a-Setup.R"), local = FALSE)
