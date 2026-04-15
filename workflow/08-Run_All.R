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

steps <- c(
  "01-Setup_And_Project_Grouping.R",
  "02-Stage1_Coverage_Filter.R",
  "04-Stage3_Group_Meta.R",
  "05-Stage4_TCGA_Translation_Filter.R",
  "06-Stage5_Preclinical_Merge.R",
  "07-Stage6_CDRTB_Clinical_Validation.R"
)

for (step in steps) {
  source(file.path(script_dir, step), local = FALSE)
}
