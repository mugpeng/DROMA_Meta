script_dir <- if (basename(getwd()) == "workflow") "." else "workflow"

steps <- c(
  "01-Setup_And_Project_Grouping.R",
  "02-Stage1_Coverage_Filter.R",
  "03-Stage2_Group_Meta.R",
  "04-Stage3_TCGA_Translation_Filter.R",
  "05-Stage4_Preclinical_Merge.R",
  "06-Stage5_CDRTB_Clinical_Validation.R"
)

for (step in steps) {
  source(file.path(script_dir, step), local = FALSE)
}
