script_dir <- file.path(".", "workflow")

steps <- c(
  "01-Coverage_Filter.R",
  "02-Group_Meta.R",
  "03-TCGA_Translation_Filter.R",
  "04-Preclinical_Merge.R",
  "05-CDRTB_Clinical_Validation.R"
)

for (step in steps) {
  source(file.path(script_dir, step), local = FALSE)
}
