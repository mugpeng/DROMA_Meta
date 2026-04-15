# Batch driver for DROMA.Meta ----

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DROMA.Set))
suppressPackageStartupMessages(library(DROMA.R))

source("R/FuncHelper.R", local = FALSE)
source("R/FuncValidCheck.R", local = FALSE)
source("R/FuncTcgaAD.R", local = FALSE)
source("R/FuncMetaWorkflow.R", local = FALSE)

project_root <- normalizePath(getwd(), mustWork = TRUE)
defaults <- getMetaWorkflowDefaults(project_root = project_root)

# Driver-level input and output locations. Keep these runtime choices here so
# the files under R/ remain reusable pure function definitions.
valid_drugs_csv <- file.path(project_root, "workflow", "Output", "valid_drugs.csv")
valid_tumor_types_csv <- file.path(project_root, "workflow", "Output", "valid_tumor_types.csv")
summary_csv <- file.path(project_root, "workflow", "Output", "meta_workflow_batch_summary.csv")

summary_dt <- runMetaWorkflowBatch(
  valid_drugs_csv = valid_drugs_csv,
  valid_tumor_types_csv = valid_tumor_types_csv,
  output_base = defaults$output_base,
  summary_csv = summary_csv,
  verbose = TRUE
)

print(summary_dt)
