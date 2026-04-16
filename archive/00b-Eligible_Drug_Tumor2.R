# ============================================================================
# 00b-Eligible_Drug_Tumor2.R
# Filter drug-tumor pairs to those supported by CTRDB annotations
# ============================================================================

library(data.table)
suppressPackageStartupMessages(library(DROMA.Meta))
library(DROMA.Set)
library(DROMA.R)

db_path <- "/Users/peng/Desktop/Project/DROMA/Data/droma.sqlite"
project_root_candidates <- c(
  file.path(getwd(), "Meta_Example"),
  file.path(getwd(), "..", "Meta_Example"),
  file.path(getwd(), "..", "..", "Meta_Example")
)
project_root <- project_root_candidates[file.exists(project_root_candidates)][1]
if (is.na(project_root) || !nzchar(project_root)) {
  stop("Could not locate Meta_Example from the current working directory", call. = FALSE)
}
project_root <- normalizePath(project_root, mustWork = TRUE)
defaults <- getMetaWorkflowDefaults(project_root = project_root)
output_base <- defaults$output_base
dir.create(output_base, showWarnings = FALSE, recursive = TRUE)

input_pairs_path <- file.path(output_base, "eligible_drug_pairs_filter.csv")
output_pairs_path <- file.path(output_base, "eligible_drug_tumor_pairs_ctrdb.csv")

cat("\n=== 00b: Eligible Drug And Tumor Type In CTRDB ===\n")

connectDROMADatabase(db_path)

drug_anno <- as.data.table(getDROMAAnnotation("drug"))
sample_anno <- as.data.table(getDROMAAnnotation("sample"))
input_pairs <- fread(input_pairs_path)

required_cols <- c("drug", "tumor_type")
missing_cols <- setdiff(required_cols, names(input_pairs))
if (length(missing_cols) > 0) {
  stop(
    sprintf(
      "Missing required column(s) in %s: %s",
      input_pairs_path,
      paste(missing_cols, collapse = ", ")
    ),
    call. = FALSE
  )
}

input_pairs <- unique(input_pairs[, ..required_cols])

ctrdb_drugs <- unique(
  drug_anno[
    ProjectID == "CTRDB" & !is.na(DrugName) & nzchar(DrugName),
    .(drug = as.character(DrugName))
  ]
)

ctrdb_tumors <- unique(
  sample_anno[
    ProjectID == "CTRDB" & !is.na(TumorType) & nzchar(TumorType),
    .(tumor_type = as.character(TumorType))
  ]
)

eligible_pairs_ctrdb <- input_pairs[
  ctrdb_drugs,
  on = "drug",
  nomatch = 0L
][
  ctrdb_tumors,
  on = "tumor_type",
  nomatch = 0L
]

eligible_pairs_ctrdb <- unique(eligible_pairs_ctrdb[order(drug, tumor_type)])
if (nrow(eligible_pairs_ctrdb) == 0L) {
  stop("Refusing to write empty eligible_drug_tumor_pairs_ctrdb table", call. = FALSE)
}
fwrite(eligible_pairs_ctrdb, output_pairs_path)

cat(sprintf("  Input pairs: %d\n", nrow(input_pairs)))
cat(sprintf("  CTRDB drugs: %d\n", nrow(ctrdb_drugs)))
cat(sprintf("  CTRDB tumor types: %d\n", nrow(ctrdb_tumors)))
cat(sprintf("  Eligible CTRDB pairs: %d\n", nrow(eligible_pairs_ctrdb)))
cat(sprintf("  Wrote: %s\n", output_pairs_path))
