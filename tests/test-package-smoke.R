fallback_if_missing <- function(x, y) {
  if (is.null(x) || length(x) == 0 || !nzchar(x)) y else x
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- sub("^--file=", "", fallback_if_missing(file_arg[1], ""))
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path, mustWork = TRUE)) else getwd()
root_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

required_files <- c(
  file.path(root_dir, "DESCRIPTION"),
  file.path(root_dir, "NAMESPACE"),
  file.path(root_dir, "one.R"),
  file.path(root_dir, "two.R"),
  file.path(root_dir, "two_bk.R"),
  file.path(root_dir, "R", "FuncHelper.R"),
  file.path(root_dir, "R", "FuncValidCheck.R"),
  file.path(root_dir, "R", "FuncTcgaAD.R"),
  file.path(root_dir, "R", "FuncMetaWorkflow.R")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(sprintf("Missing required files:\n%s", paste(missing_files, collapse = "\n")))
}

desc <- read.dcf(file.path(root_dir, "DESCRIPTION"))
stopifnot(desc[1, "Package"] == "DROMA.Meta")

ns_lines <- readLines(file.path(root_dir, "NAMESPACE"), warn = FALSE)
stopifnot(any(grepl("^export\\(runMetaWorkflow\\)$", ns_lines)))
stopifnot(!any(grepl("^export\\(createEmptyMetaDf\\)$", ns_lines)))
stopifnot(!any(grepl("^export\\(createEmptyFinalBiomarkersDf\\)$", ns_lines)))
stopifnot(!any(grepl("^export\\(fwriteNonEmptyOrStop\\)$", ns_lines)))

one_lines <- readLines(file.path(root_dir, "one.R"), warn = FALSE)
stopifnot(any(grepl("library(DROMA.Meta)", one_lines, fixed = TRUE)))
stopifnot(!any(grepl("source(\"R/", one_lines, fixed = TRUE)))
stopifnot(any(grepl("eligible_pairs_csv <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("eligible_pairs <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("eligible_pair_flags <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("eligible_pairs_run <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("pdcpdx_ge_2_pair", one_lines, fixed = TRUE)))
stopifnot(any(grepl("ctrdb_status", one_lines, fixed = TRUE)))
stopifnot(any(grepl("output_base_batch <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("summary_csv <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("for (i in seq_len(nrow(eligible_pairs_run)))", one_lines, fixed = TRUE)))
stopifnot(any(grepl("runMetaWorkflow(", one_lines, fixed = TRUE)))
stopifnot(any(grepl("summary_dt <- merge(", one_lines, fixed = TRUE)))
stopifnot(any(grepl("db_path <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("ctrdb_path <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("tcga_rna_counts_dir <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("gene_probe_map_path <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("data_type <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("cores <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("cell_min_intersected_cells <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("pdcpdx_min_intersected_cells <-", one_lines, fixed = TRUE)))
stopifnot(any(grepl("pdcpdx_n_datasets_t <- 1", one_lines, fixed = TRUE)))

two_lines <- readLines(file.path(root_dir, "two.R"), warn = FALSE)
stopifnot(any(grepl("library(DROMA.Meta)", two_lines, fixed = TRUE)))
stopifnot(any(grepl("summary_csv <-", two_lines, fixed = TRUE)))
stopifnot(any(grepl("annotation_summary_csv <-", two_lines, fixed = TRUE)))
stopifnot(any(grepl("summary_dt <- data.table::fread", two_lines, fixed = TRUE)))
stopifnot(any(grepl("for \\(i in seq_len\\(nrow\\(summary_dt\\)\\)\\)", two_lines)))
stopifnot(any(grepl("final_biomarkers_annotated.csv", two_lines, fixed = TRUE)))
stopifnot(any(grepl("pair_biomarker_annotation_summary.csv", two_lines, fixed = TRUE)))
stopifnot(!any(grepl("annotatePair <- function", two_lines, fixed = TRUE)))

two_bk_lines <- readLines(file.path(root_dir, "two_bk.R"), warn = FALSE)
stopifnot(any(grepl("eligible_pairs_csv <-", two_bk_lines, fixed = TRUE)))
stopifnot(any(grepl("summary_csv <-", two_bk_lines, fixed = TRUE)))
stopifnot(any(grepl("# 01: Load driver pair table", two_bk_lines, fixed = TRUE)))
stopifnot(any(grepl("# 02: Annotate pair outputs", two_bk_lines, fixed = TRUE)))
stopifnot(any(grepl("# 03: Write backup annotation summary", two_bk_lines, fixed = TRUE)))
stopifnot(any(grepl("file.exists(summary_csv)", two_bk_lines, fixed = TRUE)))
stopifnot(any(grepl("data.table::fread(eligible_pairs_csv", two_bk_lines, fixed = TRUE)))
stopifnot(any(grepl("batch_driver_source", two_bk_lines, fixed = TRUE)))
stopifnot(any(grepl("missing_output_dir", two_bk_lines, fixed = TRUE)))
stopifnot(any(grepl("missing_final_biomarkers", two_bk_lines, fixed = TRUE)))
stopifnot(any(grepl("empty_final_biomarkers", two_bk_lines, fixed = TRUE)))
stopifnot(any(grepl("unreadable_final_biomarkers", two_bk_lines, fixed = TRUE)))

suppressPackageStartupMessages(library(DROMA.Meta))

stopifnot(exists("runMetaWorkflow", mode = "function"))
stopifnot(exists("buildDrugTumorGrid", mode = "function"))

defaults <- getMetaWorkflowDefaults(root_dir)
stopifnot(is.null(defaults$valid_drugs_csv))
stopifnot(is.null(defaults$valid_tumor_types_csv))
stopifnot(is.null(defaults$batch_summary_csv))
stopifnot(identical(defaults$droma_db_path, "/Users/peng/Desktop/Project/DROMA/Data/droma.sqlite"))

tmp_dir <- tempfile("package-smoke-")
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
valid_drugs_csv <- file.path(tmp_dir, "valid_drugs.csv")
valid_tumor_types_csv <- file.path(tmp_dir, "valid_tumor_types.csv")

writeLines(
  c(
    "drug",
    "DrugA",
    "DrugB"
  ),
  valid_drugs_csv
)
writeLines(
  c(
    "tumor_type",
    "TumorA",
    "TumorB"
  ),
  valid_tumor_types_csv
)

grid <- buildDrugTumorGrid(
  valid_drugs_csv,
  valid_tumor_types_csv
)
stopifnot(nrow(grid) > 0)
stopifnot(all(c("drug", "tumor_type") %in% colnames(grid)))

cat("package smoke test passed\n")
