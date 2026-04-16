fallback_if_missing <- function(x, y) {
  if (is.null(x) || length(x) == 0 || !nzchar(x)) y else x
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- sub("^--file=", "", fallback_if_missing(file_arg[1], ""))
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path, mustWork = TRUE)) else getwd()
root_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

source(file.path(root_dir, "visualization", "R", "FuncVisHelper.R"))
source(file.path(root_dir, "visualization", "R", "FuncVisPipeline.R"))

summary_dt <- data.table(
  drug = c("DrugA", "DrugB"),
  tumor_type = c("TumorA", "TumorB"),
  ctrdb_status = c("tumor_type_matched", "no_ctrdb_data"),
  ctrdb_fallback = c(FALSE, FALSE),
  clinical_query_tumor_type = c("TumorA", "all"),
  n_clinical_sig = c(3L, 0L),
  n_final_biomarkers = c(2L, 0L),
  output_dir = c("a", "b")
)

p <- plotPipelineFunnel(summary_dt)
stopifnot(inherits(p, "ggplot"))

cat("visualization pipeline tests passed\n")
