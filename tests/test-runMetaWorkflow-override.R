fallback_if_missing <- function(x, y) {
  if (is.null(x) || length(x) == 0 || !nzchar(x)) y else x
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- sub("^--file=", "", fallback_if_missing(file_arg[1], ""))
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path, mustWork = TRUE)) else getwd()
root_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

suppressPackageStartupMessages(library(data.table))

source(file.path(root_dir, "R", "FuncHelper.R"))
source(file.path(root_dir, "R", "FuncMetaWorkflow.R"))

local_tempdir <- function(pattern = "file") {
  path <- tempfile(pattern = pattern)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

test_skip_existing_outputs <- function() {
  output_base <- local_tempdir("meta-output-")
  output_dir <- getMetaWorkflowOutputDir("DrugA", "Tumor A", output_base)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  fwrite(data.table(name = "GENE1"), file.path(output_dir, "batch_cell_sets_mRNA.csv"))
  fwrite(data.table(name = "GENE1"), file.path(output_dir, "batch_pdcpdx_sets_mRNA.csv"))
  fwrite(data.table(name = "GENE1"), file.path(output_dir, "mRNA_cell_sig.csv"))
  fwrite(data.table(name = "GENE1"), file.path(output_dir, "mRNA_pdcpdx_sig.csv"))
  fwrite(data.table(name = "GENE1"), file.path(output_dir, "selected_genes.csv"))
  fwrite(data.table(name = "GENE1", ccle_vs_tcga_concordant = TRUE), file.path(output_dir, "selected_genes_ad_stats.csv"))
  fwrite(data.table(name = "GENE1", ccle_vs_tcga_concordant = TRUE), file.path(output_dir, "selected_genes_ad_filtered.csv"))
  fwrite(data.table(name = "GENE1"), file.path(output_dir, "clinical_sig_mRNA.csv"))
  fwrite(data.table(name = "GENE1"), file.path(output_dir, "final_biomarkers.csv"))

  blockers <- c(
    "connectDROMADatabase",
    "listDROMAProjects",
    "createMultiDromaSetFromAllProjects",
    "batchFindSignificantFeatures",
    "getSignificantFeatures",
    "getIntersectSignificantFeatures",
    "batchFindTcgaADConcordantFeatures",
    "connectCTRDatabase",
    "batchFindClinicalSigResponse"
  )
  for (fn in blockers) {
    assign(fn, function(...) stop(sprintf("%s should not be called", fn), call. = FALSE), envir = .GlobalEnv)
  }

  on.exit(rm(list = blockers, envir = .GlobalEnv), add = TRUE)

  result <- runMetaWorkflow(
    drug = "DrugA",
    tumor_type = "Tumor A",
    output_base = output_base,
    override = FALSE,
    verbose = FALSE
  )

  stopifnot(identical(result$status[[1]], "success"))
  stopifnot(identical(result$n_batch_cell[[1]], 1L))
  stopifnot(identical(result$n_batch_pdcpdx[[1]], 1L))
  stopifnot(identical(result$n_selected_genes[[1]], 1L))
  stopifnot(identical(result$n_ad_filtered[[1]], 1L))
  stopifnot(identical(result$n_clinical_sig[[1]], 1L))
  stopifnot(identical(result$n_final_biomarkers[[1]], 1L))
}

test_force_override_recomputes_outputs <- function() {
  output_base <- local_tempdir("meta-output-")

  calls <- new.env(parent = emptyenv())
  calls$batch_sig <- 0L
  calls$ad <- 0L

  assign("connectDROMADatabase", function(...) textConnection("droma", open = "r"), envir = .GlobalEnv)
  assign(
    "listDROMAProjects",
    function(...) data.frame(
      dataset_type = c("CellLine", "PDC", "PDO", "PDX", "Other"),
      project_name = c("CELL", "PDC1", "PDO1", "PDX1", "CCLE"),
      stringsAsFactors = FALSE
    ),
    envir = .GlobalEnv
  )
  assign(
    "createMultiDromaSetFromAllProjects",
    function(..., include_projects, con = NULL) {
      structure(list(projects = include_projects), class = "MultiDromaSet")
    },
    envir = .GlobalEnv
  )
  assign(
    "batchFindSignificantFeatures",
    function(...) {
      calls$batch_sig <- calls$batch_sig + 1L
      data.table(name = "GENE1", direction = "up")
    },
    envir = .GlobalEnv
  )
  assign(
    "getSignificantFeatures",
    function(x, ...) as.data.table(x),
    envir = .GlobalEnv
  )
  assign(
    "getIntersectSignificantFeatures",
    function(...) {
      args <- list(...)
      if (!is.null(args$cell) && !is.null(args$pdcpdx) && is.null(args$ctrdb)) {
        return(data.table(name = "GENE1", direction_pdcpdx = "up"))
      }
      data.frame(name = "GENE1", direction_pdcpdx = "up", direction = "up", stringsAsFactors = FALSE)
    },
    envir = .GlobalEnv
  )
  assign(
    "batchFindTcgaADConcordantFeatures",
    function(selected_features, ...) {
      calls$ad <- calls$ad + 1L
      out <- as.data.table(selected_features)
      out[, ccle_vs_tcga_concordant := TRUE]
      out
    },
    envir = .GlobalEnv
  )
  assign("connectCTRDatabase", function(...) TRUE, envir = .GlobalEnv)
  assign(
    "batchFindClinicalSigResponse",
    function(...) data.frame(name = "GENE1", direction = "up", stringsAsFactors = FALSE),
    envir = .GlobalEnv
  )

  cleanup <- c(
    "connectDROMADatabase",
    "listDROMAProjects",
    "createMultiDromaSetFromAllProjects",
    "batchFindSignificantFeatures",
    "getSignificantFeatures",
    "getIntersectSignificantFeatures",
    "batchFindTcgaADConcordantFeatures",
    "connectCTRDatabase",
    "batchFindClinicalSigResponse"
  )
  on.exit(rm(list = cleanup, envir = .GlobalEnv), add = TRUE)

  result <- runMetaWorkflow(
    drug = "DrugA",
    tumor_type = "Tumor A",
    output_base = output_base,
    override = TRUE,
    verbose = FALSE
  )

  stopifnot(identical(result$status[[1]], "success"))
  stopifnot(identical(calls$batch_sig, 2L))
  stopifnot(identical(calls$ad, 1L))
}

test_stage03_filters_concordant_rows_without_get_lookup <- function() {
  output_base <- local_tempdir("meta-output-")

  assign("connectDROMADatabase", function(...) textConnection("droma", open = "r"), envir = .GlobalEnv)
  assign(
    "listDROMAProjects",
    function(...) data.frame(
      dataset_type = c("CellLine", "PDC", "PDO", "PDX", "Other"),
      project_name = c("CELL", "PDC1", "PDO1", "PDX1", "CCLE"),
      stringsAsFactors = FALSE
    ),
    envir = .GlobalEnv
  )
  assign(
    "createMultiDromaSetFromAllProjects",
    function(..., include_projects, con = NULL) {
      structure(list(projects = include_projects), class = "MultiDromaSet")
    },
    envir = .GlobalEnv
  )
  assign("batchFindSignificantFeatures", function(...) data.table(name = "GENE1"), envir = .GlobalEnv)
  assign("getSignificantFeatures", function(x, ...) as.data.table(x), envir = .GlobalEnv)
  assign(
    "getIntersectSignificantFeatures",
    function(...) {
      args <- list(...)
      if (!is.null(args$cell) && !is.null(args$pdcpdx) && is.null(args$ctrdb)) {
        return(data.table(name = c("GENE1", "GENE2")))
      }
      data.frame(name = "GENE1", direction = "up", stringsAsFactors = FALSE)
    },
    envir = .GlobalEnv
  )
  assign(
    "batchFindTcgaADConcordantFeatures",
    function(selected_features, ...) {
      out <- as.data.frame(selected_features, stringsAsFactors = FALSE)
      out$ccle_vs_tcga_concordant <- c(TRUE, FALSE)
      out
    },
    envir = .GlobalEnv
  )
  assign("connectCTRDatabase", function(...) TRUE, envir = .GlobalEnv)
  assign(
    "batchFindClinicalSigResponse",
    function(...) data.frame(name = "GENE1", direction = "up", stringsAsFactors = FALSE),
    envir = .GlobalEnv
  )

  cleanup <- c(
    "connectDROMADatabase",
    "listDROMAProjects",
    "createMultiDromaSetFromAllProjects",
    "batchFindSignificantFeatures",
    "getSignificantFeatures",
    "getIntersectSignificantFeatures",
    "batchFindTcgaADConcordantFeatures",
    "connectCTRDatabase",
    "batchFindClinicalSigResponse"
  )
  on.exit(rm(list = cleanup, envir = .GlobalEnv), add = TRUE)

  result <- runMetaWorkflow(
    drug = "DrugA",
    tumor_type = "Tumor A",
    output_base = output_base,
    override = TRUE,
    verbose = FALSE
  )

  stopifnot(identical(result$status[[1]], "success"))
  stopifnot(identical(result$n_ad_filtered[[1]], 1L))

  output_dir <- getMetaWorkflowOutputDir("DrugA", "Tumor A", output_base)
  selected_genes_ad_filtered <- fread(file.path(output_dir, "selected_genes_ad_filtered.csv"))
  stopifnot(nrow(selected_genes_ad_filtered) == 1L)
  stopifnot(identical(selected_genes_ad_filtered$name[[1]], "GENE1"))
}

test_skip_existing_outputs()
test_force_override_recomputes_outputs()
test_stage03_filters_concordant_rows_without_get_lookup()

cat("runMetaWorkflow override tests passed\n")
