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
  fwrite(
    data.table(
      drug = "DrugA",
      tumor_type = "Tumor A",
      ctrdb_status = "tumor_type_matched",
      ctrdb_fallback = FALSE,
      clinical_query_tumor_type = "Tumor A",
      n_clinical_sig = 1L,
      n_final_biomarkers = 1L
    ),
    file.path(output_dir, "clinical_validation_info.csv")
  )

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
  stopifnot(identical(calls$batch_sig, 4L))
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

test_batch_probe_failure_returns_failed_summary <- function() {
  output_base <- local_tempdir("meta-output-")

  calls <- new.env(parent = emptyenv())
  calls$probe <- 0L
  calls$full <- 0L

  assign("connectDROMADatabase", function(...) textConnection("droma", open = "r"), envir = .GlobalEnv)
  assign(
    "listDROMAProjects",
    function(...) data.frame(
      dataset_type = c("CellLine", "PDC", "PDO", "PDX"),
      project_name = c("CELL", "PDC1", "PDO1", "PDX1"),
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
    function(..., test_top_n = NULL) {
      if (!is.null(test_top_n)) {
        calls$probe <- calls$probe + 1L
        if (calls$probe == 2L) {
          stop("No data available for the selected feature 1.", call. = FALSE)
        }
        return(data.table(name = "GENE1"))
      }

      calls$full <- calls$full + 1L
      data.table(name = "GENE1")
    },
    envir = .GlobalEnv
  )

  cleanup <- c(
    "connectDROMADatabase",
    "listDROMAProjects",
    "createMultiDromaSetFromAllProjects",
    "batchFindSignificantFeatures"
  )
  on.exit(rm(list = cleanup, envir = .GlobalEnv), add = TRUE)

  result <- runMetaWorkflow(
    drug = "DrugA",
    tumor_type = "Tumor A",
    output_base = output_base,
    override = TRUE,
    verbose = FALSE
  )

  stopifnot(identical(result$status[[1]], "failed"))
  stopifnot(grepl("batch_pdcpdx test failed", result$error_message[[1]], fixed = TRUE))
  stopifnot(identical(calls$probe, 2L))
  stopifnot(identical(calls$full, 0L))
}

test_batch_probe_is_silent <- function() {
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
  assign(
    "batchFindSignificantFeatures",
    function(..., test_top_n = NULL) {
      if (!is.null(test_top_n)) {
        cat("Analysis completed in 1 seconds\n")
        message("Found 0 significant associations out of 5 valid features.")
        return(data.table(name = "GENE1", direction = "up"))
      }
      data.table(name = "GENE1", direction = "up")
    },
    envir = .GlobalEnv
  )
  assign("getSignificantFeatures", function(x, ...) as.data.table(x), envir = .GlobalEnv)
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

  leaked_output <- capture.output(
    leaked_messages <- capture.output(
      result <- runMetaWorkflow(
        drug = "DrugA",
        tumor_type = "Tumor A",
        output_base = output_base,
        override = TRUE,
        verbose = FALSE
      ),
      type = "message"
    ),
    type = "output"
  )

  stopifnot(identical(result$status[[1]], "success"))
  stopifnot(length(leaked_output) == 0L)
  stopifnot(length(leaked_messages) == 0L)
}

test_clinical_fallback_to_all_marks_outputs <- function() {
  output_base <- local_tempdir("meta-output-")

  calls <- new.env(parent = emptyenv())
  calls$clinical_tumor_types <- character(0)

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
  assign("batchFindSignificantFeatures", function(...) data.table(name = "GENE1", direction = "up"), envir = .GlobalEnv)
  assign("getSignificantFeatures", function(x, ...) as.data.table(x), envir = .GlobalEnv)
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
      out <- as.data.table(selected_features)
      out[, ccle_vs_tcga_concordant := TRUE]
      out
    },
    envir = .GlobalEnv
  )
  assign("connectCTRDatabase", function(...) TRUE, envir = .GlobalEnv)
  assign(
    "batchFindClinicalSigResponse",
    function(..., tumor_type) {
      calls$clinical_tumor_types <- c(calls$clinical_tumor_types, tumor_type)
      if (!identical(tumor_type, "all")) {
        return(data.frame())
      }
      data.frame(name = "GENE1", direction = "up", stringsAsFactors = FALSE)
    },
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
  stopifnot(identical(result$ctrdb_status[[1]], "fallback_all"))
  stopifnot(isTRUE(result$ctrdb_fallback[[1]]))
  stopifnot(identical(calls$clinical_tumor_types, c("Tumor A", "all")))

  output_dir <- getMetaWorkflowOutputDir("DrugA", "Tumor A", output_base)
  final_biomarkers <- fread(file.path(output_dir, "final_biomarkers.csv"))
  stopifnot("ctrdb_status" %in% names(final_biomarkers))
  stopifnot("ctrdb_fallback" %in% names(final_biomarkers))
  stopifnot(identical(final_biomarkers$ctrdb_status[[1]], "fallback_all"))
  stopifnot(isTRUE(final_biomarkers$ctrdb_fallback[[1]]))
}

test_clinical_missing_everywhere_returns_empty_outputs <- function() {
  output_base <- local_tempdir("meta-output-")

  calls <- new.env(parent = emptyenv())
  calls$clinical_tumor_types <- character(0)

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
  assign("batchFindSignificantFeatures", function(...) data.table(name = "GENE1", direction = "up"), envir = .GlobalEnv)
  assign("getSignificantFeatures", function(x, ...) as.data.table(x), envir = .GlobalEnv)
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
      out <- as.data.table(selected_features)
      out[, ccle_vs_tcga_concordant := TRUE]
      out
    },
    envir = .GlobalEnv
  )
  assign("connectCTRDatabase", function(...) TRUE, envir = .GlobalEnv)
  assign(
    "batchFindClinicalSigResponse",
    function(..., tumor_type) {
      calls$clinical_tumor_types <- c(calls$clinical_tumor_types, tumor_type)
      data.frame()
    },
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
  stopifnot(identical(result$ctrdb_status[[1]], "no_ctrdb_data"))
  stopifnot(identical(result$n_clinical_sig[[1]], 0L))
  stopifnot(identical(result$n_final_biomarkers[[1]], 0L))
  stopifnot(identical(calls$clinical_tumor_types, c("Tumor A", "all")))

  output_dir <- getMetaWorkflowOutputDir("DrugA", "Tumor A", output_base)
  clinical_sig <- fread(file.path(output_dir, "clinical_sig_mRNA.csv"))
  final_biomarkers <- fread(file.path(output_dir, "final_biomarkers.csv"))
  stopifnot(nrow(clinical_sig) == 0L)
  stopifnot(nrow(final_biomarkers) == 0L)
  stopifnot("ctrdb_status" %in% names(final_biomarkers))
  stopifnot("ctrdb_fallback" %in% names(final_biomarkers))
}

test_clinical_tumor_type_error_is_marked <- function() {
  output_base <- local_tempdir("meta-output-")

  calls <- new.env(parent = emptyenv())
  calls$clinical_tumor_types <- character(0)

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
  assign("batchFindSignificantFeatures", function(...) data.table(name = "GENE1", direction = "up"), envir = .GlobalEnv)
  assign("getSignificantFeatures", function(x, ...) as.data.table(x), envir = .GlobalEnv)
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
      out <- as.data.table(selected_features)
      out[, ccle_vs_tcga_concordant := TRUE]
      out
    },
    envir = .GlobalEnv
  )
  assign("connectCTRDatabase", function(...) TRUE, envir = .GlobalEnv)
  assign(
    "batchFindClinicalSigResponse",
    function(..., tumor_type) {
      calls$clinical_tumor_types <- c(calls$clinical_tumor_types, tumor_type)
      stop("Database connection lost", call. = FALSE)
    },
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
  stopifnot(identical(result$ctrdb_status[[1]], "tumor_type_error"))
  stopifnot(identical(result$n_clinical_sig[[1]], 0L))
  stopifnot(identical(result$n_final_biomarkers[[1]], 0L))
  stopifnot(identical(calls$clinical_tumor_types, "Tumor A"))
}

test_clinical_all_error_is_marked <- function() {
  output_base <- local_tempdir("meta-output-")

  calls <- new.env(parent = emptyenv())
  calls$clinical_tumor_types <- character(0)

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
  assign("batchFindSignificantFeatures", function(...) data.table(name = "GENE1", direction = "up"), envir = .GlobalEnv)
  assign("getSignificantFeatures", function(x, ...) as.data.table(x), envir = .GlobalEnv)
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
      out <- as.data.table(selected_features)
      out[, ccle_vs_tcga_concordant := TRUE]
      out
    },
    envir = .GlobalEnv
  )
  assign("connectCTRDatabase", function(...) TRUE, envir = .GlobalEnv)
  assign(
    "batchFindClinicalSigResponse",
    function(..., tumor_type) {
      calls$clinical_tumor_types <- c(calls$clinical_tumor_types, tumor_type)
      if (!identical(tumor_type, "all")) {
        return(data.frame())
      }
      stop("CTRDB query timeout", call. = FALSE)
    },
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
  stopifnot(identical(result$ctrdb_status[[1]], "all_error"))
  stopifnot(identical(result$n_clinical_sig[[1]], 0L))
  stopifnot(identical(result$n_final_biomarkers[[1]], 0L))
  stopifnot(identical(calls$clinical_tumor_types, c("Tumor A", "all")))
}

test_skip_existing_outputs()
test_force_override_recomputes_outputs()
test_stage03_filters_concordant_rows_without_get_lookup()
test_batch_probe_failure_returns_failed_summary()
test_batch_probe_is_silent()
test_clinical_fallback_to_all_marks_outputs()
test_clinical_missing_everywhere_returns_empty_outputs()
test_clinical_tumor_type_error_is_marked()
test_clinical_all_error_is_marked()

cat("runMetaWorkflow override tests passed\n")
