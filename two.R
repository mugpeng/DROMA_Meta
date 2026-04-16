fallback_if_missing <- function(x, y) {
  if (is.null(x) || length(x) == 0 || !nzchar(as.character(x[[1]]))) y else x[[1]]
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- sub("^--file=", "", fallback_if_missing(file_arg[1], ""))
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path, mustWork = TRUE)) else getwd()
root_dir <- normalizePath(script_dir, mustWork = TRUE)

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DROMA.Meta))

source(file.path(root_dir, "R", "FuncHelper.R"))
source(file.path(root_dir, "R", "FuncMetaWorkflow.R"))

project_root <- file.path(normalizePath(getwd(), mustWork = TRUE), "Meta_Example")
defaults <- getMetaWorkflowDefaults(project_root = project_root)

pdcpdx_es_t <- 0.1
pdcpdx_P_t <- 0.05
pdcpdx_ge_2_n_datasets_t <- 2

batch_roots <- c(file.path(defaults$output_base, "meta_batch"))

listPairDirs <- function(batch_root) {
  drug_dirs <- list.dirs(batch_root, recursive = FALSE, full.names = TRUE)
  unlist(lapply(drug_dirs, function(drug_dir) {
    list.dirs(drug_dir, recursive = FALSE, full.names = TRUE)
  }), use.names = FALSE)
}

readMaybe <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  data.table::fread(path)
}

annotatePair <- function(output_dir, summary_dt = NULL) {
  final_path <- file.path(output_dir, "final_biomarkers.csv")
  annotated_path <- file.path(output_dir, "final_biomarkers_annotated.csv")
  batch_pdcpdx_path <- file.path(output_dir, "batch_pdcpdx_sets_mRNA.csv")
  clinical_info_path <- file.path(output_dir, "clinical_validation_info.csv")

  final_biomarkers <- readMaybe(final_path)
  if (is.null(final_biomarkers)) {
    return(NULL)
  }

  clinical_info <- readMaybe(clinical_info_path)
  summary_row <- NULL
  if (!is.null(summary_dt) && "output_dir" %in% names(summary_dt)) {
    summary_row <- summary_dt[summary_dt$output_dir == output_dir]
    if (!nrow(summary_row)) {
      summary_row <- NULL
    }
  }

  ctrdb_status_value <- NA_character_
  ctrdb_fallback_value <- NA

  if (!is.null(clinical_info) && nrow(clinical_info) > 0) {
    ctrdb_status_value <- fallback_if_missing(clinical_info$ctrdb_status, NA_character_)
    ctrdb_fallback_value <- fallback_if_missing(clinical_info$ctrdb_fallback, NA)
  } else if (!is.null(summary_row)) {
    ctrdb_status_value <- fallback_if_missing(summary_row$ctrdb_status, NA_character_)
    ctrdb_fallback_value <- fallback_if_missing(summary_row$ctrdb_fallback, NA)
  }

  pdcpdx_ge_2_names <- character(0)
  batch_pdcpdx <- readMaybe(batch_pdcpdx_path)
  if (!is.null(batch_pdcpdx) && nrow(batch_pdcpdx) > 0) {
    pdcpdx_sig_ge_2 <- getSignificantFeatures(
      batch_pdcpdx,
      es_t = pdcpdx_es_t,
      P_t = pdcpdx_P_t,
      use_p_value = TRUE,
      n_datasets_t = pdcpdx_ge_2_n_datasets_t
    )
    if (nrow(pdcpdx_sig_ge_2) > 0 && "name" %in% names(pdcpdx_sig_ge_2)) {
      pdcpdx_ge_2_names <- unique(as.character(pdcpdx_sig_ge_2$name))
    }
  }

  if (!"ctrdb_status" %in% names(final_biomarkers)) {
    final_biomarkers[, ctrdb_status := ctrdb_status_value]
  } else if (nrow(final_biomarkers) > 0) {
    final_biomarkers[, ctrdb_status := ctrdb_status_value]
  }

  if (!"ctrdb_fallback" %in% names(final_biomarkers)) {
    final_biomarkers[, ctrdb_fallback := ctrdb_fallback_value]
  } else if (nrow(final_biomarkers) > 0) {
    final_biomarkers[, ctrdb_fallback := ctrdb_fallback_value]
  }

  if ("name" %in% names(final_biomarkers)) {
    final_biomarkers[, pdcpdx_ge_2 := name %in% pdcpdx_ge_2_names]
  } else {
    final_biomarkers[, pdcpdx_ge_2 := logical(.N)]
  }

  data.table::fwrite(final_biomarkers, annotated_path)

  drug_value <- if ("drug" %in% names(final_biomarkers) && nrow(final_biomarkers) > 0) {
    final_biomarkers$drug[[1]]
  } else if (!is.null(clinical_info) && "drug" %in% names(clinical_info) && nrow(clinical_info) > 0) {
    clinical_info$drug[[1]]
  } else if (!is.null(summary_row)) {
    summary_row$drug[[1]]
  } else {
    basename(dirname(output_dir))
  }

  tumor_value <- if ("tumor_type" %in% names(final_biomarkers) && nrow(final_biomarkers) > 0) {
    final_biomarkers$tumor_type[[1]]
  } else if (!is.null(clinical_info) && "tumor_type" %in% names(clinical_info) && nrow(clinical_info) > 0) {
    clinical_info$tumor_type[[1]]
  } else if (!is.null(summary_row)) {
    summary_row$tumor_type[[1]]
  } else {
    basename(output_dir)
  }

  data.table(
    drug = drug_value,
    tumor_type = tumor_value,
    output_dir = output_dir,
    ctrdb_status = ctrdb_status_value,
    ctrdb_fallback = ctrdb_fallback_value,
    n_final_biomarkers = nrow(final_biomarkers),
    n_pdcpdx_ge_2_biomarkers = sum(final_biomarkers$pdcpdx_ge_2, na.rm = TRUE),
    has_pdcpdx_ge_2_biomarker = any(final_biomarkers$pdcpdx_ge_2, na.rm = TRUE)
  )
}

for (batch_root in batch_roots) {
  if (!dir.exists(batch_root)) {
    next
  }

  summary_path <- file.path(batch_root, "meta_workflow_batch_summary.csv")
  summary_dt <- readMaybe(summary_path)
  pair_dirs <- listPairDirs(batch_root)

  annotation_rows <- Filter(
    Negate(is.null),
    lapply(pair_dirs, annotatePair, summary_dt = summary_dt)
  )

  if (length(annotation_rows) == 0L) {
    next
  }

  annotation_dt <- rbindlist(annotation_rows, fill = TRUE)
  fwrite(annotation_dt, file.path(batch_root, "pair_biomarker_annotation_summary.csv"))
  print(annotation_dt)
}
