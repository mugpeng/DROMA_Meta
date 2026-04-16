suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DROMA.Meta))

# project_root <- file.path(normalizePath(getwd(), mustWork = TRUE), "Meta_Example")
project_root <- normalizePath(getwd(), mustWork = TRUE)
defaults <- getMetaWorkflowDefaults(project_root = project_root)

# Driver-level input and output locations. Keep these runtime choices here so
# the files under R/ remain reusable pure function definitions.
output_base_batch <- file.path(defaults$output_base, "meta_batch")
summary_csv <- file.path(output_base_batch, "meta_workflow_batch_summary.csv")
annotation_summary_csv <- file.path(output_base_batch, "pair_biomarker_annotation_summary.csv")

# Batch-annotation parameters. Edit these values in the script when you want to
# change how the annotation summary is built.
pdcpdx_es_t <- 0.1
pdcpdx_P_t <- 0.05
pdcpdx_ge_2_n_datasets_t <- 2

summary_dt <- data.table::fread(summary_csv, data.table = TRUE)
batch_results <- vector("list", nrow(summary_dt))

for (i in seq_len(nrow(summary_dt))) {
  summary_row <- summary_dt[i]
  output_dir <- if ("output_dir" %in% names(summary_row) &&
    !is.na(summary_row$output_dir[[1]]) &&
    nzchar(summary_row$output_dir[[1]])) {
    summary_row$output_dir[[1]]
  } else {
    getMetaWorkflowOutputDir(
      drug = summary_row$drug[[1]],
      tumor_type = summary_row$tumor_type[[1]],
      output_base = output_base_batch
    )
  }

  final_path <- file.path(output_dir, "final_biomarkers.csv")
  if (!file.exists(final_path)) {
    next
  }

  annotated_path <- file.path(output_dir, "final_biomarkers_annotated.csv")
  batch_pdcpdx_path <- file.path(output_dir, "batch_pdcpdx_sets_mRNA.csv")
  clinical_info_path <- file.path(output_dir, "clinical_validation_info.csv")

  final_biomarkers <- data.table::fread(final_path, data.table = TRUE)
  clinical_info <- if (file.exists(clinical_info_path)) {
    data.table::fread(clinical_info_path, data.table = TRUE)
  } else {
    NULL
  }
  batch_pdcpdx <- if (file.exists(batch_pdcpdx_path)) {
    data.table::fread(batch_pdcpdx_path, data.table = TRUE)
  } else {
    NULL
  }

  ctrdb_status_value <- if (!is.null(clinical_info) &&
    nrow(clinical_info) > 0 &&
    "ctrdb_status" %in% names(clinical_info) &&
    !is.na(clinical_info$ctrdb_status[[1]]) &&
    nzchar(as.character(clinical_info$ctrdb_status[[1]]))) {
    as.character(clinical_info$ctrdb_status[[1]])
  } else if ("ctrdb_status" %in% names(summary_row) &&
    !is.na(summary_row$ctrdb_status[[1]]) &&
    nzchar(as.character(summary_row$ctrdb_status[[1]]))) {
    as.character(summary_row$ctrdb_status[[1]])
  } else {
    NA_character_
  }

  ctrdb_fallback_value <- if (!is.null(clinical_info) &&
    nrow(clinical_info) > 0 &&
    "ctrdb_fallback" %in% names(clinical_info) &&
    !is.na(clinical_info$ctrdb_fallback[[1]])) {
    as.logical(clinical_info$ctrdb_fallback[[1]])
  } else if ("ctrdb_fallback" %in% names(summary_row) &&
    !is.na(summary_row$ctrdb_fallback[[1]])) {
    as.logical(summary_row$ctrdb_fallback[[1]])
  } else {
    NA
  }

  pdcpdx_ge_2_names <- character(0)
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

  final_biomarkers[, ctrdb_status := ctrdb_status_value]
  final_biomarkers[, ctrdb_fallback := ctrdb_fallback_value]
  final_biomarkers[, pdcpdx_ge_2 := "name" %in% names(final_biomarkers) & name %in% pdcpdx_ge_2_names]
  data.table::fwrite(final_biomarkers, annotated_path)

  drug_value <- if ("drug" %in% names(final_biomarkers) && nrow(final_biomarkers) > 0) {
    final_biomarkers$drug[[1]]
  } else {
    summary_row$drug[[1]]
  }
  tumor_value <- if ("tumor_type" %in% names(final_biomarkers) && nrow(final_biomarkers) > 0) {
    final_biomarkers$tumor_type[[1]]
  } else {
    summary_row$tumor_type[[1]]
  }

  batch_results[[i]] <- data.table::data.table(
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

annotation_dt <- if (length(Filter(Negate(is.null), batch_results))) {
  data.table::rbindlist(batch_results, fill = TRUE)
} else {
  data.table::data.table()
}
data.table::fwrite(annotation_dt, annotation_summary_csv)
print(annotation_dt)
