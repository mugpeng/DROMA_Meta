suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(DROMA.Meta))
suppressPackageStartupMessages(library(DROMA.R))

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "two_bk.R"
script_dir <- dirname(normalizePath(script_path, mustWork = TRUE))
project_root <- normalizePath(file.path(script_dir, "..", "Meta_Example"), mustWork = TRUE)
defaults <- getMetaWorkflowDefaults(project_root = project_root)

# Driver-level input and output locations. This backup driver can annotate a
# meta_batch directory even when meta_workflow_batch_summary.csv was not written.
eligible_pairs_csv <- file.path(defaults$output_base, "eligible_drug_tumor_pairs.csv")
output_base_batch <- file.path(defaults$output_base, "meta_batch")
summary_csv <- file.path(output_base_batch, "meta_workflow_batch_summary.csv")
annotation_summary_csv <- file.path(output_base_batch, "pair_biomarker_annotation_summary_bk.csv")

# Batch-annotation parameters. Edit these values in the script when you want to
# change how the annotation summary is built.
pdcpdx_es_t <- 0.1
pdcpdx_P_t <- 0.05
pdcpdx_ge_2_n_datasets_t <- 2

if (file.exists(summary_csv)) {
  summary_dt <- data.table::fread(summary_csv, data.table = TRUE)
  if (!"output_dir" %in% names(summary_dt)) {
    summary_dt[
      ,
      output_dir := mapply(
        getMetaWorkflowOutputDir,
        drug = drug,
        tumor_type = tumor_type,
        MoreArgs = list(output_base = output_base_batch),
        USE.NAMES = FALSE
      )
    ]
  }
} else if (file.exists(eligible_pairs_csv)) {
  summary_dt <- unique(data.table::fread(eligible_pairs_csv, data.table = TRUE))
  summary_dt[
    ,
    output_dir := mapply(
      getMetaWorkflowOutputDir,
      drug = drug,
      tumor_type = tumor_type,
      MoreArgs = list(output_base = output_base_batch),
      USE.NAMES = FALSE
    )
  ]
} else {
  drug_dirs <- list.dirs(output_base_batch, recursive = FALSE, full.names = TRUE)
  pair_dirs <- unlist(
    lapply(drug_dirs, list.dirs, recursive = FALSE, full.names = TRUE),
    use.names = FALSE
  )
  summary_dt <- data.table(
    drug = basename(dirname(pair_dirs)),
    tumor_type = gsub("_", " ", basename(pair_dirs), fixed = TRUE),
    output_dir = pair_dirs
  )
}

batch_results <- vector("list", nrow(summary_dt))

for (i in seq_len(nrow(summary_dt))) {
  summary_row <- summary_dt[i]
  output_dir <- summary_row$output_dir[[1]]
  final_path <- file.path(output_dir, "final_biomarkers.csv")
  annotated_path <- file.path(output_dir, "final_biomarkers_annotated.csv")
  batch_pdcpdx_path <- file.path(output_dir, "batch_pdcpdx_sets_mRNA.csv")
  clinical_info_path <- file.path(output_dir, "clinical_validation_info.csv")

  output_dir_exists <- dir.exists(output_dir)
  final_exists <- file.exists(final_path)
  final_readable <- FALSE
  wrote_annotated <- FALSE
  annotation_status <- NA_character_
  ctrdb_status_value <- NA_character_
  ctrdb_fallback_value <- NA
  final_biomarkers <- data.table()

  if (!output_dir_exists) {
    annotation_status <- "missing_output_dir"
  } else if (!final_exists) {
    annotation_status <- "missing_final_biomarkers"
  } else {
    final_biomarkers <- tryCatch(
      data.table::fread(final_path, data.table = TRUE),
      error = function(e) NULL
    )
    final_readable <- !is.null(final_biomarkers)

    if (!final_readable) {
      annotation_status <- "unreadable_final_biomarkers"
      final_biomarkers <- data.table()
    } else {
      clinical_info <- if (file.exists(clinical_info_path)) {
        tryCatch(data.table::fread(clinical_info_path, data.table = TRUE), error = function(e) NULL)
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
      if (file.exists(batch_pdcpdx_path)) {
        batch_pdcpdx <- tryCatch(data.table::fread(batch_pdcpdx_path, data.table = TRUE), error = function(e) NULL)
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
      }

      final_biomarkers[, ctrdb_status := ctrdb_status_value]
      final_biomarkers[, ctrdb_fallback := ctrdb_fallback_value]
      final_biomarkers[
        ,
        pdcpdx_ge_2 := if ("name" %in% names(final_biomarkers)) {
          name %in% pdcpdx_ge_2_names
        } else {
          logical(.N)
        }
      ]

      data.table::fwrite(final_biomarkers, annotated_path)
      wrote_annotated <- TRUE
      annotation_status <- if (nrow(final_biomarkers) > 0) {
        "annotated"
      } else {
        "empty_final_biomarkers"
      }
    }
  }

  batch_results[[i]] <- data.table::data.table(
    drug = summary_row$drug[[1]],
    tumor_type = summary_row$tumor_type[[1]],
    output_dir = output_dir,
    annotation_status = annotation_status,
    output_dir_exists = output_dir_exists,
    final_biomarkers_exists = final_exists,
    final_biomarkers_readable = final_readable,
    wrote_annotated = wrote_annotated,
    ctrdb_status = ctrdb_status_value,
    ctrdb_fallback = ctrdb_fallback_value,
    n_final_biomarkers = nrow(final_biomarkers),
    n_pdcpdx_ge_2_biomarkers = if ("pdcpdx_ge_2" %in% names(final_biomarkers)) {
      sum(final_biomarkers$pdcpdx_ge_2, na.rm = TRUE)
    } else {
      0L
    },
    has_pdcpdx_ge_2_biomarker = if ("pdcpdx_ge_2" %in% names(final_biomarkers)) {
      any(final_biomarkers$pdcpdx_ge_2, na.rm = TRUE)
    } else {
      FALSE
    }
  )
}

annotation_dt <- data.table::rbindlist(batch_results, fill = TRUE)
data.table::fwrite(annotation_dt, annotation_summary_csv)
print(annotation_dt)
