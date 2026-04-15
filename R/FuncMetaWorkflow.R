# Meta workflow orchestration helpers ----

#' Get Default Paths for the Meta Workflow
#'
#' @description Returns package-level default paths for database files and the
#' workflow output directory. Batch input CSV paths are intentionally excluded
#' and should be configured by driver scripts such as `one.R`.
#' @param project_root Root directory of the DROMA.Meta project.
#' @return A named list of default paths.
#' @export
getMetaWorkflowDefaults <- function(project_root = "/Users/peng/Desktop/Project/DROMA/Meta_project3") {
  list(
    project_root = project_root,
    droma_db_path = "/Users/peng/Desktop/Project/DROMA/Data/droma.sqlite",
    ctrdb_path = "/Users/peng/Desktop/Project/DROMA/Data/ctrdb.sqlite",
    tcga_rna_counts_dir = "/Users/peng/Library/CloudStorage/OneDrive-Personal/28PHD_peng/250301-DROMA_project/archive260314/251112-DROMA_align/benchmark_mini/Input/TCGA/rna_counts",
    gene_probe_map_path = "/Users/peng/Desktop/Project/DROMA/Data/gencode.human.v49.annotation.gene.probeMap",
    output_base = file.path(project_root, "workflow", "Output")
  )
}

#' Get Output Directory for One Drug and Tumor-Type Pair
#'
#' @description Builds the standard output directory path used by the workflow.
#' @param drug Drug name.
#' @param tumor_type Tumor type.
#' @param output_base Base workflow output directory.
#' @return A single path string.
#' @export
getMetaWorkflowOutputDir <- function(drug,
                                     tumor_type,
                                     output_base = getMetaWorkflowDefaults()$output_base) {
  file.path(output_base, drug, sanitizeName(tumor_type))
}

#' Read a Single-Column CSV File
#'
#' @description Reads a simple CSV file and returns one named column as a unique
#' character vector.
#' @param csv_path Path to the CSV file.
#' @param column_name Optional column name. Defaults to the first column.
#' @return A unique character vector.
#' @export
readSingleColumnCsv <- function(csv_path, column_name = NULL) {
  dt <- utils::read.csv(
    csv_path,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  if (!nrow(dt)) {
    return(character(0))
  }

  if (is.null(column_name)) {
    column_name <- names(dt)[1]
  }

  if (!column_name %in% names(dt)) {
    stop(sprintf("Column '%s' not found in %s", column_name, csv_path), call. = FALSE)
  }

  values <- as.character(dt[[column_name]])
  unique(values[!is.na(values) & nzchar(values)])
}

#' Build the Drug and Tumor-Type Combination Grid
#'
#' @description Reads the valid-drug and valid-tumor-type CSV files and returns
#' their Cartesian product as a `data.table`.
#' @param valid_drugs_csv Path to `valid_drugs.csv`.
#' @param valid_tumor_types_csv Path to `valid_tumor_types.csv`.
#' @return A `data.table` with `drug` and `tumor_type` columns.
#' @export
buildDrugTumorGrid <- function(valid_drugs_csv, valid_tumor_types_csv) {
  drugs <- readSingleColumnCsv(valid_drugs_csv, "drug")
  tumor_types <- readSingleColumnCsv(valid_tumor_types_csv, "tumor_type")

  if (!length(drugs)) {
    stop("No valid drugs found in valid_drugs.csv", call. = FALSE)
  }
  if (!length(tumor_types)) {
    stop("No valid tumor types found in valid_tumor_types.csv", call. = FALSE)
  }

  data.table::CJ(
    drug = drugs,
    tumor_type = tumor_types,
    unique = TRUE
  )
}

#' Run the Full Meta Workflow for One Drug and Tumor Type
#'
#' @description Executes the four workflow stages for a single `drug` and
#' `tumor_type` pair, writes intermediate outputs, and returns a one-row summary.
#' @param drug Drug name.
#' @param tumor_type Tumor type.
#' @param feature2_type Molecular feature type, typically `"mRNA"`.
#' @param data_type Data type forwarded into DROMA workflow calls.
#' @param cores Number of cores used by downstream functions.
#' @param cell_min_intersected_cells Minimum intersected cells threshold for cell/PDC.
#' @param pdcpdx_min_intersected_cells Minimum intersected cells threshold for PDO/PDX.
#' @param cell_es_t Cell-line effect-size threshold.
#' @param cell_P_t Cell-line p-value threshold.
#' @param cell_n_datasets_t Cell-line dataset-count threshold.
#' @param pdcpdx_es_t PDO/PDX effect-size threshold.
#' @param pdcpdx_P_t PDO/PDX p-value threshold.
#' @param pdcpdx_n_datasets_t PDO/PDX dataset-count threshold.
#' @param tcga_ad_p_t Anderson-Darling concordance threshold.
#' @param clinical_es_t Clinical effect-size threshold.
#' @param clinical_P_t Clinical p-value threshold.
#' @param clinical_n_datasets_t Clinical dataset-count threshold.
#' @param droma_db_path Path to the DROMA SQLite database.
#' @param ctrdb_path Path to the CTRDB SQLite database.
#' @param tcga_rna_counts_dir Directory containing TCGA/TARGET counts.
#' @param gene_probe_map_path Path to the probe-map file.
#' @param output_base Base directory for workflow outputs.
#' @param verbose Logical, whether to print progress messages.
#' @return A one-row `data.table` summarizing workflow status and counts.
#' @export
runMetaWorkflowOne <- function(drug,
                               tumor_type,
                               feature2_type = "mRNA",
                               data_type = "all",
                               cores = 3,
                               cell_min_intersected_cells = 20,
                               pdcpdx_min_intersected_cells = 8,
                               cell_es_t = 0.1,
                               cell_P_t = 0.05,
                               cell_n_datasets_t = 3,
                               pdcpdx_es_t = 0.1,
                               pdcpdx_P_t = 0.05,
                               pdcpdx_n_datasets_t = 2,
                               tcga_ad_p_t = 0.01,
                               clinical_es_t = 0.1,
                               clinical_P_t = 0.05,
                               clinical_n_datasets_t = NULL,
                               droma_db_path = getMetaWorkflowDefaults()$droma_db_path,
                               ctrdb_path = getMetaWorkflowDefaults()$ctrdb_path,
                               tcga_rna_counts_dir = getMetaWorkflowDefaults()$tcga_rna_counts_dir,
                               gene_probe_map_path = getMetaWorkflowDefaults()$gene_probe_map_path,
                               output_base = getMetaWorkflowDefaults()$output_base,
                               verbose = TRUE) {
  output_dir <- getMetaWorkflowOutputDir(
    drug = drug,
    tumor_type = tumor_type,
    output_base = output_base
  )
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  if (verbose) {
    cat("\n============================================================\n")
    cat("Running meta workflow\n")
    cat("drug:", drug, "\n")
    cat("tumor_type:", tumor_type, "\n")
    cat("output_dir:", output_dir, "\n")
  }

  result <- tryCatch(
    {
      if (verbose) {
        cat("\n=== 01: Batch Preclinical ===\n")
      }

      con <- connectDROMADatabase(droma_db_path)
      on.exit(try(close(con), silent = TRUE), add = TRUE)

      project_anno <- listDROMAProjects()

      cell_names_all <- project_anno[
        project_anno$dataset_type %in% c("CellLine", "PDC"),
        "project_name"
      ]
      cell_sets <- createMultiDromaSetFromAllProjects(
        db_path = droma_db_path,
        include_projects = cell_names_all,
        con = con
      )

      pdcpdx_names_all <- project_anno[
        project_anno$dataset_type %in% c("PDO", "PDX"),
        "project_name"
      ]
      pdcpdx_sets <- createMultiDromaSetFromAllProjects(
        db_path = droma_db_path,
        include_projects = pdcpdx_names_all,
        con = con
      )

      batch_cell <- batchFindSignificantFeatures(
        cell_sets,
        feature1_type = "drug",
        feature1_name = drug,
        feature2_type = feature2_type,
        data_type = data_type,
        tumor_type = tumor_type,
        overlap_only = FALSE,
        cores = cores,
        min_intersected_cells = cell_min_intersected_cells
      )
      data.table::fwrite(batch_cell, file.path(output_dir, "batch_cell_sets_mRNA.csv"))

      batch_pdcpdx <- batchFindSignificantFeatures(
        pdcpdx_sets,
        feature1_type = "drug",
        feature1_name = drug,
        feature2_type = feature2_type,
        data_type = data_type,
        tumor_type = tumor_type,
        overlap_only = FALSE,
        cores = cores,
        min_intersected_cells = pdcpdx_min_intersected_cells
      )
      data.table::fwrite(batch_pdcpdx, file.path(output_dir, "batch_pdcpdx_sets_mRNA.csv"))

      if (verbose) {
        cat(sprintf("  OK batch_cell: %d rows\n", nrow(batch_cell)))
        cat(sprintf("  OK batch_pdcpdx: %d rows\n", nrow(batch_pdcpdx)))
        cat("\n=== 02: Select Preclinical ===\n")
      }

      mRNA_cell_sig <- getSignificantFeatures(
        batch_cell,
        es_t = cell_es_t,
        P_t = cell_P_t,
        use_p_value = TRUE,
        n_datasets_t = cell_n_datasets_t
      )
      mRNA_pdcpdx_sig <- getSignificantFeatures(
        batch_pdcpdx,
        es_t = pdcpdx_es_t,
        P_t = pdcpdx_P_t,
        use_p_value = TRUE,
        n_datasets_t = pdcpdx_n_datasets_t
      )
      selected_genes <- getIntersectSignificantFeatures(
        cell = mRNA_cell_sig,
        pdcpdx = mRNA_pdcpdx_sig
      )

      data.table::fwrite(mRNA_cell_sig, file.path(output_dir, "mRNA_cell_sig.csv"))
      data.table::fwrite(mRNA_pdcpdx_sig, file.path(output_dir, "mRNA_pdcpdx_sig.csv"))
      data.table::fwrite(selected_genes, file.path(output_dir, "selected_genes.csv"))

      if (verbose) {
        cat(sprintf("  OK mRNA_cell_sig: %d biomarkers\n", nrow(mRNA_cell_sig)))
        cat(sprintf("  OK mRNA_pdcpdx_sig: %d biomarkers\n", nrow(mRNA_pdcpdx_sig)))
        cat(sprintf("  OK selected_genes: %d biomarkers\n", nrow(selected_genes)))
        cat("\n=== 03: TCGA AD Filter ===\n")
      }

      ccle_set <- createMultiDromaSetFromAllProjects(
        db_path = droma_db_path,
        include_projects = "CCLE",
        con = con
      )

      ad_stats <- batchFindTcgaADConcordantFeatures(
        selected_features = selected_genes,
        preclinical_set = ccle_set,
        tumor_type = tumor_type,
        tcga_rna_counts_dir = tcga_rna_counts_dir,
        gene_probe_map_path = gene_probe_map_path,
        feature_type = feature2_type,
        data_type = data_type,
        p_t = tcga_ad_p_t,
        preclinical_label = "ccle"
      )
      data.table::fwrite(ad_stats, file.path(output_dir, "selected_genes_ad_stats.csv"))

      selected_genes_ad_filtered <- ad_stats[get("ccle_vs_tcga_concordant") == TRUE]
      data.table::fwrite(
        selected_genes_ad_filtered,
        file.path(output_dir, "selected_genes_ad_filtered.csv")
      )

      if (verbose) {
        cat(sprintf("  OK selected_genes_ad_stats: %d biomarkers\n", nrow(ad_stats)))
        cat(sprintf("  OK selected_genes_ad_filtered: %d biomarkers\n", nrow(selected_genes_ad_filtered)))
        cat("\n=== 04: Clinical Validation ===\n")
      }

      connectCTRDatabase(ctrdb_path)

      clinical_batch <- tryCatch(
        batchFindClinicalSigResponse(
          select_omics = unique(selected_genes_ad_filtered$name),
          select_drugs = drug,
          data_type = data_type,
          tumor_type = tumor_type,
          cores = cores
        ),
        error = function(e) {
          warning("Clinical validation returned no usable result: ", e$message, call. = FALSE)
          createEmptyMetaDf()
        }
      )

      if (nrow(clinical_batch) > 0) {
        clinical_sig <- getSignificantFeatures(
          clinical_batch,
          es_t = clinical_es_t,
          P_t = clinical_P_t,
          use_p_value = TRUE,
          n_datasets_t = clinical_n_datasets_t
        )
      } else {
        clinical_sig <- data.frame(name = character(0), stringsAsFactors = FALSE)
      }

      if (nrow(selected_genes_ad_filtered) > 0 && nrow(clinical_sig) > 0) {
        final_biomarkers <- getIntersectSignificantFeatures(
          pdcpdx = as.data.frame(selected_genes_ad_filtered),
          ctrdb = as.data.frame(clinical_sig),
          direction_cols = c(pdcpdx = "direction_pdcpdx", ctrdb = "direction")
        )

        if (nrow(final_biomarkers) > 0) {
          final_biomarkers$drug <- drug
          final_biomarkers$tumor_type <- tumor_type

          dir_cols <- grep("direction", names(final_biomarkers), value = TRUE)
          if (length(dir_cols) > 0) {
            final_biomarkers$direction <- final_biomarkers[[tail(dir_cols, 1)]]
            cols_to_remove <- setdiff(dir_cols, "direction")
            if (length(cols_to_remove) > 0) {
              final_biomarkers[, cols_to_remove] <- NULL
            }
          }
        }
      } else {
        final_biomarkers <- data.frame(name = character(0), stringsAsFactors = FALSE)
      }

      data.table::fwrite(clinical_sig, file.path(output_dir, "clinical_sig_mRNA.csv"))
      data.table::fwrite(final_biomarkers, file.path(output_dir, "final_biomarkers.csv"))

      if (verbose) {
        cat(sprintf("  OK clinical_sig: %d biomarkers\n", nrow(clinical_sig)))
        cat(sprintf("  OK final_biomarkers: %d biomarkers\n", nrow(final_biomarkers)))
      }

      data.table::data.table(
        drug = drug,
        tumor_type = tumor_type,
        output_dir = output_dir,
        status = "success",
        error_message = NA_character_,
        n_batch_cell = nrow(batch_cell),
        n_batch_pdcpdx = nrow(batch_pdcpdx),
        n_selected_genes = nrow(selected_genes),
        n_ad_filtered = nrow(selected_genes_ad_filtered),
        n_clinical_sig = nrow(clinical_sig),
        n_final_biomarkers = nrow(final_biomarkers)
      )
    },
    error = function(e) {
      if (verbose) {
        cat("  FAILED:", conditionMessage(e), "\n")
      }

      data.table::data.table(
        drug = drug,
        tumor_type = tumor_type,
        output_dir = output_dir,
        status = "failed",
        error_message = conditionMessage(e),
        n_batch_cell = NA_integer_,
        n_batch_pdcpdx = NA_integer_,
        n_selected_genes = NA_integer_,
        n_ad_filtered = NA_integer_,
        n_clinical_sig = NA_integer_,
        n_final_biomarkers = NA_integer_
      )
    }
  )

  result
}

#' Run the Meta Workflow Over All Valid Combinations
#'
#' @description Builds the Cartesian product of valid drugs and tumor types, runs
#' the full workflow for each pair, and writes an optional batch summary table.
#' @param valid_drugs_csv Path to `valid_drugs.csv`.
#' @param valid_tumor_types_csv Path to `valid_tumor_types.csv`.
#' @param output_base Base directory for workflow outputs.
#' @param summary_csv Optional output path for the batch summary table.
#' @param verbose Logical, whether to print progress messages.
#' @param ... Additional arguments forwarded to `runMetaWorkflowOne()`.
#' @return A `data.table` summarizing all combinations.
#' @export
runMetaWorkflowBatch <- function(valid_drugs_csv,
                                 valid_tumor_types_csv,
                                 output_base = getMetaWorkflowDefaults()$output_base,
                                 summary_csv = NULL,
                                 verbose = TRUE,
                                 ...) {
  grid <- buildDrugTumorGrid(
    valid_drugs_csv = valid_drugs_csv,
    valid_tumor_types_csv = valid_tumor_types_csv
  )

  if (verbose) {
    cat("\n============================================================\n")
    cat("Batch meta workflow\n")
    cat("valid_drugs_csv:", valid_drugs_csv, "\n")
    cat("valid_tumor_types_csv:", valid_tumor_types_csv, "\n")
    cat("total combinations:", nrow(grid), "\n")
  }

  results <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    drug <- as.character(grid$drug[[i]])
    tumor_type <- as.character(grid$tumor_type[[i]])

    if (verbose) {
      cat("\n------------------------------------------------------------\n")
      cat(sprintf("[%d/%d] %s | %s\n", i, nrow(grid), drug, tumor_type))
    }

    results[[i]] <- runMetaWorkflowOne(
      drug = drug,
      tumor_type = tumor_type,
      output_base = output_base,
      verbose = verbose,
      ...
    )
  }

  summary_dt <- data.table::rbindlist(results, fill = TRUE)

  if (!is.null(summary_csv) && nzchar(summary_csv)) {
    dir.create(dirname(summary_csv), showWarnings = FALSE, recursive = TRUE)
    data.table::fwrite(summary_dt, summary_csv)
  }

  summary_dt
}
