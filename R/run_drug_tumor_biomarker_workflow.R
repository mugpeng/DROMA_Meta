`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || identical(x, "") || (is.character(x) && !nzchar(x[1]))) y else x
}

.sanitize_name <- function(x) {
  out <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  out <- gsub("^_+|_+$", "", out)
  if (!nzchar(out)) "unknown" else out
}

.find_droma_repo_root <- function(start_dir) {
  current <- normalizePath(start_dir, mustWork = TRUE)

  repeat {
    if (file.exists(file.path(current, "Data", "droma.sqlite"))) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate DROMA repository root from: ", start_dir, call. = FALSE)
    }
    current <- parent
  }
}

.load_biomarker_runtime <- function(repo_root) {
  suppressPackageStartupMessages({
    library(data.table)
    library(DROMA.Set)
    library(DROMA.R)
  })

  droma_r_dir <- file.path(repo_root, "DB_project", "DROMA_R", "R")
  for (path in c(
    "FuncPairDataLoading.R",
    "FuncPairDataPairing.R",
    "FuncPairMetaAnalysis.R",
    "FuncPairBatchFeature.R",
    "FuncClinical.R"
  )) {
    source(file.path(droma_r_dir, path), local = FALSE)
  }
}

.resolve_biomarker_workflow_config <- function(drug,
                                               tumor_type,
                                               db_path = NULL,
                                               ctrdb_path = NULL,
                                               feature2_type = "mRNA",
                                               data_type = "all",
                                               cores = 3L,
                                               output_base = "Output",
                                               cell_min_intersected_cells = 20L,
                                               pdcpdx_min_intersected_cells = 8L,
                                               cell_es_t = 0.1,
                                               cell_P_t = 0.05,
                                               cell_n_datasets_t = 3L,
                                               pdcpdx_es_t = 0.1,
                                               pdcpdx_P_t = 0.05,
                                               pdcpdx_n_datasets_t = 2L,
                                               clinical_es_t = 0.05,
                                               clinical_P_t = 0.1,
                                               clinical_n_datasets_t = NULL,
                                               use_p_value = TRUE,
                                               workflow_dir = file.path(getwd(), "workflow")) {
  if (missing(drug) || !nzchar(drug)) {
    stop("drug must be a non-empty string", call. = FALSE)
  }
  if (missing(tumor_type) || !nzchar(tumor_type)) {
    stop("tumor_type must be a non-empty string", call. = FALSE)
  }
  if (!identical(feature2_type, "mRNA")) {
    stop("This workflow currently supports feature2_type = 'mRNA' only", call. = FALSE)
  }

  workflow_dir <- normalizePath(workflow_dir, mustWork = FALSE)
  project_root <- normalizePath(file.path(workflow_dir, ".."), mustWork = FALSE)
  repo_root <- .find_droma_repo_root(project_root)

  db_path <- normalizePath(db_path %||% file.path(repo_root, "Data", "droma.sqlite"), mustWork = FALSE)
  ctrdb_path <- normalizePath(ctrdb_path %||% file.path(repo_root, "Data", "ctrdb.sqlite"), mustWork = FALSE)

  if (!file.exists(db_path)) {
    stop("DROMA database not found: ", db_path, call. = FALSE)
  }
  if (!file.exists(ctrdb_path)) {
    stop("CTRDB database not found: ", ctrdb_path, call. = FALSE)
  }

  output_base <- if (grepl("^(/|~)", output_base)) {
    normalizePath(output_base, mustWork = FALSE)
  } else {
    normalizePath(file.path(workflow_dir, output_base), mustWork = FALSE)
  }

  tumor_slug <- .sanitize_name(tumor_type)
  drug_slug <- .sanitize_name(drug)
  output_dir <- file.path(output_base, drug, tumor_slug)

  list(
    drug = drug,
    drug_slug = drug_slug,
    tumor_type = tumor_type,
    tumor_slug = tumor_slug,
    db_path = db_path,
    ctrdb_path = ctrdb_path,
    feature2_type = feature2_type,
    data_type = data_type,
    cores = as.integer(cores),
    output_base = output_base,
    workflow_dir = workflow_dir,
    project_root = project_root,
    repo_root = repo_root,
    output_dir = output_dir,
    batch_cell_csv = file.path(output_dir, "batch_cell_sets_mRNA.csv"),
    batch_pdcpdx_csv = file.path(output_dir, "batch_pdcpdx_sets_mRNA.csv"),
    cell_sig_csv = file.path(output_dir, "mRNA_cell_sig.csv"),
    pdcpdx_sig_csv = file.path(output_dir, "mRNA_pdcpdx_sig.csv"),
    clinical_batch_csv = file.path(output_dir, "clinical_batch_mRNA.csv"),
    clinical_sig_csv = file.path(output_dir, "clinical_sig_mRNA.csv"),
    final_biomarkers_csv = file.path(output_dir, "final_biomarkers.csv"),
    batch_cell_rds = file.path(output_dir, "batch_cell_sets_mRNA.rds"),
    batch_pdcpdx_rds = file.path(output_dir, "batch_pdcpdx_sets_mRNA.rds"),
    cell_sig_rds = file.path(output_dir, "mRNA_cell_sig.rds"),
    pdcpdx_sig_rds = file.path(output_dir, "mRNA_pdcpdx_sig.rds"),
    clinical_batch_rds = file.path(output_dir, "clinical_batch_mRNA.rds"),
    clinical_sig_rds = file.path(output_dir, "clinical_sig_mRNA.rds"),
    final_biomarkers_rds = file.path(output_dir, "final_biomarkers.rds"),
    cell_min_intersected_cells = as.integer(cell_min_intersected_cells),
    pdcpdx_min_intersected_cells = as.integer(pdcpdx_min_intersected_cells),
    cell_es_t = cell_es_t,
    cell_P_t = cell_P_t,
    cell_n_datasets_t = as.integer(cell_n_datasets_t),
    pdcpdx_es_t = pdcpdx_es_t,
    pdcpdx_P_t = pdcpdx_P_t,
    pdcpdx_n_datasets_t = as.integer(pdcpdx_n_datasets_t),
    clinical_es_t = clinical_es_t,
    clinical_P_t = clinical_P_t,
    clinical_n_datasets_t = clinical_n_datasets_t,
    use_p_value = use_p_value
  )
}

.empty_meta_df <- function() {
  data.frame(
    name = character(0),
    effect_size = numeric(0),
    p_value = numeric(0),
    q_value = numeric(0),
    n_datasets = integer(0),
    stringsAsFactors = FALSE
  )
}

.read_or_empty <- function(path) {
  if (file.exists(path)) readRDS(path) else .empty_meta_df()
}

.run_batch_preclinical <- function(config) {
  .load_biomarker_runtime(config$repo_root)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  con <- DROMA.Set::connectDROMADatabase(config$db_path)
  on.exit(try(close(con), silent = TRUE), add = TRUE)

  project_anno <- DROMA.Set::listDROMAProjects()
  cell_names <- project_anno[project_anno$dataset_type %in% c("CellLine", "PDC"), "project_name"]
  pdcpdx_names <- project_anno[project_anno$dataset_type %in% c("PDO", "PDX"), "project_name"]

  if (length(cell_names) == 0) {
    stop("No CellLine/PDC projects found for cell_sets", call. = FALSE)
  }
  if (length(pdcpdx_names) == 0) {
    stop("No PDO/PDX projects found for pdcpdx_sets", call. = FALSE)
  }

  cell_sets <- DROMA.Set::createMultiDromaSetFromAllProjects(
    db_path = config$db_path,
    include_projects = cell_names,
    con = con
  )
  pdcpdx_sets <- DROMA.Set::createMultiDromaSetFromAllProjects(
    db_path = config$db_path,
    include_projects = pdcpdx_names,
    con = con
  )

  batch_cell <- DROMA.R::batchFindSignificantFeatures(
    dromaset_object = cell_sets,
    feature1_type = "drug",
    feature1_name = config$drug,
    feature2_type = config$feature2_type,
    data_type = config$data_type,
    tumor_type = config$tumor_type,
    overlap_only = FALSE,
    cores = config$cores,
    min_intersected_cells = config$cell_min_intersected_cells
  )

  batch_pdcpdx <- DROMA.R::batchFindSignificantFeatures(
    dromaset_object = pdcpdx_sets,
    feature1_type = "drug",
    feature1_name = config$drug,
    feature2_type = config$feature2_type,
    data_type = config$data_type,
    tumor_type = config$tumor_type,
    overlap_only = FALSE,
    cores = config$cores,
    min_intersected_cells = config$pdcpdx_min_intersected_cells
  )

  data.table::fwrite(batch_cell, config$batch_cell_csv)
  data.table::fwrite(batch_pdcpdx, config$batch_pdcpdx_csv)
  saveRDS(batch_cell, config$batch_cell_rds)
  saveRDS(batch_pdcpdx, config$batch_pdcpdx_rds)

  list(batch_cell = batch_cell, batch_pdcpdx = batch_pdcpdx)
}

.run_select_preclinical <- function(config,
                                    batch_cell = NULL,
                                    batch_pdcpdx = NULL) {
  .load_biomarker_runtime(config$repo_root)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  batch_cell <- batch_cell %||% .read_or_empty(config$batch_cell_rds)
  batch_pdcpdx <- batch_pdcpdx %||% .read_or_empty(config$batch_pdcpdx_rds)

  if (!nrow(batch_cell)) {
    warning("batch_cell is empty; writing empty significant result", call. = FALSE)
  }
  if (!nrow(batch_pdcpdx)) {
    warning("batch_pdcpdx is empty; writing empty significant result", call. = FALSE)
  }

  cell_sig <- if (nrow(batch_cell) > 0) {
    DROMA.R::getSignificantFeatures(
      meta_df = batch_cell,
      es_t = config$cell_es_t,
      P_t = config$cell_P_t,
      use_p_value = config$use_p_value,
      n_datasets_t = config$cell_n_datasets_t
    )
  } else {
    data.frame(name = character(0), stringsAsFactors = FALSE)
  }

  pdcpdx_sig <- if (nrow(batch_pdcpdx) > 0) {
    DROMA.R::getSignificantFeatures(
      meta_df = batch_pdcpdx,
      es_t = config$pdcpdx_es_t,
      P_t = config$pdcpdx_P_t,
      use_p_value = config$use_p_value,
      n_datasets_t = config$pdcpdx_n_datasets_t
    )
  } else {
    data.frame(name = character(0), stringsAsFactors = FALSE)
  }

  data.table::fwrite(cell_sig, config$cell_sig_csv)
  data.table::fwrite(pdcpdx_sig, config$pdcpdx_sig_csv)
  saveRDS(cell_sig, config$cell_sig_rds)
  saveRDS(pdcpdx_sig, config$pdcpdx_sig_rds)

  list(cell_sig = cell_sig, pdcpdx_sig = pdcpdx_sig)
}

.run_clinical_validation <- function(config,
                                     pdcpdx_sig = NULL,
                                     cell_sig = NULL) {
  .load_biomarker_runtime(config$repo_root)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  pdcpdx_sig <- pdcpdx_sig %||% .read_or_empty(config$pdcpdx_sig_rds)
  cell_sig <- cell_sig %||% .read_or_empty(config$cell_sig_rds)

  if (!nrow(pdcpdx_sig)) {
    warning("mRNA_pdcpdx_sig is empty; clinical validation will return empty outputs", call. = FALSE)
    clinical_batch <- .empty_meta_df()
    clinical_sig <- data.frame(name = character(0), stringsAsFactors = FALSE)
    final_biomarkers <- data.frame(name = character(0), stringsAsFactors = FALSE)
  } else {
    DROMA.Set::connectCTRDatabase(config$ctrdb_path)
    on.exit(
      {
        if (exists("ctrdb_connection", envir = .GlobalEnv)) {
          rm("ctrdb_connection", envir = .GlobalEnv)
        }
      },
      add = TRUE
    )

    clinical_batch <- tryCatch(
      DROMA.R::batchFindClinicalSigResponse(
        select_omics = unique(as.character(pdcpdx_sig$name)),
        select_drugs = config$drug,
        data_type = config$data_type,
        tumor_type = config$tumor_type,
        cores = config$cores
      ),
      error = function(e) {
        warning("Clinical validation returned no usable result: ", e$message, call. = FALSE)
        .empty_meta_df()
      }
    )

    if (!nrow(clinical_batch)) {
      warning("Clinical batch result is empty", call. = FALSE)
      clinical_sig <- data.frame(name = character(0), stringsAsFactors = FALSE)
      final_biomarkers <- data.frame(name = character(0), stringsAsFactors = FALSE)
    } else {
      clinical_sig <- DROMA.R::getSignificantFeatures(
        meta_df = clinical_batch,
        es_t = config$clinical_es_t,
        P_t = config$clinical_P_t,
        use_p_value = config$use_p_value,
        n_datasets_t = config$clinical_n_datasets_t
      )

      if (!nrow(clinical_sig)) {
        warning("No significant CTRDB biomarkers passed the configured thresholds", call. = FALSE)
        final_biomarkers <- data.frame(name = character(0), stringsAsFactors = FALSE)
      } else {
        pdcpdx_for_intersect <- as.data.frame(pdcpdx_sig)
        if (!"direction" %in% names(pdcpdx_for_intersect)) {
          stop("pdcpdx_sig must contain a direction column", call. = FALSE)
        }
        names(pdcpdx_for_intersect)[names(pdcpdx_for_intersect) == "direction"] <- "direction_pdcpdx"

        final_biomarkers <- DROMA.R::getIntersectSignificantFeatures(
          pdcpdx = pdcpdx_for_intersect,
          ctrdb = as.data.frame(clinical_sig),
          direction_cols = c(pdcpdx = "direction_pdcpdx", ctrdb = "direction")
        )

        if (nrow(final_biomarkers) > 0) {
          final_biomarkers$drug <- config$drug
          final_biomarkers$tumor_type <- config$tumor_type
          final_biomarkers$cell_supported <- final_biomarkers$name %in% cell_sig$name
        }
      }
    }
  }

  data.table::fwrite(clinical_batch, config$clinical_batch_csv)
  data.table::fwrite(clinical_sig, config$clinical_sig_csv)
  data.table::fwrite(final_biomarkers, config$final_biomarkers_csv)
  saveRDS(clinical_batch, config$clinical_batch_rds)
  saveRDS(clinical_sig, config$clinical_sig_rds)
  saveRDS(final_biomarkers, config$final_biomarkers_rds)

  list(
    clinical_batch = clinical_batch,
    clinical_sig = clinical_sig,
    final_biomarkers = final_biomarkers
  )
}

run_drug_tumor_biomarker_workflow <- function(drug,
                                              tumor_type,
                                              db_path = NULL,
                                              ctrdb_path = NULL,
                                              feature2_type = "mRNA",
                                              data_type = "all",
                                              cores = 3L,
                                              output_base = "Output",
                                              cell_min_intersected_cells = 20L,
                                              pdcpdx_min_intersected_cells = 8L,
                                              cell_es_t = 0.1,
                                              cell_P_t = 0.05,
                                              cell_n_datasets_t = 3L,
                                              pdcpdx_es_t = 0.1,
                                              pdcpdx_P_t = 0.05,
                                              pdcpdx_n_datasets_t = 2L,
                                              clinical_es_t = 0.05,
                                              clinical_P_t = 0.1,
                                              clinical_n_datasets_t = NULL,
                                              use_p_value = TRUE,
                                              workflow_dir = file.path(getwd(), "workflow")) {
  config <- .resolve_biomarker_workflow_config(
    drug = drug,
    tumor_type = tumor_type,
    db_path = db_path,
    ctrdb_path = ctrdb_path,
    feature2_type = feature2_type,
    data_type = data_type,
    cores = cores,
    output_base = output_base,
    cell_min_intersected_cells = cell_min_intersected_cells,
    pdcpdx_min_intersected_cells = pdcpdx_min_intersected_cells,
    cell_es_t = cell_es_t,
    cell_P_t = cell_P_t,
    cell_n_datasets_t = cell_n_datasets_t,
    pdcpdx_es_t = pdcpdx_es_t,
    pdcpdx_P_t = pdcpdx_P_t,
    pdcpdx_n_datasets_t = pdcpdx_n_datasets_t,
    clinical_es_t = clinical_es_t,
    clinical_P_t = clinical_P_t,
    clinical_n_datasets_t = clinical_n_datasets_t,
    use_p_value = use_p_value,
    workflow_dir = workflow_dir
  )

  batch_results <- .run_batch_preclinical(config)
  select_results <- .run_select_preclinical(
    config,
    batch_cell = batch_results$batch_cell,
    batch_pdcpdx = batch_results$batch_pdcpdx
  )
  clinical_results <- .run_clinical_validation(
    config,
    pdcpdx_sig = select_results$pdcpdx_sig,
    cell_sig = select_results$cell_sig
  )

  c(batch_results, select_results, clinical_results, list(config = config))
}
