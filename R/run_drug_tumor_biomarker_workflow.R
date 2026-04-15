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

  source(file.path(repo_root, "Meta_project3", "R", "FuncTcgaAD.R"), local = FALSE)
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
                                               tcga_ad_p_t = 0.05,
                                               tcga_rna_counts_dir = file.path(
                                                 path.expand("~"),
                                                 "Library/CloudStorage/OneDrive-Personal/28PHD_peng/250301-DROMA_project/archive260314/251112-DROMA_align/benchmark_mini/Input/TCGA/rna_counts"
                                               ),
                                               gene_probe_map_path = NULL,
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
  gene_probe_map_path <- normalizePath(
    gene_probe_map_path %||% file.path(repo_root, "Data", "gencode.human.v49.annotation.gene.probeMap"),
    mustWork = FALSE
  )

  if (!file.exists(db_path)) {
    stop("DROMA database not found: ", db_path, call. = FALSE)
  }
  if (!file.exists(ctrdb_path)) {
    stop("CTRDB database not found: ", ctrdb_path, call. = FALSE)
  }

  tcga_rna_counts_dir <- normalizePath(tcga_rna_counts_dir, mustWork = FALSE)
  gene_probe_map_path <- normalizePath(gene_probe_map_path, mustWork = FALSE)
  if (!dir.exists(tcga_rna_counts_dir)) {
    stop("TCGA RNA counts directory not found: ", tcga_rna_counts_dir, call. = FALSE)
  }
  if (!file.exists(gene_probe_map_path)) {
    stop("Gene probe map not found: ", gene_probe_map_path, call. = FALSE)
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
    selected_genes_csv = file.path(output_dir, "selected_genes.csv"),
    ad_stats_csv = file.path(output_dir, "selected_genes_ad_stats.csv"),
    ad_filtered_csv = file.path(output_dir, "selected_genes_ad_filtered.csv"),
    clinical_batch_csv = file.path(output_dir, "clinical_batch_mRNA.csv"),
    clinical_sig_csv = file.path(output_dir, "clinical_sig_mRNA.csv"),
    final_biomarkers_csv = file.path(output_dir, "final_biomarkers.csv"),
    batch_cell_rds = file.path(output_dir, "batch_cell_sets_mRNA.rds"),
    batch_pdcpdx_rds = file.path(output_dir, "batch_pdcpdx_sets_mRNA.rds"),
    cell_sig_rds = file.path(output_dir, "mRNA_cell_sig.rds"),
    pdcpdx_sig_rds = file.path(output_dir, "mRNA_pdcpdx_sig.rds"),
    selected_genes_rds = file.path(output_dir, "selected_genes.rds"),
    ad_stats_rds = file.path(output_dir, "selected_genes_ad_stats.rds"),
    ad_filtered_rds = file.path(output_dir, "selected_genes_ad_filtered.rds"),
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
    tcga_ad_p_t = tcga_ad_p_t,
    tcga_rna_counts_dir = tcga_rna_counts_dir,
    gene_probe_map_path = gene_probe_map_path,
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

.get_valid_values <- function(project_anno,
                             anno_df,
                             dataset_types,
                             value_col,
                             min_project_count,
                             exclude_values = NULL) {
  group_projects <- unique(as.character(
    project_anno[project_anno$dataset_type %in% dataset_types, "project_name"]
  ))
  group_projects <- group_projects[!is.na(group_projects) & nzchar(group_projects)]

  anno_filtered <- anno_df[anno_df$ProjectID %in% group_projects &
                            !is.na(anno_df[[value_col]]) &
                            nzchar(anno_df[[value_col]]), ]

  counts <- table(anno_filtered[[value_col]])
  keep_values <- names(counts)[counts >= min_project_count]

  if (!is.null(exclude_values)) {
    keep_values <- setdiff(keep_values, exclude_values)
  }

  sort(keep_values)
}

.filter_projects_for_drug_tumor <- function(project_names,
                                            project_anno,
                                            drug_anno,
                                            sample_anno,
                                            dataset_types,
                                            drug,
                                            tumor_type,
                                            min_project_count) {
  group_projects <- unique(as.character(
    project_anno[project_anno$dataset_type %in% dataset_types, "project_name"]
  ))
  group_projects <- group_projects[!is.na(group_projects) & nzchar(group_projects)]

  drug_projects <- unique(as.character(
    drug_anno$ProjectID[!is.na(drug_anno$DrugName) & drug_anno$DrugName == drug]
  ))
  tumor_projects <- unique(as.character(
    sample_anno$ProjectID[!is.na(sample_anno$TumorType) & sample_anno$TumorType == tumor_type]
  ))

  keep_projects <- intersect(group_projects, intersect(drug_projects, tumor_projects))
  keep_projects <- keep_projects[!is.na(keep_projects) & nzchar(keep_projects)]

  if (length(keep_projects) < min_project_count) {
    stop(
      sprintf(
        paste0(
          "Not enough eligible projects for drug '%s' and tumor_type '%s' in [%s]. ",
          "Need at least %d project(s) to satisfy n_datasets_t, found %d."
        ),
        drug,
        tumor_type,
        paste(dataset_types, collapse = ","),
        min_project_count,
        length(keep_projects)
      ),
      call. = FALSE
    )
  }

  keep_projects
}

.get_valid_drugs_and_tumor_types <- function(project_anno,
                                             drug_anno,
                                             sample_anno,
                                             cell_n_datasets_t,
                                             pdcpdx_n_datasets_t) {
  valid_cell_drugs <- .get_valid_values(
    project_anno = project_anno,
    anno_df = drug_anno,
    dataset_types = c("CellLine", "PDC"),
    value_col = "DrugName",
    min_project_count = cell_n_datasets_t
  )
  valid_pdcpdx_drugs <- .get_valid_values(
    project_anno = project_anno,
    anno_df = drug_anno,
    dataset_types = c("PDO", "PDX"),
    value_col = "DrugName",
    min_project_count = pdcpdx_n_datasets_t
  )

  valid_cell_tumor_types <- .get_valid_values(
    project_anno = project_anno,
    anno_df = sample_anno,
    dataset_types = c("CellLine", "PDC"),
    value_col = "TumorType",
    min_project_count = cell_n_datasets_t,
    exclude_values = "non-cancer"
  )
  valid_pdcpdx_tumor_types <- .get_valid_values(
    project_anno = project_anno,
    anno_df = sample_anno,
    dataset_types = c("PDO", "PDX"),
    value_col = "TumorType",
    min_project_count = pdcpdx_n_datasets_t,
    exclude_values = "non-cancer"
  )

  list(
    valid_drugs = intersect(valid_cell_drugs, valid_pdcpdx_drugs),
    valid_tumor_types = intersect(valid_cell_tumor_types, valid_pdcpdx_tumor_types)
  )
}

.run_batch_preclinical <- function(config) {
  .load_biomarker_runtime(config$repo_root)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  con <- DROMA.Set::connectDROMADatabase(config$db_path)
  on.exit(try(close(con), silent = TRUE), add = TRUE)

  project_anno <- DROMA.Set::listDROMAProjects()
  drug_anno <- DROMA.Set::getDROMAAnnotation("drug")
  sample_anno <- DROMA.Set::getDROMAAnnotation("sample")

  valid_inputs <- .get_valid_drugs_and_tumor_types(
    project_anno = project_anno,
    drug_anno = drug_anno,
    sample_anno = sample_anno,
    cell_n_datasets_t = config$cell_n_datasets_t,
    pdcpdx_n_datasets_t = config$pdcpdx_n_datasets_t
  )
  valid_drugs <- valid_inputs$valid_drugs
  valid_tumor_types <- valid_inputs$valid_tumor_types

  if (!config$drug %in% valid_drugs) {
    stop("drug does not satisfy both cell_sets and pdcpdx_sets project-count requirements: ", config$drug, call. = FALSE)
  }
  if (!config$tumor_type %in% valid_tumor_types) {
    stop("tumor_type does not satisfy both cell_sets and pdcpdx_sets project-count requirements: ", config$tumor_type, call. = FALSE)
  }

  cell_names <- .filter_projects_for_drug_tumor(
    project_names = project_anno$project_name,
    project_anno = project_anno,
    drug_anno = drug_anno,
    sample_anno = sample_anno,
    dataset_types = c("CellLine", "PDC"),
    drug = config$drug,
    tumor_type = config$tumor_type,
    min_project_count = config$cell_n_datasets_t
  )
  pdcpdx_names <- .filter_projects_for_drug_tumor(
    project_names = project_anno$project_name,
    project_anno = project_anno,
    drug_anno = drug_anno,
    sample_anno = sample_anno,
    dataset_types = c("PDO", "PDX"),
    drug = config$drug,
    tumor_type = config$tumor_type,
    min_project_count = config$pdcpdx_n_datasets_t
  )

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

  writeLines(cell_names, file.path(config$output_dir, "cell_sets_projects.txt"))
  writeLines(pdcpdx_names, file.path(config$output_dir, "pdcpdx_sets_projects.txt"))

  list(
    batch_cell = batch_cell,
    batch_pdcpdx = batch_pdcpdx,
    cell_projects = cell_names,
    pdcpdx_projects = pdcpdx_names
  )
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

  selected_genes <- if (nrow(cell_sig) > 0 && nrow(pdcpdx_sig) > 0) {
    DROMA.R::getIntersectSignificantFeatures(
      cell = cell_sig,
      pdcpdx = pdcpdx_sig
    )
  } else {
    data.frame(name = character(0), stringsAsFactors = FALSE)
  }

  data.table::fwrite(cell_sig, config$cell_sig_csv)
  data.table::fwrite(pdcpdx_sig, config$pdcpdx_sig_csv)
  data.table::fwrite(selected_genes, config$selected_genes_csv)
  saveRDS(cell_sig, config$cell_sig_rds)
  saveRDS(pdcpdx_sig, config$pdcpdx_sig_rds)
  saveRDS(selected_genes, config$selected_genes_rds)

  list(cell_sig = cell_sig, pdcpdx_sig = pdcpdx_sig, selected_genes = selected_genes)
}

.run_tcga_ad <- function(config,
                         selected_genes = NULL) {
  .load_biomarker_runtime(config$repo_root)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  selected_genes <- selected_genes %||% .read_or_empty(config$selected_genes_rds)
  if (!nrow(selected_genes)) {
    warning("selected_genes is empty; writing empty TCGA AD outputs", call. = FALSE)
    ad_stats <- data.frame(name = character(0), stringsAsFactors = FALSE)
    ad_filtered <- data.frame(name = character(0), stringsAsFactors = FALSE)
  } else {
    con <- DROMA.Set::connectDROMADatabase(config$db_path)
    on.exit(try(close(con), silent = TRUE), add = TRUE)

    project_anno <- DROMA.Set::listDROMAProjects()
    drug_anno <- DROMA.Set::getDROMAAnnotation("drug")
    sample_anno <- DROMA.Set::getDROMAAnnotation("sample")

    cell_names <- .filter_projects_for_drug_tumor(
      project_names = project_anno$project_name,
      project_anno = project_anno,
      drug_anno = drug_anno,
      sample_anno = sample_anno,
      dataset_types = c("CellLine", "PDC"),
      drug = config$drug,
      tumor_type = config$tumor_type,
      min_project_count = config$cell_n_datasets_t
    )
    pdcpdx_names <- .filter_projects_for_drug_tumor(
      project_names = project_anno$project_name,
      project_anno = project_anno,
      drug_anno = drug_anno,
      sample_anno = sample_anno,
      dataset_types = c("PDO", "PDX"),
      drug = config$drug,
      tumor_type = config$tumor_type,
      min_project_count = config$pdcpdx_n_datasets_t
    )

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

    ad_results <- batchFindTcgaADConcordantFeatures(
      selected_features = selected_genes,
      cell_set = cell_sets,
      pdcpdx_set = pdcpdx_sets,
      tumor_type = config$tumor_type,
      tcga_rna_counts_dir = config$tcga_rna_counts_dir,
      gene_probe_map_path = config$gene_probe_map_path,
      feature_type = config$feature2_type,
      data_type = config$data_type,
      p_t = config$tcga_ad_p_t
    )
    ad_stats <- ad_results$ad_stats
    ad_filtered <- ad_results$ad_filtered
  }

  data.table::fwrite(ad_stats, config$ad_stats_csv)
  data.table::fwrite(ad_filtered, config$ad_filtered_csv)
  saveRDS(ad_stats, config$ad_stats_rds)
  saveRDS(ad_filtered, config$ad_filtered_rds)

  list(
    ad_stats = ad_stats,
    ad_filtered = ad_filtered
  )
}

.run_clinical_validation <- function(config,
                                     preclinical_sig = NULL,
                                     cell_sig = NULL) {
  .load_biomarker_runtime(config$repo_root)
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)

  preclinical_sig <- preclinical_sig %||% .read_or_empty(config$ad_filtered_rds)
  cell_sig <- cell_sig %||% .read_or_empty(config$cell_sig_rds)

  if (!nrow(preclinical_sig)) {
    warning("AD-filtered preclinical genes are empty; clinical validation will return empty outputs", call. = FALSE)
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
        select_omics = unique(as.character(preclinical_sig$name)),
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
        preclinical_for_intersect <- as.data.frame(preclinical_sig)
        if (!"direction_pdcpdx" %in% names(preclinical_for_intersect)) {
          stop("AD-filtered preclinical data must contain direction_pdcpdx column", call. = FALSE)
        }

        final_biomarkers <- DROMA.R::getIntersectSignificantFeatures(
          pdcpdx = preclinical_for_intersect,
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
                                              tcga_ad_p_t = 0.05,
                                              tcga_rna_counts_dir = file.path(
                                                path.expand("~"),
                                                "Library/CloudStorage/OneDrive-Personal/28PHD_peng/250301-DROMA_project/archive260314/251112-DROMA_align/benchmark_mini/Input/TCGA/rna_counts"
                                              ),
                                              gene_probe_map_path = NULL,
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
    tcga_ad_p_t = tcga_ad_p_t,
    tcga_rna_counts_dir = tcga_rna_counts_dir,
    gene_probe_map_path = gene_probe_map_path,
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
  ad_results <- .run_tcga_ad(
    config,
    selected_genes = select_results$selected_genes
  )
  clinical_results <- .run_clinical_validation(
    config,
    preclinical_sig = ad_results$ad_filtered,
    cell_sig = select_results$cell_sig
  )

  c(batch_results, select_results, ad_results, clinical_results, list(config = config))
}
