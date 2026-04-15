buildWorkflowConfig <- function(repo_root,
                                db_path = file.path(repo_root, "..", "Data", "droma.sqlite"),
                                ctrdb_path = file.path(repo_root, "..", "Data", "ctrdb.sqlite"),
                                tcga_dir = file.path(repo_root, "..", "Data", "TCGA", "rna_counts"),
                                gene_probe_map = file.path(repo_root, "..", "Data", "gencode.v22.annotation.gene.probeMap"),
                                output_base,
                                drug_names = "Paclitaxel",
                                tumor_types = "breast cancer",
                                feature_type = "mRNA",
                                es_t = 0.1,
                                fdr_t = 0.1,
                                pair_cor_t = 0.2,
                                pair_p_t = 0.05,
                                tcga_fdr_t = 0.01,
                                requested_cores = 3L) {
  list(
    repo_root = normalizePath(repo_root, mustWork = TRUE),
    db_path = normalizePath(db_path, mustWork = TRUE),
    ctrdb_path = normalizePath(ctrdb_path, mustWork = TRUE),
    tcga_dir = normalizePath(tcga_dir, mustWork = TRUE),
    gene_probe_map = normalizePath(gene_probe_map, mustWork = FALSE),
    output_base = normalizePath(output_base, mustWork = FALSE),
    drug_names = unique(as.character(drug_names)),
    tumor_types = unique(as.character(tumor_types)),
    feature_type = feature_type,
    es_t = as.numeric(es_t),
    fdr_t = as.numeric(fdr_t),
    pair_cor_t = as.numeric(pair_cor_t),
    pair_p_t = as.numeric(pair_p_t),
    tcga_fdr_t = as.numeric(tcga_fdr_t),
    requested_cores = as.integer(requested_cores),
    groups = list(
      cellline = list(
        dataset_types = c("CellLine", "PDC"),
        min_overlap_per_pair = 20L,
        min_pairs = 3L
      ),
      pdcpdx = list(
        dataset_types = c("PDO", "PDX"),
        min_overlap_per_pair = 10L,
        min_pairs = 2L
      )
    )
  )
}

createWorkflowProjectGroups <- function(project_anno, groups = NULL) {
  if (!is.data.frame(project_anno)) {
    stop("project_anno must be a data.frame")
  }
  if (!all(c("project_name", "dataset_type") %in% colnames(project_anno))) {
    stop("project_anno must contain project_name and dataset_type")
  }

  if (is.null(groups)) {
    groups <- list(
      cellline = c("CellLine", "PDC"),
      pdcpdx = c("PDO", "PDX")
    )
  }

  lapply(groups, function(dataset_types) {
    sort(unique(project_anno$project_name[project_anno$dataset_type %in% dataset_types]))
  })
}

getSharedGroupFeatures <- function(group_sets, feature_type = "mRNA") {
  group_feature_sets <- lapply(group_sets, function(group_set) {
    project_features <- lapply(names(group_set@DromaSets), function(project_name) {
      tryCatch(
        DROMA.Set::listDROMAFeatures(project_name, feature_type),
        error = function(e) character(0)
      )
    })
    sort(unique(unlist(project_features, use.names = FALSE)))
  })

  shared <- Reduce(intersect, group_feature_sets)
  shared[nzchar(shared)]
}

.sanitize_name <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

.get_safe_cores <- function(requested_cores = 3L) {
  detected_cores <- if (requireNamespace("parallel", quietly = TRUE)) {
    parallel::detectCores()
  } else {
    1L
  }
  if (is.na(detected_cores) || detected_cores < 2L) {
    return(1L)
  }
  max(1L, min(as.integer(requested_cores), detected_cores - 1L))
}

.runParallelRows <- function(X, FUN, cores = 1L) {
  X <- as.list(X)
  if (length(X) == 0) {
    return(list())
  }
  lapply(X, FUN)
}

.empty_coverage_candidates <- function() {
  data.table::data.table(
    drug = character(),
    tumor_type = character(),
    eligible_pair_count = integer(),
    eligible_pairs = list(),
    eligible_drug_projects = list(),
    eligible_expr_projects = list(),
    filtered_projects = list()
  )
}

.empty_meta <- function() {
  data.table::data.table(
    drug = character(),
    tumor_type = character(),
    name = character(),
    p_value = numeric(),
    q_value = numeric(),
    effect_size = numeric(),
    n_datasets = integer(),
    model_group = character(),
    heterogeneity_p = numeric(),
    i2 = numeric(),
    direction = character()
  )
}
