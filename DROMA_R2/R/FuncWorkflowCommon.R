buildWorkflowConfig <- function(repo_root,
                                db_path = file.path(repo_root, "Data", "droma.sqlite"),
                                ctrdb_path = file.path(repo_root, "Data", "ctrdb.sqlite"),
                                tcga_dir,
                                output_base,
                                drug_names,
                                tumor_types,
                                feature_type = "mRNA",
                                es_t = 0.1,
                                fdr_t = 0.1,
                                requested_cores = 3L) {
  list(
    repo_root = normalizePath(repo_root, mustWork = TRUE),
    db_path = normalizePath(db_path, mustWork = TRUE),
    ctrdb_path = ctrdb_path,
    tcga_dir = normalizePath(tcga_dir, mustWork = TRUE),
    output_base = output_base,
    drug_names = unique(as.character(drug_names)),
    tumor_types = unique(as.character(tumor_types)),
    feature_type = feature_type,
    es_t = es_t,
    fdr_t = fdr_t,
    requested_cores = as.integer(requested_cores),
    groups = list(
      invitro = list(
        dataset_types = c("CellLine", "PDC"),
        min_overlap_per_study = 20L,
        min_studies = 3L
      ),
      invivo = list(
        dataset_types = c("PDO", "PDX"),
        min_overlap_per_study = 10L,
        min_studies = 2L
      )
    )
  )
}

#' Subset a workflow MultiDromaSet to selected projects
#'
#' @description Creates a new \code{MultiDromaSet} containing only the selected
#'   workflow projects. This keeps the workflow implementation aligned with
#'   \code{DROMA.Set} objects and allows downstream code to reuse
#'   \code{DROMA.R::batchFindSignificantFeatures()} directly.
#' @param group_set A \code{MultiDromaSet} object.
#' @param project_names Character vector of project names to keep.
#' @return A \code{MultiDromaSet} restricted to the requested projects.
#' @export
subsetWorkflowGroupSet <- function(group_set, project_names) {
  if (!inherits(group_set, "MultiDromaSet")) {
    stop("group_set must be a MultiDromaSet")
  }

  project_names <- unique(as.character(project_names))
  keep_projects <- intersect(project_names, names(group_set@DromaSets))
  if (length(keep_projects) == 0) {
    stop("No matching projects found in group_set")
  }

  DROMA.Set::createMultiDromaSetFromObjects(
    group_set@DromaSets[keep_projects],
    project_names = keep_projects
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
      invitro = c("CellLine", "PDC"),
      invivo = c("PDO", "PDX")
    )
  }

  lapply(groups, function(dataset_types) {
    sort(unique(project_anno$project_name[project_anno$dataset_type %in% dataset_types]))
  })
}

createWorkflowGroupSet <- function(db_path, con, project_names) {
  createMultiDromaSetFromAllProjects(
    db_path = db_path,
    include_projects = sort(unique(project_names)),
    con = con
  )
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
    detected_cores <- 2L
  }
  min(as.integer(requested_cores), max(1L, detected_cores - 1L))
}

.empty_within_study <- function() {
  data.table::data.table(
    project_scope = character(),
    drug = character(),
    tumor_type = character(),
    gene = character(),
    n_total = integer(),
    p_value = numeric(),
    effect_size = numeric(),
    q_value = numeric(),
    n_datasets = integer(),
    eligible_project_count = integer(),
    eligible_projects = character()
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
    i2 = numeric()
  )
}
