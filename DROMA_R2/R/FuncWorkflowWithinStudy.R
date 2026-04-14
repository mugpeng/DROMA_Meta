# Workflow Within-Study Screening Functions ----

.formatWorkflowWithinStudyResults <- function(meta_dt,
                                              project_scope,
                                              drug_name,
                                              tumor_type,
                                              eligible_projects) {
  meta_dt <- data.table::as.data.table(meta_dt)
  if (nrow(meta_dt) == 0) {
    return(.empty_within_study())
  }

  meta_dt <- data.table::copy(meta_dt)
  meta_dt[, `:=`(
    project_scope = project_scope,
    drug = drug_name,
    tumor_type = tumor_type,
    gene = name,
    n_total = NA_integer_,
    eligible_project_count = length(eligible_projects),
    eligible_projects = paste(sort(unique(eligible_projects)), collapse = ";")
  )]

  needed_cols <- c(
    "project_scope", "drug", "tumor_type", "gene", "n_total",
    "p_value", "effect_size", "q_value", "n_datasets",
    "eligible_project_count", "eligible_projects"
  )
  for (col in needed_cols) {
    if (!col %in% colnames(meta_dt)) {
      meta_dt[, (col) := NA]
    }
  }

  meta_dt[, ..needed_cols]
}

#' Run grouped biomarker screening for a workflow candidate
#'
#' @description Reuses \code{DROMA.R::batchFindSignificantFeatures()} to screen
#'   expression features against a drug across a subset of workflow projects.
#'   In contrast to the original within-project implementation, the expression
#'   side is allowed to come from any project inside the same workflow group as
#'   long as sample/model identifiers match, mirroring the \code{DROMA.R}
#'   batch-analysis logic.
#' @param group_set A \code{MultiDromaSet} object covering one workflow group.
#' @param drug_name Character scalar giving the drug name.
#' @param tumor_type Character scalar giving the tumor type.
#' @param feature_type Molecular feature type. Defaults to \code{"mRNA"}.
#' @param feature_names Optional feature names passed to
#'   \code{batchFindSignificantFeatures()}.
#' @param projects Optional character vector of eligible projects. When
#'   provided, the group is subset before analysis.
#' @param cores Number of cores for batch screening.
#' @param preloaded Whether to preload feature matrices before the batch run.
#' @param verbose Logical, whether to show verbose output from the underlying
#'   batch screening function.
#' @param batch_fun Function used to run the batch screen. Defaults to
#'   \code{DROMA.R::batchFindSignificantFeatures()} and is exposed mainly for
#'   testing.
#' @return A data table with workflow-oriented columns derived from the batch
#'   meta-analysis output.
#' @export
runWithinStudyScreen <- function(group_set,
                                 drug_name,
                                 tumor_type,
                                 feature_type = "mRNA",
                                 feature_names = NULL,
                                 projects = NULL,
                                 cores = 1L,
                                 preloaded = TRUE,
                                 verbose = FALSE,
                                 batch_fun = DROMA.R::batchFindSignificantFeatures) {
  if (!inherits(group_set, "MultiDromaSet")) {
    stop("group_set must be a MultiDromaSet")
  }

  selected_projects <- if (is.null(projects)) {
    names(group_set@DromaSets)
  } else {
    intersect(unique(as.character(projects)), names(group_set@DromaSets))
  }
  if (length(selected_projects) == 0) {
    return(.empty_within_study())
  }

  working_set <- if (length(selected_projects) == length(group_set@DromaSets)) {
    group_set
  } else {
    subsetWorkflowGroupSet(group_set, selected_projects)
  }

  meta_dt <- tryCatch(
    batch_fun(
      dromaset_object = working_set,
      feature1_type = "drug",
      feature1_name = drug_name,
      feature2_type = feature_type,
      feature2_name = feature_names,
      data_type = "all",
      tumor_type = tumor_type,
      overlap_only = FALSE,
      cores = as.integer(cores),
      show_progress = verbose,
      preloaded = preloaded,
      verbose = verbose
    ),
    error = function(e) NULL
  )
  if (is.null(meta_dt) || !is.data.frame(meta_dt) || nrow(meta_dt) == 0) {
    return(.empty_within_study())
  }

  .formatWorkflowWithinStudyResults(
    meta_dt = meta_dt,
    project_scope = paste(selected_projects, collapse = "|"),
    drug_name = drug_name,
    tumor_type = tumor_type,
    eligible_projects = selected_projects
  )
}
