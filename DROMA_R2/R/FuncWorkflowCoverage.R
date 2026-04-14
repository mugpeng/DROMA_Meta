# Workflow Coverage Functions ----

#' Compute Stage-1 coverage for a workflow group
#'
#' @description Evaluates whether each drug-bearing project in a workflow group
#'   has enough matched models to support downstream biomarker screening.
#'   Coverage follows the same cross-project matching intent as
#'   \code{DROMA.R::batchFindSignificantFeatures()}: drug response is taken from
#'   the focal project, while expression/model availability is traversed across
#'   all projects in the same workflow group and matched by shared sample/model
#'   identifiers.
#' @param project_names Character vector of workflow projects in one model group.
#' @param drug_names Character vector of drugs to evaluate.
#' @param tumor_types Character vector of tumor types to evaluate.
#' @param db_path Path to the DROMA SQLite database.
#' @param con Open database connection.
#' @param feature_type Molecular feature type used for screening.
#' @param dataset_label Workflow group label.
#' @param min_overlap_per_study Minimum matched sample count required per focal
#'   drug-response project.
#' @param min_studies Minimum number of focal projects required for a
#'   drug/tumor combination to advance.
#' @param feature_sample_resolver Optional function for resolving expression
#'   samples. Exposed mainly for testing.
#' @param drug_sample_resolver Optional function for resolving drug-response
#'   samples. Exposed mainly for testing.
#' @return A list containing the detailed coverage table plus runtime and export
#'   candidate tables.
#' @export
computeCoverageForGroup <- function(project_names,
                                    drug_names,
                                    tumor_types,
                                    db_path,
                                    con,
                                    feature_type = "mRNA",
                                    dataset_label = "all",
                                    min_overlap_per_study = 20L,
                                    min_studies = 3L,
                                    feature_sample_resolver = NULL,
                                    drug_sample_resolver = NULL) {
  dromaset_cache <- new.env(parent = emptyenv())
  feature_sample_cache <- new.env(parent = emptyenv())
  drug_sample_cache <- new.env(parent = emptyenv())

  get_project_set <- function(project_name) {
    if (!exists(project_name, envir = dromaset_cache, inherits = FALSE)) {
      assign(
        project_name,
        createDromaSetFromDatabase(
          projects = project_name,
          db_path = db_path,
          dataset_type = "all",
          con = con
        ),
        envir = dromaset_cache
      )
    }
    get(project_name, envir = dromaset_cache, inherits = FALSE)
  }

  if (is.null(feature_sample_resolver)) {
    feature_sample_resolver <- function(project_name, tumor_type) {
      cache_key <- paste(project_name, tumor_type, sep = "::")
      if (!exists(cache_key, envir = feature_sample_cache, inherits = FALSE)) {
        samples <- tryCatch(
          listDROMASamples(
            projects = project_name,
            feature_type = feature_type,
            tumor_type = tumor_type,
            connection = con
          ),
          error = function(e) character(0)
        )
        assign(cache_key, unique(samples), envir = feature_sample_cache)
      }
      get(cache_key, envir = feature_sample_cache, inherits = FALSE)
    }
  }

  if (is.null(drug_sample_resolver)) {
    drug_sample_resolver <- function(project_name, drug_name, tumor_type) {
      cache_key <- paste(project_name, drug_name, tumor_type, sep = "::")
      if (!exists(cache_key, envir = drug_sample_cache, inherits = FALSE)) {
        project_set <- get_project_set(project_name)
        drug_mat <- tryCatch(
          loadTreatmentResponse(
            object = project_set,
            select_drugs = drug_name,
            return_data = TRUE,
            data_type = "all",
            tumor_type = tumor_type,
            zscore = FALSE
          ),
          error = function(e) NULL
        )
        if (!is.matrix(drug_mat) || nrow(drug_mat) == 0 || ncol(drug_mat) == 0) {
          samples <- character(0)
        } else {
          drug_values <- as.numeric(drug_mat[1, , drop = TRUE])
          samples <- unique(colnames(drug_mat)[!is.na(drug_values)])
        }
        assign(cache_key, samples, envir = drug_sample_cache)
      }
      get(cache_key, envir = drug_sample_cache, inherits = FALSE)
    }
  }

  get_group_feature_samples <- function(tumor_type) {
    project_samples <- lapply(project_names, feature_sample_resolver, tumor_type = tumor_type)
    names(project_samples) <- project_names
    project_samples <- lapply(project_samples, unique)
    non_empty <- project_samples[vapply(project_samples, length, integer(1)) > 0]
    pooled_samples <- if (length(non_empty) > 0) {
      sort(unique(unlist(non_empty, use.names = FALSE)))
    } else {
      character(0)
    }
    list(
      pooled_samples = pooled_samples,
      project_samples = project_samples
    )
  }

  rows <- vector("list", length(project_names) * length(drug_names) * length(tumor_types))
  idx <- 1L
  for (drug_name in drug_names) {
    for (tumor_type in tumor_types) {
      feature_lookup <- get_group_feature_samples(tumor_type)
      expr_samples <- feature_lookup$pooled_samples

      for (project_name in project_names) {
        drug_samples <- unique(drug_sample_resolver(project_name, drug_name, tumor_type))
        overlap_samples <- intersect(expr_samples, drug_samples)
        matched_expr_projects <- names(feature_lookup$project_samples)[vapply(
          feature_lookup$project_samples,
          function(samples) length(intersect(samples, drug_samples)),
          integer(1)
        ) > 0]

        rows[[idx]] <- data.table::data.table(
          project_name = project_name,
          drug = drug_name,
          tumor_type = tumor_type,
          model_group = dataset_label,
          feature_type = feature_type,
          expr_sample_n = length(expr_samples),
          expr_project_count = length(matched_expr_projects),
          expr_projects = paste(sort(unique(matched_expr_projects)), collapse = ";"),
          drug_sample_n = length(drug_samples),
          overlap_n = length(overlap_samples),
          eligible = length(overlap_samples) >= min_overlap_per_study
        )
        idx <- idx + 1L
      }
    }
  }

  coverage_dt <- data.table::rbindlist(rows, fill = TRUE)
  coverage_dt <- coverage_dt[!is.na(coverage_dt[["project_name"]]), ]

  eligible_dt <- coverage_dt[coverage_dt[["eligible"]] == TRUE, ]
  if (nrow(eligible_dt) > 0) {
    split_keys <- interaction(eligible_dt[["drug"]], eligible_dt[["tumor_type"]], drop = TRUE)
    candidate_list <- lapply(split(eligible_dt, split_keys), function(df) {
      data.table::data.table(
        drug = df[["drug"]][1],
        tumor_type = df[["tumor_type"]][1],
        eligible_project_count = length(unique(df[["project_name"]])),
        eligible_projects = list(sort(unique(df[["project_name"]])))
      )
    })
    candidates <- data.table::rbindlist(candidate_list, fill = TRUE)
    candidates <- candidates[candidates[["eligible_project_count"]] >= min_studies, ]
  } else {
    candidates <- data.table::data.table(
      drug = character(),
      tumor_type = character(),
      eligible_project_count = integer(),
      eligible_projects = list()
    )
  }

  candidate_save <- data.table::copy(candidates)
  if (nrow(candidate_save) > 0) {
    candidate_save[["eligible_projects"]] <- vapply(
      candidate_save[["eligible_projects"]],
      function(x) paste(x, collapse = ";"),
      FUN.VALUE = character(1)
    )
    candidate_save[["model_group"]] <- dataset_label
    candidate_save[["feature_type"]] <- feature_type
    candidate_save[["min_overlap_per_study"]] <- min_overlap_per_study
    candidate_save[["min_studies"]] <- min_studies
  }

  list(
    coverage = coverage_dt,
    candidates_runtime = candidates,
    candidates_save = candidate_save
  )
}
