computeCoverageForGroup <- function(project_names,
                                    drug_names,
                                    tumor_types,
                                    db_path,
                                    con,
                                    feature_type = "mRNA",
                                    dataset_label = "all",
                                    min_overlap_per_study = 20L,
                                    min_studies = 3L) {
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

  get_feature_samples <- function(project_name, tumor_type) {
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

  get_drug_samples <- function(project_name, drug_name, tumor_type) {
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

  rows <- vector("list", length(project_names) * length(drug_names) * length(tumor_types))
  idx <- 1L
  for (drug_name in drug_names) {
    for (tumor_type in tumor_types) {
      for (project_name in project_names) {
        expr_samples <- get_feature_samples(project_name, tumor_type)
        drug_samples <- get_drug_samples(project_name, drug_name, tumor_type)
        overlap_n <- length(intersect(expr_samples, drug_samples))
        rows[[idx]] <- data.table::data.table(
          project_name = project_name,
          drug = drug_name,
          tumor_type = tumor_type,
          model_group = dataset_label,
          feature_type = feature_type,
          expr_sample_n = length(expr_samples),
          drug_sample_n = length(drug_samples),
          overlap_n = overlap_n,
          eligible = overlap_n >= min_overlap_per_study
        )
        idx <- idx + 1L
      }
    }
  }

  coverage_dt <- data.table::rbindlist(rows, fill = TRUE)
  coverage_dt <- coverage_dt[!is.na(project_name)]

  candidates <- coverage_dt[
    eligible == TRUE,
    .(eligible_project_count = .N, eligible_projects = list(sort(project_name))),
    by = .(drug, tumor_type)
  ][eligible_project_count >= min_studies]

  candidate_save <- data.table::copy(candidates)
  if (nrow(candidate_save) > 0) {
    candidate_save[, eligible_projects := vapply(
      eligible_projects,
      function(x) paste(x, collapse = ";"),
      FUN.VALUE = character(1)
    )]
    candidate_save[, `:=`(
      model_group = dataset_label,
      feature_type = feature_type,
      min_overlap_per_study = min_overlap_per_study,
      min_studies = min_studies
    )]
  }

  list(
    coverage = coverage_dt,
    candidates_runtime = candidates,
    candidates_save = candidate_save
  )
}
