# Workflow Coverage Functions ----

computeCoverageForGroup <- function(project_names,
                                    drug_names,
                                    tumor_types,
                                    db_path,
                                    con,
                                    feature_type = "mRNA",
                                    dataset_label = "all",
                                    min_overlap_per_pair = 20L,
                                    min_pairs = 3L,
                                    feature_sample_resolver = NULL,
                                    drug_sample_resolver = NULL) {
  dromaset_cache <- new.env(parent = emptyenv())
  feature_sample_cache <- new.env(parent = emptyenv())
  drug_sample_cache <- new.env(parent = emptyenv())

  get_project_set <- function(project_name) {
    if (!exists(project_name, envir = dromaset_cache, inherits = FALSE)) {
      assign(
        project_name,
        DROMA.Set::createDromaSetFromDatabase(
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
          DROMA.Set::listDROMASamples(
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

  rows <- list()
  idx <- 1L
  for (drug_name in drug_names) {
    for (tumor_type in tumor_types) {
      for (drug_project in project_names) {
        drug_samples <- unique(drug_sample_resolver(drug_project, drug_name, tumor_type))
        for (expr_project in project_names) {
          expr_samples <- unique(feature_sample_resolver(expr_project, tumor_type))
          overlap_samples <- intersect(drug_samples, expr_samples)
          rows[[idx]] <- data.table::data.table(
            drug = drug_name,
            tumor_type = tumor_type,
            model_group = dataset_label,
            drug_project = drug_project,
            expr_project = expr_project,
            drug_sample_n = length(drug_samples),
            expr_sample_n = length(expr_samples),
            overlap_n = length(overlap_samples),
            eligible = length(overlap_samples) >= min_overlap_per_pair
          )
          idx <- idx + 1L
        }
      }
    }
  }

  coverage_dt <- data.table::rbindlist(rows, fill = TRUE)
  eligible_dt <- coverage_dt[eligible == TRUE]
  if (nrow(eligible_dt) == 0) {
    return(list(
      coverage = coverage_dt,
      candidates_runtime = .empty_coverage_candidates(),
      candidates_save = .empty_coverage_candidates()
    ))
  }

  split_keys <- interaction(eligible_dt$drug, eligible_dt$tumor_type, drop = TRUE)
  candidate_list <- lapply(split(eligible_dt, split_keys), function(df) {
    pair_dt <- unique(df[, c("drug_project", "expr_project"), with = FALSE])
    data.table::data.table(
      drug = df$drug[1],
      tumor_type = df$tumor_type[1],
      eligible_pair_count = nrow(pair_dt),
      eligible_pairs = list(pair_dt),
      eligible_drug_projects = list(sort(unique(pair_dt$drug_project))),
      eligible_expr_projects = list(sort(unique(pair_dt$expr_project)))
    )
  })
  candidates <- data.table::rbindlist(candidate_list, fill = TRUE)
  candidates <- candidates[eligible_pair_count >= min_pairs]

  candidates_save <- data.table::copy(candidates)
  if (nrow(candidates_save) > 0) {
    candidates_save[, eligible_pairs := vapply(
      eligible_pairs,
      function(x) paste(paste(x$drug_project, x$expr_project, sep = "->"), collapse = ";"),
      FUN.VALUE = character(1)
    )]
    candidates_save[, eligible_drug_projects := vapply(
      eligible_drug_projects,
      function(x) paste(x, collapse = ";"),
      FUN.VALUE = character(1)
    )]
    candidates_save[, eligible_expr_projects := vapply(
      eligible_expr_projects,
      function(x) paste(x, collapse = ";"),
      FUN.VALUE = character(1)
    )]
    candidates_save[, model_group := dataset_label]
    candidates_save[, feature_type := feature_type]
    candidates_save[, min_overlap_per_pair := min_overlap_per_pair]
    candidates_save[, min_pairs := min_pairs]
  }

  list(
    coverage = coverage_dt,
    candidates_runtime = candidates,
    candidates_save = candidates_save
  )
}
