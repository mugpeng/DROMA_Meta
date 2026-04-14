.extract_feature_vector <- function(project_set, feature_type, feature_name, tumor_type) {
  prof <- tryCatch(
    loadMolecularProfiles(
      object = project_set,
      feature_type = feature_type,
      select_features = feature_name,
      return_data = TRUE,
      data_type = "all",
      tumor_type = tumor_type,
      zscore = FALSE,
      format = "wide"
    ),
    error = function(e) NULL
  )
  if (is.null(prof) || !is.matrix(prof) || nrow(prof) == 0 || ncol(prof) == 0) {
    return(NULL)
  }
  values <- as.numeric(prof[1, , drop = TRUE])
  names(values) <- colnames(prof)
  values[!is.na(values)]
}

.extract_drug_vector <- function(project_set, drug_name, tumor_type) {
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
  if (is.null(drug_mat) || !is.matrix(drug_mat) || nrow(drug_mat) == 0 || ncol(drug_mat) == 0) {
    return(NULL)
  }
  values <- as.numeric(drug_mat[1, , drop = TRUE])
  names(values) <- colnames(drug_mat)
  values[!is.na(values)]
}

.run_single_gene_wilcox <- function(expr_values, drug_values) {
  common_samples <- intersect(names(expr_values), names(drug_values))
  if (length(common_samples) < 6) {
    return(NULL)
  }

  expr_common <- expr_values[common_samples]
  drug_common <- drug_values[common_samples]
  split_cut <- stats::median(expr_common, na.rm = TRUE)
  high_idx <- expr_common > split_cut
  low_idx <- expr_common <= split_cut
  if (sum(high_idx) < 3 || sum(low_idx) < 3) {
    return(NULL)
  }

  high_drug <- as.numeric(drug_common[high_idx])
  low_drug <- as.numeric(drug_common[low_idx])
  wilcox_re <- tryCatch(
    suppressWarnings(stats::wilcox.test(high_drug, low_drug)),
    error = function(e) NULL
  )
  cliff_re <- tryCatch(
    suppressWarnings(effsize::cliff.delta(high_drug, low_drug)),
    error = function(e) NULL
  )
  if (is.null(wilcox_re) || is.null(cliff_re)) {
    return(NULL)
  }

  data.table::data.table(
    n_total = length(common_samples),
    n_high = length(high_drug),
    n_low = length(low_drug),
    p_value = wilcox_re$p.value,
    effect_size = as.numeric(cliff_re$estimate)
  )
}

runWithinStudyScreen <- function(group_set,
                                 drug_name,
                                 tumor_type,
                                 feature_type = "mRNA",
                                 feature_names = NULL,
                                 projects = NULL) {
  if (is.null(projects)) {
    projects <- names(group_set@DromaSets)
  }
  rows <- list()

  for (project_name in projects) {
    project_set <- group_set@DromaSets[[project_name]]
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
    if (is.null(drug_mat) || !is.matrix(drug_mat) || nrow(drug_mat) == 0 || ncol(drug_mat) == 0) {
      next
    }
    drug_values <- as.numeric(drug_mat[1, , drop = TRUE])
    names(drug_values) <- colnames(drug_mat)
    drug_values <- drug_values[!is.na(drug_values)]
    if (length(drug_values) < 6) {
      next
    }

    expr_mat <- tryCatch(
      loadMolecularProfiles(
        object = project_set,
        feature_type = feature_type,
        select_features = feature_names,
        return_data = TRUE,
        data_type = "all",
        tumor_type = tumor_type,
        zscore = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(expr_mat) || !is.matrix(expr_mat) || nrow(expr_mat) == 0 || ncol(expr_mat) == 0) {
      next
    }

    common_samples <- intersect(colnames(expr_mat), names(drug_values))
    if (length(common_samples) < 6) {
      next
    }
    expr_mat <- expr_mat[, common_samples, drop = FALSE]
    drug_common <- as.numeric(drug_values[common_samples])
    names(drug_common) <- common_samples

    study_rows <- lapply(seq_len(nrow(expr_mat)), function(idx) {
      expr_values <- as.numeric(expr_mat[idx, , drop = TRUE])
      if (all(is.na(expr_values))) {
        return(NULL)
      }
      names(expr_values) <- colnames(expr_mat)
      expr_values <- expr_values[!is.na(expr_values)]
      stat_row <- .run_single_gene_wilcox(expr_values, drug_common)
      if (is.null(stat_row)) {
        return(NULL)
      }
      stat_row[, `:=`(
        project_name = project_name,
        drug = drug_name,
        tumor_type = tumor_type,
        gene = rownames(expr_mat)[idx]
      )]
      stat_row
    })

    study_dt <- data.table::rbindlist(study_rows, fill = TRUE)
    if (nrow(study_dt) == 0) {
      next
    }
    study_dt[, q_value := stats::p.adjust(p_value, method = "BH")]
    data.table::setcolorder(
      study_dt,
      c("project_name", "drug", "tumor_type", "gene", "n_total", "n_high", "n_low", "p_value", "effect_size", "q_value")
    )
    rows[[project_name]] <- study_dt
  }

  if (length(rows) == 0) {
    return(.empty_within_study())
  }
  data.table::rbindlist(rows, fill = TRUE)
}
