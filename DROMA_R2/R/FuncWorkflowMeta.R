# Workflow Meta Functions ----

#' Convert workflow screening output into stage-3 meta results
#'
#' @description The workflow now reuses
#'   \code{DROMA.R::batchFindSignificantFeatures()} for grouped Spearman
#'   screening, so stage 3 mainly standardizes and annotates the resulting
#'   meta-analysis table for downstream workflow stages.
#' @param within_study_dt Data frame returned by \code{runWithinStudyScreen()}.
#' @param model_group Workflow model group label, typically \code{"invitro"}
#'   or \code{"invivo"}.
#' @param fdr_t FDR threshold retained for interface compatibility.
#' @param es_t Effect-size threshold retained for interface compatibility.
#' @param min_studies Minimum dataset count required to keep a feature.
#' @return A standardized meta-analysis result table.
#' @export
runGroupedMetaAnalysis <- function(within_study_dt,
                                   model_group,
                                   fdr_t = 0.1,
                                   es_t = 0.1,
                                   min_studies = 2L) {
  if (!is.data.frame(within_study_dt) || nrow(within_study_dt) == 0) {
    return(.empty_meta())
  }

  meta_dt <- data.table::as.data.table(within_study_dt)
  keep_idx <- is.finite(meta_dt[["effect_size"]]) &
    !is.na(meta_dt[["effect_size"]]) &
    is.finite(meta_dt[["p_value"]]) &
    !is.na(meta_dt[["p_value"]]) &
    !is.na(meta_dt[["n_datasets"]]) &
    meta_dt[["n_datasets"]] >= min_studies
  meta_dt <- meta_dt[keep_idx, ]
  if (nrow(meta_dt) == 0) {
    return(.empty_meta())
  }

  meta_dt <- unique(meta_dt[, c(
    "drug", "tumor_type", "gene", "p_value",
    "q_value", "effect_size", "n_datasets"
  ), with = FALSE])
  data.table::setnames(meta_dt, "gene", "name")

  if (!"q_value" %in% colnames(meta_dt) || all(is.na(meta_dt$q_value))) {
    meta_dt[, q_value := stats::p.adjust(p_value, method = "BH")]
  }

  meta_dt <- meta_dt[abs(meta_dt[["effect_size"]]) >= es_t, ]
  if (nrow(meta_dt) == 0) {
    return(.empty_meta())
  }

  meta_dt[["model_group"]] <- model_group
  meta_dt[["heterogeneity_p"]] <- NA_real_
  meta_dt[["i2"]] <- NA_real_

  data.table::setcolorder(
    meta_dt,
    c(
      "drug", "tumor_type", "name", "p_value", "q_value",
      "effect_size", "n_datasets", "model_group",
      "heterogeneity_p", "i2"
    )
  )
  meta_dt[]
}

#' Merge preclinical candidate tables across workflow groups
#'
#' @description Filters significant meta-analysis hits for the in vitro and in
#'   vivo groups, merges them by feature name, and records directional
#'   concordance across evidence sources.
#' @param invitro_meta In vitro meta-analysis result table.
#' @param invivo_meta In vivo meta-analysis result table.
#' @param fdr_t FDR threshold.
#' @param es_t Effect-size threshold.
#' @return A merged candidate table.
#' @export
mergePreclinicalCandidates <- function(invitro_meta, invivo_meta, fdr_t = 0.1, es_t = 0.1) {
  to_sig <- function(dt, suffix) {
    dt <- data.table::as.data.table(dt)
    keep_idx <- dt[["q_value"]] < fdr_t & abs(dt[["effect_size"]]) >= es_t
    keep_idx[is.na(keep_idx)] <- FALSE
    dt <- dt[keep_idx, ]
    if (nrow(dt) == 0) {
      return(data.table::data.table(name = character()))
    }
    out <- dt[, c("name", "effect_size", "q_value", "n_datasets"), with = FALSE]
    data.table::setnames(
      out,
      c("effect_size", "q_value", "n_datasets"),
      paste0(c("effect_size", "q_value", "n_datasets"), "_", suffix)
    )
    direction_col <- paste0("direction_", suffix)
    effect_col <- paste0("effect_size_", suffix)
    out[[direction_col]] <- ifelse(out[[effect_col]] >= 0, "Up", "Down")
    out
  }

  cell_sig <- to_sig(invitro_meta, "cell")
  invivo_sig <- to_sig(invivo_meta, "invivo")
  merged <- merge(cell_sig, invivo_sig, by = "name", all = TRUE)
  if (nrow(merged) == 0) {
    return(merged)
  }

  needed_cols <- c(
    "effect_size_cell", "q_value_cell", "n_datasets_cell", "direction_cell",
    "effect_size_invivo", "q_value_invivo", "n_datasets_invivo", "direction_invivo"
  )
  for (col in needed_cols) {
    if (!col %in% colnames(merged)) {
      merged[[col]] <- NA
    }
  }

  merged[["invitro_supported"]] <- !is.na(merged[["effect_size_cell"]])
  merged[["invivo_supported"]] <- !is.na(merged[["effect_size_invivo"]])
  merged[["direction_concordant"]] <- merged[["invitro_supported"]] &
    merged[["invivo_supported"]] &
    merged[["direction_cell"]] == merged[["direction_invivo"]]
  merged[["effect_size"]] <- rowMeans(
    cbind(merged[["effect_size_cell"]], merged[["effect_size_invivo"]]),
    na.rm = TRUE
  )
  merged[]
}
