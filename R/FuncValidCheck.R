# Validity checks for drug and tumor-type grids ----

#' Get Valid Annotation Values by Dataset Group
#'
#' @description Filters annotation values by dataset type and keeps only entries
#' that appear in at least `min_project_count` projects.
#' @param project_anno Project annotation table.
#' @param anno_df Annotation table containing `ProjectID` and `value_col`.
#' @param dataset_types Dataset types to retain.
#' @param value_col Column name holding the candidate values.
#' @param min_project_count Minimum number of projects required.
#' @param exclude_values Optional values to remove after counting.
#' @return A sorted character vector of retained values.
#' @export
getValidValues <- function(project_anno,
                           anno_df,
                           dataset_types,
                           value_col,
                           min_project_count,
                           exclude_values = NULL) {
  group_projects <- unique(as.character(
    project_anno[project_anno$dataset_type %in% dataset_types, "project_name"]
  ))
  group_projects <- group_projects[!is.na(group_projects) & nzchar(group_projects)]

  anno_filtered <- anno_df[
    anno_df$ProjectID %in% group_projects &
      !is.na(anno_df[[value_col]]) &
      nzchar(anno_df[[value_col]]),
  ]

  counts <- table(anno_filtered[[value_col]])
  keep_values <- names(counts)[counts >= min_project_count]

  if (!is.null(exclude_values)) {
    keep_values <- setdiff(keep_values, exclude_values)
  }

  sort(keep_values)
}

#' Get Valid Drugs and Tumor Types
#'
#' @description Computes the intersection of drug and tumor-type values that
#' satisfy minimum project-count thresholds in both cell/PDC and PDO/PDX data.
#' @param project_anno Project annotation table.
#' @param drug_anno Drug annotation table.
#' @param sample_anno Sample annotation table.
#' @param cell_n_datasets_t Minimum project count for CellLine/PDC.
#' @param pdcpdx_n_datasets_t Minimum project count for PDO/PDX.
#' @return A list with `valid_drugs` and `valid_tumor_types`.
#' @export
getValidDrugsAndTumorTypes <- function(project_anno,
                                       drug_anno,
                                       sample_anno,
                                       cell_n_datasets_t,
                                       pdcpdx_n_datasets_t) {
  valid_cell_drugs <- getValidValues(
    project_anno = project_anno,
    anno_df = drug_anno,
    dataset_types = c("CellLine", "PDC"),
    value_col = "DrugName",
    min_project_count = cell_n_datasets_t
  )
  valid_pdcpdx_drugs <- getValidValues(
    project_anno = project_anno,
    anno_df = drug_anno,
    dataset_types = c("PDO", "PDX"),
    value_col = "DrugName",
    min_project_count = pdcpdx_n_datasets_t
  )

  valid_cell_tumor_types <- getValidValues(
    project_anno = project_anno,
    anno_df = sample_anno,
    dataset_types = c("CellLine", "PDC"),
    value_col = "TumorType",
    min_project_count = cell_n_datasets_t,
    exclude_values = "non-cancer"
  )
  valid_pdcpdx_tumor_types <- getValidValues(
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

#' Filter Projects for a Drug and Tumor-Type Pair
#'
#' @description Identifies projects from a dataset group that contain both the
#' requested drug and tumor type, and enforces a minimum project-count threshold.
#' @param project_anno Project annotation table.
#' @param drug_anno Drug annotation table.
#' @param sample_anno Sample annotation table.
#' @param dataset_types Dataset types to search.
#' @param drug Drug name to retain.
#' @param tumor_type Tumor type to retain.
#' @param min_project_count Minimum number of projects required.
#' @return A character vector of eligible project IDs.
#' @export
filterProjectsForDrugTumor <- function(project_anno,
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
          "Need at least %d project(s), found %d."
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
