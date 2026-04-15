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
