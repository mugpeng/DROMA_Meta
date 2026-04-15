# Meta-analysis and merge helpers ----

runGroupedMetaAnalysis <- function(group_set,
                                   group_name,
                                   drug_name,
                                   tumor_type,
                                   feature_names,
                                   feature_type = "mRNA",
                                   cores = 1L,
                                   preloaded = TRUE) {
  if (length(feature_names) == 0) {
    return(.empty_meta())
  }

  meta_dt <- batchFindSignificantFeatures(
    dromaset_object = group_set,
    feature1_type = "drug",
    feature1_name = drug_name,
    feature2_type = feature_type,
    feature2_name = feature_names,
    data_type = "all",
    tumor_type = tumor_type,
    overlap_only = FALSE,
    cores = as.integer(cores),
    show_progress = FALSE,
    preloaded = preloaded,
    verbose = FALSE
  )
  if (!is.data.frame(meta_dt) || nrow(meta_dt) == 0) {
    return(.empty_meta())
  }

  meta_dt <- data.table::as.data.table(meta_dt)
  data.table::setnames(meta_dt, "name", "name", skip_absent = TRUE)
  if (!"q_value" %in% colnames(meta_dt)) {
    meta_dt[, q_value := stats::p.adjust(p_value, method = "BH")]
  }
  meta_dt[, `:=`(
    drug = drug_name,
    tumor_type = tumor_type,
    model_group = group_name,
    heterogeneity_p = if ("pval.Q" %in% colnames(meta_dt)) pval.Q else NA_real_,
    i2 = if ("I2" %in% colnames(meta_dt)) I2 else NA_real_,
    direction = ifelse(effect_size >= 0, "Up", "Down")
  )]

  keep_cols <- c(
    "drug", "tumor_type", "name", "p_value", "q_value", "effect_size",
    "n_datasets", "model_group", "heterogeneity_p", "i2", "direction"
  )
  for (col in keep_cols) {
    if (!col %in% colnames(meta_dt)) {
      meta_dt[, (col) := NA]
    }
  }
  unique(meta_dt[, ..keep_cols])
}

mergePreclinicalCandidates <- function(cellline_meta,
                                       pdcpdx_meta,
                                       tcga_results = NULL,
                                       fdr_t = 0.1,
                                       es_t = 0.1) {
  to_sig <- function(dt, suffix) {
    dt <- data.table::as.data.table(dt)
    if (nrow(dt) == 0) {
      return(data.table::data.table(name = character(), drug = character(), tumor_type = character()))
    }
    keep <- dt$q_value < fdr_t & abs(dt$effect_size) >= es_t
    keep[is.na(keep)] <- FALSE
    dt <- dt[keep]
    if (nrow(dt) == 0) {
      return(data.table::data.table(name = character(), drug = character(), tumor_type = character()))
    }
    out <- unique(dt[, .(drug, tumor_type, name, effect_size, q_value, n_datasets, direction)])
    data.table::setnames(
      out,
      c("effect_size", "q_value", "n_datasets", "direction"),
      paste0(c("effect_size", "q_value", "n_datasets", "direction"), "_", suffix)
    )
    out
  }

  cell_sig <- to_sig(cellline_meta, "cellline")
  pdcpdx_sig <- to_sig(pdcpdx_meta, "pdcpdx")
  merged <- merge(cell_sig, pdcpdx_sig, by = c("drug", "tumor_type", "name"))
  if (nrow(merged) == 0) {
    return(merged)
  }

  merged[, `:=`(
    invitro_supported = TRUE,
    pdcpdx_supported = TRUE,
    direction_concordant = direction_cellline == direction_pdcpdx,
    preclinical_direction = ifelse(effect_size_cellline + effect_size_pdcpdx >= 0, "Up", "Down"),
    effect_size = rowMeans(cbind(effect_size_cellline, effect_size_pdcpdx), na.rm = TRUE)
  )]

  if (!is.null(tcga_results) && nrow(tcga_results) > 0) {
    merged <- merge(merged, tcga_results, by = c("drug", "tumor_type", "name"), all.x = TRUE)
  }

  merged[]
}
