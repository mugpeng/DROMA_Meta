# Preclinical merge helpers ----

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
    out <- unique(dt[, c("drug", "tumor_type", "name", "effect_size", "q_value", "n_datasets", "direction"), with = FALSE])
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

  merged$invitro_supported <- TRUE
  merged$pdcpdx_supported <- TRUE
  merged$direction_concordant <- merged$direction_cellline == merged$direction_pdcpdx
  merged$preclinical_direction <- ifelse(
    merged$effect_size_cellline + merged$effect_size_pdcpdx >= 0,
    "Up",
    "Down"
  )
  merged$effect_size <- rowMeans(
    cbind(merged$effect_size_cellline, merged$effect_size_pdcpdx),
    na.rm = TRUE
  )

  if (!is.null(tcga_results) && nrow(tcga_results) > 0) {
    merged <- merge(merged, tcga_results, by = c("drug", "tumor_type", "name"), all.x = TRUE)
  }

  merged[]
}
