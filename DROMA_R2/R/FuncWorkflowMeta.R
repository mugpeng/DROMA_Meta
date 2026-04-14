runGroupedMetaAnalysis <- function(within_study_dt,
                                   model_group,
                                   fdr_t = 0.1,
                                   es_t = 0.1,
                                   min_studies = 2L) {
  if (!is.data.frame(within_study_dt) || nrow(within_study_dt) == 0) {
    return(.empty_meta())
  }

  sig_dt <- data.table::as.data.table(within_study_dt)
  sig_dt <- sig_dt[is.finite(effect_size) & !is.na(effect_size) & is.finite(p_value) & !is.na(p_value)]
  if (nrow(sig_dt) == 0) {
    return(.empty_meta())
  }

  grouped <- split(sig_dt, interaction(sig_dt$drug, sig_dt$tumor_type, sig_dt$gene, drop = TRUE))
  meta_rows <- lapply(grouped, function(df) {
    if (nrow(df) < min_studies) {
      return(NULL)
    }
    df <- data.table::copy(df)
    df[, se := sqrt((1 - effect_size^2) * (n_high + n_low + 1) / (12 * n_high * n_low))]
    df <- df[is.finite(se) & !is.na(se) & se > 0]
    if (nrow(df) < min_studies) {
      return(NULL)
    }
    w_fixed <- 1 / (df$se ^ 2)
    effect_fixed <- sum(w_fixed * df$effect_size) / sum(w_fixed)
    se_fixed <- sqrt(1 / sum(w_fixed))
    q_stat <- sum(w_fixed * (df$effect_size - effect_fixed) ^ 2)
    df_q <- nrow(df) - 1
    heterogeneity_p <- stats::pchisq(q_stat, df = df_q, lower.tail = FALSE)
    c_val <- sum(w_fixed) - (sum(w_fixed ^ 2) / sum(w_fixed))
    tau2 <- max(0, (q_stat - df_q) / c_val)
    w_random <- 1 / (df$se ^ 2 + tau2)
    effect_random <- sum(w_random * df$effect_size) / sum(w_random)
    se_random <- sqrt(1 / sum(w_random))
    use_random <- is.finite(heterogeneity_p) && !is.na(heterogeneity_p) && heterogeneity_p < 0.1
    effect <- if (use_random) effect_random else effect_fixed
    se_meta <- if (use_random) se_random else se_fixed
    z_val <- effect / se_meta
    p_value <- 2 * stats::pnorm(abs(z_val), lower.tail = FALSE)
    data.table::data.table(
      drug = df$drug[[1]],
      tumor_type = df$tumor_type[[1]],
      name = df$gene[[1]],
      p_value = ifelse(is.na(p_value), 1, p_value),
      effect_size = ifelse(is.na(effect), 0, effect),
      n_datasets = nrow(df),
      model_group = model_group,
      heterogeneity_p = heterogeneity_p,
      i2 = if (q_stat > 0) max(0, (q_stat - df_q) / q_stat) * 100 else 0
    )
  })

  meta_dt <- data.table::rbindlist(meta_rows, fill = TRUE)
  if (nrow(meta_dt) == 0) {
    return(.empty_meta())
  }
  meta_dt[, q_value := stats::p.adjust(p_value, method = "BH")]
  meta_dt <- meta_dt[abs(effect_size) >= es_t]
  data.table::setcolorder(meta_dt, c("drug", "tumor_type", "name", "p_value", "q_value", "effect_size", "n_datasets", "model_group", "heterogeneity_p", "i2"))
  meta_dt[]
}

mergePreclinicalCandidates <- function(invitro_meta, invivo_meta, fdr_t = 0.1, es_t = 0.1) {
  to_sig <- function(dt, suffix) {
    dt <- data.table::as.data.table(dt)
    dt <- dt[q_value < fdr_t & abs(effect_size) >= es_t]
    if (nrow(dt) == 0) {
      return(data.table::data.table(name = character()))
    }
    out <- dt[, .(name, effect_size, q_value, n_datasets)]
    data.table::setnames(out, c("effect_size", "q_value", "n_datasets"),
                         paste0(c("effect_size", "q_value", "n_datasets"), "_", suffix))
    out[, paste0("direction_", suffix) := ifelse(get(paste0("effect_size_", suffix)) >= 0, "Up", "Down")]
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
      merged[, (col) := NA]
    }
  }
  merged[, invitro_supported := !is.na(effect_size_cell)]
  merged[, invivo_supported := !is.na(effect_size_invivo)]
  merged[, direction_concordant := invitro_supported & invivo_supported & direction_cell == direction_invivo]
  merged[, effect_size := rowMeans(cbind(effect_size_cell, effect_size_invivo), na.rm = TRUE)]
  merged[]
}
