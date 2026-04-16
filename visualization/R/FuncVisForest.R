# Forest plot visualization ----

#' Plot Stage-wise Forest Plot for Top Biomarkers
#'
#' @description Creates a horizontal forest-style plot showing the effect size
#' at each validation stage (cell line, PDC/PDX, clinical) for the top final
#' biomarkers. Each gene is a row, with points and horizontal lines representing
#' the effect size and a nominal confidence interval.
#' @param final_biomarkers A \code{data.table} from \code{collectAllFinalBiomarkers()}.
#' @param top_n Number of top biomarkers to show, ranked by clinical effect size.
#' @param es_cols Named character vector mapping stage labels to effect-size columns.
#'   Defaults to \code{c("Cell" = "effect_size_cell_invitro", "PDC/PDX" = "effect_size_pdcpdx_invitro",
#'   "Clinical" = "effect_size_ctrdb")}.
#' @param title Plot title.
#' @return A ggplot2 object.
#' @export
plotStageForest <- function(final_biomarkers,
                            top_n = 10,
                            es_cols = c(
                              "Cell"     = "effect_size_cell_invitro",
                              "PDC/PDX" = "effect_size_pdcpdx_invitro",
                              "Clinical" = "effect_size_ctrdb"
                            ),
                            title = "Top biomarkers across validation stages") {
  final_biomarkers <- data.table::as.data.table(final_biomarkers)
  if (!all(es_cols %in% names(final_biomarkers))) {
    stop("final_biomarkers must contain all effect-size columns specified in es_cols", call. = FALSE)
  }
  if (!"name" %in% names(final_biomarkers)) {
    stop("final_biomarkers must contain 'name' column", call. = FALSE)
  }

  # Rank by absolute clinical effect size
  clinical_col <- es_cols["Clinical"]
  final_biomarkers <- final_biomarkers[order(-abs(get(clinical_col)))]
  top_genes <- final_biomarkers$name[seq_len(min(top_n, nrow(final_biomarkers)))]
  sub <- final_biomarkers[name %in% top_genes]

  # Melt to long format
  melt_list <- lapply(names(es_cols), function(stage) {
    col <- es_cols[stage]
    dt <- data.table::data.table(
      name  = sub$name,
      stage = stage,
      es    = sub[[col]]
    )
    dt
  })
  long_dt <- data.table::rbindlist(melt_list)
  long_dt <- long_dt[!is.na(es)]

  # Order genes by clinical effect size
  gene_order <- sub$name
  long_dt[, name := factor(name, levels = rev(gene_order))]
  long_dt[, stage := factor(stage, levels = c("Cell", "PDC/PDX", "Clinical"))]

  ev_colors <- getMetaVisColors("evidence")
  stage_colors <- ev_colors[c("Cell line", "PDC/PDX", "Clinical")]
  names(stage_colors) <- c("Cell", "PDC/PDX", "Clinical")

  ggplot2::ggplot(long_dt, ggplot2::aes(x = es, y = name, color = stage)) +
    # Zero reference
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey50", linewidth = 0.3) +
    # Connect stages per gene with light line
    ggplot2::geom_path(
      ggplot2::aes(group = name),
      color = "grey70", linewidth = 0.4,
      position = ggplot2::position_dodge(width = 0)
    ) +
    # Points
    ggplot2::geom_point(size = 3, position = ggplot2::position_dodge(width = 0)) +
    ggplot2::scale_color_manual(values = stage_colors, name = "Stage") +
    ggplot2::labs(
      title = title,
      x = "Effect size",
      y = NULL
    ) +
    themeMetaPaper() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_line(color = "grey92", linewidth = 0.3),
      legend.position = "bottom"
    )
}
