# Summary bar-plot visualizations ----

#' Plot Biomarker Direction Distribution
#'
#' @description Horizontal stacked bar chart showing the number of Up and Down
#' biomarkers per drug-tumor pair with percentage labels.
#' @param final_biomarkers A \code{data.table} from \code{collectAllFinalBiomarkers()}.
#' @param title Plot title.
#' @return A ggplot2 object.
#' @export
plotDirectionSummary <- function(final_biomarkers,
                                 title = "Biomarker direction by drug-tumor pair") {
  final_biomarkers <- data.table::as.data.table(final_biomarkers)
  if (!all(c("direction", "drug", "tumor_type") %in% names(final_biomarkers))) {
    stop("final_biomarkers must contain 'direction', 'drug', 'tumor_type'", call. = FALSE)
  }

  final_biomarkers[, pair := paste(drug, tumor_type, sep = "\n")]
  dir_colors <- getMetaVisColors("direction")

  count_dt <- final_biomarkers[, .N, by = .(pair, direction)]
  count_dt[, direction := factor(direction, levels = c("Up", "Down"))]

  # Add percentage within each pair
  count_dt[, total := sum(N), by = pair]
  count_dt[, pct := sprintf("%.0f%%", N / total * 100)]

  ggplot2::ggplot(count_dt, ggplot2::aes(x = reorder(pair, total, sum), y = N, fill = direction)) +
    ggplot2::geom_col(width = 0.65, alpha = 0.92) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0(N, " (", pct, ")")),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 3.2, color = "white", fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = dir_colors, name = "Direction") +
    ggplot2::coord_flip() +
    ggplot2::labs(title = title, x = NULL, y = "Number of biomarkers") +
    themeMetaPaper() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}

#' Plot Workflow Stage Statistics as Dumbbell Chart
#'
#' @description Dumbbell chart where each row is a drug-tumor pair. The left end
#' marks the preclinical intersection count and the right end marks the final
#' validated count. The connecting line length represents the attrition.
#' @param summary_dt A \code{data.table} from \code{collectWorkflowResults()}.
#' @param title Plot title.
#' @return A ggplot2 object.
#' @export
plotStageCountComparison <- function(summary_dt,
                                     title = "Biomarker attrition across validation stages") {
  summary_dt <- data.table::as.data.table(summary_dt)
  if (!nrow(summary_dt)) {
    stop("summary_dt has no rows", call. = FALSE)
  }

  required <- c("n_selected_genes", "n_final_biomarkers")
  if (!all(required %in% names(summary_dt))) {
    stop("summary_dt must contain n_selected_genes and n_final_biomarkers", call. = FALSE)
  }

  summary_dt[, pair := paste(drug, tumor_type, sep = "\n")]

  # Dumbbell: wide format
  dumb_dt <- data.table::data.table(
    pair      = summary_dt$pair,
    start     = summary_dt$n_selected_genes,
    end       = summary_dt$n_final_biomarkers,
    attrition = summary_dt$n_selected_genes - summary_dt$n_final_biomarkers
  )
  dumb_dt[, pct_retained := {
    ifelse(start == 0, "N/A", sprintf("%.0f%%", end / start * 100))
  }]
  dumb_dt <- dumb_dt[order(start), ]
  dumb_dt[, pair := factor(pair, levels = pair)]

  ev_colors <- getMetaVisColors("evidence")

  # Build a long-format table for proper legend
  legend_dt <- data.table::rbindlist(list(
    data.table::data.table(pair = dumb_dt$pair, x = dumb_dt$start, Stage = "Preclinical intersection"),
    data.table::data.table(pair = dumb_dt$pair, x = dumb_dt$end, Stage = "Final validated")
  ))
  legend_dt[, pair := factor(pair, levels = levels(dumb_dt$pair))]
  legend_dt[, Stage := factor(Stage, levels = c("Preclinical intersection", "Final validated"))]

  stage_colors <- c(
    "Preclinical intersection" = ev_colors[["Cell line"]],
    "Final validated" = ev_colors[["Clinical"]]
  )

  ggplot2::ggplot() +
    # Connecting segment
    ggplot2::geom_segment(
      data = dumb_dt,
      ggplot2::aes(x = start, xend = end, y = pair, yend = pair),
      color = "grey60", linewidth = 1.2
    ) +
    # Points with proper legend
    ggplot2::geom_point(
      data = legend_dt,
      ggplot2::aes(x = x, y = pair, fill = Stage),
      shape = 21, size = 4, stroke = 0.8, color = "grey30"
    ) +
    ggplot2::scale_fill_manual(values = stage_colors, name = NULL) +
    # Labels
    ggplot2::geom_text(
      data = dumb_dt,
      ggplot2::aes(x = start, y = pair, label = start),
      hjust = 1.2, size = 3, color = "grey30"
    ) +
    ggplot2::geom_text(
      data = dumb_dt,
      ggplot2::aes(x = end, y = pair, label = paste0(end, " (", pct_retained, ")")),
      hjust = -0.2, size = 3, color = "grey30"
    ) +
    ggplot2::labs(
      title = title,
      x = "Number of features",
      y = NULL
    ) +
    themeMetaPaper() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}
