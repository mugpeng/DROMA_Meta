# Pipeline overview / funnel visualization ----

#' Plot Pipeline Funnel Summary
#'
#' @description Creates a horizontal funnel chart showing the number of features
#' remaining after each workflow stage. Bars are centred and coloured by stage.
#' Each bar annotates the count and the attrition percentage from the previous stage.
#' @param summary_dt A \code{data.table} produced by \code{collectWorkflowResults()}
#'   or from \code{meta_workflow_batch_summary.csv}. Uses whichever pipeline
#'   count columns are available among \code{n_batch_cell},
#'   \code{n_selected_genes}, \code{n_ad_filtered}, \code{n_clinical_sig}, and
#'   \code{n_final_biomarkers}.
#' @param title Plot title.
#' @return A ggplot2 object.
#' @export
plotPipelineFunnel <- function(summary_dt,
                               title = "Biomarker discovery pipeline") {
  summary_dt <- data.table::as.data.table(summary_dt)
  if (!nrow(summary_dt)) {
    stop("summary_dt has no rows", call. = FALSE)
  }

  stages <- c(
    "All genes tested"  = "n_batch_cell",
    "Preclinical intersection" = "n_selected_genes",
    "TCGA-AD filtered"  = "n_ad_filtered",
    "Clinical significant" = "n_clinical_sig",
    "Clinical validated" = "n_final_biomarkers"
  )
  available_stages <- stages[stages %in% names(summary_dt)]
  if (!length(available_stages)) {
    stop(
      "summary_dt must contain at least one pipeline count column: ",
      paste(unname(stages), collapse = ", "),
      call. = FALSE
    )
  }

  # Aggregate across all pairs
  agg <- summary_dt[, lapply(.SD, sum, na.rm = TRUE),
    .SDcols = unname(available_stages)
  ]

  counts <- vapply(available_stages, function(k) {
    as.numeric(agg[[k]])
  }, numeric(1))

  stage_colors <- getMetaVisColors("stage")
  funnel_colors <- c(
    "All genes tested" = stage_colors[["Cell batch"]],
    "Preclinical intersection" = stage_colors[["Intersection"]],
    "TCGA-AD filtered" = stage_colors[["TCGA-AD filtered"]],
    "Clinical significant" = stage_colors[["Clinical validated"]],
    "Clinical validated" = stage_colors[["Clinical validated"]]
  )

  plot_dt <- data.table::data.table(
    stage = factor(names(available_stages), levels = rev(names(available_stages))),
    count = counts,
    stage_key = unname(available_stages)
  )

  # Compute attrition percentage
  plot_dt[, pct := {
    prev <- data.table::shift(count, type = "lag")
    ifelse(is.na(prev), NA_character_, sprintf("%.0f%%", count / prev * 100))
  }]

  ggplot2::ggplot(plot_dt, ggplot2::aes(x = count, y = stage, fill = stage)) +
    ggplot2::geom_col(width = 0.65, alpha = 0.92) +
    ggplot2::geom_text(
      ggplot2::aes(label = format(count, big.mark = ",")),
      hjust = -0.15, size = 3.8, fontface = "bold", color = "black"
    ) +
    ggplot2::geom_text(
      data = plot_dt[!is.na(pct)],
      ggplot2::aes(x = count * 0.5, label = pct),
      hjust = 0.5, size = 3.2, color = "white", fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = funnel_colors[names(available_stages)], guide = "none") +
    ggplot2::labs(
      title = title,
      x = "Number of features",
      y = NULL
    ) +
    themeMetaPaper() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 11, face = "bold"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
    ) +
    ggplot2::expand_limits(x = max(1, max(counts, na.rm = TRUE) * 1.18))
}

#' Plot Number of Final Biomarkers per Drug-Tumor Pair
#'
#' @description Horizontal grouped bar chart showing how many biomarkers each
#' drug-tumor pair yielded, coloured by direction (Up/Down).
#' @param final_biomarkers A \code{data.table} from \code{collectAllFinalBiomarkers()}.
#' @param title Plot title.
#' @return A ggplot2 object.
#' @export
plotBiomarkerCountByPair <- function(final_biomarkers,
                                     title = "Final biomarkers per drug-tumor pair") {
  final_biomarkers <- data.table::as.data.table(final_biomarkers)
  if (!nrow(final_biomarkers) || !all(c("drug", "tumor_type") %in% names(final_biomarkers))) {
    stop("final_biomarkers must contain 'drug' and 'tumor_type' columns", call. = FALSE)
  }

  final_biomarkers[, pair := paste(drug, tumor_type, sep = "\n")]
  dir_colors <- getMetaVisColors("direction")

  if ("direction" %in% names(final_biomarkers)) {
    count_dt <- final_biomarkers[, .N, by = .(pair, direction)]
    count_dt[, direction := factor(direction, levels = c("Up", "Down"))]

    p <- ggplot2::ggplot(count_dt, ggplot2::aes(x = stats::reorder(pair, N, sum), y = N, fill = direction)) +
      ggplot2::geom_col(width = 0.7) +
      ggplot2::geom_text(
        ggplot2::aes(label = N),
        position = ggplot2::position_stack(vjust = 0.5),
        size = 3.5, color = "white", fontface = "bold"
      ) +
      ggplot2::scale_fill_manual(values = dir_colors, name = "Direction")
  } else {
    count_dt <- final_biomarkers[, .N, by = .(pair)]
    p <- ggplot2::ggplot(count_dt, ggplot2::aes(x = stats::reorder(pair, N), y = N)) +
      ggplot2::geom_col(fill = "#5B9BD5", width = 0.7) +
      ggplot2::geom_text(ggplot2::aes(label = N), hjust = -0.2, size = 3.5)
  }

  p +
    ggplot2::coord_flip() +
    ggplot2::labs(title = title, x = NULL, y = "Number of biomarkers") +
    themeMetaPaper() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}
