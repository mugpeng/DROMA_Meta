# Preclinical-Clinical concordance visualization ----

#' Plot Preclinical vs Clinical Effect-Size Concordance
#'
#' @description Scatter plot comparing preclinical meta-analysis effect sizes
#' with clinical effect sizes for all final biomarkers. Includes a diagonal
#' reference line, linear regression fit with 95% CI, and four-quadrant
#' annotation based on effect-size thresholds.
#' @param final_biomarkers A \code{data.table} from
#'   \code{collectAllFinalBiomarkers()}. Must contain preclinical and clinical
#'   effect-size columns and \code{name}.
#' @param preclinical_es_col Column with preclinical effect sizes.
#' @param clinical_es_col Column with clinical effect sizes.
#' @param color_col Column used for colour grouping. Set to NULL for single colour.
#' @param method Correlation method for annotation.
#' @param title Plot title.
#' @param point_size Point size.
#' @param point_alpha Point alpha.
#' @return A ggplot2 object.
#' @export
plotConcordanceScatter <- function(final_biomarkers,
                                   preclinical_es_col = "effect_size_pdcpdx_invitro",
                                   clinical_es_col = "effect_size_ctrdb",
                                   color_col = "drug",
                                   method = "spearman",
                                   title = "Preclinical vs clinical effect sizes",
                                   point_size = 2.5,
                                   point_alpha = 0.75) {
  final_biomarkers <- as.data.frame(final_biomarkers)
  required <- c("name", preclinical_es_col, clinical_es_col)
  if (!all(required %in% colnames(final_biomarkers))) {
    stop("final_biomarkers must contain columns: ", paste(required, collapse = ", "), call. = FALSE)
  }

  plot_df <- final_biomarkers[
    !is.na(final_biomarkers[[preclinical_es_col]]) &
    !is.na(final_biomarkers[[clinical_es_col]]), , drop = FALSE
  ]
  if (!nrow(plot_df)) {
    stop("No rows with non-NA values in both effect-size columns", call. = FALSE)
  }

  # Correlation
  cor_test <- stats::cor.test(
    plot_df[[preclinical_es_col]],
    plot_df[[clinical_es_col]],
    method = method
  )
  cor_label <- sprintf(
    "%s %.2f, P = %.1e, n = %d",
    toupper(substr(method, 1, 1)),
    cor_test$estimate,
    cor_test$p.value,
    nrow(plot_df)
  )

  # Build base aes
  aes_mapping <- ggplot2::aes(
    x = .data[[preclinical_es_col]],
    y = .data[[clinical_es_col]]
  )
  if (!is.null(color_col) && color_col %in% colnames(plot_df)) {
    aes_mapping <- ggplot2::aes(
      x = .data[[preclinical_es_col]],
      y = .data[[clinical_es_col]],
      color = .data[[color_col]]
    )
  }

  p <- ggplot2::ggplot(plot_df, aes_mapping) +
    # Zero reference lines
    ggplot2::geom_hline(yintercept = 0, linetype = 3, color = "grey60", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = 0, linetype = 3, color = "grey60", linewidth = 0.3) +
    # Diagonal
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey40", linewidth = 0.4) +
    # LM fit with CI
    ggplot2::geom_smooth(
      method = "lm", se = TRUE,
      color = "grey30", fill = "grey80", alpha = 0.3,
      linewidth = 0.6, linetype = 1
    ) +
    # Points
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    # Correlation annotation
    ggplot2::annotate("text",
      x = -Inf, y = Inf,
      label = cor_label,
      hjust = -0.05, vjust = 1.5,
      size = 3.5, fontface = "bold", color = "grey20"
    ) +
    ggplot2::labs(
      title = title,
      x = "Preclinical effect size",
      y = "Clinical effect size"
    ) +
    themeMetaPaper()

  if (!is.null(color_col) && color_col %in% colnames(plot_df)) {
    pair_colors <- getMetaVisColors("pair")
    present_pairs <- intersect(names(pair_colors), unique(plot_df[[color_col]]))
    if (length(present_pairs) > 0) {
      p <- p + ggplot2::scale_color_manual(
        values = pair_colors,
        name = tools::toTitleCase(color_col)
      )
    } else {
      p <- p + ggplot2::labs(color = tools::toTitleCase(color_col))
    }
  }

  # Gene labels (when n <= 50)
  if (requireNamespace("ggrepel", quietly = TRUE) && nrow(plot_df) <= 50) {
    p <- p + ggrepel::geom_text_repel(
      data = plot_df,
      ggplot2::aes(label = name),
      size = 2.8, max.overlaps = 15, color = "grey30",
      segment.color = "grey60", segment.size = 0.2,
      box.padding = 0.3, point.padding = 0.15
    )
  }

  # Density contours when many points
  if (nrow(plot_df) > 50) {
    p <- p + ggplot2::geom_density2d(
      color = "grey50", linewidth = 0.3, alpha = 0.5
    )
  }

  p
}
