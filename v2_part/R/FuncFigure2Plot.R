suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

plotFigure2ADistanceDistribution <- function(distance_dt) {
  plot_dt <- copy(distance_dt)[!is.na(distance)]
  plot_dt[, stage_label := orderFigure2StageLabels(labelFigure2Stage(stage))]

  ggplot(plot_dt, aes(x = distance, fill = stage_label)) +
    geom_histogram(aes(y = after_stat(count / sum(count))), binwidth = 1, boundary = -0.5, color = "white") +
    facet_wrap(~stage_label, nrow = 1) +
    scale_x_continuous(breaks = seq(0, max(plot_dt$distance, na.rm = TRUE), by = 1)) +
    labs(
      title = "Figure 2A: Biomarker-to-target distance distributions",
      x = "Distance to drug target",
      y = "Proportion",
      fill = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )
}

plotFigure2BDistanceVsEffect <- function(distance_dt) {
  plot_dt <- copy(distance_dt)[!is.na(distance) & !is.na(effect_size)]
  plot_dt[, distance_bin := factor(distance, levels = sort(unique(distance)))]
  plot_dt[, stage_label := orderFigure2StageLabels(labelFigure2Stage(stage))]

  ggplot(plot_dt, aes(x = distance_bin, y = effect_size, fill = distance_bin)) +
    geom_boxplot(outlier.alpha = 0.2) +
    geom_jitter(width = 0.15, alpha = 0.25, size = 0.9) +
    facet_wrap(~stage_label, scales = "free_y") +
    labs(
      title = "Figure 2B: Target distance versus effect size",
      x = "Shortest path distance",
      y = "Effect size"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )
}

plotFigure2CDistanceVsStage <- function(distance_dt) {
  plot_dt <- copy(distance_dt)[!is.na(distance)]

  ggplot(plot_dt, aes(x = stage, y = distance, fill = stage)) +
    geom_violin(alpha = 0.35, scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.8) +
    labs(
      title = "Figure 2C: Target distance across validation stages",
      x = NULL,
      y = "Shortest path distance"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )
}
