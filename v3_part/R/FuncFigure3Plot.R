# Figure 3 breast-state plots ----

plotFigure3ProgramGeneOverlay <- function(marker_gene_overlay) {
  dt <- data.table::as.data.table(marker_gene_overlay)
  required <- c("category", "program", "name", "drug", "effect_size")
  if (!all(required %in% names(dt))) {
    stop("marker_gene_overlay must contain columns: ", paste(required, collapse = ", "), call. = FALSE)
  }
  if (!nrow(dt)) stop("marker_gene_overlay is empty", call. = FALSE)

  max_abs <- max(abs(dt$effect_size), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1

  dt[, label := paste(program, name, sep = ": ")]

  ggplot2::ggplot(dt, ggplot2::aes(x = drug, y = label, fill = effect_size)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::facet_grid(category ~ program, scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_gradient2(
      low = getMetaVisColors("heatmap")[["low"]],
      mid = getMetaVisColors("heatmap")[["mid"]],
      high = getMetaVisColors("heatmap")[["high"]],
      limits = c(-max_abs, max_abs),
      name = "Effect size"
    ) +
    ggplot2::labs(
      title = "PAM50 and breast-state marker overlay",
      x = NULL,
      y = NULL
    ) +
    themeMetaPaper(base_size = 9) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      strip.text.y = ggplot2::element_text(angle = 0)
    )
}

plotFigure3ClusterProgramAnnotation <- function(cluster_program_annotation, top_n = 3L) {
  dt <- data.table::as.data.table(cluster_program_annotation)
  required <- c("gene_cluster", "program", "overlap_n", "overlap_frac", "rank_within_cluster")
  if (!all(required %in% names(dt))) {
    stop("cluster_program_annotation must contain columns: ", paste(required, collapse = ", "), call. = FALSE)
  }
  dt <- dt[rank_within_cluster <= top_n]
  if (!nrow(dt)) stop("cluster_program_annotation is empty", call. = FALSE)

  dt[, y_label := stats::reorder(program, overlap_frac)]
  ggplot2::ggplot(dt, ggplot2::aes(x = overlap_frac, y = y_label, fill = gene_cluster)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::facet_wrap(~ gene_cluster, scales = "free_y") +
    ggplot2::labs(
      title = "Breast-state programs annotating each shared gene cluster",
      x = "Marker overlap fraction",
      y = NULL,
      fill = "Gene cluster"
    ) +
    themeMetaPaper(base_size = 10)
}

plotFigure3PseudoAtlasAnnotation <- function(pseudo_atlas_annotation, top_n = 3L) {
  dt <- data.table::as.data.table(pseudo_atlas_annotation)
  required <- c("gene_cluster", "cell_state", "overlap_n", "overlap_frac", "rank_within_cluster")
  if (!all(required %in% names(dt))) {
    stop("pseudo_atlas_annotation must contain columns: ", paste(required, collapse = ", "), call. = FALSE)
  }
  dt <- dt[rank_within_cluster <= top_n]
  if (!nrow(dt)) stop("pseudo_atlas_annotation is empty", call. = FALSE)

  dt[, y_label := stats::reorder(cell_state, overlap_frac)]
  ggplot2::ggplot(dt, ggplot2::aes(x = overlap_frac, y = y_label, fill = gene_cluster)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = paste0("n=", overlap_n)),
      hjust = -0.1,
      size = 3
    ) +
    ggplot2::facet_wrap(~ gene_cluster, scales = "free_y") +
    ggplot2::coord_cartesian(xlim = c(0, max(dt$overlap_frac, na.rm = TRUE) * 1.15)) +
    ggplot2::labs(
      title = "Marker-guided pseudo-atlas annotation of shared gene clusters",
      x = "Marker overlap fraction",
      y = NULL,
      fill = "Gene cluster"
    ) +
    themeMetaPaper(base_size = 10)
}
