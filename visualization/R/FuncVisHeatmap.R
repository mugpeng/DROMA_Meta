# Heatmap visualization ----

#' Plot Biomarker Effect-Size Heatmap
#'
#' @description Creates a publication-quality heatmap where rows are gene names
#' and columns are drug-tumor pairs, coloured by a diverging effect-size scale.
#' Includes a top annotation bar for biomarker direction and optional cell text
#' for effect-size values.
#' @param final_biomarkers A \code{data.table} from \code{collectAllFinalBiomarkers()}.
#' @param es_col Column name for the effect size to plot. Defaults to
#'   \code{"effect_size_ctrdb"}.
#' @param direction_col Column name for direction annotation. NULL to skip.
#' @param title Plot title.
#' @param show_row_names Logical; show gene names on y-axis.
#' @param show_cell_text Logical; show effect-size values inside cells.
#' @param row_font_size Font size for row names.
#' @return A \code{ComplexHeatmap} object.
#' @export
plotBiomarkerHeatmap <- function(final_biomarkers,
                                 es_col = "effect_size_ctrdb",
                                 direction_col = "direction",
                                 title = "Biomarker effect sizes",
                                 show_row_names = TRUE,
                                 show_cell_text = FALSE,
                                 row_font_size = 8) {
  final_biomarkers <- data.table::as.data.table(final_biomarkers)
  required <- c("name", es_col, "drug", "tumor_type")
  if (!all(required %in% names(final_biomarkers))) {
    stop("final_biomarkers must contain columns: ", paste(required, collapse = ", "), call. = FALSE)
  }

  final_biomarkers[, pair := paste(drug, tumor_type, sep = "\n")]

  # Wide matrix: gene x pair
  mat <- data.table::dcast(final_biomarkers, name ~ pair,
    value.var = es_col, fun.aggregate = mean
  )
  rn <- mat$name
  mat[, name := NULL]
  mat <- as.matrix(mat)
  rownames(mat) <- rn

  # Direction annotation
  dir_colors <- getMetaVisColors("direction")
  ha_top <- NULL
  if (!is.null(direction_col) && direction_col %in% names(final_biomarkers)) {
    dir_map <- final_biomarkers[, .(direction = unique(get(direction_col))), by = pair]
    dir_vec <- setNames(dir_map$direction, dir_map$pair)
    ha_top <- ComplexHeatmap::HeatmapAnnotation(
      Direction = dir_vec[colnames(mat)],
      col = list(Direction = dir_colors),
      annotation_legend_param = list(
        Direction = list(
          title = "Direction",
          title_gp = grid::gpar(fontface = "bold", fontsize = 10),
          labels_gp = grid::gpar(fontsize = 9)
        )
      )
    )
  }

  # Diverging colour scale
  max_abs <- max(abs(mat), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  hm_colors <- getMetaVisColors("heatmap")
  col_fun <- circlize::colorRamp2(
    c(-max_abs, 0, max_abs),
    c(hm_colors[["low"]], hm_colors[["mid"]], hm_colors[["high"]])
  )

  # Cell text layer
  cell_fun <- NULL
  if (show_cell_text && nrow(mat) <= 40) {
    cell_fun <- function(j, i, x, y, w, h, fill) {
      val <- mat[i, j]
      if (is.na(val)) return()
      text_col <- if (abs(val) > max_abs * 0.6) "white" else "black"
      grid::grid.text(sprintf("%.2f", val), x, y,
        gp = grid::gpar(fontsize = 7, col = text_col)
      )
    }
  }

  ComplexHeatmap::Heatmap(
    mat,
    name             = "Effect\nsize",
    col              = col_fun,
    top_annotation   = ha_top,
    show_row_names   = show_row_names,
    show_column_names = TRUE,
    row_names_gp     = grid::gpar(fontsize = row_font_size),
    column_names_gp  = grid::gpar(fontsize = 9),
    column_title     = title,
    column_title_gp  = grid::gpar(fontface = "bold", fontsize = 12),
    na_col           = "grey92",
    cluster_rows     = TRUE,
    cluster_columns  = TRUE,
    show_row_dend    = FALSE,
    show_column_dend = FALSE,
    rect_gp          = grid::gpar(col = "white", lwd = 0.5),
    cell_fun         = cell_fun,
    heatmap_legend_param = list(
      title = "Effect size",
      title_gp = grid::gpar(fontface = "bold", fontsize = 10),
      labels_gp = grid::gpar(fontsize = 9),
      direction = "horizontal",
      legend_width = unit(4, "cm")
    )
  )
}
