# ============================================================================
# 13-Heatmap_Biomarkers.R
# Generate biomarker effect-size heatmap
# ============================================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(DROMA.Meta.visualization))

project_root <- getVisProjectRoot()
defaults   <- getVisWorkflowDefaults(project_root = project_root)
vis_output <- defaults$vis_output
dir.create(vis_output, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 13: Heatmap Biomarkers ===\n")

final_path <- file.path(vis_output, "collected_final_biomarkers.csv")
if (!file.exists(final_path)) {
  stop("Run 10-Collect_Results.R first", call. = FALSE)
}

final_biomarkers <- fread(final_path)

buildFallbackHeatmap <- function(final_biomarkers,
                                 es_col = "effect_size_ctrdb",
                                 direction_col = "direction",
                                 title = "Final biomarker effect sizes") {
  final_biomarkers <- as.data.table(copy(final_biomarkers))
  final_biomarkers[, pair := paste(drug, tumor_type, sep = "\n")]

  mat_dt <- dcast(final_biomarkers, name ~ pair, value.var = es_col, fun.aggregate = mean)
  rn <- mat_dt$name
  mat_dt[, name := NULL]
  mat <- as.matrix(mat_dt)
  rownames(mat) <- rn
  mat[is.nan(mat)] <- NA_real_

  dir_colors <- getMetaVisColors("direction")
  ha_top <- NULL
  if (direction_col %in% names(final_biomarkers)) {
    dir_map <- final_biomarkers[, .(direction = unique(get(direction_col))[1]), by = pair]
    dir_vec <- stats::setNames(dir_map$direction, dir_map$pair)
    ha_top <- ComplexHeatmap::HeatmapAnnotation(
      Direction = dir_vec[colnames(mat)],
      col = list(Direction = dir_colors)
    )
  }

  max_abs <- max(abs(mat), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  hm_colors <- getMetaVisColors("heatmap")
  col_fun <- circlize::colorRamp2(
    c(-max_abs, 0, max_abs),
    c(hm_colors[["low"]], hm_colors[["mid"]], hm_colors[["high"]])
  )

  ComplexHeatmap::Heatmap(
    mat,
    name = "Effect\nsize",
    col = col_fun,
    top_annotation = ha_top,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 8),
    column_names_gp = grid::gpar(fontsize = 9),
    column_title = paste0(title, " (fallback: no clustering)"),
    column_title_gp = grid::gpar(fontface = "bold", fontsize = 12),
    na_col = "grey92",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    rect_gp = grid::gpar(col = "white", lwd = 0.5)
  )
}

if (nrow(final_biomarkers) == 0) {
  cat("  SKIP heatmap: no final biomarkers\n")
} else {
  ht <- tryCatch(
    plotBiomarkerHeatmap(final_biomarkers,
      es_col = "effect_size_ctrdb",
      title = "Final biomarker effect sizes",
      show_cell_text = nrow(final_biomarkers) <= 40
    ),
    error = function(e) {
      cat("  WARN heatmap clustering failed during build; using fallback heatmap\n")
      buildFallbackHeatmap(final_biomarkers,
        es_col = "effect_size_ctrdb",
        title = "Final biomarker effect sizes"
      )
    }
  )

  pdf_path <- file.path(vis_output, "biomarker_heatmap.pdf")
  tryCatch(
    saveComplexHeatmapPdf(ht, pdf_path,
      width = max(7.2, length(unique(paste(final_biomarkers$drug, final_biomarkers$tumor_type))) * 2 + 3),
      height = max(5, nrow(final_biomarkers) * 0.3 + 3)
    ),
    error = function(e) {
      cat("  WARN heatmap clustering failed during save; retrying without clustering\n")
      ht_fallback <- buildFallbackHeatmap(final_biomarkers,
        es_col = "effect_size_ctrdb",
        title = "Final biomarker effect sizes"
      )
      saveComplexHeatmapPdf(ht_fallback, pdf_path,
        width = max(7.2, length(unique(paste(final_biomarkers$drug, final_biomarkers$tumor_type))) * 2 + 3),
        height = max(5, nrow(final_biomarkers) * 0.3 + 3)
      )
    }
  )
  cat(sprintf("  OK biomarker_heatmap.pdf\n"))
}

cat("\nDone.\n")
