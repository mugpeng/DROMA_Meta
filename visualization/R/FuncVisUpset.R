# Upset / intersection visualization ----

#' Plot Multi-Stage Intersection as Upset Plot
#'
#' @description Produces a publication-quality Upset plot showing the overlap
#' between cell-significant, PDC/PDX-significant, and clinical-significant gene
#' sets. Includes combination size bars and left annotation for set sizes.
#' @param cell_sig_genes Character vector of gene names significant in cell line data.
#' @param pdcpdx_sig_genes Character vector of gene names significant in PDC/PDX data.
#' @param clinical_sig_genes Character vector of gene names significant in clinical data.
#' @param final_genes Optional character vector of final validated biomarkers.
#' @param title Plot title.
#' @return A \code{ComplexHeatmap} object.
#' @export
plotStageUpset <- function(cell_sig_genes,
                           pdcpdx_sig_genes,
                           clinical_sig_genes,
                           final_genes = NULL,
                           title = "Biomarker overlap across stages") {
  gene_sets <- list(
    "Cell line" = unique(cell_sig_genes),
    "PDC/PDX"  = unique(pdcpdx_sig_genes),
    "Clinical" = unique(clinical_sig_genes)
  )

  # Remove empty sets
  gene_sets <- gene_sets[vapply(gene_sets, length, integer(1)) > 0]
  if (length(gene_sets) < 2) {
    stop("Need at least two non-empty gene sets for upset plot", call. = FALSE)
  }

  combination_mat <- ComplexHeatmap::make_comb_mat(gene_sets)

  ev_colors <- getMetaVisColors("evidence")
  set_colors <- ev_colors[names(gene_sets)]

  ht <- ComplexHeatmap::UpSet(
    combination_mat,
    top_annotation = ComplexHeatmap::upset_top_annotation(
      combination_mat,
      add_numbers = TRUE,
      gp = grid::gpar(fill = "#5B9BD5", col = NA),
      annotation_legend_param = list(
        labels_gp = grid::gpar(fontsize = 9),
        title_gp = grid::gpar(fontface = "bold", fontsize = 10)
      )
    ),
    left_annotation = ComplexHeatmap::rowAnnotation(
      "Set size" = ComplexHeatmap::anno_barplot(
        ComplexHeatmap::set_size(combination_mat),
        border = FALSE,
        gp = grid::gpar(fill = set_colors, col = NA),
        bar_width = 0.7
      )
    ),
    column_title = title,
    column_title_gp = grid::gpar(fontface = "bold", fontsize = 12),
    row_names_gp = grid::gpar(fontsize = 10, fontface = "bold"),
    set_order = names(gene_sets)
  )

  ht
}

#' Build Upset Plot from One Drug-Tumor Output Directory
#'
#' @description Reads the intermediate CSVs from a single pair output directory
#' and calls \code{plotStageUpset()}.
#' @param pair_dir Path to one drug-tumor output directory.
#' @param title Plot title. NULL to auto-generate from directory path.
#' @return A \code{ComplexHeatmap} object.
#' @export
plotUpsetFromPairDir <- function(pair_dir, title = NULL) {
  readGenes <- function(path, col = "name") {
    if (!file.exists(path)) return(character(0))
    dt <- data.table::fread(path)
    if (!col %in% names(dt) || !nrow(dt)) return(character(0))
    unique(as.character(dt[[col]]))
  }

  cell     <- readGenes(file.path(pair_dir, "mRNA_cell_sig.csv"))
  pdcpdx   <- readGenes(file.path(pair_dir, "mRNA_pdcpdx_sig.csv"))
  clinical <- readGenes(file.path(pair_dir, "clinical_sig_mRNA.csv"))
  final    <- readGenes(file.path(pair_dir, "final_biomarkers.csv"))

  if (is.null(title)) {
    drug_name  <- basename(dirname(pair_dir))
    tumor_name <- basename(pair_dir)
    title <- paste(drug_name, "-", gsub("_", " ", tumor_name))
  }

  plotStageUpset(cell, pdcpdx, clinical, final, title = title)
}
