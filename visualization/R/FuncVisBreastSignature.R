# Shared selected-gene signature visualizations ----

#' Collect AD-Filtered Selected Genes Across All Drug-Tumor Pairs
#'
#' @description Reads every \code{selected_genes_ad_filtered.csv} under
#' \code{output_base} and annotates each row with drug and tumor type from the
#' workflow directory structure.
#' @param output_base Root output directory produced by the meta workflow.
#' @return A \code{data.table} containing merged AD-filtered selected genes.
#' @export
collectAllSelectedAdFiltered <- function(output_base) {
  if (!dir.exists(output_base)) {
    stop("output_base directory does not exist: ", output_base, call. = FALSE)
  }

  all_files <- list.files(
    output_base,
    pattern = "^selected_genes_ad_filtered\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (!length(all_files)) {
    warning("No selected_genes_ad_filtered.csv files found under: ", output_base, call. = FALSE)
    return(data.table::data.table())
  }

  data.table::rbindlist(lapply(all_files, function(f) {
    dt <- data.table::fread(f)
    if (!nrow(dt)) return(NULL)
    dt[, drug := basename(dirname(dirname(f)))]
    dt[, tumor_type := basename(dirname(f))]
    dt
  }), fill = TRUE)
}

#' Prepare a Shared Drug-Response Signature Matrix
#'
#' @description Filters AD-selected genes to one tumor type, retains genes
#' shared by at least \code{min_drugs_per_gene} drugs, and returns a gene-by-drug
#' effect-size matrix for heatmap plotting.
#' @param selected_ad_filtered Merged AD-filtered selected-gene table.
#' @param tumor_type Tumor type to analyze. Directory-safe labels such as
#'   \code{"breast_cancer"} and display labels such as \code{"breast cancer"}
#'   are both accepted.
#' @param effect_col Column containing effect size values.
#' @param min_drugs_per_gene Minimum number of drugs per retained gene.
#' @param top_n_genes Maximum number of genes retained by cross-drug variance.
#' @param min_genes_per_drug Minimum AD-filtered selected genes required for a
#'   drug to enter the heatmap.
#' @return A list with \code{matrix}, \code{genes}, and \code{drugs}.
#' @export
prepareSharedSignatureMatrix <- function(selected_ad_filtered,
                                         tumor_type = "breast_cancer",
                                         effect_col = NULL,
                                         min_drugs_per_gene = 4L,
                                         top_n_genes = 200L,
                                         min_genes_per_drug = 1L) {
  dt <- data.table::as.data.table(selected_ad_filtered)
  if (is.null(effect_col)) {
    effect_candidates <- c("effect_size_cell", "effect_size_cell_invitro", "effect_size")
    effect_col <- effect_candidates[effect_candidates %in% names(dt)][1]
    if (is.na(effect_col)) {
      stop(
        "No supported effect-size column found. Expected one of: ",
        paste(effect_candidates, collapse = ", "),
        call. = FALSE
      )
    }
  }
  required <- c("name", effect_col, "drug", "tumor_type")
  if (!all(required %in% names(dt))) {
    stop("selected_ad_filtered must contain columns: ", paste(required, collapse = ", "), call. = FALSE)
  }

  tumor_keys <- unique(c(tumor_type, gsub(" ", "_", tumor_type), gsub("_", " ", tumor_type)))
  dt <- dt[tumor_type %in% tumor_keys]
  if (!nrow(dt)) {
    stop("No selected AD-filtered genes found for tumor_type: ", tumor_type, call. = FALSE)
  }

  dt <- dt[is.finite(get(effect_col))]
  drug_counts <- dt[, .N, by = drug][N >= min_genes_per_drug]
  dt <- dt[drug %in% drug_counts$drug]
  if (!nrow(dt)) {
    stop("No drugs meet min_genes_per_drug for tumor_type: ", tumor_type, call. = FALSE)
  }

  gene_counts <- dt[, .(n_drugs = data.table::uniqueN(drug)), by = name]
  shared_genes <- gene_counts[n_drugs >= min_drugs_per_gene]
  if (!nrow(shared_genes)) {
    stop(
      "No genes are shared by at least ", min_drugs_per_gene,
      " drugs for tumor_type: ", tumor_type,
      call. = FALSE
    )
  }

  wide <- data.table::dcast(
    dt[name %in% shared_genes$name],
    name ~ drug,
    value.var = effect_col,
    fun.aggregate = mean
  )
  rn <- wide$name
  wide[, name := NULL]
  mat <- as.matrix(wide)
  rownames(mat) <- rn

  row_var <- apply(mat, 1L, function(x) stats::var(x, na.rm = TRUE))
  row_var[!is.finite(row_var)] <- 0
  keep <- names(sort(row_var, decreasing = TRUE))[seq_len(min(length(row_var), top_n_genes))]
  mat <- mat[keep, , drop = FALSE]

  genes <- shared_genes[name %in% rownames(mat)][order(-n_drugs, name)]
  drugs <- data.table::data.table(
    drug = colnames(mat),
    n_selected_ad_filtered = colSums(is.finite(mat))
  )[order(-n_selected_ad_filtered, drug)]

  list(matrix = mat, genes = genes, drugs = drugs)
}

#' Plot a Shared Selected-Gene Signature Heatmap
#'
#' @description Creates a ComplexHeatmap from the matrix returned by
#' \code{prepareSharedSignatureMatrix()}.
#' @param signature Output from \code{prepareSharedSignatureMatrix()}.
#' @param title Heatmap title.
#' @param row_font_size Font size for gene labels.
#' @return A \code{ComplexHeatmap} object.
#' @export
plotSharedSignatureHeatmap <- function(signature,
                                       title = "Shared AD-filtered response signature",
                                       row_font_size = 5) {
  mat <- signature$matrix
  if (!is.matrix(mat) || !nrow(mat) || !ncol(mat)) {
    stop("signature$matrix must be a non-empty matrix", call. = FALSE)
  }

  max_abs <- max(abs(mat), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1
  hm_colors <- getMetaVisColors("heatmap")
  col_fun <- circlize::colorRamp2(
    c(-max_abs, 0, max_abs),
    c(hm_colors[["low"]], hm_colors[["mid"]], hm_colors[["high"]])
  )

  cluster_dist <- function(x) {
    x[!is.finite(x)] <- 0
    stats::dist(x)
  }

  ComplexHeatmap::Heatmap(
    mat,
    name = "Effect\nsize",
    col = col_fun,
    na_col = "grey92",
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = row_font_size),
    column_names_gp = grid::gpar(fontsize = 9),
    column_title = title,
    column_title_gp = grid::gpar(fontface = "bold", fontsize = 12),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    clustering_distance_rows = cluster_dist,
    clustering_distance_columns = cluster_dist,
    show_row_dend = FALSE,
    show_column_dend = TRUE,
    rect_gp = grid::gpar(col = "white", lwd = 0.25),
    heatmap_legend_param = list(
      title = "Effect size",
      title_gp = grid::gpar(fontface = "bold", fontsize = 10),
      labels_gp = grid::gpar(fontsize = 9),
      direction = "horizontal",
      legend_width = grid::unit(4, "cm")
    )
  )
}

#' Assign Row and Column Clusters for a Shared Signature Matrix
#'
#' @description Clusters rows and columns after replacing missing values with
#' zero, returning tidy assignments for downstream enrichment and summaries.
#' @param mat Numeric gene-by-drug matrix.
#' @param k Number of gene clusters.
#' @return A \code{data.table} with gene and cluster assignments.
#' @export
assignSharedSignatureClusters <- function(mat, k = 3L) {
  if (!is.matrix(mat) || !nrow(mat)) {
    stop("mat must be a non-empty matrix", call. = FALSE)
  }
  k <- as.integer(k)
  if (!is.finite(k) || k < 1L) {
    stop("k must be a positive integer", call. = FALSE)
  }
  k <- min(k, nrow(mat))

  mat0 <- mat
  mat0[!is.finite(mat0)] <- 0
  if (k == 1L || nrow(mat0) == 1L) {
    clusters <- rep(1L, nrow(mat0))
  } else {
    clusters <- stats::cutree(stats::hclust(stats::dist(mat0)), k = k)
  }

  data.table::data.table(
    name = rownames(mat),
    gene_cluster = paste0("Cluster ", as.integer(clusters))
  )
}

#' Assign Drug Clusters for a Shared Signature Matrix
#'
#' @description Clusters drug columns using the shared gene signature matrix.
#' @param mat Numeric gene-by-drug matrix.
#' @param k Number of drug clusters.
#' @return A \code{data.table} with drug and cluster assignments.
#' @export
assignSharedSignatureDrugClusters <- function(mat, k = 3L) {
  if (!is.matrix(mat) || !ncol(mat)) {
    stop("mat must be a non-empty matrix", call. = FALSE)
  }
  k <- as.integer(k)
  if (!is.finite(k) || k < 1L) {
    stop("k must be a positive integer", call. = FALSE)
  }
  k <- min(k, ncol(mat))

  mat0 <- mat
  mat0[!is.finite(mat0)] <- 0
  if (k == 1L || ncol(mat0) == 1L) {
    clusters <- rep(1L, ncol(mat0))
  } else {
    clusters <- stats::cutree(stats::hclust(stats::dist(t(mat0))), k = k)
  }

  data.table::data.table(
    drug = colnames(mat),
    drug_cluster = paste0("Drug cluster ", as.integer(clusters))
  )
}

#' Run GO Biological Process Enrichment by Gene Cluster
#'
#' @description Converts gene symbols to Entrez identifiers and runs
#' \code{clusterProfiler::enrichGO()} for each gene cluster.
#' @param gene_clusters Output from \code{assignSharedSignatureClusters()}.
#' @param universe Optional background gene symbols.
#' @param pvalue_cutoff P-value cutoff passed to \code{enrichGO()}.
#' @param qvalue_cutoff Q-value cutoff passed to \code{enrichGO()}.
#' @return A \code{data.table} of enrichment results.
#' @export
runClusterGoEnrichment <- function(gene_clusters,
                                   universe = NULL,
                                   pvalue_cutoff = 0.05,
                                   qvalue_cutoff = 0.2) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler is required for GO enrichment", call. = FALSE)
  }
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    stop("org.Hs.eg.db is required for GO enrichment", call. = FALSE)
  }

  dt <- data.table::as.data.table(gene_clusters)
  required <- c("name", "gene_cluster")
  if (!all(required %in% names(dt))) {
    stop("gene_clusters must contain columns: ", paste(required, collapse = ", "), call. = FALSE)
  }

  symbols <- unique(dt$name)
  gene_map <- suppressMessages(clusterProfiler::bitr(
    symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
  ))
  if (!nrow(gene_map)) {
    warning("No gene symbols could be mapped to Entrez IDs", call. = FALSE)
    return(data.table::data.table())
  }

  universe_ids <- NULL
  if (!is.null(universe)) {
    universe_map <- suppressMessages(clusterProfiler::bitr(
      unique(as.character(universe)),
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db::org.Hs.eg.db
    ))
    universe_ids <- unique(as.character(universe_map[["ENTREZID"]]))
  }

  out <- list()
  for (cluster_id in unique(dt$gene_cluster)) {
    cluster_symbols <- unique(dt[gene_cluster == cluster_id, name])
    cluster_ids <- unique(as.character(gene_map[gene_map[["SYMBOL"]] %in% cluster_symbols, "ENTREZID"]))
    if (length(cluster_ids) < 3L) next

    ego <- suppressMessages(clusterProfiler::enrichGO(
      gene = cluster_ids,
      universe = universe_ids,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      pvalueCutoff = pvalue_cutoff,
      qvalueCutoff = qvalue_cutoff,
      readable = TRUE
    ))
    ego_dt <- data.table::as.data.table(as.data.frame(ego))
    if (!nrow(ego_dt)) next
    ego_dt[, gene_cluster := cluster_id]
    out[[length(out) + 1L]] <- ego_dt
  }

  if (!length(out)) {
    return(data.table::data.table())
  }
  data.table::rbindlist(out, fill = TRUE)
}

#' Plot Cluster GO Enrichment Terms
#'
#' @description Dot plot of top GO terms per gene cluster.
#' @param go_dt Result from \code{runClusterGoEnrichment()}.
#' @param top_n Number of terms per cluster.
#' @return A ggplot object.
#' @export
plotClusterGoEnrichment <- function(go_dt, top_n = 8L) {
  go_dt <- data.table::as.data.table(go_dt)
  required <- c("Description", "pvalue", "p.adjust", "Count", "gene_cluster")
  if (!all(required %in% names(go_dt))) {
    stop("go_dt must contain columns: ", paste(required, collapse = ", "), call. = FALSE)
  }
  if (!nrow(go_dt)) {
    stop("go_dt is empty", call. = FALSE)
  }

  plot_dt <- go_dt[order(p.adjust), head(.SD, top_n), by = gene_cluster]
  plot_dt[, term := factor(Description, levels = rev(unique(Description)))]
  plot_dt[, enrichment_score := -log10(pmax(pvalue, .Machine$double.xmin))]

  ggplot2::ggplot(plot_dt, ggplot2::aes(x = gene_cluster, y = term)) +
    ggplot2::geom_point(ggplot2::aes(size = Count, color = enrichment_score)) +
    ggplot2::scale_color_gradient(low = "#6BAED6", high = "#B2182B", name = "-log10(P)") +
    ggplot2::labs(
      title = "GO Biological Process enrichment by breast signature cluster",
      x = NULL,
      y = NULL,
      size = "Genes"
    ) +
    themeMetaPaper(base_size = 10) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
}

#' Build Breast Marker Overlay Matrix
#'
#' @description Extracts known breast biology markers from AD-filtered selected
#' genes and returns a marker-by-drug effect-size table.
#' @param selected_ad_filtered Merged AD-filtered selected-gene table.
#' @param markers Gene symbols to overlay.
#' @param tumor_type Tumor type to analyze.
#' @param effect_col Effect-size column. Defaults to detected cell-line effect.
#' @return A \code{data.table}.
#' @export
buildBreastMarkerOverlay <- function(selected_ad_filtered,
                                     markers = c(
                                       "ESR1", "PGR", "ERBB2", "GRB7", "FOXA1", "GATA3",
                                       "KRT5", "KRT14", "KRT17", "EPCAM", "VIM", "CDH1",
                                       "MKI67", "TOP2A", "AURKA"
                                     ),
                                     tumor_type = "breast_cancer",
                                     effect_col = NULL) {
  dt <- data.table::as.data.table(selected_ad_filtered)
  if (is.null(effect_col)) {
    effect_candidates <- c("effect_size_cell", "effect_size_cell_invitro", "effect_size")
    effect_col <- effect_candidates[effect_candidates %in% names(dt)][1]
  }
  required <- c("name", effect_col, "drug", "tumor_type")
  if (!all(required %in% names(dt))) {
    stop("selected_ad_filtered must contain columns: ", paste(required, collapse = ", "), call. = FALSE)
  }
  tumor_keys <- unique(c(tumor_type, gsub(" ", "_", tumor_type), gsub("_", " ", tumor_type)))
  dt <- dt[tumor_type %in% tumor_keys & name %in% markers]
  if (!nrow(dt)) {
    return(data.table::data.table(
      name = character(),
      drug = character(),
      effect_size = numeric()
    ))
  }
  dt[, .(effect_size = mean(get(effect_col), na.rm = TRUE)), by = .(name, drug)]
}

#' Plot Breast Marker Overlay
#'
#' @description Tile plot of known breast marker effect sizes across drugs.
#' @param marker_dt Output from \code{buildBreastMarkerOverlay()}.
#' @return A ggplot object.
#' @export
plotBreastMarkerOverlay <- function(marker_dt) {
  marker_dt <- data.table::as.data.table(marker_dt)
  required <- c("name", "drug", "effect_size")
  if (!all(required %in% names(marker_dt))) {
    stop("marker_dt must contain columns: ", paste(required, collapse = ", "), call. = FALSE)
  }
  if (!nrow(marker_dt)) {
    stop("marker_dt is empty", call. = FALSE)
  }
  max_abs <- max(abs(marker_dt$effect_size), na.rm = TRUE)
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1

  ggplot2::ggplot(marker_dt, ggplot2::aes(x = drug, y = name, fill = effect_size)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_gradient2(
      low = getMetaVisColors("heatmap")[["low"]],
      mid = getMetaVisColors("heatmap")[["mid"]],
      high = getMetaVisColors("heatmap")[["high"]],
      limits = c(-max_abs, max_abs),
      name = "Effect size"
    ) +
    ggplot2::labs(
      title = "Breast subtype and response marker overlay",
      x = NULL,
      y = NULL
    ) +
    themeMetaPaper(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Summarize Final Biomarkers by Shared Signature Cluster
#'
#' @description Joins final biomarkers to gene clusters and counts sparse
#' clinical-validation hits per cluster and drug.
#' @param final_biomarkers Merged final biomarker table.
#' @param gene_clusters Output from \code{assignSharedSignatureClusters()}.
#' @return A \code{data.table} summary.
#' @export
buildFinalBiomarkerClusterOverlay <- function(final_biomarkers, gene_clusters) {
  final_dt <- data.table::as.data.table(final_biomarkers)
  clusters <- data.table::as.data.table(gene_clusters)
  if (!all(c("name", "drug", "tumor_type") %in% names(final_dt))) {
    stop("final_biomarkers must contain columns: name, drug, tumor_type", call. = FALSE)
  }
  if (!all(c("name", "gene_cluster") %in% names(clusters))) {
    stop("gene_clusters must contain columns: name, gene_cluster", call. = FALSE)
  }

  breast_keys <- c("breast cancer", "breast_cancer")
  joined <- merge(final_dt[tumor_type %in% breast_keys], clusters, by = "name", all.x = FALSE, all.y = FALSE)
  if (!nrow(joined)) {
    return(data.table::data.table())
  }
  joined[, .(n_final_biomarkers = data.table::uniqueN(name)), by = .(gene_cluster, drug)][
    order(gene_cluster, -n_final_biomarkers, drug)
  ]
}

#' Plot Final Biomarker Overlay by Cluster
#'
#' @description Bar plot of final biomarker counts assigned to shared signature
#' clusters.
#' @param overlay_dt Output from \code{buildFinalBiomarkerClusterOverlay()}.
#' @return A ggplot object.
#' @export
plotFinalBiomarkerClusterOverlay <- function(overlay_dt) {
  overlay_dt <- data.table::as.data.table(overlay_dt)
  required <- c("gene_cluster", "drug", "n_final_biomarkers")
  if (!all(required %in% names(overlay_dt))) {
    stop("overlay_dt must contain columns: ", paste(required, collapse = ", "), call. = FALSE)
  }
  if (!nrow(overlay_dt)) {
    stop("overlay_dt is empty", call. = FALSE)
  }

  ggplot2::ggplot(
    overlay_dt,
    ggplot2::aes(x = gene_cluster, y = n_final_biomarkers, fill = drug)
  ) +
    ggplot2::geom_col(width = 0.7, position = "stack") +
    ggplot2::labs(
      title = "Final biomarkers mapped onto shared AD-filtered clusters",
      x = NULL,
      y = "Final biomarkers"
    ) +
    themeMetaPaper(base_size = 10)
}
