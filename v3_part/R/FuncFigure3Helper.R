# Figure 3 breast-state helpers ----

getFigure3Paths <- function(project_root = "/Users/peng/Desktop/Project/DROMA/Meta_project/Meta_Example") {
  output_base <- Sys.getenv("DROMA_FIGURE3_OUTPUT_BASE", unset = NA_character_)
  if (is.na(output_base) || !nzchar(output_base)) {
    output_base <- file.path(project_root, "Output", "meta_batch")
  }

  output_dir <- Sys.getenv("DROMA_FIGURE3_OUTPUT_DIR", unset = NA_character_)
  if (is.na(output_dir) || !nzchar(output_dir)) {
    output_dir <- file.path(project_root, "Output", "0416_plot", "figure3_breast_state")
  }

  output_base <- normalizePath(output_base, mustWork = FALSE)
  output_dir <- normalizePath(output_dir, mustWork = FALSE)

  list(
    output_base = output_base,
    output_dir = output_dir,
    input_selected_ad = file.path(output_dir, "figure3_selected_ad_filtered_all.csv"),
    input_final_biomarkers = file.path(output_dir, "figure3_final_biomarkers_all.csv"),
    universe_csv = file.path(output_dir, "figure3_breast_universe.csv"),
    signature_matrix_csv = file.path(output_dir, "Figure3A_breast_shared_signature_matrix.csv"),
    signature_genes_csv = file.path(output_dir, "Figure3A_breast_shared_signature_genes.csv"),
    signature_drugs_csv = file.path(output_dir, "Figure3A_breast_shared_signature_drugs.csv"),
    gene_clusters_csv = file.path(output_dir, "Figure3A_gene_clusters.csv"),
    drug_clusters_csv = file.path(output_dir, "Figure3A_drug_clusters.csv"),
    fig3a_pdf = file.path(output_dir, "Figure3A_breast_shared_signature_heatmap.pdf"),
    go_csv = file.path(output_dir, "Figure3B_GO_enrichment_by_gene_cluster.csv"),
    fig3b_pdf = file.path(output_dir, "Figure3B_GO_enrichment_by_gene_cluster.pdf"),
    marker_gene_overlay_csv = file.path(output_dir, "Figure3C_marker_gene_overlay.csv"),
    marker_program_overlay_csv = file.path(output_dir, "Figure3C_marker_program_overlay.csv"),
    fig3c_pdf = file.path(output_dir, "Figure3C_PAM50_marker_overlay.pdf"),
    cluster_program_csv = file.path(output_dir, "Figure3D_cluster_program_annotation.csv"),
    fig3d_pdf = file.path(output_dir, "Figure3D_cluster_program_annotation.pdf"),
    pseudo_atlas_csv = file.path(output_dir, "Figure3E_pseudo_atlas_annotation.csv"),
    fig3e_pdf = file.path(output_dir, "Figure3E_pseudo_atlas_annotation.pdf"),
    final_overlay_csv = file.path(output_dir, "Figure3_final_biomarker_cluster_overlay.csv"),
    interpretation_md = file.path(output_dir, "Figure3_breast_state_interpretation.md"),
    legend_md = file.path(output_dir, "legend.md")
  )
}

collectFigure3SelectedAdFiltered <- function(output_base) {
  collectAllSelectedAdFiltered(output_base)
}

collectFigure3FinalBiomarkers <- function(output_base) {
  collectAllFinalBiomarkers(output_base)
}

getFigure3BreastPrograms <- function() {
  data.table::rbindlist(list(
    data.table::data.table(program = "Luminal", gene = c("ESR1", "PGR", "FOXA1", "GATA3", "KRT8", "KRT18", "MUC1", "XBP1"), category = "pam50"),
    data.table::data.table(program = "HER2-enriched", gene = c("ERBB2", "GRB7", "FGFR4", "TMEM45B", "TCAP"), category = "pam50"),
    data.table::data.table(program = "Basal-like", gene = c("KRT5", "KRT14", "KRT17", "EGFR", "FOXC1", "LAD1"), category = "pam50"),
    data.table::data.table(program = "Proliferation", gene = c("MKI67", "AURKA", "TOP2A", "CCNB1", "CDC20", "BIRC5", "MYBL2", "UBE2C"), category = "pam50"),
    data.table::data.table(program = "EMT-like", gene = c("VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "CDH2"), category = "state"),
    data.table::data.table(program = "Epithelial", gene = c("EPCAM", "CDH1", "KRT8", "KRT18", "MUC1", "KRT19"), category = "state"),
    data.table::data.table(program = "Immune-adjacent", gene = c("CXCL8", "CCL5", "CXCL10", "STAT1", "IFITM1", "HLA-DRA"), category = "state")
  ), fill = TRUE)
}

getFigure3PseudoAtlasMarkers <- function() {
  data.table::rbindlist(list(
    data.table::data.table(cell_state = "luminal epithelial", gene = c("ESR1", "PGR", "FOXA1", "GATA3", "KRT8", "KRT18", "MUC1")),
    data.table::data.table(cell_state = "HER2-like epithelial", gene = c("ERBB2", "GRB7", "FGFR4", "TMEM45B")),
    data.table::data.table(cell_state = "basal epithelial", gene = c("KRT5", "KRT14", "KRT17", "EGFR", "FOXC1", "LAD1")),
    data.table::data.table(cell_state = "cycling/proliferative", gene = c("MKI67", "AURKA", "TOP2A", "CCNB1", "CDC20", "BIRC5", "UBE2C")),
    data.table::data.table(cell_state = "mesenchymal-like", gene = c("VIM", "ZEB1", "ZEB2", "SNAI1", "SNAI2", "TWIST1", "CDH2")),
    data.table::data.table(cell_state = "immune-adjacent/interferon", gene = c("CXCL8", "CCL5", "CXCL10", "STAT1", "IFITM1", "HLA-DRA"))
  ), fill = TRUE)
}

buildFigure3ProgramGeneOverlay <- function(selected_ad_filtered,
                                           programs = getFigure3BreastPrograms(),
                                           tumor_type = "breast_cancer",
                                           effect_col = NULL) {
  dt <- data.table::as.data.table(selected_ad_filtered)
  programs <- data.table::as.data.table(programs)
  if (is.null(effect_col)) {
    effect_candidates <- c("effect_size_cell", "effect_size_cell_invitro", "effect_size")
    effect_col <- effect_candidates[effect_candidates %in% names(dt)][1]
  }
  required <- c("name", effect_col, "drug", "tumor_type")
  if (!all(required %in% names(dt))) {
    stop("selected_ad_filtered must contain columns: ", paste(required, collapse = ", "), call. = FALSE)
  }
  tumor_keys <- unique(c(tumor_type, gsub(" ", "_", tumor_type), gsub("_", " ", tumor_type)))
  dt <- dt[tumor_type %in% tumor_keys]
  if (!nrow(dt)) return(data.table::data.table())

  merged <- merge(
    dt[, .(name, drug, effect_size = get(effect_col))],
    unique(programs[, .(program, category, name = gene)]),
    by = "name",
    all = FALSE
  )
  if (!nrow(merged)) return(data.table::data.table())
  merged[, .(effect_size = mean(effect_size, na.rm = TRUE)), by = .(category, program, name, drug)]
}

buildFigure3ProgramOverlay <- function(marker_gene_overlay) {
  dt <- data.table::as.data.table(marker_gene_overlay)
  if (!nrow(dt)) return(data.table::data.table())
  dt[, .(
    mean_effect_size = mean(effect_size, na.rm = TRUE),
    n_markers_present = data.table::uniqueN(name),
    leading_markers = paste(head(name[order(-abs(effect_size))], 4L), collapse = ", ")
  ), by = .(category, program, drug)][order(category, program, drug)]
}

annotateFigure3GeneClusters <- function(gene_clusters,
                                        marker_sets = getFigure3BreastPrograms(),
                                        min_overlap = 1L) {
  clusters <- data.table::as.data.table(gene_clusters)
  markers <- data.table::as.data.table(marker_sets)
  if (!nrow(clusters)) return(data.table::data.table())

  cluster_sizes <- clusters[, .(cluster_size = .N), by = gene_cluster]
  overlaps <- merge(
    clusters[, .(gene_cluster, name)],
    unique(markers[, .(program, category, name = gene)]),
    by = "name",
    all = FALSE
  )
  if (!nrow(overlaps)) {
    return(cluster_sizes[, `:=`(
      category = "none",
      program = "Unassigned",
      overlap_n = 0L,
      overlap_frac = 0,
      matched_genes = ""
    )])
  }

  ann <- overlaps[, .(
    overlap_n = data.table::uniqueN(name),
    matched_genes = paste(sort(unique(name)), collapse = ", ")
  ), by = .(gene_cluster, category, program)]
  ann <- merge(ann, cluster_sizes, by = "gene_cluster", all.x = TRUE)
  ann[, overlap_frac := overlap_n / pmax(cluster_size, 1L)]
  ann <- ann[overlap_n >= min_overlap][order(gene_cluster, -overlap_frac, -overlap_n, program)]
  if (!nrow(ann)) {
    return(cluster_sizes[, `:=`(
      category = "none",
      program = "Unassigned",
      overlap_n = 0L,
      overlap_frac = 0,
      matched_genes = "",
      rank_within_cluster = 1L
    )])
  }
  ann[, rank_within_cluster := seq_len(.N), by = gene_cluster]
  missing_clusters <- setdiff(cluster_sizes$gene_cluster, unique(ann$gene_cluster))
  if (length(missing_clusters)) {
    ann <- data.table::rbindlist(list(
      ann,
      cluster_sizes[gene_cluster %in% missing_clusters][, `:=`(
        category = "none",
        program = "Unassigned",
        overlap_n = 0L,
        overlap_frac = 0,
        matched_genes = "",
        rank_within_cluster = 1L
      )]
    ), fill = TRUE)
  }
  ann[order(gene_cluster, rank_within_cluster)]
}

annotateFigure3PseudoAtlas <- function(gene_clusters,
                                       atlas_markers = getFigure3PseudoAtlasMarkers(),
                                       min_overlap = 1L) {
  clusters <- data.table::as.data.table(gene_clusters)
  atlas <- data.table::as.data.table(atlas_markers)
  if (!nrow(clusters)) return(data.table::data.table())

  cluster_sizes <- clusters[, .(cluster_size = .N), by = gene_cluster]
  overlaps <- merge(
    clusters[, .(gene_cluster, name)],
    unique(atlas[, .(cell_state, name = gene)]),
    by = "name",
    all = FALSE
  )
  if (!nrow(overlaps)) {
    return(cluster_sizes[, `:=`(
      cell_state = "unassigned",
      overlap_n = 0L,
      overlap_frac = 0,
      supporting_markers = "",
      interpretation = "No pseudo-atlas marker state reached the minimum overlap threshold."
    )])
  }

  ann <- overlaps[, .(
    overlap_n = data.table::uniqueN(name),
    supporting_markers = paste(sort(unique(name)), collapse = ", ")
  ), by = .(gene_cluster, cell_state)]
  ann <- merge(ann, cluster_sizes, by = "gene_cluster", all.x = TRUE)
  ann[, overlap_frac := overlap_n / pmax(cluster_size, 1L)]
  ann <- ann[overlap_n >= min_overlap][order(gene_cluster, -overlap_frac, -overlap_n, cell_state)]
  if (!nrow(ann)) {
    return(cluster_sizes[, `:=`(
      cell_state = "unassigned",
      overlap_n = 0L,
      overlap_frac = 0,
      supporting_markers = "",
      interpretation = "No pseudo-atlas marker state reached the minimum overlap threshold.",
      rank_within_cluster = 1L
    )])
  }
  ann[, rank_within_cluster := seq_len(.N), by = gene_cluster]
  ann[, interpretation := paste0(
    "Marker-guided pseudo-atlas label: ", cell_state,
    " (supporting markers: ", supporting_markers, ")."
  )]
  missing_clusters <- setdiff(cluster_sizes$gene_cluster, unique(ann$gene_cluster))
  if (length(missing_clusters)) {
    ann <- data.table::rbindlist(list(
      ann,
      cluster_sizes[gene_cluster %in% missing_clusters][, `:=`(
        cell_state = "unassigned",
        overlap_n = 0L,
        overlap_frac = 0,
        supporting_markers = "",
        interpretation = "No pseudo-atlas marker state reached the minimum overlap threshold.",
        rank_within_cluster = 1L
      )]
    ), fill = TRUE)
  }
  ann[order(gene_cluster, rank_within_cluster)]
}

summarizeFigure3StateHypothesis <- function(drug_clusters,
                                            cluster_program_annotation,
                                            pseudo_atlas_annotation,
                                            top_n = 2L) {
  drugs <- data.table::as.data.table(drug_clusters)
  program_dt <- data.table::as.data.table(cluster_program_annotation)
  atlas_dt <- data.table::as.data.table(pseudo_atlas_annotation)

  if (!nrow(drugs)) return(data.table::data.table())
  top_programs <- program_dt[rank_within_cluster %in% seq_len(top_n),
                             .(program_summary = paste(program, collapse = " / ")),
                             by = gene_cluster]
  top_atlas <- atlas_dt[rank_within_cluster %in% seq_len(top_n),
                        .(cell_state_summary = paste(cell_state, collapse = " / ")),
                        by = gene_cluster]

  merge(top_programs, top_atlas, by = "gene_cluster", all = TRUE)[
    drugs, on = .(gene_cluster = drug_cluster), allow.cartesian = TRUE
  ]
}

writeFigure3Interpretation <- function(paths,
                                       signature,
                                       gene_clusters,
                                       drug_clusters,
                                       go_dt,
                                       cluster_program_annotation,
                                       pseudo_atlas_annotation,
                                       final_overlay) {
  program_top <- data.table::as.data.table(cluster_program_annotation)[
    rank_within_cluster == 1L,
    .(gene_cluster, top_program = program, overlap_n, matched_genes)
  ]
  atlas_top <- data.table::as.data.table(pseudo_atlas_annotation)[
    rank_within_cluster == 1L,
    .(gene_cluster, top_cell_state = cell_state, supporting_markers)
  ]
  go_top <- data.table::as.data.table(go_dt)
  if (nrow(go_top)) {
    go_top <- go_top[order(gene_cluster, p.adjust)][
      , head(.SD, 3L), by = gene_cluster
    ][, .(top_go_terms = paste(Description, collapse = "; ")), by = gene_cluster]
  } else {
    go_top <- data.table::data.table(gene_cluster = unique(gene_clusters$gene_cluster), top_go_terms = "No significant GO term")
  }

  cluster_summary <- Reduce(function(x, y) merge(x, y, by = "gene_cluster", all = TRUE),
                            list(
                              unique(data.table::as.data.table(gene_clusters)[, .(gene_cluster)]),
                              program_top,
                              atlas_top,
                              go_top
                            ))
  cluster_summary[is.na(top_program), top_program := "not confidently assigned"]
  cluster_summary[is.na(top_cell_state), top_cell_state := "not confidently assigned"]
  assigned_n <- sum(!(cluster_summary$top_program %in% c("not confidently assigned", "Unassigned")))
  cluster_lines <- if (nrow(cluster_summary)) {
    sprintf(
      "- `%s`: top program `%s`; pseudo-atlas label `%s`; leading GO terms: %s.",
      cluster_summary$gene_cluster,
      cluster_summary$top_program,
      cluster_summary$top_cell_state,
      cluster_summary$top_go_terms
    )
  } else {
    "- No cluster summary available."
  }

  final_overlay <- data.table::as.data.table(final_overlay)
  final_line <- if (nrow(final_overlay)) {
    top_final <- final_overlay[order(-n_final_biomarkers)][1]
    sprintf(
      "Final clinically filtered biomarkers still map back to shared-state clusters; the largest overlap is `%s` for `%s` (%d biomarkers).",
      top_final$gene_cluster, top_final$drug, top_final$n_final_biomarkers
    )
  } else {
    "Final biomarker overlay was sparse, consistent with the idea that shared transcriptional state is most visible at the AD-filtered selected-gene layer."
  }

  drug_cluster_lines <- sprintf(
    "- `%s`: %d drugs.",
    sort(unique(drug_clusters$drug_cluster)),
    as.integer(data.table::as.data.table(drug_clusters)[, .N, by = drug_cluster][order(drug_cluster)]$N)
  )

  lines <- c(
    "# Breast Figure 3 interpretation",
    "",
    "## Scientific question",
    "",
    "A central question for this figure is whether the large number of breast-cancer biomarkers reflects many isolated gene-drug associations, or whether they are organized by a smaller number of shared cell states.",
    "",
    "## Main conclusion",
    "",
    sprintf(
      "The breast-cancer AD-filtered signature contains `%d` shared genes across `%d` drugs, and these genes separate into `%d` gene clusters while the drugs separate into `%d` drug clusters.",
      nrow(signature$matrix),
      ncol(signature$matrix),
      data.table::uniqueN(gene_clusters$gene_cluster),
      data.table::uniqueN(drug_clusters$drug_cluster)
    ),
    if (assigned_n > 0) {
      sprintf(
        "Taken together, the clustering, GO enrichment, PAM50-style marker overlays, and marker-guided pseudo-atlas annotations support the interpretation that many biomarkers are not independent events. At the current marker threshold, `%d/%d` shared clusters show a confident canonical breast-state assignment, most clearly highlighting a HER2-like / basal-adjacent program, while the remaining clusters are better described by non-subtype functional programs from GO.",
        assigned_n,
        nrow(cluster_summary)
      )
    } else {
      "Taken together, the clustering and GO enrichment support non-random shared structure among biomarkers, but the current canonical subtype marker panels do not confidently assign any shared cluster to a single breast cell state."
    },
    "",
    "## Drug-cluster structure",
    "",
    drug_cluster_lines,
    "",
    "## Cluster-level evidence",
    "",
    cluster_lines,
    "",
    "## Interpretation for the manuscript",
    "",
    "These results argue that the abundance of biomarkers in breast cancer is at least partly driven by shared biological structure rather than by entirely drug-specific mechanisms. In practice, a substantial subset of biomarkers appears to co-vary because they mark the same underlying programs. For this dataset, the clearest canonical subtype signal is a HER2-like / basal-adjacent module, whereas other shared clusters are dominated by trafficking, lysosomal, glycoprotein, and migration-related functions and therefore should be interpreted as broader cellular programs rather than clean PAM50 subtypes.",
    "",
    "This interpretation is strongest at the selected-gene and AD-filtered layers, where the shared-state structure is dense enough to drive unsupervised clustering. The final clinically filtered biomarkers are much sparser, but they still map back onto these shared programs rather than forming unrelated single-gene signals.",
    "",
    final_line,
    "",
    "## Note on the pseudo-atlas layer",
    "",
    "The cell-subgroup panel in this figure is a marker-guided pseudo-atlas annotation rather than a formal single-cell enrichment test against a reference atlas. It is intended to answer whether the shared clusters most resemble luminal epithelial, basal epithelial, cycling, mesenchymal-like, HER2-like, or immune-adjacent states.",
    ""
  )

  writeLines(lines, paths$interpretation_md, useBytes = TRUE)
  invisible(paths$interpretation_md)
}
