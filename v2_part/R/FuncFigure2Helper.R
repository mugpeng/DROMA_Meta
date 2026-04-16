suppressPackageStartupMessages(library(data.table))

getFigure2MetaBatchDir <- function() {
  normalizePath(
    "/Users/peng/Desktop/Project/DROMA/Meta_project/Meta_Example/Output/meta_batch",
    mustWork = TRUE
  )
}

getFigure2Root <- function() {
  normalizePath(
    "/Users/peng/Desktop/Project/DROMA/Meta_project/DROMA_Meta/v2_part",
    mustWork = TRUE
  )
}

getFigure2OutputDir <- function() {
  output_dir <- Sys.getenv("DROMA_FIGURE2_OUTPUT_DIR", unset = "")
  if (!nzchar(output_dir)) {
    output_dir <- "/Users/peng/Desktop/Project/DROMA/Meta_project/Meta_Example/Output/visualization/figure2"
  }
  normalizePath(output_dir, mustWork = FALSE)
}

getFigure2Paths <- function() {
  root <- getFigure2Root()
  data_dir <- file.path(root, "data")
  output_dir <- getFigure2OutputDir()
  tables_dir <- file.path(output_dir, "tables")
  figures_dir <- output_dir
  workflow_dir <- file.path(root, "workflow")
  r_dir <- file.path(root, "R")

  dirs <- c(data_dir, output_dir, tables_dir, figures_dir, workflow_dir, r_dir)
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

  output_dir <- normalizePath(output_dir, mustWork = TRUE)
  tables_dir <- normalizePath(tables_dir, mustWork = TRUE)
  figures_dir <- normalizePath(figures_dir, mustWork = TRUE)
  workflow_dir <- normalizePath(workflow_dir, mustWork = TRUE)
  r_dir <- normalizePath(r_dir, mustWork = TRUE)

  list(
    root = root,
    data_dir = data_dir,
    output_dir = output_dir,
    tables_dir = tables_dir,
    figures_dir = figures_dir,
    workflow_dir = workflow_dir,
    r_dir = r_dir,
    input_features_csv = file.path(tables_dir, "figure2_input_features.csv"),
    stage_summary_csv = file.path(tables_dir, "figure2_stage_summary.csv"),
    drug_target_csv = file.path(data_dir, "drug_target_map.csv"),
    reactome_gmt_zip = file.path(data_dir, "ReactomePathways.gmt.zip"),
    reactome_gmt = file.path(data_dir, "ReactomePathways.gmt"),
    reactome_edges_csv = file.path(tables_dir, "reactome_edges.csv"),
    reactome_nodes_csv = file.path(tables_dir, "reactome_nodes.csv"),
    distance_results_csv = file.path(tables_dir, "distance_results.csv"),
    background_distance_results_csv = file.path(tables_dir, "background_distance_results.csv"),
    distance_summary_csv = file.path(tables_dir, "distance_summary_by_stage.csv"),
    fig2a_pdf = file.path(figures_dir, "figure2A_distance_distribution.pdf"),
    fig2b_pdf = file.path(figures_dir, "figure2B_distance_vs_effect.pdf"),
    fig2c_pdf = file.path(figures_dir, "figure2C_distance_vs_validation_stage.pdf")
  )
}

chooseFigure2EffectSize <- function(dt, stage) {
  dt <- as.data.table(copy(dt))

  candidates <- switch(stage,
    selected_genes = c("effect_size_pdcpdx", "effect_size_cell"),
    selected_genes_ad_filtered = c("effect_size_pdcpdx", "effect_size_cell"),
    clinical_sig_mRNA = c("effect_size"),
    final_biomarkers = c("effect_size_ctrdb", "effect_size_pdcpdx_invitro", "effect_size_cell_invitro"),
    stop("Unsupported stage: ", stage, call. = FALSE)
  )

  chosen <- rep(NA_real_, nrow(dt))
  source_col <- rep(NA_character_, nrow(dt))
  for (col in candidates) {
    if (!col %in% names(dt)) next
    idx <- is.na(chosen) & !is.na(dt[[col]])
    chosen[idx] <- dt[[col]][idx]
    source_col[idx] <- col
  }

  dt[, effect_size := chosen]
  dt[, effect_source := source_col]
  dt[]
}

standardizeFigure2Stage <- function(stage_path, stage) {
  dt <- fread(stage_path)
  dt <- chooseFigure2EffectSize(dt, stage)

  tumor_dir <- basename(dirname(stage_path))
  drug_dir <- basename(dirname(dirname(stage_path)))

  if (!"name" %in% names(dt)) {
    stop("Stage file missing name column: ", stage_path, call. = FALSE)
  }

  dt[, drug := drug_dir]
  dt[, tumor_type := gsub("_", " ", tumor_dir, fixed = TRUE)]
  dt[, stage := stage]
  dt[, source_file := basename(stage_path)]

  keep <- intersect(
    c(
      "drug", "tumor_type", "stage", "name", "effect_size", "effect_source",
      "direction", "direction_cell", "direction_pdcpdx", "direction_cell_invitro",
      "direction_pdcpdx_invitro", "ctrdb_status", "ctrdb_fallback",
      "p_value_ctrdb", "q_value_ctrdb", "effect_size_ctrdb",
      "effect_size_pdcpdx", "effect_size_cell",
      "effect_size_pdcpdx_invitro", "effect_size_cell_invitro",
      "source_file"
    ),
    names(dt)
  )
  unique(dt[, ..keep])
}

collectFigure2Inputs <- function(meta_batch_dir = getFigure2MetaBatchDir()) {
  patterns <- list(
    selected_genes_ad_filtered = "^selected_genes_ad_filtered\\.csv$",
    clinical_sig_mRNA = "^clinical_sig_mRNA\\.csv$",
    final_biomarkers = "^final_biomarkers\\.csv$"
  )

  rows <- list()
  for (stage in names(patterns)) {
    files <- list.files(
      meta_batch_dir,
      pattern = patterns[[stage]],
      recursive = TRUE,
      full.names = TRUE
    )
    if (!length(files)) next
    rows[[stage]] <- rbindlist(lapply(files, standardizeFigure2Stage, stage = stage), fill = TRUE)
  }

  if (!length(rows)) {
    stop("No figure2 input files found under: ", meta_batch_dir, call. = FALSE)
  }

  rbindlist(rows, fill = TRUE)
}

labelFigure2Stage <- function(stage) {
  fifelse(stage == "selected_genes_ad_filtered", "TCGA-AD filtered",
    fifelse(stage == "clinical_sig_mRNA", "Clinical significant",
      fifelse(stage == "final_biomarkers", "Clinical validated", stage)
    )
  )
}

orderFigure2StageLabels <- function(labels) {
  factor(
    labels,
    levels = c("TCGA-AD filtered", "Clinical significant", "Clinical validated")
  )
}

readGmtLines <- function(gmt_path) {
  lines <- readLines(gmt_path, warn = FALSE)
  lines[nzchar(lines)]
}

parseReactomeGmt <- function(gmt_path) {
  lines <- readGmtLines(gmt_path)
  pathway_rows <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 4) {
      return(NULL)
    }
    pathway_name <- fields[[1]]
    genes <- unique(fields[-c(1, 2)])
    genes <- genes[!is.na(genes) & nzchar(genes) & genes != "NA"]
    if (!length(genes)) {
      return(NULL)
    }
    data.table(pathway = pathway_name, gene = genes)
  })
  rbindlist(pathway_rows, fill = TRUE)
}

buildReactomeEdges <- function(pathway_gene_dt) {
  pathway_gene_dt <- as.data.table(pathway_gene_dt)
  pathway_gene_dt <- unique(pathway_gene_dt[, .(pathway, gene)])
  edge_rows <- lapply(split(pathway_gene_dt$gene, pathway_gene_dt$pathway), function(genes) {
    genes <- sort(unique(genes[!is.na(genes) & nzchar(genes) & genes != "NA"]))
    if (length(genes) < 2) {
      return(NULL)
    }
    combos <- t(utils::combn(genes, 2))
    data.table(from = combos[, 1], to = combos[, 2])
  })

  edge_dt <- unique(rbindlist(edge_rows, fill = TRUE))
  if (!nrow(edge_dt)) {
    return(data.table(from = character(0), to = character(0)))
  }
  edge_dt <- edge_dt[
    !is.na(from) & !is.na(to) & nzchar(from) & nzchar(to) &
      from != "NA" & to != "NA"
  ]
  edge_dt[order(from, to)]
}

buildReactomeGraph <- function(edge_dt) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' is required for figure2 distance calculations", call. = FALSE)
  }
  edge_dt <- as.data.table(copy(edge_dt))
  edge_dt <- edge_dt[
    !is.na(from) & !is.na(to) & nzchar(from) & nzchar(to) &
      !from %in% c("NA", "N/A") & !to %in% c("NA", "N/A")
  ]
  igraph::graph_from_data_frame(as.data.frame(edge_dt), directed = FALSE)
}

computeMinTargetDistance <- function(gene, target_genes, edge_dt = NULL, graph = NULL) {
  if (is.null(graph)) {
    if (is.null(edge_dt)) {
      stop("Provide edge_dt or graph", call. = FALSE)
    }
    graph <- buildReactomeGraph(edge_dt)
  }

  target_genes <- unique(target_genes[!is.na(target_genes) & nzchar(target_genes)])
  vertex_names <- igraph::V(graph)$name

  has_gene_mapping <- gene %in% vertex_names
  mapped_targets <- intersect(target_genes, vertex_names)
  has_target_mapping <- length(mapped_targets) > 0

  if (!has_gene_mapping || !has_target_mapping) {
    return(list(
      distance = NA_integer_,
      target_used = NA_character_,
      has_gene_mapping = has_gene_mapping,
      has_target_mapping = has_target_mapping
    ))
  }

  dists <- igraph::distances(graph, v = gene, to = mapped_targets, mode = "all")
  dists <- as.numeric(dists[1, ])
  if (!length(dists) || all(!is.finite(dists))) {
    return(list(
      distance = NA_integer_,
      target_used = NA_character_,
      has_gene_mapping = has_gene_mapping,
      has_target_mapping = has_target_mapping
    ))
  }

  best_idx <- which.min(dists)
  list(
    distance = as.integer(dists[[best_idx]]),
    target_used = mapped_targets[[best_idx]],
    has_gene_mapping = has_gene_mapping,
    has_target_mapping = has_target_mapping
  )
}

computeFigure2Distances <- function(feature_dt, target_dt, edge_dt) {
  graph <- buildReactomeGraph(edge_dt)
  target_split <- split(target_dt$target_gene, target_dt$drug)
  vertex_names <- igraph::V(graph)$name

  per_drug <- lapply(split(feature_dt, feature_dt$drug), function(dt) {
    target_genes <- unique(target_split[[dt$drug[[1]]]])
    mapped_targets <- intersect(target_genes, vertex_names)
    has_target_mapping <- length(mapped_targets) > 0
    genes <- unique(dt$name)
    has_gene_mapping <- genes %in% vertex_names

    gene_map <- data.table(name = genes, has_gene_mapping = has_gene_mapping)
    if (has_target_mapping) {
      mapped_genes <- genes[has_gene_mapping]
      if (length(mapped_genes)) {
        dist_mat <- igraph::distances(graph, v = mapped_genes, to = mapped_targets, mode = "all")
        dist_dt <- data.table(
          name = mapped_genes,
          distance = apply(dist_mat, 1, function(x) {
            if (all(!is.finite(x))) NA_integer_ else as.integer(min(x[is.finite(x)]))
          }),
          target_used = apply(dist_mat, 1, function(x) {
            if (all(!is.finite(x))) NA_character_ else mapped_targets[[which.min(x)]]
          })
        )
      } else {
        dist_dt <- data.table(name = character(0), distance = integer(0), target_used = character(0))
      }
    } else {
      dist_dt <- data.table(name = genes, distance = NA_integer_, target_used = NA_character_)
    }

    merged <- merge(gene_map, dist_dt, by = "name", all.x = TRUE)
    merged[, has_target_mapping := has_target_mapping]
    merge(dt, merged, by = "name", all.x = TRUE)
  })

  rbindlist(per_drug, fill = TRUE)
}

computePairwiseGraphDistances <- function(from_genes, to_genes, graph) {
  stopifnot(length(from_genes) == length(to_genes))

  vertex_names <- igraph::V(graph)$name
  from_idx <- match(from_genes, vertex_names)
  to_idx <- match(to_genes, vertex_names)

  out <- data.table(
    distance = rep(NA_integer_, length(from_genes)),
    has_gene_mapping = !is.na(from_idx),
    has_target_mapping = !is.na(to_idx)
  )

  valid <- which(!is.na(from_idx) & !is.na(to_idx))
  if (!length(valid)) {
    return(out)
  }

  unique_from <- unique(from_genes[valid])
  unique_to <- unique(to_genes[valid])
  dist_mat <- igraph::distances(graph, v = unique_from, to = unique_to, mode = "all")

  row_idx <- match(from_genes[valid], unique_from)
  col_idx <- match(to_genes[valid], unique_to)
  dists <- dist_mat[cbind(row_idx, col_idx)]
  out$distance[valid] <- ifelse(is.finite(dists), as.integer(dists), NA_integer_)
  out
}

computeMinDistancesToTargets <- function(genes, mapped_targets, graph) {
  genes <- as.character(genes)
  out <- data.table(
    random_gene = genes,
    distance = NA_integer_,
    has_gene_mapping = FALSE,
    has_target_mapping = length(mapped_targets) > 0
  )

  if (!length(genes)) {
    return(out)
  }

  vertex_names <- igraph::V(graph)$name
  out[, has_gene_mapping := random_gene %in% vertex_names]
  valid_genes <- unique(out[has_gene_mapping == TRUE, random_gene])

  if (!length(valid_genes) || !length(mapped_targets)) {
    return(out)
  }

  dist_mat <- igraph::distances(graph, v = valid_genes, to = mapped_targets, mode = "all")
  min_dist <- apply(dist_mat, 1, function(x) {
    finite_x <- x[is.finite(x)]
    if (!length(finite_x)) NA_integer_ else as.integer(min(finite_x))
  })

  dist_dt <- data.table(random_gene = valid_genes, distance = as.integer(min_dist))
  out[dist_dt, distance := i.distance, on = "random_gene"]
  out
}

computeFigure2Background <- function(feature_dt, target_dt, edge_dt, seed = 20260416L) {
  set.seed(seed)
  graph <- buildReactomeGraph(edge_dt)
  graph_genes <- igraph::V(graph)$name
  target_split <- split(target_dt$target_gene, target_dt$drug)
  bg_rows <- lapply(split(feature_dt, paste(feature_dt$drug, feature_dt$stage, sep = "||")), function(dt) {
    drug_value <- dt$drug[[1]]
    stage_value <- dt$stage[[1]]
    target_genes <- unique(target_split[[drug_value]])
    mapped_targets <- intersect(target_genes, graph_genes)
    sample_n <- nrow(dt)

    sampled_genes <- sample(graph_genes, size = sample_n, replace = sample_n > length(graph_genes))
    if (length(mapped_targets)) {
      gene_target_dt <- computeMinDistancesToTargets(sampled_genes, mapped_targets, graph)
      gene_target_dt[, `:=`(
        drug = drug_value,
        stage = stage_value,
        background_type = "random_gene_to_target"
      )]
      data.table::setcolorder(
        gene_target_dt,
        c("drug", "stage", "background_type", "random_gene", "distance", "has_gene_mapping", "has_target_mapping")
      )
    } else {
      gene_target_dt <- data.table(
        drug = drug_value,
        stage = stage_value,
        background_type = "random_gene_to_target",
        random_gene = sampled_genes,
        distance = NA_integer_,
        has_gene_mapping = TRUE,
        has_target_mapping = FALSE
      )
    }

    random_targets <- sample(graph_genes, size = sample_n, replace = sample_n > length(graph_genes))
    pair_dist_dt <- computePairwiseGraphDistances(sampled_genes, random_targets, graph)
    pair_dt <- cbind(
      data.table(
        drug = drug_value,
        stage = stage_value,
        background_type = "random_gene_to_random_gene",
        random_gene = sampled_genes,
        random_target = random_targets
      ),
      pair_dist_dt
    )
    data.table::setcolorder(
      pair_dt,
      c("drug", "stage", "background_type", "random_gene", "random_target", "distance", "has_gene_mapping", "has_target_mapping")
    )

    rbindlist(list(gene_target_dt, pair_dt), fill = TRUE)
  })

  rbindlist(bg_rows, fill = TRUE)
}

summarizeFigure2Distances <- function(distance_dt) {
  distance_dt[
    ,
    .(
      n = .N,
      n_finite_distance = sum(!is.na(distance)),
      median_distance = as.numeric(suppressWarnings(median(distance, na.rm = TRUE))),
      mean_distance = as.numeric(suppressWarnings(mean(distance, na.rm = TRUE))),
      mapped_gene_fraction = mean(has_gene_mapping, na.rm = TRUE),
      mapped_target_fraction = mean(has_target_mapping, na.rm = TRUE)
    ),
    by = .(stage)
  ]
}

downloadFigure2File <- function(url, destfile) {
  utils::download.file(url, destfile = destfile, mode = "wb", quiet = FALSE)
  invisible(destfile)
}

maybeUnzipReactomeGmt <- function(zip_path, out_path) {
  if (file.exists(out_path)) {
    return(out_path)
  }
  utils::unzip(zip_path, exdir = dirname(out_path))
  if (!file.exists(out_path)) {
    stop("Expected Reactome GMT after unzip: ", out_path, call. = FALSE)
  }
  out_path
}

fetchChemblTargets <- function(drugs) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required for ChEMBL target fetch", call. = FALSE)
  }
  old_timeout <- getOption("timeout")
  options(timeout = max(60, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  fetch_one <- function(drug_name) {
    encoded <- utils::URLencode(drug_name, reserved = TRUE)
    mol_url <- sprintf(
      "https://www.ebi.ac.uk/chembl/api/data/molecule/search.json?q=%s",
      encoded
    )
    mol_json <- jsonlite::fromJSON(mol_url)
    mols <- mol_json$molecules
    if (is.null(mols) || !nrow(mols)) {
      return(NULL)
    }
    mols <- as.data.table(mols)
    if ("score" %in% names(mols)) {
      data.table::setorderv(mols, cols = "score", order = -1)
    }
    mols <- head(mols, 3)

    extract_target_rows <- function(chembl_id) {
      mech_url <- sprintf(
        "https://www.ebi.ac.uk/chembl/api/data/mechanism.json?molecule_chembl_id=%s&limit=1000",
        chembl_id
      )
      mech_json <- jsonlite::fromJSON(mech_url)
      mechs <- mech_json$mechanisms
      if (!is.data.frame(mechs) || !nrow(mechs)) {
        return(NULL)
      }

      target_rows <- lapply(seq_len(nrow(mechs)), function(i) {
        target_id <- mechs$target_chembl_id[[i]]
        if (is.null(target_id) || is.na(target_id) || !nzchar(target_id)) {
          return(NULL)
        }
        target_url <- sprintf(
          "https://www.ebi.ac.uk/chembl/api/data/target/%s.json",
          target_id
        )
        target_json <- jsonlite::fromJSON(target_url)
        comps <- target_json$target_components
        if (!is.data.frame(comps) || !nrow(comps)) {
          return(NULL)
        }

        comp_rows <- lapply(seq_len(nrow(comps)), function(j) {
          synonyms <- comps$target_component_synonyms[[j]]
          gene_symbol <- NA_character_
          if (is.data.frame(synonyms) && nrow(synonyms) > 0 && "component_synonym" %in% names(synonyms)) {
            symbol_rows <- synonyms[synonyms$syn_type == "GENE_SYMBOL", , drop = FALSE]
            if (nrow(symbol_rows) > 0) {
              gene_symbol <- symbol_rows$component_synonym[[1]]
            } else {
              gene_symbol <- synonyms$component_synonym[[1]]
            }
          }
          data.table(
            drug = drug_name,
            target_gene = gene_symbol,
            target_source = "ChEMBL",
            evidence_level = "mechanism",
            molecule_chembl_id = chembl_id,
            target_chembl_id = target_id,
            mechanism_of_action = if ("mechanism_of_action" %in% names(mechs)) mechs$mechanism_of_action[[i]] else NA_character_
          )
        })

        rbindlist(comp_rows, fill = TRUE)
      })

      rows <- rbindlist(target_rows, fill = TRUE)
      if (!nrow(rows)) {
        return(NULL)
      }
      rows[!is.na(target_gene) & nzchar(target_gene)]
    }

    chembl_ids <- unique(mols$molecule_chembl_id[!is.na(mols$molecule_chembl_id) & nzchar(mols$molecule_chembl_id)])
    collected <- rbindlist(lapply(chembl_ids, extract_target_rows), fill = TRUE)
    if (!nrow(collected)) {
      return(NULL)
    }
    unique(collected)
  }

  rbindlist(lapply(unique(drugs), function(drug_name) {
    cat(sprintf("  Fetching targets for %s\n", drug_name))
    tryCatch(
      fetch_one(drug_name),
      error = function(e) {
        message("    target fetch failed for ", drug_name, ": ", conditionMessage(e))
        NULL
      }
    )
  }), fill = TRUE)
}
