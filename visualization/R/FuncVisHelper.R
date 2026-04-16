# Visualization helpers for DROMA.Meta ----

# Color / Theme Constants ----

#' Get Default DROMA.Meta Visualization Colors
#'
#' @description Returns the shared color palette used across all DROMA.Meta
#' visualization functions. Palettes use Tableau 10 colorblind-safe values.
#' @param palette One of "direction", "evidence", "heatmap", "stage", "pair".
#' @return A named character vector of hex colors.
#' @export
getMetaVisColors <- function(palette = c("direction", "evidence", "heatmap",
                                         "stage", "pair")) {
  palette <- match.arg(palette)
  switch(palette,
    # Tableau 10 colorblind-safe
    direction = c(
      "Up"   = "#D62728",
      "Down" = "#1F77B4",
      "NS"   = "#C7C7C7"
    ),
    evidence = c(
      "Cell line" = "#D6A834",
      "PDC/PDX"   = "#1F77B4",
      "TCGA-AD"   = "#2CA02C",
      "Clinical"  = "#D62728"
    ),
    heatmap = c(
      "low"  = "#2166AC",
      "mid"  = "#F7F7F7",
      "high" = "#B2182B"
    ),
    stage = c(
      "Cell batch"        = "#AEC7E8",
      "Cell significant"  = "#6BAED6",
      "PDC/PDX significant" = "#3182BD",
      "Intersection"      = "#FDAE6B",
      "TCGA-AD filtered"  = "#E6550D",
      "Clinical validated" = "#D62728"
    ),
    pair = c(
      "Paclitaxel"       = "#1F77B4",
      "5-Fluorouracil"   = "#FF7F0E",
      "4'-Epiadriamycin" = "#2CA02C",
      "Afatinib"         = "#D62728"
    )
  )
}

#' Get Publication-Quality ggplot2 Theme
#'
#' @description A \code{theme_bw()} base tuned for journal figures: black text,
#' no minor grid, clean margins. Matches Nature/Cell formatting conventions.
#' @param base_size Base font size in pt. Default 11.
#' @return A \code{theme} object.
#' @export
themeMetaPaper <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", size = base_size + 1,
                                               margin = ggplot2::margin(b = 4)),
      plot.subtitle    = ggplot2::element_text(size = base_size - 1, color = "#6F685E",
                                               margin = ggplot2::margin(b = 3)),
      axis.title       = ggplot2::element_text(size = base_size, colour = "black"),
      axis.text        = ggplot2::element_text(size = base_size - 1, colour = "black"),
      axis.text.x      = ggplot2::element_text(angle = 0, hjust = 0.5),
      axis.line        = ggplot2::element_line(linewidth = 0.5, colour = "black"),
      legend.title     = ggplot2::element_text(size = base_size - 1, face = "bold"),
      legend.text      = ggplot2::element_text(size = base_size - 2),
      legend.position  = "right",
      legend.background = ggplot2::element_rect(fill = NA, colour = NA),
      legend.key       = ggplot2::element_rect(fill = NA, colour = NA),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.3),
      panel.border     = ggplot2::element_rect(colour = "black", linewidth = 0.5),
      plot.margin      = ggplot2::margin(10, 10, 10, 10),
      strip.background = ggplot2::element_rect(fill = "grey95", colour = "black",
                                               linewidth = 0.4),
      strip.text       = ggplot2::element_text(size = base_size - 1, face = "bold"),
      text             = ggplot2::element_text(colour = "black", family = "sans")
    )
}

# Path Management ----

#' Get Default Paths for the Visualization Workflow
#'
#' @description Returns paths for visualization output, mirroring
#' \code{getMetaWorkflowDefaults()} from the parent package.
#' @param project_root Root directory of the project.
#' @return A named list with \code{project_root}, \code{output_base},
#'   \code{vis_output}, \code{droma_db_path}, and \code{ctrdb_path}.
#' @export
getVisWorkflowDefaults <- function(project_root = normalizePath(getwd(), mustWork = TRUE)) {
  output_base <- file.path(project_root, "Output", "meta_batch")
  list(
    project_root  = project_root,
    output_base   = output_base,
    vis_output    = file.path(project_root, "Output", "visualization"),
    droma_db_path = "/Users/peng/Desktop/Project/DROMA/Data/droma.sqlite",
    ctrdb_path    = "/Users/peng/Desktop/Project/DROMA/Data/ctrdb.sqlite"
  )
}

# Data Aggregation ----

#' Collect Workflow Results Across All Drug-Tumor Pairs
#'
#' @description Scans the output directory produced by \code{runMetaWorkflow()} and
#' builds a master summary \code{data.table} with one row per drug-tumor pair.
#' @param output_base Root output directory (the \code{output_base} used by the
#'   workflow, typically ending in \code{"meta_batch"}).
#' @return A \code{data.table} with columns \code{drug}, \code{tumor_type},
#'   \code{output_dir}, and all count columns from
#'   \code{clinical_validation_info.csv}.
#' @export
collectWorkflowResults <- function(output_base) {
  if (!dir.exists(output_base)) {
    stop("output_base directory does not exist: ", output_base, call. = FALSE)
  }

  drug_dirs <- list.dirs(output_base, recursive = FALSE, full.names = TRUE)
  if (!length(drug_dirs)) {
    stop("No drug subdirectories found under: ", output_base, call. = FALSE)
  }

  rows <- list()
  for (drug_dir in drug_dirs) {
    tumor_dirs <- list.dirs(drug_dir, recursive = FALSE, full.names = TRUE)
    for (pair_dir in tumor_dirs) {
      info_path <- file.path(pair_dir, "clinical_validation_info.csv")
      if (!file.exists(info_path)) next

      info <- data.table::fread(info_path)
      info$output_dir <- pair_dir

      # Recover drug/tumor_type from directory names as canonical source
      info$drug <- basename(drug_dir)
      info$tumor_type <- basename(pair_dir)

      rows[[length(rows) + 1L]] <- info
    }
  }

  if (!length(rows)) {
    stop("No clinical_validation_info.csv files found under: ", output_base, call. = FALSE)
  }

  data.table::rbindlist(rows, fill = TRUE)
}

#' Collect Final Biomarkers Across All Drug-Tumor Pairs
#'
#' @description Reads every \code{final_biomarkers.csv} under \code{output_base}
#' and returns a single merged \code{data.table}.
#' @param output_base Root output directory produced by the workflow.
#' @return A \code{data.table} of all final biomarkers.
#' @export
collectAllFinalBiomarkers <- function(output_base) {
  if (!dir.exists(output_base)) {
    stop("output_base directory does not exist: ", output_base, call. = FALSE)
  }

  all_files <- list.files(
    output_base,
    pattern = "^final_biomarkers\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (!length(all_files)) {
    warning("No final_biomarkers.csv files found under: ", output_base, call. = FALSE)
    return(data.table::data.table())
  }

  data.table::rbindlist(lapply(all_files, data.table::fread), fill = TRUE)
}

#' Collect Batch Results Across All Drug-Tumor Pairs
#'
#' @description Reads every \code{batch_cell_sets_mRNA.csv} or
#' \code{batch_pdcpdx_sets_mRNA.csv} and annotates with drug and tumor_type
#' extracted from the directory structure.
#' @param output_base Root output directory produced by the workflow.
#' @param type One of "cell" or "pdcpdx".
#' @return A \code{data.table} of all batch results.
#' @export
collectAllBatchResults <- function(output_base,
                                   type = c("cell", "pdcpdx")) {
  type <- match.arg(type)
  pattern <- switch(type,
    cell   = "^batch_cell_sets_mRNA\\.csv$",
    pdcpdx = "^batch_pdcpdx_sets_mRNA\\.csv$"
  )

  all_files <- list.files(
    output_base,
    pattern = pattern,
    recursive = TRUE,
    full.names = TRUE
  )

  if (!length(all_files)) {
    warning("No ", type, " batch files found under: ", output_base, call. = FALSE)
    return(data.table::data.table())
  }

  data.table::rbindlist(lapply(all_files, function(f) {
    dt <- data.table::fread(f)
    # Extract drug/tumor_type from path: .../meta_batch/<drug>/<tumor>/
    # Normalize then split; tumor_type is the immediate parent, drug is one above
    tumor_dir  <- basename(dirname(f))
    drug_dir   <- basename(dirname(dirname(f)))
    dt[, drug := drug_dir]
    dt[, tumor_type := tumor_dir]
    dt
  }), fill = TRUE)
}

#' Collect Selected Genes AD Stats Across All Drug-Tumor Pairs
#'
#' @description Reads every \code{selected_genes_ad_stats.csv} and merges them.
#' @param output_base Root output directory produced by the workflow.
#' @return A \code{data.table} of all AD statistics.
#' @export
collectAllAdStats <- function(output_base) {
  all_files <- list.files(
    output_base,
    pattern = "^selected_genes_ad_stats\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (!length(all_files)) {
    warning("No selected_genes_ad_stats.csv files found", call. = FALSE)
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

#' Collect Significant Features Across All Drug-Tumor Pairs
#'
#' @description Reads mRNA_cell_sig or mRNA_pdcpdx_sig CSVs and merges them.
#' @param output_base Root output directory produced by the workflow.
#' @param type One of "cell" or "pdcpdx".
#' @return A \code{data.table} of all significant features.
#' @export
collectAllSigFeatures <- function(output_base,
                                  type = c("cell", "pdcpdx")) {
  type <- match.arg(type)
  pattern <- switch(type,
    cell   = "^mRNA_cell_sig\\.csv$",
    pdcpdx = "^mRNA_pdcpdx_sig\\.csv$"
  )

  all_files <- list.files(
    output_base,
    pattern = pattern,
    recursive = TRUE,
    full.names = TRUE
  )

  if (!length(all_files)) {
    warning("No ", type, " sig files found", call. = FALSE)
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

# PDF Export ----

#' Save a ggplot Object to PDF
#'
#' @description A thin wrapper around \code{ggplot2::ggsave()} with journal defaults.
#' Default dimensions target a Nature double-column figure (183 mm wide).
#' @param plot A ggplot object.
#' @param filename Output file path (should end in ".pdf").
#' @param width Plot width in inches. Default 7.2 (183 mm).
#' @param height Plot height in inches. Default 5.6.
#' @return The \code{filename}, invisibly.
#' @export
saveMetaVisPdf <- function(plot,
                           filename,
                           width = 7.2,
                           height = 5.6) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(
    filename = filename,
    plot     = plot,
    width    = width,
    height   = height,
    device   = "pdf",
    useDingbats = FALSE
  )
  invisible(filename)
}

#' Save a ComplexHeatmap to PDF
#'
#' @description Opens a PDF device, draws a ComplexHeatmap object, and closes
#' the device.
#' @param ht A \code{ComplexHeatmap} or \code{HeatmapList} object.
#' @param filename Output file path.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @return The \code{filename}, invisibly.
#' @export
saveComplexHeatmapPdf <- function(ht, filename,
                                  width = 7.2, height = 5.6) {
  dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
  pdf(filename, width = width, height = height)
  ComplexHeatmap::draw(ht)
  dev.off()
  invisible(filename)
}
