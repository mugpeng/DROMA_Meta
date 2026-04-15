# TCGA Anderson-Darling concordance helpers ----

#' Get Tumor-Type Mapping to TCGA or TARGET Cohorts
#'
#' @description Returns the internal lookup used to map DROMA tumor types onto
#' TCGA or TARGET cohort identifiers for concordance filtering.
#' @return A named character vector from tumor type to cohort code.
#' @export
getTcgaTumorTypeMapping <- function() {
  c(
    "haematopoietic/lymphoid cancer" = "TCGA-LAML",
    "nervous system cancer" = "TCGA-GBM",
    "sarcoma" = "TCGA-SARC",
    "lung cancer" = "TCGA-LUAD",
    "breast cancer" = "TCGA-BRCA",
    "prostate cancer" = "TCGA-PRAD",
    "stomach cancer" = "TCGA-STAD",
    "bladder cancer" = "TCGA-BLCA",
    "skin cancer" = "TCGA-SKCM",
    "ovarian cancer" = "TCGA-OV",
    "kidney cancer" = "TCGA-KIRC",
    "thyroid cancer" = "TCGA-THCA",
    "aerodigestive tract cancer" = "TCGA-HNSC",
    "vulvar cancer" = NA_character_,
    "endometrial cancer" = "TCGA-UCEC",
    "non-cancer" = NA_character_,
    "colon cancer" = "TCGA-COAD",
    "pancreatic cancer" = "TCGA-PAAD",
    "cervical cancer" = "TCGA-CESC",
    "liver cancer" = "TCGA-LIHC",
    "choriocarcinoma" = NA_character_,
    "uterine cancer" = "TCGA-UCS",
    "testicular cancer" = "TCGA-TGCT",
    "nasopharyngeal cancer" = NA_character_,
    "retinoblastoma" = "TARGET-RT"
  )
}

#' Match a Tumor Type to Its TCGA or TARGET Cohort
#'
#' @description Resolves one DROMA tumor type into the corresponding cohort code
#' used by the TCGA Anderson-Darling filtering step.
#' @param tumor_type Tumor type string.
#' @return A single TCGA or TARGET cohort label.
#' @export
getMatchedTcgaTumorType <- function(tumor_type) {
  mapping <- getTcgaTumorTypeMapping()
  if (is.null(tumor_type) || !nzchar(tumor_type)) {
    stop("tumor_type must be a non-empty string", call. = FALSE)
  }

  if (!tumor_type %in% names(mapping)) {
    stop("No TCGA/TARGET mapping found for tumor_type: ", tumor_type, call. = FALSE)
  }

  matched <- unname(mapping[[tumor_type]])
  if (is.na(matched) || !nzchar(matched)) {
    stop("tumor_type is not supported for TCGA/TARGET AD filtering: ", tumor_type, call. = FALSE)
  }

  matched
}

#' Load TCGA Feature Data for Selected Genes
#'
#' @description Reads TCGA/TARGET count matrices and the probe-map file, caches
#' them in-memory, and returns per-gene vectors ready for concordance testing.
#' @param select_features Character vector of gene names.
#' @param tcga_tumor_type TCGA/TARGET cohort code.
#' @param tcga_rna_counts_dir Directory containing count matrices.
#' @param gene_probe_map_path Path to the probe-map file.
#' @return A named list of per-gene TCGA data summaries.
#' @export
loadTcgaFeatureData <- local({
  tcga_cache <- new.env(parent = emptyenv())
  gene_map_cache <- new.env(parent = emptyenv())

  readProbeMap <- function(path) {
    header_line <- tryCatch(readLines(path, n = 1L, warn = FALSE), error = function(e) "")
    if (!length(header_line)) {
      return(data.frame())
    }

    if (grepl(",", header_line[[1]], fixed = TRUE)) {
      utils::read.csv(path, stringsAsFactors = FALSE)
    } else {
      utils::read.delim(path, stringsAsFactors = FALSE)
    }
  }

  function(select_features,
           tcga_tumor_type,
           tcga_rna_counts_dir,
           gene_probe_map_path) {
    if (missing(select_features) || !length(select_features)) {
      stop("select_features must contain at least one gene", call. = FALSE)
    }

    tcga_file <- file.path(tcga_rna_counts_dir, paste0(tcga_tumor_type, ".htseq_counts.tsv.gz"))
    if (!file.exists(tcga_file)) {
      stop("TCGA/TARGET counts file not found: ", tcga_file, call. = FALSE)
    }
    if (!file.exists(gene_probe_map_path)) {
      stop("gene_probe_map_path not found: ", gene_probe_map_path, call. = FALSE)
    }

    file_key <- normalizePath(tcga_file, mustWork = TRUE)
    if (!exists(file_key, envir = tcga_cache, inherits = FALSE)) {
      tcga_dt <- utils::read.delim(gzfile(tcga_file), check.names = FALSE, stringsAsFactors = FALSE)
      gene_col <- colnames(tcga_dt)[1]
      clean_ids <- sub("\\.[0-9]+$", "", tcga_dt[[gene_col]])
      gene_index <- stats::setNames(seq_len(nrow(tcga_dt)), clean_ids)
      assign(file_key, list(data = tcga_dt, gene_col = gene_col, gene_index = gene_index), envir = tcga_cache)
    }
    tcga_obj <- get(file_key, envir = tcga_cache, inherits = FALSE)

    map_key <- normalizePath(gene_probe_map_path, mustWork = TRUE)
    if (!exists(map_key, envir = gene_map_cache, inherits = FALSE)) {
      gene_map_dt <- readProbeMap(gene_probe_map_path)
      if (!all(c("id", "gene") %in% colnames(gene_map_dt))) {
        stop("gene_probe_map_path must contain 'id' and 'gene' columns", call. = FALSE)
      }
      gene_ids <- sub("\\.[0-9]+$", "", as.character(gene_map_dt$id))
      gene_map_vec <- stats::setNames(gene_ids, as.character(gene_map_dt$gene))
      assign(map_key, gene_map_vec, envir = gene_map_cache)
    }
    gene_map_vec <- get(map_key, envir = gene_map_cache, inherits = FALSE)

    selected_features <- unique(as.character(select_features))
    rows <- vector("list", length(selected_features))
    names(rows) <- selected_features

    for (gene in selected_features) {
      status <- "ok"
      lookup_gene <- sub("\\.[0-9]+$", "", gene)
      row_idx <- unname(tcga_obj$gene_index[lookup_gene])

      if (length(row_idx) == 0 || is.na(row_idx)) {
        mapped_gene <- unname(gene_map_vec[gene])
        mapped_gene <- sub("\\.[0-9]+$", "", mapped_gene)
        if (!is.na(mapped_gene) && nzchar(mapped_gene)) {
          row_idx <- unname(tcga_obj$gene_index[mapped_gene])
        }
      }

      if (length(row_idx) == 0 || is.na(row_idx)) {
        status <- "gene_not_found_in_tcga"
        values <- numeric()
      } else {
        row <- tcga_obj$data[row_idx, , drop = FALSE]
        values <- suppressWarnings(as.numeric(row[1, -1, drop = TRUE]))
        values <- values[is.finite(values)]
        values <- log2(values + 1)
      }

      rows[[gene]] <- list(
        gene = gene,
        tcga_tumor_type = tcga_tumor_type,
        tcga_values = values,
        tcga_n = length(values),
        tcga_status = status
      )
    }

    rows
  }
})

#' Calculate Anderson-Darling Concordance Between Two Feature Vectors
#'
#' @description Standardizes two numeric vectors, runs `kSamples::ad.test()`,
#' and reports the test statistic, p-value, and concordance call.
#' @param x First numeric vector.
#' @param y Second numeric vector.
#' @param p_t P-value threshold used to call concordance.
#' @return A list with AD statistic, p-value, concordance flag, and status.
#' @export
calcFeatureADConcordance <- function(x, y, p_t = 0.05) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]

  result <- list(
    ad_stat = NA_real_,
    ad_p_value = NA_real_,
    concordant = FALSE,
    status = "ok"
  )

  if (length(x) < 3) {
    result$status <- "insufficient_x_data"
    return(result)
  }
  if (length(y) < 3) {
    result$status <- "insufficient_y_data"
    return(result)
  }

  x_sd <- stats::sd(x)
  y_sd <- stats::sd(y)
  if (!is.finite(x_sd) || x_sd == 0) {
    result$status <- "zero_variance_x"
    return(result)
  }
  if (!is.finite(y_sd) || y_sd == 0) {
    result$status <- "zero_variance_y"
    return(result)
  }

  x_z <- as.numeric(scale(x))
  y_z <- as.numeric(scale(y))

  ad_re <- tryCatch(kSamples::ad.test(x_z, y_z), error = function(e) e)
  if (inherits(ad_re, "error")) {
    result$status <- paste0("ad_test_error: ", conditionMessage(ad_re))
    return(result)
  }

  if (!("ad" %in% names(ad_re))) {
    result$status <- "ad_test_missing_output"
    return(result)
  }

  stat_col <- intersect(c("T.AD", "ad"), colnames(ad_re$ad))
  p_col <- intersect(c(" asympt. P-value", "asympt. P-value", " asymptotic P-value"), colnames(ad_re$ad))

  if (!length(stat_col) || !length(p_col)) {
    result$status <- "ad_test_missing_columns"
    return(result)
  }

  result$ad_stat <- as.numeric(ad_re$ad[1, stat_col[1]])
  result$ad_p_value <- as.numeric(ad_re$ad[1, p_col[1]])
  result$concordant <- is.finite(result$ad_p_value) && result$ad_p_value >= p_t
  result
}

#' Batch TCGA Concordance Filtering for Selected Features
#'
#' @description Compares selected preclinical features against the matched
#' TCGA/TARGET cohort and appends Anderson-Darling concordance statistics.
#' @param selected_features Feature table containing at least a `name` column.
#' @param preclinical_set `MultiDromaSet` object used as the preclinical source.
#' @param tumor_type Tumor type string.
#' @param tcga_rna_counts_dir Directory containing TCGA/TARGET count matrices.
#' @param gene_probe_map_path Path to the probe-map file.
#' @param feature_type Feature type to load, typically `"mRNA"`.
#' @param data_type Data type passed through to `loadFeatureData()`.
#' @param p_t P-value threshold used to call concordance.
#' @param preclinical_label Label used as a prefix in output column names.
#' @return A `data.table` with appended concordance statistics.
#' @export
batchFindTcgaADConcordantFeatures <- function(selected_features,
                                              preclinical_set,
                                              tumor_type,
                                              tcga_rna_counts_dir,
                                              gene_probe_map_path,
                                              feature_type = "mRNA",
                                              data_type = "all",
                                              p_t = 0.05,
                                              preclinical_label = "ccle") {
  selected_dt <- data.table::as.data.table(selected_features)
  if (!nrow(selected_dt)) {
    return(data.table::as.data.table(selected_dt))
  }

  tcga_tumor_type <- getMatchedTcgaTumorType(tumor_type)
  genes <- unique(as.character(selected_dt$name))
  if (!inherits(preclinical_set, "MultiDromaSet")) {
    stop("preclinical_set must be a MultiDromaSet object", call. = FALSE)
  }

  preclinical_data <- loadFeatureData(
    dromaset_object = preclinical_set,
    feature_type = feature_type,
    select_features = genes,
    data_type = data_type,
    tumor_type = tumor_type,
    overlap_only = FALSE,
    is_continuous = TRUE,
    zscore = FALSE
  ) %||% list()
  tcga_data <- loadTcgaFeatureData(
    select_features = genes,
    tcga_tumor_type = tcga_tumor_type,
    tcga_rna_counts_dir = tcga_rna_counts_dir,
    gene_probe_map_path = gene_probe_map_path
  )

  collectValues <- function(preloaded_data, gene) {
    values <- numeric()
    if (!length(preloaded_data)) {
      return(values)
    }

    for (dataset in preloaded_data) {
      if (!is.matrix(dataset) || !gene %in% rownames(dataset)) {
        next
      }
      gene_values <- suppressWarnings(as.numeric(dataset[gene, , drop = TRUE]))
      gene_values <- gene_values[is.finite(gene_values)]
      if (length(gene_values)) {
        values <- c(values, gene_values)
      }
    }
    values
  }

  stat_col <- paste0(preclinical_label, "_vs_tcga_ad_stat")
  p_col <- paste0(preclinical_label, "_vs_tcga_ad_p")
  concordant_col <- paste0(preclinical_label, "_vs_tcga_concordant")
  n_col <- paste0(preclinical_label, "_n")

  ad_rows <- lapply(seq_len(nrow(selected_dt)), function(i) {
    row <- selected_dt[i]
    gene <- as.character(row$name[[1]])

    preclinical_values <- collectValues(preclinical_data, gene)
    tcga_row <- tcga_data[[gene]]
    tcga_values <- tcga_row$tcga_values %||% numeric()
    ad_result <- calcFeatureADConcordance(preclinical_values, tcga_values, p_t = p_t)

    ad_status <- c(
      if (!identical(tcga_row$tcga_status, "ok")) paste0("tcga=", tcga_row$tcga_status),
      if (!identical(ad_result$status, "ok")) paste0(preclinical_label, "=", ad_result$status)
    )
    if (!length(ad_status)) {
      ad_status <- "ok"
    }

    out <- data.table::data.table(
      row,
      tcga_tumor_type = tcga_tumor_type,
      tcga_n = length(tcga_values),
      ad_status = paste(ad_status, collapse = ";")
    )
    out[[n_col]] <- length(preclinical_values)
    out[[stat_col]] <- ad_result$ad_stat
    out[[p_col]] <- ad_result$ad_p_value
    out[[concordant_col]] <- ad_result$concordant
    out
  })

  data.table::rbindlist(ad_rows, fill = TRUE)
}
