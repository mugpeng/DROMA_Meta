# TCGA and clinical validation helpers ----

mapTumorTypeToTcgaCode <- function(tumor_type) {
  cancer_to_project_simple <- c(
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
  key <- tolower(tumor_type)
  if (!key %in% names(cancer_to_project_simple)) {
    return(NULL)
  }
  project_code <- unname(cancer_to_project_simple[[key]])
  if (is.na(project_code) || !nzchar(project_code)) {
    return(NULL)
  }
  project_code
}

loadTcgaExpressionVector <- local({
  tcga_cache <- new.env(parent = emptyenv())
  gene_map_cache <- new.env(parent = emptyenv())
  counts_file_cache <- new.env(parent = emptyenv())
  read_gene_map <- function(map_path) {
    header_line <- tryCatch(readLines(map_path, n = 1L, warn = FALSE), error = function(e) "")
    if (!length(header_line)) {
      return(data.frame())
    }
    if (grepl(",", header_line[[1]], fixed = TRUE)) {
      utils::read.csv(map_path, stringsAsFactors = FALSE)
    } else {
      utils::read.delim(map_path, stringsAsFactors = FALSE)
    }
  }

  function(tcga_dir, tumor_type, gene, gene_probe_map = NULL) {
    tcga_code <- mapTumorTypeToTcgaCode(tumor_type)
    if (is.null(tcga_code)) {
      return(NULL)
    }

    dir_key <- normalizePath(tcga_dir, mustWork = FALSE)
    cache_key <- paste(dir_key, tcga_code, sep = "::")
    if (!exists(cache_key, envir = counts_file_cache, inherits = FALSE)) {
      candidate_files <- c(
        file.path(tcga_dir, paste0(tcga_code, ".htseq_counts.tsv.gz")),
        file.path(tcga_dir, "rna_counts", paste0(tcga_code, ".htseq_counts.tsv.gz"))
      )
      matched <- candidate_files[file.exists(candidate_files)]
      assign(
        cache_key,
        if (length(matched) > 0) matched[[1]] else NA_character_,
        envir = counts_file_cache
      )
    }
    tcga_file <- get(cache_key, envir = counts_file_cache, inherits = FALSE)

    if (is.na(tcga_file) || !nzchar(tcga_file)) {
      return(NULL)
    }
    if (!file.exists(tcga_file)) {
      return(NULL)
    }

    tcga_data_key <- paste(dir_key, tcga_code, sep = "::")
    if (!exists(tcga_data_key, envir = tcga_cache, inherits = FALSE)) {
      tcga_dt <- utils::read.delim(gzfile(tcga_file), check.names = FALSE, stringsAsFactors = FALSE)
      gene_col <- colnames(tcga_dt)[1]
      clean_ids <- sub("\\.[0-9]+$", "", tcga_dt[[gene_col]])
      gene_index <- stats::setNames(seq_len(nrow(tcga_dt)), clean_ids)
      assign(tcga_data_key, list(data = tcga_dt, gene_col = gene_col, gene_index = gene_index), envir = tcga_cache)
    }
    tcga_obj <- get(tcga_data_key, envir = tcga_cache, inherits = FALSE)

    lookup_gene <- gene
    candidate_map_paths <- unique(c(
      if (!is.null(gene_probe_map) && nzchar(gene_probe_map)) normalizePath(gene_probe_map, mustWork = FALSE) else character(0),
      file.path(tcga_dir, "gencode.human.v49.annotation.gene.probeMap"),
      file.path(tcga_dir, "gencode.v22.annotation.gene.probeMap"),
      file.path(dirname(tcga_dir), "gencode.human.v49.annotation.gene.probeMap"),
      file.path(dirname(tcga_dir), "gencode.v22.annotation.gene.probeMap")
    ))
    candidate_map_paths <- candidate_map_paths[file.exists(candidate_map_paths)]

    if (!lookup_gene %in% names(tcga_obj$gene_index) && length(candidate_map_paths) > 0) {
      for (map_key in candidate_map_paths) {
        if (!exists(map_key, envir = gene_map_cache, inherits = FALSE)) {
          gene_map_dt <- read_gene_map(map_key)
          if (!all(c("id", "gene") %in% colnames(gene_map_dt))) {
            assign(map_key, character(0), envir = gene_map_cache)
          } else {
            gene_ids <- sub("\\.[0-9]+$", "", gene_map_dt$id)
            assign(
              map_key,
              stats::setNames(gene_ids, gene_map_dt$gene),
              envir = gene_map_cache
            )
          }
        }
        gene_map_vec <- get(map_key, envir = gene_map_cache, inherits = FALSE)
        if (lookup_gene %in% names(gene_map_vec)) {
          mapped_gene <- unname(gene_map_vec[[lookup_gene]])
          if (nzchar(mapped_gene) && mapped_gene %in% names(tcga_obj$gene_index)) {
            lookup_gene <- mapped_gene
            break
          }
        }
      }
    }
    lookup_gene <- sub("\\.[0-9]+$", "", lookup_gene)
    row_idx <- unname(tcga_obj$gene_index[lookup_gene])
    if (length(row_idx) == 0 || is.na(row_idx)) {
      return(NULL)
    }
    row <- tcga_obj$data[row_idx, , drop = FALSE]
    values <- as.numeric(row[1, -1, drop = TRUE])
    values <- log2(values + 1)
    values[is.finite(values)]
  }
})

.ad_two_sample_test <- function(x, y) {
  if (length(x) < 3 || length(y) < 3) {
    return(list(p_value = NA_real_, method = NA_character_))
  }
  if (requireNamespace("kSamples", quietly = TRUE)) {
    ad_re <- tryCatch(kSamples::ad.test(x, y), error = function(e) NULL)
    if (!is.null(ad_re) && "ad" %in% names(ad_re)) {
      return(list(
        p_value = as.numeric(ad_re$ad[1, " asympt. P-value"]),
        method = "anderson-darling"
      ))
    }
  }
  ks_re <- tryCatch(stats::ks.test(x, y), error = function(e) NULL)
  if (is.null(ks_re)) {
    return(list(p_value = NA_real_, method = NA_character_))
  }
  list(
    p_value = as.numeric(ks_re$p.value),
    method = "kolmogorov-smirnov"
  )
}

.collectGroupExpression <- function(group_set, tumor_type, gene, feature_type = "mRNA") {
  values <- numeric()
  for (project_name in names(group_set@DromaSets)) {
    expr_mat <- tryCatch(withCallingHandlers(
      loadMolecularProfiles(
        object = group_set@DromaSets[[project_name]],
        feature_type = feature_type,
        select_features = gene,
        return_data = TRUE,
        data_type = "all",
        tumor_type = tumor_type,
        zscore = FALSE
      ),
      warning = function(w) {
        if (grepl("^The following samples do not exist:", conditionMessage(w)) ||
            grepl("^The following features do not exist in the .* data:", conditionMessage(w)) ||
            grepl("^None of the specified features exist", conditionMessage(w))) {
          invokeRestart("muffleWarning")
        }
      }
    ), error = function(e) NULL)
    if (is.matrix(expr_mat) && gene %in% rownames(expr_mat)) {
      expr_values <- as.numeric(expr_mat[gene, , drop = TRUE])
      expr_values <- expr_values[is.finite(expr_values)]
      if (length(expr_values) > 0) {
        values <- c(values, expr_values)
      }
    }
  }
  values
}

.createTranslationProgressCallback <- function(show_progress = TRUE, update_interval = 10) {
  if (!isTRUE(show_progress)) {
    return(function(done, total, elapsed) invisible(NULL))
  }

  last_update_time <- 0

  function(done, total, elapsed) {
    if ((elapsed - last_update_time) < update_interval && done < total) {
      return(invisible(NULL))
    }

    last_update_time <<- elapsed
    progress_pct <- (done / total) * 100

    if (done < total) {
      message(sprintf(
        "TCGA translation filter: %d/%d (%.1f%%) | Elapsed: %.1fs",
        done, total, progress_pct, elapsed
      ))
    } else {
      message(sprintf(
        "TCGA translation filter: %d/%d (100%%) | Total time: %.1fs",
        done, total, elapsed
      ))
    }

    invisible(NULL)
  }
}

.runTcgaTranslationOne <- function(row, cellline_set, pdcpdx_set, tcga_dir, feature_type = "mRNA",
                                   gene_probe_map = NULL) {
  tcga_values <- loadTcgaExpressionVector(
    tcga_dir, row$tumor_type[[1]], row$name[[1]],
    gene_probe_map = gene_probe_map
  )
  cellline_values <- .collectGroupExpression(cellline_set, row$tumor_type[[1]], row$name[[1]], feature_type)
  pdcpdx_values <- .collectGroupExpression(pdcpdx_set, row$tumor_type[[1]], row$name[[1]], feature_type)

  tcga_values <- if (is.null(tcga_values)) numeric() else tcga_values
  tcga_values <- if (length(tcga_values) >= 3) as.numeric(scale(tcga_values)) else tcga_values
  cellline_values <- if (length(cellline_values) >= 3) as.numeric(scale(cellline_values)) else cellline_values
  pdcpdx_values <- if (length(pdcpdx_values) >= 3) as.numeric(scale(pdcpdx_values)) else pdcpdx_values

  cellline_test <- if (length(tcga_values) >= 3 && length(cellline_values) >= 3) {
    .ad_two_sample_test(cellline_values, tcga_values)
  } else {
    list(p_value = NA_real_, method = NA_character_)
  }
  pdcpdx_test <- if (length(tcga_values) >= 3 && length(pdcpdx_values) >= 3) {
    .ad_two_sample_test(pdcpdx_values, tcga_values)
  } else {
    list(p_value = NA_real_, method = NA_character_)
  }

  data.table::data.table(
    drug = row$drug[[1]],
    tumor_type = row$tumor_type[[1]],
    name = row$name[[1]],
    tcga_ad_cellline_p = cellline_test$p_value,
    tcga_ad_pdcpdx_p = pdcpdx_test$p_value,
    tcga_ad_cellline_method = cellline_test$method,
    tcga_ad_pdcpdx_method = pdcpdx_test$method
  )
}

runTcgaTranslationFilter <- function(preclinical_candidates,
                                     cellline_set,
                                     pdcpdx_set,
                                     tcga_dir,
                                     feature_type = "mRNA",
                                     p_t = 0.01,
                                     cores = 1L,
                                     show_progress = TRUE,
                                     test_top_n = NULL,
                                     gene_probe_map = NULL) {
  candidates <- data.table::as.data.table(preclinical_candidates)
  if (nrow(candidates) == 0) {
    return(candidates)
  }

  if (!is.logical(show_progress) || length(show_progress) != 1) {
    stop("show_progress must be TRUE or FALSE")
  }

  if (!is.null(test_top_n)) {
    if (!is.numeric(test_top_n) || length(test_top_n) != 1 || test_top_n < 1 || test_top_n != as.integer(test_top_n)) {
      stop("test_top_n must be a positive integer or NULL")
    }
    if (nrow(candidates) > test_top_n) {
      candidates <- candidates[seq_len(test_top_n), ]
    }
  }

  max_cores <- if (requireNamespace("parallel", quietly = TRUE)) {
    detected_cores <- tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
    if (is.na(detected_cores) || detected_cores < 2L) 1L else detected_cores - 1L
  } else {
    1L
  }
  if (!is.numeric(cores) || length(cores) != 1 || cores < 1 || cores != as.integer(cores) || cores > max_cores) {
    stop(sprintf("cores must be a positive integer between 1 and %d", max_cores))
  }
  cores <- as.integer(cores)

  start_time <- Sys.time()
  progress_callback <- if (cores == 1L) {
    .createTranslationProgressCallback(show_progress = show_progress, update_interval = 10)
  } else {
    NULL
  }

  worker_function <- function(i) {
    result <- .runTcgaTranslationOne(
      row = candidates[i, ],
      cellline_set = cellline_set,
      pdcpdx_set = pdcpdx_set,
      tcga_dir = tcga_dir,
      feature_type = feature_type,
      gene_probe_map = gene_probe_map
    )
    if (!is.null(progress_callback)) {
      progress_callback(
        i,
        nrow(candidates),
        as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      )
    }
    result
  }

  if (show_progress) {
    message(sprintf("Running TCGA translation filter on %d candidate(s)...", nrow(candidates)))
  }

  if (cores > 1L && requireNamespace("future", quietly = TRUE) && requireNamespace("furrr", quietly = TRUE)) {
    old_size <- getOption("future.globals.maxSize")
    options(future.globals.maxSize = 2 * 1024^3)

    if (.Platform$OS.type == "unix") {
      future::plan(future::multicore, workers = cores)
    } else {
      future::plan(future::multisession, workers = cores)
    }
    on.exit({
      future::plan(future::sequential)
      options(future.globals.maxSize = old_size)
    }, add = TRUE)

    n_candidates <- nrow(candidates)
    chunk_size <- if (n_candidates < 100L) {
      max(5L, ceiling(n_candidates / cores))
    } else if (n_candidates < 1000L) {
      ceiling(n_candidates / (cores * 4L))
    } else {
      max(25L, min(100L, ceiling(n_candidates / (cores * 8L))))
    }
    indices <- seq_len(n_candidates)
    chunks <- split(indices, ceiling(indices / chunk_size))

    if (show_progress && requireNamespace("progressr", quietly = TRUE)) {
      rows <- progressr::with_progress({
        p <- progressr::progressor(steps = n_candidates)
        chunk_results <- furrr::future_map(chunks, function(chunk_indices) {
          lapply(chunk_indices, function(i) {
            result <- .runTcgaTranslationOne(
              row = candidates[i, ],
              cellline_set = cellline_set,
              pdcpdx_set = pdcpdx_set,
              tcga_dir = tcga_dir,
              feature_type = feature_type,
              gene_probe_map = gene_probe_map
            )
            p()
            result
          })
        }, .options = furrr::furrr_options(seed = TRUE))
        unlist(chunk_results, recursive = FALSE)
      })
    } else {
      if (show_progress && !requireNamespace("progressr", quietly = TRUE)) {
        message("Note: Install 'progressr' for parallel progress tracking: install.packages('progressr')")
      }
      chunk_results <- furrr::future_map(chunks, function(chunk_indices) {
        lapply(chunk_indices, function(i) {
          .runTcgaTranslationOne(
            row = candidates[i, ],
            cellline_set = cellline_set,
            pdcpdx_set = pdcpdx_set,
            tcga_dir = tcga_dir,
            feature_type = feature_type,
            gene_probe_map = gene_probe_map
          )
        })
      }, .options = furrr::furrr_options(seed = TRUE))
      rows <- unlist(chunk_results, recursive = FALSE)
    }
  } else {
    if (cores > 1L && show_progress) {
      message("Note: Install 'future' and 'furrr' to enable parallel TCGA translation filtering. Falling back to serial mode.")
    }
    rows <- lapply(seq_len(nrow(candidates)), worker_function)
  }

  tcga_dt <- data.table::rbindlist(rows, fill = TRUE)
  tcga_dt$tcga_ad_cellline_fdr <- stats::p.adjust(tcga_dt$tcga_ad_cellline_p, method = "BH")
  tcga_dt$tcga_ad_pdcpdx_fdr <- stats::p.adjust(tcga_dt$tcga_ad_pdcpdx_p, method = "BH")
  tcga_dt$tcga_supported <- (
    (!is.na(tcga_dt$tcga_ad_cellline_p) & tcga_dt$tcga_ad_cellline_p > p_t) |
      (!is.na(tcga_dt$tcga_ad_pdcpdx_p) & tcga_dt$tcga_ad_pdcpdx_p > p_t)
  )
  if (show_progress) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    message(sprintf("TCGA translation filter completed in %.1fs", elapsed))
  }
  tcga_dt[]
}
