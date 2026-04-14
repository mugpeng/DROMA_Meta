mapTumorTypeToTcgaCode <- function(tumor_type) {
  tumor_map <- c(
    "breast cancer" = "TCGA-BRCA",
    "lung adenocarcinoma" = "TCGA-LUAD",
    "lung cancer" = "TCGA-LUAD",
    "lung squamous cell carcinoma" = "TCGA-LUSC",
    "colon cancer" = "TCGA-COAD",
    "rectal cancer" = "TCGA-READ",
    "pancreatic cancer" = "TCGA-PAAD",
    "ovarian cancer" = "TCGA-OV",
    "prostate cancer" = "TCGA-PRAD",
    "liver cancer" = "TCGA-LIHC",
    "stomach cancer" = "TCGA-STAD",
    "melanoma" = "TCGA-SKCM",
    "glioma" = "TCGA-LGG",
    "glioblastoma" = "TCGA-GBM",
    "head and neck cancer" = "TCGA-HNSC"
  )
  key <- tolower(tumor_type)
  if (!key %in% names(tumor_map)) {
    return(NULL)
  }
  unname(tumor_map[[key]])
}

loadTcgaExpressionVector <- local({
  tcga_cache <- new.env(parent = emptyenv())
  gene_map_cache <- NULL

  function(tcga_dir, tumor_type, gene) {
  tcga_code <- mapTumorTypeToTcgaCode(tumor_type)
  if (is.null(tcga_code) || !nzchar(tcga_code)) {
    return(NULL)
  }
  tcga_file <- file.path(tcga_dir, paste0(tcga_code, ".htseq_counts.tsv.gz"))
  if (!file.exists(tcga_file)) {
    return(NULL)
  }

  if (!exists(tcga_code, envir = tcga_cache, inherits = FALSE)) {
    tcga_dt <- utils::read.delim(
      gzfile(tcga_file),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    gene_col <- colnames(tcga_dt)[1]
    clean_ids <- sub("\\.[0-9]+$", "", tcga_dt[[gene_col]])
    gene_index <- stats::setNames(seq_len(nrow(tcga_dt)), clean_ids)
    assign(tcga_code, list(data = tcga_dt, gene_col = gene_col, gene_index = gene_index), envir = tcga_cache)
  }
  tcga_obj <- get(tcga_code, envir = tcga_cache, inherits = FALSE)
  if (is.null(gene_map_cache)) {
    map_file <- "/Users/peng/Desktop/Project/DROMA/others/archive/Align/Input/bulkformer/bulkformer_gene_info.csv"
    if (file.exists(map_file)) {
      gene_map_dt <- utils::read.csv(map_file, stringsAsFactors = FALSE)
      gene_map_cache <<- stats::setNames(gene_map_dt$ensg_id, gene_map_dt$gene_symbol)
    } else {
      gene_map_cache <<- character(0)
    }
  }
  lookup_gene <- gene
  if (!lookup_gene %in% names(tcga_obj$gene_index) && lookup_gene %in% names(gene_map_cache)) {
    lookup_gene <- unname(gene_map_cache[[lookup_gene]])
  }
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

.ad_two_sample_p <- function(x, y) {
  if (length(x) < 3 || length(y) < 3) {
    return(NA_real_)
  }
  if (requireNamespace("kSamples", quietly = TRUE)) {
    ad_re <- tryCatch(kSamples::ad.test(x, y), error = function(e) NULL)
    if (!is.null(ad_re) && "ad" %in% names(ad_re)) {
      return(as.numeric(ad_re$ad[1, " asympt. P-value"]))
    }
  }
  ks_re <- tryCatch(stats::ks.test(x, y), error = function(e) NULL)
  if (is.null(ks_re)) {
    return(NA_real_)
  }
  as.numeric(ks_re$p.value)
}

runTcgaTranslationFilter <- function(meta_candidates,
                                     group_set,
                                     tcga_dir,
                                     feature_type = "mRNA",
                                     fdr_t = 0.1,
                                     es_t = 0.1) {
  meta_candidates <- data.table::as.data.table(meta_candidates)
  if (nrow(meta_candidates) == 0) {
    return(meta_candidates[, `:=`(tcga_p_value = numeric(), tcga_q_value = numeric(), tcga_supported = logical())])
  }
  meta_candidates <- meta_candidates[q_value < fdr_t & abs(effect_size) >= es_t]
  if (nrow(meta_candidates) == 0) {
    return(meta_candidates[, `:=`(tcga_p_value = numeric(), tcga_q_value = numeric(), tcga_supported = logical())])
  }

  tumor_gene_map <- split(meta_candidates$name, meta_candidates$tumor_type)
  expr_cache <- lapply(names(group_set@DromaSets), function(project_name) {
    project_set <- group_set@DromaSets[[project_name]]
    tumor_cache <- lapply(names(tumor_gene_map), function(tumor_type) {
      genes <- unique(tumor_gene_map[[tumor_type]])
      expr_mat <- tryCatch(
        loadMolecularProfiles(
          object = project_set,
          feature_type = feature_type,
          select_features = genes,
          return_data = TRUE,
          data_type = "all",
          tumor_type = tumor_type,
          zscore = FALSE
        ),
        error = function(e) NULL
      )
      expr_mat
    })
    names(tumor_cache) <- names(tumor_gene_map)
    tumor_cache
  })
  names(expr_cache) <- names(group_set@DromaSets)

  rows <- lapply(seq_len(nrow(meta_candidates)), function(i) {
    row <- meta_candidates[i]
    model_values <- numeric()
    for (project_name in names(expr_cache)) {
      expr_mat <- expr_cache[[project_name]][[row$tumor_type[[1]]]]
      if (is.matrix(expr_mat) && row$name[[1]] %in% rownames(expr_mat)) {
        expr_vec <- as.numeric(expr_mat[row$name[[1]], , drop = TRUE])
        expr_vec <- expr_vec[!is.na(expr_vec)]
        if (length(expr_vec) > 0) {
          model_values <- c(model_values, expr_vec)
        }
      }
    }
    tcga_values <- loadTcgaExpressionVector(tcga_dir, row$tumor_type[[1]], row$name[[1]])
    p_value <- if (length(model_values) >= 3 && length(tcga_values) >= 3) {
      .ad_two_sample_p(scale(model_values), scale(tcga_values))
    } else {
      NA_real_
    }
    data.table::data.table(
      drug = row$drug[[1]],
      tumor_type = row$tumor_type[[1]],
      name = row$name[[1]],
      tcga_p_value = p_value
    )
  })

  tcga_dt <- data.table::rbindlist(rows, fill = TRUE)
  tcga_dt[, tcga_q_value := stats::p.adjust(tcga_p_value, method = "BH")]
  tcga_dt[, tcga_supported := is.na(tcga_q_value) | tcga_q_value >= 0.1]
  merge(meta_candidates, tcga_dt, by = c("drug", "tumor_type", "name"), all.x = TRUE)
}

runClinicalValidation <- function(candidate_df,
                                  drug_name,
                                  tumor_type = "all",
                                  db_path = NULL,
                                  ctrdb_path = NULL,
                                  cores = 3L,
                                  es_t = 0.1,
                                  fdr_t = 0.1) {
  candidate_df <- data.table::as.data.table(candidate_df)
  if (nrow(candidate_df) == 0) {
    return(candidate_df)
  }
  if (!is.null(db_path) && file.exists(db_path)) {
    con <- connectDROMADatabase(db_path)
    on.exit(closeDROMADatabase(con), add = TRUE)
  }
  if (!is.null(ctrdb_path) && file.exists(ctrdb_path)) {
    connectCTRDatabase(ctrdb_path)
  }
  cli_df <- batchFindClinicalSigResponse(
    select_omics = candidate_df$name,
    select_drugs = drug_name,
    data_type = "all",
    tumor_type = tumor_type,
    cores = cores
  )
  cli_sig <- getClinicalCandidateFeatures(
    clinical_results = cli_df,
    es_t = es_t,
    P_t = fdr_t,
    use_p_value = FALSE
  )
  merged <- merge(candidate_df, cli_df, by = "name", all.x = TRUE, suffixes = c("_preclinical", "_clinical"))
  merged[, clinical_supported := name %in% cli_sig$name]
  if ("effect_size_clinical" %in% colnames(merged)) {
    merged[, direction_preclinical := ifelse(effect_size_preclinical >= 0, "Up", "Down")]
    merged[, direction_clinical := ifelse(effect_size_clinical >= 0, "Up", "Down")]
    merged[, direction_concordant := clinical_supported & direction_preclinical == direction_clinical]
    merged[, retained := clinical_supported & direction_concordant]
  }
  merged[]
}
