# TCGA and clinical validation helpers ----

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
    if (is.null(tcga_code)) {
      return(NULL)
    }
    tcga_file <- file.path(tcga_dir, paste0(tcga_code, ".htseq_counts.tsv.gz"))
    if (!file.exists(tcga_file)) {
      return(NULL)
    }

    if (!exists(tcga_code, envir = tcga_cache, inherits = FALSE)) {
      tcga_dt <- utils::read.delim(gzfile(tcga_file), check.names = FALSE, stringsAsFactors = FALSE)
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

.collectGroupExpression <- function(group_set, tumor_type, gene, feature_type = "mRNA") {
  values <- numeric()
  for (project_name in names(group_set@DromaSets)) {
    expr_mat <- tryCatch(
      loadMolecularProfiles(
        object = group_set@DromaSets[[project_name]],
        feature_type = feature_type,
        select_features = gene,
        return_data = TRUE,
        data_type = "all",
        tumor_type = tumor_type,
        zscore = FALSE
      ),
      error = function(e) NULL
    )
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

runTcgaTranslationFilter <- function(preclinical_candidates,
                                     cellline_set,
                                     pdcpdx_set,
                                     tcga_dir,
                                     feature_type = "mRNA",
                                     fdr_t = 0.01,
                                     cores = 1L) {
  candidates <- data.table::as.data.table(preclinical_candidates)
  if (nrow(candidates) == 0) {
    return(candidates)
  }

  rows <- .runParallelRows(seq_len(nrow(candidates)), function(i) {
    row <- candidates[i]
    tcga_values <- loadTcgaExpressionVector(tcga_dir, row$tumor_type[[1]], row$name[[1]])
    cellline_values <- .collectGroupExpression(cellline_set, row$tumor_type[[1]], row$name[[1]], feature_type)
    pdcpdx_values <- .collectGroupExpression(pdcpdx_set, row$tumor_type[[1]], row$name[[1]], feature_type)

    data.table::data.table(
      drug = row$drug[[1]],
      tumor_type = row$tumor_type[[1]],
      name = row$name[[1]],
      tcga_ad_cellline_p = if (!is.null(tcga_values)) .ad_two_sample_p(scale(cellline_values), scale(tcga_values)) else NA_real_,
      tcga_ad_pdcpdx_p = if (!is.null(tcga_values)) .ad_two_sample_p(scale(pdcpdx_values), scale(tcga_values)) else NA_real_
    )
  }, cores = cores)

  tcga_dt <- data.table::rbindlist(rows, fill = TRUE)
  tcga_dt[, tcga_ad_cellline_fdr := stats::p.adjust(tcga_ad_cellline_p, method = "BH")]
  tcga_dt[, tcga_ad_pdcpdx_fdr := stats::p.adjust(tcga_ad_pdcpdx_p, method = "BH")]
  tcga_dt[, tcga_supported := (tcga_ad_cellline_fdr < fdr_t) | (tcga_ad_pdcpdx_fdr < fdr_t)]
  tcga_dt[]
}

runClinicalValidationForCandidates <- function(candidate_df,
                                               drug_name,
                                               tumor_type,
                                               ctrdb_path,
                                               cores = 1L,
                                               es_t = 0.1,
                                               fdr_t = 0.1) {
  candidate_df <- data.table::as.data.table(candidate_df)
  if (nrow(candidate_df) == 0) {
    return(candidate_df)
  }

  DROMA.Set::connectCTRDatabase(ctrdb_path)
  on.exit({
    if (exists("ctrdb_connection", envir = .GlobalEnv)) {
      rm("ctrdb_connection", envir = .GlobalEnv)
    }
  }, add = TRUE)

  run_one_scope <- function(scope) {
    tryCatch(
      DROMA.R::batchFindClinicalSigResponse(
        select_omics = candidate_df$name,
        select_drugs = drug_name,
        data_type = "all",
        tumor_type = scope,
        cores = cores
      ),
      error = function(e) NULL
    )
  }

  clinical_df <- run_one_scope(tumor_type)
  clinical_scope <- tumor_type
  fallback_to_all <- FALSE
  if (is.null(clinical_df) || !is.data.frame(clinical_df) || nrow(clinical_df) == 0) {
    clinical_df <- run_one_scope("all")
    clinical_scope <- "all"
    fallback_to_all <- TRUE
  }

  if (is.null(clinical_df) || !is.data.frame(clinical_df) || nrow(clinical_df) == 0) {
    out <- data.table::copy(candidate_df)
    out[, `:=`(
      clinical_scope = clinical_scope,
      fallback_to_all = fallback_to_all,
      clinical_supported = FALSE,
      direction_clinical = NA_character_,
      direction_concordant = FALSE,
      retained = FALSE
    )]
    return(out)
  }

  clinical_candidates <- getClinicalCandidateFeatures(
        select_omics = candidate_df$name,
    clinical_results = clinical_df,
    es_t = es_t,
    P_t = fdr_t,
    use_p_value = FALSE
  )

  if (!"direction" %in% colnames(clinical_df) && "effect_size" %in% colnames(clinical_df)) {
    clinical_df$direction <- ifelse(clinical_df$effect_size >= 0, "Up", "Down")
  }

  merged <- merge(
    candidate_df,
    clinical_df[, c("name", "p_value", "q_value", "effect_size", "n_datasets", "direction"), with = FALSE],
    by = "name",
    all.x = TRUE,
    suffixes = c("_preclinical", "_clinical")
  )
  merged[, `:=`(
    clinical_scope = clinical_scope,
    fallback_to_all = fallback_to_all,
    clinical_supported = name %in% clinical_candidates$name,
    direction_clinical = direction,
    preclinical_direction_concordant = direction_concordant
  )]
  merged[, direction_concordant := clinical_supported & !is.na(direction_clinical) & preclinical_direction == direction_clinical]
  merged[, retained := clinical_supported & direction_concordant]
  merged[]
}
