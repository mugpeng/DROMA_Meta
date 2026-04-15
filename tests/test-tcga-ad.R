fallback_if_missing <- function(x, y) {
  if (is.null(x) || length(x) == 0 || !nzchar(x)) y else x
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- sub("^--file=", "", fallback_if_missing(file_arg[1], ""))
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path, mustWork = TRUE)) else getwd()
root_dir <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

suppressPackageStartupMessages(library(data.table))

source(file.path(root_dir, "R", "FuncHelper.R"))
source(file.path(root_dir, "R", "FuncTcgaAD.R"))

test_batch_find_tcga_ad_handles_non_matrix_entries <- function() {
  selected_features <- data.table(name = "GENE1")

  droma_r_ns <- asNamespace("DROMA.R")
  original_load_feature_data <- get("loadFeatureData", envir = droma_r_ns)
  unlockBinding("loadFeatureData", droma_r_ns)
  assign(
    "loadFeatureData",
    function(...) {
      list(
        bad_entry = structure(list(), class = "not_a_matrix"),
        good_entry = structure(
          c(1, 2, 3),
          dim = c(1, 3),
          dimnames = list("GENE1", c("S1", "S2", "S3"))
        )
      )
    },
    envir = droma_r_ns
  )
  assign(
    "loadTcgaFeatureData",
    function(...) {
      list(
        GENE1 = list(
          tcga_values = c(1, 2, 3),
          tcga_n = 3L,
          tcga_status = "ok"
        )
      )
    },
    envir = .GlobalEnv
  )

  on.exit({
    assign("loadFeatureData", original_load_feature_data, envir = droma_r_ns)
    lockBinding("loadFeatureData", droma_r_ns)
    rm("loadTcgaFeatureData", envir = .GlobalEnv)
  }, add = TRUE)

  preclinical_set <- structure(list(), class = "MultiDromaSet")
  result <- batchFindTcgaADConcordantFeatures(
    selected_features = selected_features,
    preclinical_set = preclinical_set,
    tumor_type = "breast cancer",
    tcga_rna_counts_dir = tempfile("tcga-dir-"),
    gene_probe_map_path = tempfile("probe-map-"),
    feature_type = "mRNA",
    data_type = "all",
    p_t = 0.05,
    preclinical_label = "ccle"
  )

  stopifnot(nrow(result) == 1L)
  stopifnot(identical(result$name[[1]], "GENE1"))
  stopifnot(identical(result$ccle_n[[1]], 3L))
  stopifnot(identical(result$tcga_n[[1]], 3L))
}

test_batch_find_tcga_ad_handles_multiple_rows <- function() {
  selected_features <- data.table(name = c("GENE1", "GENE2"))

  droma_r_ns <- asNamespace("DROMA.R")
  original_load_feature_data <- get("loadFeatureData", envir = droma_r_ns)
  unlockBinding("loadFeatureData", droma_r_ns)
  assign(
    "loadFeatureData",
    function(...) {
      structure(
        list(
          good_entry = structure(
            c(1, 2, 3, 3, 4, 5),
            dim = c(2, 3),
            dimnames = list(c("GENE1", "GENE2"), c("S1", "S2", "S3"))
          )
        ),
        names = "CCLE"
      )
    },
    envir = droma_r_ns
  )
  assign(
    "loadTcgaFeatureData",
    function(...) {
      list(
        GENE1 = list(tcga_values = c(1, 2, 3), tcga_n = 3L, tcga_status = "ok"),
        GENE2 = list(tcga_values = c(3, 4, 5), tcga_n = 3L, tcga_status = "ok")
      )
    },
    envir = .GlobalEnv
  )

  on.exit({
    assign("loadFeatureData", original_load_feature_data, envir = droma_r_ns)
    lockBinding("loadFeatureData", droma_r_ns)
    rm("loadTcgaFeatureData", envir = .GlobalEnv)
  }, add = TRUE)

  preclinical_set <- structure(list(), class = "MultiDromaSet")
  result <- batchFindTcgaADConcordantFeatures(
    selected_features = selected_features,
    preclinical_set = preclinical_set,
    tumor_type = "breast cancer",
    tcga_rna_counts_dir = tempfile("tcga-dir-"),
    gene_probe_map_path = tempfile("probe-map-"),
    feature_type = "mRNA",
    data_type = "all",
    p_t = 0.05,
    preclinical_label = "ccle"
  )

  stopifnot(nrow(result) == 2L)
  stopifnot(identical(result$name, c("GENE1", "GENE2")))
  stopifnot(identical(result$ccle_n, c(3L, 3L)))
  stopifnot(identical(result$tcga_n, c(3L, 3L)))
}

test_batch_find_tcga_ad_handles_non_matrix_entries()
test_batch_find_tcga_ad_handles_multiple_rows()

cat("tcga ad tests passed\n")
