test_that("createWorkflowProjectGroups splits projects by dataset type", {
  project_anno <- data.frame(
    project_name = c("CCLE", "PDC1", "PDO1", "PDX1"),
    dataset_type = c("CellLine", "PDC", "PDO", "PDX"),
    stringsAsFactors = FALSE
  )

  groups <- createWorkflowProjectGroups(project_anno)

  expect_equal(groups$cellline, c("CCLE", "PDC1"))
  expect_equal(groups$pdcpdx, c("PDO1", "PDX1"))
})

test_that("mapTumorTypeToTcgaCode maps common tumor names", {
  expect_equal(mapTumorTypeToTcgaCode("breast cancer"), "TCGA-BRCA")
  expect_equal(mapTumorTypeToTcgaCode("lung cancer"), "TCGA-LUAD")
  expect_null(mapTumorTypeToTcgaCode("unknown tumor"))
})

test_that("loadTcgaExpressionVector resolves rna_counts subdirectory", {
  td <- tempdir()
  dir.create(file.path(td, "rna_counts"), showWarnings = FALSE)
  tcga_file <- file.path(td, "rna_counts", "TCGA-BRCA.htseq_counts.tsv.gz")

  con <- gzfile(tcga_file, open = "wt")
  writeLines(c(
    "gene_id\tS1\tS2\tS3",
    "ENSG_TEST\t1\t3\t7"
  ), con = con)
  close(con)

  vals <- loadTcgaExpressionVector(td, "breast cancer", "ENSG_TEST")
  expect_equal(length(vals), 3L)
  expect_true(all(vals > 0))
})

test_that("runTcgaTranslationFilter supports test_top_n and validates show_progress", {
  candidates <- data.frame(
    drug = c("D1", "D2", "D3"),
    tumor_type = c("breast cancer", "breast cancer", "breast cancer"),
    name = c("G1", "G2", "G3"),
    stringsAsFactors = FALSE
  )

  testthat::local_mocked_bindings(
    loadTcgaExpressionVector = function(tcga_dir, tumor_type, gene) c(1, 2, 3, 4),
    .collectGroupExpression = function(group_set, tumor_type, gene, feature_type = "mRNA") c(1, 2, 3, 4),
    .ad_two_sample_test = function(x, y) list(p_value = 0.01, method = "mock-ad"),
    .env = asNamespace("DROMA.R2")
  )

  expect_error(
    runTcgaTranslationFilter(
      preclinical_candidates = candidates,
      cellline_set = NULL,
      pdcpdx_set = NULL,
      tcga_dir = tempdir(),
      show_progress = "yes"
    ),
    "show_progress must be TRUE or FALSE"
  )

  out <- runTcgaTranslationFilter(
    preclinical_candidates = candidates,
    cellline_set = NULL,
    pdcpdx_set = NULL,
    tcga_dir = tempdir(),
    cores = 1L,
    show_progress = FALSE,
    test_top_n = 2L
  )

  expect_equal(nrow(out), 2L)
  expect_equal(out$name, c("G1", "G2"))
})

test_that("runTcgaTranslationFilter keeps results consistent in parallel mode", {
  skip_if_not_installed("future")
  skip_if_not_installed("furrr")
  skip_if(parallel::detectCores() < 2, "parallel verification requires at least 2 detected cores")

  candidates <- data.frame(
    drug = c("D1", "D2", "D3", "D4"),
    tumor_type = rep("breast cancer", 4),
    name = c("G1", "G2", "G3", "G4"),
    stringsAsFactors = FALSE
  )

  testthat::local_mocked_bindings(
    loadTcgaExpressionVector = function(tcga_dir, tumor_type, gene) c(1, 2, 3, 4),
    .collectGroupExpression = function(group_set, tumor_type, gene, feature_type = "mRNA") {
      idx <- match(gene, candidates$name)
      c(idx, idx + 1, idx + 2, idx + 3)
    },
    .ad_two_sample_test = function(x, y) list(p_value = sum(x, na.rm = TRUE) / 100, method = "mock-ad"),
    .env = asNamespace("DROMA.R2")
  )

  serial <- runTcgaTranslationFilter(
    preclinical_candidates = candidates,
    cellline_set = NULL,
    pdcpdx_set = NULL,
    tcga_dir = tempdir(),
    cores = 1L,
    show_progress = FALSE
  )

  parallel <- runTcgaTranslationFilter(
    preclinical_candidates = candidates,
    cellline_set = NULL,
    pdcpdx_set = NULL,
    tcga_dir = tempdir(),
    cores = 2L,
    show_progress = FALSE
  )

  expect_equal(parallel[, names(serial), with = FALSE], serial)
})

test_that("runTcgaTranslationFilter treats non-significant distribution difference as TCGA support", {
  candidates <- data.frame(
    drug = c("D1", "D2"),
    tumor_type = c("breast cancer", "breast cancer"),
    name = c("G_same", "G_diff"),
    stringsAsFactors = FALSE
  )

  testthat::local_mocked_bindings(
    loadTcgaExpressionVector = function(tcga_dir, tumor_type, gene) c(1, 2, 3, 4),
    .collectGroupExpression = function(group_set, tumor_type, gene, feature_type = "mRNA") {
      if (identical(gene, "G_same")) c(1, 2, 3, 4) else c(10, 11, 12, 13)
    },
    .ad_two_sample_test = function(x, y) {
      list(
        p_value = if (isTRUE(all.equal(x, y))) 0.8 else 0.001,
        method = "mock-ad"
      )
    },
    .env = asNamespace("DROMA.R2")
  )

  out <- runTcgaTranslationFilter(
    preclinical_candidates = candidates,
    cellline_set = NULL,
    pdcpdx_set = NULL,
    tcga_dir = tempdir(),
    cores = 1L,
    show_progress = FALSE,
    fdr_t = 0.05
  )

  expect_true(out[name == "G_same", tcga_supported][[1]])
  expect_false(out[name == "G_diff", tcga_supported][[1]])
})

test_that("runTcgaTranslationFilter compares translation on unscaled expression values", {
  candidates <- data.frame(
    drug = "D1",
    tumor_type = "breast cancer",
    name = "G_shifted",
    stringsAsFactors = FALSE
  )

  captured_inputs <- list()
  testthat::local_mocked_bindings(
    loadTcgaExpressionVector = function(tcga_dir, tumor_type, gene) c(1, 2, 3, 4),
    .collectGroupExpression = function(group_set, tumor_type, gene, feature_type = "mRNA") c(10, 11, 12, 13),
    .ad_two_sample_test = function(x, y) {
      captured_inputs[[length(captured_inputs) + 1L]] <<- list(x = x, y = y)
      list(p_value = 0.001, method = "mock-ad")
    },
    .env = asNamespace("DROMA.R2")
  )

  runTcgaTranslationFilter(
    preclinical_candidates = candidates,
    cellline_set = NULL,
    pdcpdx_set = NULL,
    tcga_dir = tempdir(),
    cores = 1L,
    show_progress = FALSE,
    fdr_t = 0.05
  )

  expect_length(captured_inputs, 2L)
  expect_equal(captured_inputs[[1]]$x, c(10, 11, 12, 13))
  expect_equal(captured_inputs[[1]]$y, c(1, 2, 3, 4))
})
