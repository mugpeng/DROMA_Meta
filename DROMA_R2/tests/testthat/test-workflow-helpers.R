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

test_that("computeCoverageForGroup counts eligible project pairs", {
  feature_map <- list(
    CCLE = c("M1", "M2", "M3"),
    GDSC = c("M4", "M5"),
    PDC1 = c("M6")
  )
  drug_map <- list(
    CCLE = c("M4", "M5"),
    GDSC = c("M2", "M3"),
    PDC1 = c("Z1")
  )

  out <- computeCoverageForGroup(
    project_names = c("CCLE", "GDSC", "PDC1"),
    drug_names = "Paclitaxel",
    tumor_types = "breast cancer",
    db_path = tempfile(fileext = ".sqlite"),
    con = NULL,
    min_overlap_per_pair = 2L,
    min_pairs = 2L,
    feature_sample_resolver = function(project_name, tumor_type) feature_map[[project_name]],
    drug_sample_resolver = function(project_name, drug_name, tumor_type) drug_map[[project_name]]
  )

  eligible <- out$coverage[out$coverage$eligible, c("drug_project", "expr_project")]
  expect_equal(nrow(eligible), 2L)
  expect_equal(sort(out$candidates_runtime$eligible_pairs[[1]]$drug_project), c("CCLE", "GDSC"))
  expect_equal(sort(out$candidates_runtime$eligible_pairs[[1]]$expr_project), c("CCLE", "GDSC"))
})

test_that("mergePreclinicalCandidates keeps only intersected support", {
  cellline <- data.frame(
    drug = "Paclitaxel",
    tumor_type = "breast cancer",
    name = "CACNA1H",
    p_value = 0.001,
    q_value = 0.01,
    effect_size = 0.2,
    n_datasets = 3,
    model_group = "cellline",
    direction = "Up"
  )
  pdcpdx <- data.frame(
    drug = "Paclitaxel",
    tumor_type = "breast cancer",
    name = "CACNA1H",
    p_value = 0.003,
    q_value = 0.02,
    effect_size = 0.15,
    n_datasets = 2,
    model_group = "pdcpdx",
    direction = "Up"
  )

  out <- mergePreclinicalCandidates(cellline, pdcpdx)
  expect_equal(nrow(out), 1L)
  expect_true(out$invitro_supported)
  expect_true(out$pdcpdx_supported)
  expect_true(out$direction_concordant)
})

test_that("mapTumorTypeToTcgaCode maps common tumor names", {
  expect_equal(mapTumorTypeToTcgaCode("breast cancer"), "TCGA-BRCA")
  expect_equal(mapTumorTypeToTcgaCode("lung cancer"), "TCGA-LUAD")
  expect_null(mapTumorTypeToTcgaCode("unknown tumor"))
})
