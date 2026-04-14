test_that("createWorkflowProjectGroups splits projects by dataset type", {
  project_anno <- data.frame(
    project_name = c("CCLE", "PDC1", "PDO1", "PDX1"),
    dataset_type = c("CellLine", "PDC", "PDO", "PDX"),
    stringsAsFactors = FALSE
  )

  groups <- createWorkflowProjectGroups(project_anno)

  expect_equal(groups$invitro, c("CCLE", "PDC1"))
  expect_equal(groups$invivo, c("PDO1", "PDX1"))
})

test_that("mapTumorTypeToTcgaCode maps common tumor names", {
  expect_equal(mapTumorTypeToTcgaCode("breast cancer"), "TCGA-BRCA")
  expect_equal(mapTumorTypeToTcgaCode("lung cancer"), "TCGA-LUAD")
  expect_null(mapTumorTypeToTcgaCode("unknown tumor"))
})

test_that("mergePreclinicalCandidates preserves support flags", {
  invitro <- data.frame(
    drug = "Paclitaxel",
    tumor_type = "breast cancer",
    name = "CACNA1H",
    p_value = 0.001,
    q_value = 0.01,
    effect_size = 0.2,
    n_datasets = 3,
    model_group = "invitro",
    heterogeneity_p = 0.2,
    i2 = 0
  )
  invivo <- data.frame(
    drug = "Paclitaxel",
    tumor_type = "breast cancer",
    name = "CACNA1H",
    p_value = 0.003,
    q_value = 0.02,
    effect_size = 0.15,
    n_datasets = 2,
    model_group = "invivo",
    heterogeneity_p = 0.4,
    i2 = 0
  )

  out <- mergePreclinicalCandidates(invitro, invivo)

  expect_equal(nrow(out), 1)
  expect_true(out$invitro_supported)
  expect_true(out$invivo_supported)
  expect_true(out$direction_concordant)
})
