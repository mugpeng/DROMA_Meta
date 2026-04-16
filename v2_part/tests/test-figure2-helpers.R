fallback_if_missing <- function(x, y) {
  if (is.null(x) || length(x) == 0 || !nzchar(x)) y else x
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- sub("^--file=", "", fallback_if_missing(file_arg[1], ""))
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path, mustWork = TRUE)) else getwd()
v2_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

suppressPackageStartupMessages(library(data.table))

helper_path <- file.path(v2_root, "R", "FuncFigure2Helper.R")
if (!file.exists(helper_path)) {
  stop(sprintf("Missing helper file: %s", helper_path))
}

source(helper_path)

stopifnot(exists("chooseFigure2EffectSize", mode = "function"))
stopifnot(exists("computeMinTargetDistance", mode = "function"))
stopifnot(exists("getFigure2Paths", mode = "function"))
stopifnot(exists("labelFigure2Stage", mode = "function"))

ad_dt <- data.table(effect_size_pdcpdx = 0.4, effect_size_cell = 0.2)
ad_res <- chooseFigure2EffectSize(ad_dt, "selected_genes_ad_filtered")
stopifnot(identical(ad_res$effect_source[[1]], "effect_size_pdcpdx"))
stopifnot(isTRUE(all.equal(ad_res$effect_size[[1]], 0.4)))

clinical_dt <- data.table(effect_size = -0.25)
clinical_res <- chooseFigure2EffectSize(clinical_dt, "clinical_sig_mRNA")
stopifnot(identical(clinical_res$effect_source[[1]], "effect_size"))
stopifnot(isTRUE(all.equal(clinical_res$effect_size[[1]], -0.25)))

final_dt <- data.table(
  effect_size_ctrdb = c(NA_real_),
  effect_size_pdcpdx_invitro = 0.3,
  effect_size_cell_invitro = 0.1
)
final_res <- chooseFigure2EffectSize(final_dt, "final_biomarkers")
stopifnot(identical(final_res$effect_source[[1]], "effect_size_pdcpdx_invitro"))
stopifnot(isTRUE(all.equal(final_res$effect_size[[1]], 0.3)))

stage_labels <- labelFigure2Stage(c("selected_genes_ad_filtered", "clinical_sig_mRNA", "final_biomarkers"))
stopifnot(identical(stage_labels, c("TCGA-AD filtered", "Clinical significant", "Clinical validated")))

edge_dt <- data.table(from = c("A", "B", "C"), to = c("B", "C", "D"))
dist_one <- computeMinTargetDistance(gene = "A", target_genes = c("D", "X"), edge_dt = edge_dt)
stopifnot(isTRUE(all.equal(dist_one$distance, 3)))
stopifnot(identical(dist_one$target_used, "D"))
stopifnot(isTRUE(dist_one$has_gene_mapping))
stopifnot(isTRUE(dist_one$has_target_mapping))

dist_missing_gene <- computeMinTargetDistance(gene = "Z", target_genes = c("D"), edge_dt = edge_dt)
stopifnot(is.na(dist_missing_gene$distance))
stopifnot(!dist_missing_gene$has_gene_mapping)

dist_missing_target <- computeMinTargetDistance(gene = "A", target_genes = c("Q"), edge_dt = edge_dt)
stopifnot(is.na(dist_missing_target$distance))
stopifnot(!dist_missing_target$has_target_mapping)

paths <- getFigure2Paths()
expected_output_dir <- normalizePath(
  "/Users/peng/Desktop/Project/DROMA/Meta_project/Meta_Example/Output/visualization/figure2",
  mustWork = FALSE
)
stopifnot(identical(paths$output_dir, expected_output_dir))
stopifnot(identical(paths$fig2a_pdf, file.path(expected_output_dir, "figure2A_distance_distribution.pdf")))

tmp_output_dir <- file.path(tempdir(), "figure2-test-output")
Sys.setenv(DROMA_FIGURE2_OUTPUT_DIR = tmp_output_dir)
override_paths <- getFigure2Paths()
stopifnot(identical(override_paths$output_dir, normalizePath(tmp_output_dir, mustWork = FALSE)))
Sys.unsetenv("DROMA_FIGURE2_OUTPUT_DIR")

drug_target_dt <- fread(file.path(v2_root, "data", "drug_target_map.csv"))
epi_rows <- drug_target_dt[drug == "4'-Epiadriamycin"]
stopifnot(nrow(epi_rows) > 0L)
stopifnot(all(c("TOP2A", "TOP2B") %in% epi_rows$target_gene))

cat("v2 figure2 helper tests passed\n")
