# ============================================================================
# 20-Organize_Article_Figures.R
# Reorganize visualization outputs into manuscript-oriented figure folders.
# ============================================================================

suppressPackageStartupMessages(library(data.table))

resolve_env_path <- function(...) {
  candidates <- c(...)
  for (name in candidates) {
    value <- Sys.getenv(name, unset = NA_character_)
    if (!is.na(value) && nzchar(value)) {
      return(normalizePath(value, mustWork = FALSE))
    }
  }
  NA_character_
}

project_root <- resolve_env_path("DROMA_META_PROJECT_ROOT")
vis_input <- resolve_env_path("DROMA_META_VIS_INPUT", "DROMA_META_VIS_OUTPUT", "VIS_OUTPUT")
vis_output <- resolve_env_path("VIS_ARTICLE_OUTPUT", "DROMA_META_VIS_ARTICLE_OUTPUT")

if (is.na(vis_input)) {
  if (!is.na(project_root)) {
    vis_input <- file.path(project_root, "Output", "visualization")
  } else {
    stop("Set VIS_OUTPUT (or DROMA_META_VIS_INPUT) before running 20-Organize_Article_Figures.R", call. = FALSE)
  }
}
if (is.na(vis_output)) {
  vis_output <- file.path(dirname(vis_input), paste0(basename(vis_input), "2"))
}

if (!dir.exists(vis_input)) {
  stop("Visualization input directory does not exist: ", vis_input, call. = FALSE)
}

if (dir.exists(vis_output)) {
  unlink(vis_output, recursive = TRUE)
}
dir.create(vis_output, recursive = TRUE, showWarnings = FALSE)

copy_figure <- function(src_rel, dest_rel, figure, panel, title) {
  src <- file.path(vis_input, src_rel)
  dest <- file.path(vis_output, dest_rel)
  if (!file.exists(src)) {
    warning("Missing source figure: ", src, call. = FALSE)
    return(NULL)
  }
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(src, dest, overwrite = TRUE)
  if (!ok) stop("Failed to copy figure: ", src, call. = FALSE)
  data.table(
    figure = figure,
    panel = panel,
    title = title,
    file = dest_rel,
    source = src_rel
  )
}

write_md <- function(path, lines) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(lines, path, useBytes = TRUE)
}

manifest <- list()

manifest[[length(manifest) + 1L]] <- copy_figure(
  "pipeline_funnel.pdf", "figure1/Figure1A_pipeline_funnel.pdf",
  "Figure 1", "A", "Pipeline-level attrition across evidence layers"
)
manifest[[length(manifest) + 1L]] <- copy_figure(
  "stage_attrition_dumbbell.pdf", "figure1/Figure1B_pairwise_attrition.pdf",
  "Figure 1", "B", "Pair-level attrition from preclinical intersection to final validation"
)
manifest[[length(manifest) + 1L]] <- copy_figure(
  "figure1_panel.pdf", "figure1/Figure1C_integrated_workflow_panel.pdf",
  "Figure 1", "C", "Integrated workflow overview panel"
)
write_md(file.path(vis_output, "figure1", "legend.md"), c(
  "# Figure 1. Overview of the DROMA.Meta biomarker discovery workflow",
  "",
  "**Suggested manuscript legend.** DROMA.Meta integrates multi-model preclinical drug response evidence with tumor expression concordance and clinical validation to prioritize drug-specific biomarkers. **A**, aggregate feature attrition across the analysis workflow, from all tested genes through preclinical intersection, TCGA Anderson-Darling filtering, clinical significance, and final biomarker nomination. **B**, pair-level attrition comparing the number of genes entering the validation stage with the number retained as final biomarkers. **C**, integrated overview panel summarizing workflow scale, an example preclinical association landscape, preclinical-clinical effect-size concordance, and final biomarker direction.",
  "",
  "**Files.**",
  "- `Figure1A_pipeline_funnel.pdf`",
  "- `Figure1B_pairwise_attrition.pdf`",
  "- `Figure1C_integrated_workflow_panel.pdf`",
  "",
  "**Article role.** Use this figure to introduce the pipeline before presenting biological results."
))

manifest[[length(manifest) + 1L]] <- copy_figure(
  "biomarker_count_by_pair.pdf", "figure2/Figure2A_final_biomarkers_by_pair.pdf",
  "Figure 2", "A", "Final biomarker counts by drug-tumor pair"
)
manifest[[length(manifest) + 1L]] <- copy_figure(
  "direction_summary.pdf", "figure2/Figure2B_biomarker_direction_summary.pdf",
  "Figure 2", "B", "Direction distribution of final biomarkers"
)
manifest[[length(manifest) + 1L]] <- copy_figure(
  "biomarker_heatmap.pdf", "figure2/Figure2C_final_biomarker_effect_heatmap.pdf",
  "Figure 2", "C", "Clinical effect-size heatmap for final biomarkers"
)
write_md(file.path(vis_output, "figure2", "legend.md"), c(
  "# Figure 2. Global landscape of final biomarkers across drug-tumor pairs",
  "",
  "**Suggested manuscript legend.** Final biomarkers were summarized across all analyzed drug-tumor pairs to identify where the workflow produced clinically supported candidates. **A**, number of final biomarkers detected for each drug-tumor pair. **B**, distribution of biomarker direction, with up- and down-associated markers shown separately. **C**, heatmap of clinical effect sizes for final biomarkers across drug-tumor contexts; grey cells indicate drug-tumor combinations without a measured final biomarker entry for that gene.",
  "",
  "**Files.**",
  "- `Figure2A_final_biomarkers_by_pair.pdf`",
  "- `Figure2B_biomarker_direction_summary.pdf`",
  "- `Figure2C_final_biomarker_effect_heatmap.pdf`",
  "",
  "**Article role.** Use this figure to move from method overview to the main discovery landscape."
))

manifest[[length(manifest) + 1L]] <- copy_figure(
  "concordance_scatter.pdf", "figure3/Figure3A_preclinical_clinical_concordance.pdf",
  "Figure 3", "A", "Preclinical-clinical effect-size concordance"
)
manifest[[length(manifest) + 1L]] <- copy_figure(
  "biomarker_forest.pdf", "figure3/Figure3B_stagewise_effect_forest.pdf",
  "Figure 3", "B", "Stage-wise effect sizes for top biomarkers"
)
write_md(file.path(vis_output, "figure3", "legend.md"), c(
  "# Figure 3. Cross-evidence support for prioritized biomarkers",
  "",
  "**Suggested manuscript legend.** Final biomarkers were evaluated for consistency between preclinical models and clinical response datasets. **A**, scatter plot comparing preclinical and clinical effect sizes for final biomarkers; points are colored by drug and labeled where space permits. The fitted line summarizes the global trend across biomarkers. **B**, forest-style summary of stage-wise effect sizes for top-ranked biomarkers, showing how evidence is distributed across cell-line, PDC/PDX, and clinical layers.",
  "",
  "**Files.**",
  "- `Figure3A_preclinical_clinical_concordance.pdf`",
  "- `Figure3B_stagewise_effect_forest.pdf`",
  "",
  "**Article role.** Use this figure to support the claim that final biomarkers are not single-dataset hits, but are supported across evidence layers."
))

manifest[[length(manifest) + 1L]] <- copy_figure(
  "tcga_ad_pvalue_histogram.pdf", "figure4/Figure4A_tcga_ad_pvalue_histogram.pdf",
  "Figure 4", "A", "TCGA Anderson-Darling p-value distribution"
)
manifest[[length(manifest) + 1L]] <- copy_figure(
  "tcga_ad_density.pdf", "figure4/Figure4B_tcga_expression_concordance_examples.pdf",
  "Figure 4", "B", "TCGA concordant and non-concordant expression examples"
)
write_md(file.path(vis_output, "figure4", "legend.md"), c(
  "# Figure 4. TCGA expression concordance as an orthogonal filter",
  "",
  "**Suggested manuscript legend.** TCGA tumor expression distributions were used as an orthogonal filter to retain preclinical biomarkers whose expression behavior is compatible with patient tumors. **A**, distribution of Anderson-Darling test p-values across evaluated candidate features; the threshold marks the concordance criterion used in the workflow. **B**, representative density plots contrasting a TCGA-concordant gene with a non-concordant gene, illustrating how the filter separates candidates with tumor-compatible and tumor-incompatible expression distributions.",
  "",
  "**Files.**",
  "- `Figure4A_tcga_ad_pvalue_histogram.pdf`",
  "- `Figure4B_tcga_expression_concordance_examples.pdf`",
  "",
  "**Article role.** Use this figure to justify the TCGA-AD step as more than a technical filter."
))

case_sources <- data.table(
  case_id = c(
    "afatinib_breast_cancer",
    "4_epiadriamycin_breast_cancer",
    "bortezomib_breast_cancer",
    "azacitidine_breast_cancer",
    "afatinib_liver_cancer"
  ),
  label = c(
    "Afatinib-breast cancer",
    "4'-Epiadriamycin-breast cancer",
    "Bortezomib-breast cancer",
    "Azacitidine-breast cancer",
    "Afatinib-liver cancer"
  ),
  volcano = c(
    "volcano/volcano_Afatinib_breast_cancer.pdf",
    "volcano/volcano_4_Epiadriamycin_breast_cancer.pdf",
    "volcano/volcano_Bortezomib_breast_cancer.pdf",
    "volcano/volcano_Azacitidine_breast_cancer.pdf",
    "volcano/volcano_Afatinib_liver_cancer.pdf"
  ),
  upset = c(
    "upset/upset_Afatinib_breast_cancer.pdf",
    "upset/upset_4_Epiadriamycin_breast_cancer.pdf",
    "upset/upset_Bortezomib_breast_cancer.pdf",
    "upset/upset_Azacitidine_breast_cancer.pdf",
    "upset/upset_Afatinib_liver_cancer.pdf"
  )
)

for (i in seq_len(nrow(case_sources))) {
  row <- case_sources[i]
  case_dir <- file.path("figure5_case_studies", row$case_id)
  manifest[[length(manifest) + 1L]] <- copy_figure(
    row$volcano, file.path(case_dir, "Figure5_volcano.pdf"),
    "Figure 5", paste0(LETTERS[i], "1"), paste(row$label, "volcano plot")
  )
  manifest[[length(manifest) + 1L]] <- copy_figure(
    row$upset, file.path(case_dir, "Figure5_upset.pdf"),
    "Figure 5", paste0(LETTERS[i], "2"), paste(row$label, "stage-wise overlap")
  )
  write_md(file.path(vis_output, case_dir, "legend.md"), c(
    paste0("# Figure 5 case study. ", row$label),
    "",
    "**Suggested manuscript legend.** This case study shows the feature-level evidence supporting final biomarker nomination for the selected drug-tumor context. The volcano plot summarizes preclinical association strength across tested genes, with effect size on the x-axis and statistical evidence on the y-axis. The upset plot summarizes overlap among workflow stages, highlighting how candidates pass through preclinical significance, PDC/PDX support, TCGA expression concordance, and clinical validation.",
    "",
    "**Files.**",
    "- `Figure5_volcano.pdf`",
    "- `Figure5_upset.pdf`",
    "",
    "**Article role.** Use this case as a detailed example after the global figures. It can be promoted to the main text or moved to supplementary figures depending on space."
  ))
}

supp_dir <- file.path(vis_output, "supplementary")
dir.create(file.path(supp_dir, "volcano_all"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(supp_dir, "upset_all"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(supp_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

invisible(file.copy(
  list.files(file.path(vis_input, "volcano"), pattern = "[.]pdf$", full.names = TRUE),
  file.path(supp_dir, "volcano_all"),
  overwrite = TRUE
))
invisible(file.copy(
  list.files(file.path(vis_input, "upset"), pattern = "[.]pdf$", full.names = TRUE),
  file.path(supp_dir, "upset_all"),
  overwrite = TRUE
))
invisible(file.copy(
  list.files(vis_input, pattern = "[.]csv$", full.names = TRUE),
  file.path(supp_dir, "tables"),
  overwrite = TRUE
))
write_md(file.path(supp_dir, "legend.md"), c(
  "# Supplementary figures and tables",
  "",
  "**Suggested manuscript legend.** Supplementary figures provide the full batch-level visualization output for all analyzed drug-tumor pairs. Volcano plots show preclinical association results for each pair with sufficient data. Upset plots show stage-wise feature-set overlap for each eligible pair. Supplementary tables provide collected workflow summaries, merged final biomarkers, batch-level association results, TCGA-AD statistics, and significant feature tables.",
  "",
  "**Folders.**",
  "- `volcano_all/`: all generated volcano plots.",
  "- `upset_all/`: all generated upset plots.",
  "- `tables/`: collected CSV tables used by the article figures."
))

manifest_dt <- rbindlist(Filter(Negate(is.null), manifest), fill = TRUE)
fwrite(manifest_dt, file.path(vis_output, "figure_manifest.csv"))

write_md(file.path(vis_output, "index.md"), c(
  "# Manuscript-oriented figure package",
  "",
  "Storyline: pipeline-first.",
  "",
  "Recommended order:",
  "",
  "1. `figure1/`: workflow overview and evidence-layer attrition.",
  "2. `figure2/`: global final biomarker landscape.",
  "3. `figure3/`: cross-evidence support for final biomarkers.",
  "4. `figure4/`: TCGA expression concordance filter.",
  "5. `figure5_case_studies/`: detailed drug-tumor case studies.",
  "6. `supplementary/`: complete batch plots and source tables.",
  "",
  "Each figure folder includes a `legend.md` with manuscript-ready legend text and a short note on how to use the figure in the article.",
  "",
  "See `figure_manifest.csv` for the source-to-output mapping."
))

cat("Article figure package written to:\n", vis_output, "\n", sep = "")
cat("Files organized: ", nrow(manifest_dt), " main/case figure files\n", sep = "")
