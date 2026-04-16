# ============================================================================
# 14-Concordance_Scatter.R
# Generate preclinical vs clinical effect-size concordance scatter plot
# ============================================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DROMA.Meta.visualization))

project_root <- getVisProjectRoot()
defaults   <- getVisWorkflowDefaults(project_root = project_root)
vis_output <- defaults$vis_output
dir.create(vis_output, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 14: Concordance Scatter ===\n")

final_path <- file.path(vis_output, "collected_final_biomarkers.csv")
if (!file.exists(final_path)) {
  stop("Run 10-Collect_Results.R first", call. = FALSE)
}

final_biomarkers <- fread(final_path)

# Detect available effect-size columns
preclinical_cols <- grep("effect_size.*(invitro|cell|pdcpdx)", names(final_biomarkers), value = TRUE)
clinical_col <- "effect_size_ctrdb"

if (nrow(final_biomarkers) == 0) {
  cat("  SKIP concordance: no final biomarkers\n")
} else if (length(preclinical_cols) == 0 || !clinical_col %in% names(final_biomarkers)) {
  cat("  SKIP concordance: required effect-size columns not found\n")
} else {
  # Prefer pdcpdx_invitro, then cell_invitro
  preferred_order <- c("effect_size_pdcpdx_invitro", "effect_size_cell_invitro")
  preclinical_col <- intersect(preferred_order, preclinical_cols)
  if (!length(preclinical_col)) preclinical_col <- preclinical_cols[1]
  preclinical_col <- preclinical_col[1]
  cat(sprintf("  Using preclinical column: %s\n", preclinical_col))

  p <- plotConcordanceScatter(
    final_biomarkers,
    preclinical_es_col = preclinical_col,
    clinical_es_col = clinical_col
  )

  saveMetaVisPdf(p, file.path(vis_output, "concordance_scatter.pdf"),
    width = 7.2, height = 6
  )
  cat("  OK concordance_scatter.pdf\n")
}

cat("\nDone.\n")
