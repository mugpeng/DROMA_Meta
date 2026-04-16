# ============================================================================
# 19-Figure1_Panel.R
# Assemble multi-panel composite figure (Figure 1)
# ============================================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(DROMA.Meta))

defaults   <- getVisWorkflowDefaults()
vis_output <- defaults$vis_output
dir.create(vis_output, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 19: Figure 1 Multi-Panel ===\n")

summary_path <- file.path(vis_output, "collected_summary.csv")
final_path   <- file.path(vis_output, "collected_final_biomarkers.csv")
batch_path   <- file.path(vis_output, "collected_batch_cell.csv")
ad_path      <- file.path(vis_output, "collected_ad_stats.csv")

if (!all(file.exists(c(summary_path, final_path, batch_path)))) {
  stop("Run 10-Collect_Results.R first", call. = FALSE)
}

summary_dt      <- fread(summary_path)
final_biomarkers <- fread(final_path)
batch_cell      <- fread(batch_path)

panels <- list()

# Panel A: Pipeline funnel
panels$A <- plotPipelineFunnel(summary_dt, title = NULL) +
  ggplot2::labs(tag = "A")

# Panel B: Volcano (first drug-tumor pair)
if (nrow(batch_cell) > 0) {
  pairs <- unique(batch_cell[, .(drug, tumor_type)])
  pair <- pairs[1]
  sub <- batch_cell[drug == pair$drug & tumor_type == pair$tumor_type]
  panels$B <- plotMetaVolcano(sub,
    title = paste(pair$drug, "-", pair$tumor_type),
    label = TRUE, top_label_each = 3
  ) +
    ggplot2::labs(tag = "B")
}

# Panel C: Concordance scatter
preclinical_cols <- grep("effect_size.*(invitro|cell|pdcpdx)", names(final_biomarkers), value = TRUE)
if (nrow(final_biomarkers) > 0 && length(preclinical_cols) > 0 && "effect_size_ctrdb" %in% names(final_biomarkers)) {
  preclinical_col <- intersect(c("effect_size_pdcpdx_invitro", "effect_size_cell_invitro"), preclinical_cols)
  if (!length(preclinical_col)) preclinical_col <- preclinical_cols[1]
  panels$C <- plotConcordanceScatter(final_biomarkers,
    preclinical_es_col = preclinical_col,
    clinical_es_col = "effect_size_ctrdb",
    title = NULL
  ) +
    ggplot2::labs(tag = "C")
}

# Panel D: Direction summary
if (nrow(final_biomarkers) > 0 && "direction" %in% names(final_biomarkers)) {
  panels$D <- plotDirectionSummary(final_biomarkers, title = NULL) +
    ggplot2::labs(tag = "D")
}

# Remove NULL panels
panels <- panels[!vapply(panels, is.null, logical(1))]

if (length(panels) > 0) {
  design <- if (length(panels) == 4) "AB\nCD" else NULL
  composite <- wrap_plots(panels, ncol = 2, byrow = TRUE) +
    plot_annotation(
      title = "DROMA.Meta biomarker discovery pipeline",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
      )
    )

  saveMetaVisPdf(composite, file.path(vis_output, "figure1_panel.pdf"),
    width = 14, height = 10
  )
  cat("  OK figure1_panel.pdf\n")
} else {
  cat("  SKIP figure1_panel: no panels could be generated\n")
}

cat("\nDone.\n")
