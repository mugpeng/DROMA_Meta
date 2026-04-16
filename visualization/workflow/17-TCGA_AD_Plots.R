# ============================================================================
# 17-TCGA_AD_Plots.R
# Generate TCGA Anderson-Darling concordance visualizations
# ============================================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(DROMA.Meta))

project_root <- file.path(normalizePath(getwd(), mustWork = TRUE), "Meta_Example")
defaults   <- getVisWorkflowDefaults(project_root = project_root)
vis_output <- defaults$vis_output
dir.create(vis_output, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 17: TCGA AD Plots ===\n")

ad_path <- file.path(vis_output, "collected_ad_stats.csv")
if (!file.exists(ad_path)) {
  stop("Run 10-Collect_Results.R first", call. = FALSE)
}

ad_stats <- fread(ad_path)

if (nrow(ad_stats) == 0) {
  cat("  SKIP TCGA AD: no AD statistics collected\n")
} else {
  # 1) AD p-value histogram
  p_hist <- plotTcgaAdPvalueHistogram(ad_stats, p_t = 0.01)
  saveMetaVisPdf(p_hist, file.path(vis_output, "tcga_ad_pvalue_histogram.pdf"),
    width = 7.2, height = 5
  )
  cat("  OK tcga_ad_pvalue_histogram.pdf\n")

  # 2) AD statistic density (concordant vs non-concordant)
  p_density <- plotTcgaAdDensity(ad_stats, output_base = defaults$output_base)
  saveMetaVisPdf(p_density, file.path(vis_output, "tcga_ad_density.pdf"),
    width = 7.2, height = 5
  )
  cat("  OK tcga_ad_density.pdf\n")
}

cat("\nDone.\n")
