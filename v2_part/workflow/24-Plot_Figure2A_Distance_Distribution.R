# ============================================================================
# 24-Plot_Figure2A_Distance_Distribution.R
# Plot Figure 2A distance distribution panel
# ============================================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

root_dir <- normalizePath(file.path(getwd(), "v2_part"), mustWork = TRUE)
source(file.path(root_dir, "R", "FuncFigure2Helper.R"))
source(file.path(root_dir, "R", "FuncFigure2Plot.R"))

paths <- getFigure2Paths()

cat("\n=== 24: Plot Figure 2A ===\n")

distance_dt <- fread(paths$distance_results_csv)

p <- plotFigure2ADistanceDistribution(distance_dt)
ggsave(paths$fig2a_pdf, plot = p, width = 10, height = 6.5, device = "pdf", useDingbats = FALSE)

cat(sprintf("  OK %s\n", basename(paths$fig2a_pdf)))
cat("\nDone.\n")
