# ============================================================================
# 26-Plot_Figure2C_Distance_vs_ValidationStage.R
# Plot Figure 2C distance across validation stages
# ============================================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ggplot2))

root_dir <- normalizePath(file.path(getwd(), "v2_part"), mustWork = TRUE)
source(file.path(root_dir, "R", "FuncFigure2Helper.R"))
source(file.path(root_dir, "R", "FuncFigure2Plot.R"))

paths <- getFigure2Paths()

cat("\n=== 26: Plot Figure 2C ===\n")

distance_dt <- fread(paths$distance_results_csv)
p <- plotFigure2CDistanceVsStage(distance_dt)
ggsave(paths$fig2c_pdf, plot = p, width = 8, height = 6, device = "pdf", useDingbats = FALSE)

cat(sprintf("  OK %s\n", basename(paths$fig2c_pdf)))
cat("\nDone.\n")
