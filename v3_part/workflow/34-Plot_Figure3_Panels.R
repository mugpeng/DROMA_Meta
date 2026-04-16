# ============================================================================
# 34-Plot_Figure3_Panels.R
# Assemble legend and article-ready Figure 3 interpretation markdown
# ============================================================================

fallback_if_missing <- function(x, y) {
  if (is.null(x) || length(x) == 0 || !nzchar(x)) y else x
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
script_path <- sub("^--file=", "", fallback_if_missing(file_arg[1], ""))
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path, mustWork = TRUE)) else getwd()
root_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)

suppressPackageStartupMessages(library(data.table))

source(file.path(root_dir, "visualization", "R", "FuncVisHelper.R"))
source(file.path(root_dir, "visualization", "R", "FuncVisBreastSignature.R"))
source(file.path(root_dir, "v3_part", "R", "FuncFigure3Helper.R"))

paths <- getFigure3Paths()
dir.create(paths$output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 34: Plot Figure 3 Panels / Write Markdown ===\n")

signature_matrix_dt <- fread(paths$signature_matrix_csv)
signature_mat <- as.matrix(signature_matrix_dt[, setdiff(names(signature_matrix_dt), "name"), with = FALSE])
rownames(signature_mat) <- signature_matrix_dt$name
signature <- list(matrix = signature_mat)

gene_clusters <- fread(paths$gene_clusters_csv)
drug_clusters <- fread(paths$drug_clusters_csv)
go_dt <- if (file.exists(paths$go_csv)) fread(paths$go_csv) else data.table()
cluster_program_annotation <- fread(paths$cluster_program_csv)
pseudo_atlas_annotation <- fread(paths$pseudo_atlas_csv)
final_overlay <- fread(paths$final_overlay_csv)

writeFigure3Interpretation(
  paths = paths,
  signature = signature,
  gene_clusters = gene_clusters,
  drug_clusters = drug_clusters,
  go_dt = go_dt,
  cluster_program_annotation = cluster_program_annotation,
  pseudo_atlas_annotation = pseudo_atlas_annotation,
  final_overlay = final_overlay
)

writeLines(c(
  "# Figure 3. Shared breast-cancer cell states underlying many biomarkers",
  "",
  "- `Figure3A_breast_shared_signature_heatmap.pdf`: shared AD-filtered gene-by-drug signature in breast cancer.",
  "- `Figure3B_GO_enrichment_by_gene_cluster.pdf`: GO Biological Process terms enriched in each shared cluster.",
  "- `Figure3C_PAM50_marker_overlay.pdf`: PAM50 and related breast-state marker overlay across drugs.",
  "- `Figure3D_cluster_program_annotation.pdf`: cluster-level mapping to luminal, basal-like, HER2-like, proliferative, and EMT-like programs.",
  "- `Figure3E_pseudo_atlas_annotation.pdf`: marker-guided pseudo-atlas annotation of shared clusters to putative breast cell subgroups.",
  "",
  "Scientific use note: Figure 3 is designed to answer whether many biomarkers are organized by a smaller number of shared cell states rather than representing isolated gene-drug effects."
), paths$legend_md, useBytes = TRUE)

cat(sprintf("  OK interpretation markdown: %s\n", paths$interpretation_md))
cat(sprintf("  OK legend markdown: %s\n", paths$legend_md))
