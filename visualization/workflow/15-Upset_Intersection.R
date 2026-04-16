# ============================================================================
# 15-Upset_Intersection.R
# Generate upset plots showing stage-wise gene overlap per drug-tumor pair
# ============================================================================

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(DROMA.Meta))

project_root <- file.path(normalizePath(getwd(), mustWork = TRUE), "Meta_Example")
defaults   <- getVisWorkflowDefaults(project_root = project_root)
output_base <- defaults$output_base
vis_output <- defaults$vis_output
upset_dir  <- file.path(vis_output, "upset")
dir.create(upset_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n=== 15: Upset Intersection ===\n")

if (!dir.exists(output_base)) {
  stop("Output directory does not exist: ", output_base, call. = FALSE)
}

drug_dirs <- list.dirs(output_base, recursive = FALSE, full.names = TRUE)
for (drug_dir in drug_dirs) {
  tumor_dirs <- list.dirs(drug_dir, recursive = FALSE, full.names = TRUE)
  for (pair_dir in tumor_dirs) {
    drug_name  <- basename(drug_dir)
    tumor_name <- basename(pair_dir)
    slug <- paste0(sanitizeName(drug_name), "_", sanitizeName(tumor_name))

    sig_cell   <- file.path(pair_dir, "mRNA_cell_sig.csv")
    sig_pdcpdx <- file.path(pair_dir, "mRNA_pdcpdx_sig.csv")
    if (!file.exists(sig_cell) || !file.exists(sig_pdcpdx)) next

    tryCatch({
      ht <- plotUpsetFromPairDir(pair_dir)
      saveComplexHeatmapPdf(ht, file.path(upset_dir, paste0("upset_", slug, ".pdf")),
        width = 7.2, height = 4.5
      )
      cat(sprintf("  OK upset_%s.pdf\n", slug))
    }, error = function(e) {
      cat(sprintf("  SKIP %s: %s\n", slug, conditionMessage(e)))
    })
  }
}

cat("\nDone.\n")
