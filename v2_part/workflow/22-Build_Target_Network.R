# ============================================================================
# 22-Build_Target_Network.R
# Prepare drug-target mappings and Reactome-derived network resources
# ============================================================================

suppressPackageStartupMessages(library(data.table))

root_dir <- normalizePath(file.path(getwd(), "v2_part"), mustWork = TRUE)
source(file.path(root_dir, "R", "FuncFigure2Helper.R"))

paths <- getFigure2Paths()

cat("\n=== 22: Build Target Network ===\n")

feature_dt <- if (file.exists(paths$input_features_csv)) {
  fread(paths$input_features_csv)
} else {
  collectFigure2Inputs()
}

if (!file.exists(paths$drug_target_csv)) {
  cat("  drug_target_map.csv not found; attempting ChEMBL fetch\n")
  target_dt <- tryCatch(
    fetchChemblTargets(unique(feature_dt$drug)),
    error = function(e) {
      message("  ChEMBL fetch failed: ", conditionMessage(e))
      data.table()
    }
  )
  if (!nrow(target_dt)) {
    stop(
      "No drug target map available. Add v2_part/data/drug_target_map.csv or rerun with network access.",
      call. = FALSE
    )
  }
  fwrite(unique(target_dt), paths$drug_target_csv)
} else {
  target_dt <- fread(paths$drug_target_csv)
}

if (!all(c("drug", "target_gene") %in% names(target_dt))) {
  stop("drug_target_map.csv must contain drug and target_gene columns", call. = FALSE)
}
target_dt <- unique(target_dt[!is.na(target_gene) & nzchar(target_gene)])
fwrite(target_dt, paths$drug_target_csv)

if (!file.exists(paths$reactome_gmt)) {
  if (!file.exists(paths$reactome_gmt_zip)) {
    cat("  Reactome GMT not found; attempting download\n")
    downloadFigure2File(
      "https://reactome.org/download/current/ReactomePathways.gmt.zip",
      paths$reactome_gmt_zip
    )
  }
  maybeUnzipReactomeGmt(paths$reactome_gmt_zip, paths$reactome_gmt)
}

pathway_gene_dt <- parseReactomeGmt(paths$reactome_gmt)
edge_dt <- buildReactomeEdges(pathway_gene_dt)
node_dt <- unique(rbindlist(list(
  data.table(gene = edge_dt$from),
  data.table(gene = edge_dt$to)
)))

fwrite(edge_dt, paths$reactome_edges_csv)
fwrite(node_dt, paths$reactome_nodes_csv)

cat(sprintf("  OK drug_target_map: %d rows\n", nrow(target_dt)))
cat(sprintf("  OK reactome_edges: %d edges\n", nrow(edge_dt)))
cat(sprintf("  OK reactome_nodes: %d genes\n", nrow(node_dt)))
cat("\nDone.\n")
