# Volcano plot visualization ----

#' Plot Volcano for One Drug-Tumor Batch Result
#'
#' @description Creates a publication-quality volcano plot (effect size vs
#' \code{-log10(p-value)}) from a single batch result table. Points are
#' coloured by direction (Up / Down / NS) using Tableau 10 colorblind-safe
#' palette. Significant regions have a subtle shaded background.
#' @param batch_df A \code{data.table} or \code{data.frame} with columns
#'   \code{effect_size}, \code{p_value} (or \code{q_value}), \code{name},
#'   and optionally \code{n_datasets}.
#' @param es_t Effect-size threshold. Default 0.1.
#' @param P_t P-value threshold. Default 0.05.
#' @param use_p_value Logical; use \code{p_value} instead of \code{q_value}.
#' @param n_datasets_t Minimum dataset count threshold. NULL to skip.
#' @param label Logical; label top genes with ggrepel.
#' @param top_label_each Number of top genes to label per direction.
#' @param custom_labels Optional character vector of gene names to force-label.
#' @param title Plot title.
#' @param point_size Base point size.
#' @param point_alpha Point alpha.
#' @return A ggplot2 object.
#' @export
plotMetaVolcano <- function(batch_df,
                            es_t = 0.1,
                            P_t = 0.05,
                            use_p_value = FALSE,
                            n_datasets_t = NULL,
                            label = TRUE,
                            top_label_each = 5,
                            custom_labels = NULL,
                            title = NULL,
                            point_size = 2,
                            point_alpha = 0.55) {
  batch_df <- as.data.frame(batch_df)

  p_val_col <- if (use_p_value) "p_value" else "q_value"
  required <- c("effect_size", "name", p_val_col)
  if (!all(required %in% colnames(batch_df))) {
    stop("batch_df must contain columns: ", paste(required, collapse = ", "), call. = FALSE)
  }

  dir_colors <- getMetaVisColors("direction")

  # Classify points
  batch_df$group <- ifelse(
    batch_df$effect_size > es_t & batch_df[[p_val_col]] < P_t, "Up",
    ifelse(
      batch_df$effect_size < -es_t & batch_df[[p_val_col]] < P_t, "Down", "NS"
    )
  )
  if (!is.null(n_datasets_t) && "n_datasets" %in% colnames(batch_df)) {
    batch_df$group <- ifelse(
      batch_df$group != "NS" & batch_df$n_datasets >= n_datasets_t,
      batch_df$group, "NS"
    )
  )

  # Summary text
  n_up   <- sum(batch_df$group == "Up", na.rm = TRUE)
  n_down <- sum(batch_df$group == "Down", na.rm = TRUE)
  n_total <- nrow(batch_df)
  sig_text <- sprintf("Up: %d   Down: %d   Total: %d", n_up, n_down, n_total)

  y_max <- max(-log10(batch_df[[p_val_col]]), na.rm = TRUE)
  if (!is.finite(y_max)) y_max <- 1

  # Build plot
  p <- ggplot2::ggplot(
    batch_df,
    ggplot2::aes(x = effect_size, y = -log10(.data[[p_val_col]]))
  ) +
    # Shaded significant quadrants
    ggplot2::annotate("rect",
      xmin = es_t, xmax = Inf, ymin = -log10(P_t), ymax = y_max * 1.05,
      alpha = 0.06, fill = dir_colors[["Up"]]
    ) +
    ggplot2::annotate("rect",
      xmin = -Inf, xmax = -es_t, ymin = -log10(P_t), ymax = y_max * 1.05,
      alpha = 0.06, fill = dir_colors[["Down"]]
    ) +
    # Points
    ggplot2::geom_point(
      ggplot2::aes(color = group),
      alpha = point_alpha, size = point_size
    ) +
    ggplot2::scale_color_manual(
      values = dir_colors,
      breaks = c("Up", "Down"),
      drop = FALSE,
      name = "Direction"
    ) +
    # Threshold lines
    ggplot2::geom_vline(xintercept = c(-es_t, es_t),
      linetype = "dashed", linewidth = 0.4, color = "grey40"
    ) +
    ggplot2::geom_hline(yintercept = -log10(P_t),
      linetype = "dashed", linewidth = 0.4, color = "grey40"
    ) +
    # Summary annotation
    ggplot2::annotate("text",
      x = 0, y = y_max * 0.95,
      label = sig_text, hjust = 0.5, size = 3.5, fontface = "bold", color = "grey30"
    ) +
    # Direction labels on sides
    ggplot2::annotate("text",
      x = max(batch_df$effect_size * 0.7, es_t * 2),
      y = y_max * 0.88,
      label = sprintf("%d up", n_up), hjust = 0.5, size = 3.2,
      color = dir_colors[["Up"]], fontface = "italic"
    ) +
    ggplot2::annotate("text",
      x = min(batch_df$effect_size * 0.7, -es_t * 2),
      y = y_max * 0.88,
      label = sprintf("%d down", n_down), hjust = 0.5, size = 3.2,
      color = dir_colors[["Down"]], fontface = "italic"
    ) +
    ggplot2::labs(
      x = "Effect size",
      y = if (use_p_value) expression(-log[10](italic(P))) else expression(-log[10](italic(Q))),
      title = title
    ) +
    themeMetaPaper()

  # Label top genes
  if (label && requireNamespace("ggrepel", quietly = TRUE)) {
    sig_df <- batch_df[batch_df$group != "NS", , drop = FALSE]
    if (nrow(sig_df) > 0) {
      if (!is.null(custom_labels)) {
        lab_df <- sig_df[sig_df$name %in% custom_labels, , drop = FALSE]
      } else {
        up_df <- sig_df[sig_df$group == "Up", , drop = FALSE]
        dn_df <- sig_df[sig_df$group == "Down", , drop = FALSE]
        up_df <- up_df[order(-up_df$effect_size), , drop = FALSE]
        dn_df <- dn_df[order(dn_df$effect_size), , drop = FALSE]
        up_df <- up_df[seq_len(min(top_label_each, nrow(up_df))), , drop = FALSE]
        dn_df <- dn_df[seq_len(min(top_label_each, nrow(dn_df))), , drop = FALSE]
        lab_df <- rbind(up_df, dn_df)
      }
      if (nrow(lab_df) > 0) {
        p <- p + ggrepel::geom_text_repel(
          data        = lab_df,
          ggplot2::aes(label = name),
          size        = 3.2,
          color       = "black",
          segment.color = "grey50",
          segment.size  = 0.25,
          box.padding = 0.4,
          point.padding = 0.2,
          force       = 4,
          max.overlaps = 20,
          min.segment.length = 0.2
        )
      }
    }
  }

  p
}

#' Batch Plot Volcano for All Drug-Tumor Pairs
#'
#' @description Generates one volcano PDF per drug-tumor pair from batch
#' results collected by \code{collectAllBatchResults()}.
#' @param batch_all A \code{data.table} from \code{collectAllBatchResults()}.
#' @param output_dir Directory to write PDFs.
#' @param es_t Effect-size threshold.
#' @param P_t P-value threshold.
#' @param use_p_value Logical.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @return Invisible vector of output file paths.
#' @export
batchPlotVolcano <- function(batch_all,
                             output_dir,
                             es_t = 0.1,
                             P_t = 0.05,
                             use_p_value = FALSE,
                             width = 7.2,
                             height = 5.6) {
  batch_all <- data.table::as.data.table(batch_all)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  pairs <- unique(batch_all[, .(drug, tumor_type)])
  out_files <- character(nrow(pairs))

  for (i in seq_len(nrow(pairs))) {
    pair <- pairs[i]
    sub <- batch_all[drug == pair$drug & tumor_type == pair$tumor_type]

    title <- paste(pair$drug, "-", pair$tumor_type)
    slug  <- paste0(DROMA.Meta::sanitizeName(pair$drug), "_",
                    DROMA.Meta::sanitizeName(pair$tumor_type))

    p <- plotMetaVolcano(
      batch_df    = sub,
      es_t        = es_t,
      P_t         = P_t,
      use_p_value = use_p_value,
      title       = title
    )

    fpath <- file.path(output_dir, paste0("volcano_", slug, ".pdf"))
    saveMetaVisPdf(p, fpath, width = width, height = height)
    out_files[i] <- fpath
    cat(sprintf("  OK %s\n", basename(fpath)))
  }

  invisible(out_files)
}
