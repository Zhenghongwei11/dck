source("scripts/figures/_setup.R")
source("scripts/figures/_common.R")

draw_marker_panel <- function(markers_path, title, show_legend = FALSE) {
  df <- utils::read.csv(markers_path, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- tolower(names(df))
  df$avg_log2fc <- safe_numeric(df$avg_log2fc)
  df$pct.1 <- safe_numeric(df$pct.1)
  df$pct.2 <- safe_numeric(df$pct.2)
  df <- df[is.finite(df$avg_log2fc) & nzchar(df$gene), , drop = FALSE]
  df <- df[order(df$avg_log2fc, decreasing = TRUE), , drop = FALSE]
  df <- utils::head(df, 12)

  df$gene <- sanitize_ascii(df$gene)
  df$pct.1[!is.finite(df$pct.1)] <- 0
  df$pct.2[!is.finite(df$pct.2)] <- 0
  df$delta_pct <- df$pct.1 - df$pct.2

  par(mar = c(4, 10, 3, 4))
  y <- rev(seq_len(nrow(df)))
  x <- df$avg_log2fc
  xlim <- range(c(0, x), finite = TRUE)
  plot(NA, xlim = xlim, ylim = c(0.5, length(y) + 0.5), yaxt = "n",
       xlab = "avg_log2FC", ylab = "", main = title, frame.plot = FALSE)
  axis(2, at = y, labels = df$gene, las = 1, cex.axis = 0.85)
  ezhu_light_grid(nx = NULL, ny = NA)

  size <- 0.8 + 2.2 * pmin(pmax(df$pct.1, 0), 1)
  pal <- okabe_ito_palette()
  cols <- ifelse(df$delta_pct >= 0, pal["orange"], pal["blue"])
  segments(0, y, x, y, col = "grey75", lwd = 2)
  points(x, y, pch = 21, bg = cols, col = "grey40", cex = size, lwd = 0.6)

  if (isTRUE(show_legend)) {
    legend(
      "bottomright",
      inset = 0.01,
      legend = c("Orange: pct.1 >= pct.2", "Blue: pct.1 < pct.2", "Point size ~ pct.1"),
      pt.bg = c(pal["orange"], pal["blue"], NA),
      pch = c(21, 21, NA),
      bty = "n",
      cex = 0.75,
      text.col = c("black", "black", "grey40")
    )
  }
}

draw_feature_umap <- function(umap_cells, expr, title, legend_title = "Expr") {
  par(mar = c(4, 10, 3, 4))
  expr <- safe_numeric(expr)
  ok <- is.finite(umap_cells$UMAP_1) & is.finite(umap_cells$UMAP_2) & is.finite(expr)
  if (sum(ok) < 10) {
    ezhu_placeholder_plot(title, c("Not enough cells with feature values."))
    return(invisible(NULL))
  }

  x <- umap_cells$UMAP_1[ok]
  y <- umap_cells$UMAP_2[ok]
  z <- expr[ok]

  q <- stats::quantile(z, probs = c(0.01, 0.99), na.rm = TRUE)
  z <- pmax(pmin(z, q[2]), q[1])
  pal <- ezhu_seq_pal(64)
  idx <- as.integer(1 + 63 * (z - min(z)) / max(diff(range(z)), 1e-9))
  idx <- pmax(1, pmin(64, idx))

  plot(
    x, y,
    pch = 16,
    cex = 0.4,
    col = grDevices::adjustcolor(pal[idx], alpha.f = 0.7),
    xlab = "UMAP_1",
    ylab = "UMAP_2",
    main = title
  )

  # Color key
  usr <- par("usr")
  xleft <- usr[2] + 0.02 * (usr[2] - usr[1])
  xright <- usr[2] + 0.06 * (usr[2] - usr[1])
  ybottom <- usr[3]
  ytop <- usr[4]
  yy <- seq(ybottom, ytop, length.out = 64)
  rect(xleft, yy[-length(yy)], xright, yy[-1], col = pal, border = NA, xpd = TRUE)
  mtext(legend_title, side = 4, line = 2.2, cex = 0.8)
}

draw_figure4_markers <- function(ann_path, markers_dir = "results/annotation/markers") {
  layout(matrix(1:8, nrow = 4, byrow = TRUE))
  # Shift panel labels left to avoid overlaps with long titles.
  label_shift <- 2.0

  # Select clusters dynamically from the repeat-run scPagwas ranking, so figure panels
  # remain valid after reclustering / integration changes.
  pv_path <- Sys.getenv("EZHU_PRIMARY_REPEAT_PVALUE", unset = "results/scpagwas/repeat_run/merged_celltype_pvalue.csv")
  pv <- read_pvalue_table(pv_path)
  pv <- pv[order(pv$pvalue, decreasing = FALSE), , drop = FALSE]
  clusters <- as.character(utils::head(pv$celltype, 6L))
  clusters <- clusters[nzchar(clusters)]
  if (length(clusters) == 0) clusters <- as.character(utils::head(pv$celltype, 3L))

  # A–C: top markers (selected clusters)
  marker_clusters <- utils::head(clusters, 3L)

  ann <- NULL
  if (file.exists(ann_path)) {
    ann <- utils::read.delim(ann_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    ann$cluster <- as.character(ann$cluster)
    ann$proposed_cell_type <- as.character(ann$proposed_cell_type)
    ann$proposed_state <- as.character(ann$proposed_state)
  }

  for (i in seq_along(marker_clusters)) {
    cl <- marker_clusters[i]
    path <- file.path(markers_dir, paste0("markers_cluster_", cl, ".csv"))
    label <- paste0("Cluster ", cl)
    if (!is.null(ann) && nrow(ann) > 0) {
      row <- ann[ann$cluster == cl, , drop = FALSE]
      if (nrow(row) > 0) {
        ctx <- trimws(paste0(row$proposed_cell_type[1] %||% "", " / ", row$proposed_state[1] %||% ""))
        if (nzchar(ctx)) label <- paste0(label, " (", ctx, ")")
      }
    }
    if (file.exists(path)) {
      draw_marker_panel(path, paste0("Top markers (", label, ")"), show_legend = (i == 1))
    } else {
      ezhu_placeholder_plot(paste0("Top markers (", label, ")"), c(paste0("Missing: ", path)))
    }
    panel_label(LETTERS[i], x_shift_lines = label_shift)
  }

  # D–F: UMAP feature plots from anchor tables (no Seurat dependency)
  umap_cells <- NULL
  umap_feat <- NULL
  if (file.exists("results/figures/gse131882_umap_cells.tsv") && file.exists("results/figures/gse131882_umap_features.tsv")) {
    umap_cells0 <- read_umap_cells("results/figures/gse131882_umap_cells.tsv")
    umap_feat0 <- read_umap_features("results/figures/gse131882_umap_features.tsv")
    # Filter to the scPagwas cell set to avoid mixing in excluded samples.
    cs <- read_cell_scores("results/scpagwas/repeat_run/cell_scores.tsv")
    keep <- cs$cell
    umap_cells <- umap_cells0[umap_cells0$cell %in% keep, , drop = FALSE]
    umap_feat <- umap_feat0[umap_feat0$cell %in% keep, , drop = FALSE]
    umap_cells <- merge(umap_cells, umap_feat, by = "cell", all.x = TRUE)
  }

  feat_panels <- list(
    list(letter = "D", gene = "TNFRSF12A", title = "UMAP feature: TNFRSF12A"),
    list(letter = "E", gene = "SLC12A1", title = "UMAP feature: SLC12A1"),
    list(letter = "F", gene = "AQP2", title = "UMAP feature: AQP2")
  )
  for (k in seq_along(feat_panels)) {
    p <- feat_panels[[k]]
    if (!is.null(umap_cells) && nrow(umap_cells) > 0 && p$gene %in% names(umap_cells)) {
      draw_feature_umap(umap_cells, umap_cells[[p$gene]], p$title, legend_title = p$gene)
    } else {
      ezhu_placeholder_plot(p$title, c("Missing UMAP anchor tables or feature column.", paste0("Expected: ", p$gene)))
    }
    panel_label(p$letter, x_shift_lines = label_shift)
  }

  # G: canonical marker support heatmap (presence in top markers)
  par(mar = c(6, 10, 3, 6)) # Increased right margin for legend
  marker_sets <- list(
    `PT` = c("SLC34A1", "LRP2", "AQP1", "SLC5A2"),
    `TAL` = c("SLC12A1", "UMOD", "CLDN16", "KCNJ1"),
    `CD-PC` = c("AQP2", "AQP3", "SCNN1A", "FXYD4", "AVPR2"),
    `CD-IC` = c("ATP6V1B1", "ATP6V0A4", "SLC4A1", "FOXI1", "AQP6"),
    `Stress` = c("HAVCR1", "KRT8", "KRT18", "VIM", "JUN", "FOS", "GDF15", "HSPA1A"),
    `Endo` = c("PECAM1", "KDR", "EMCN", "VWF")
  )

  cats <- names(marker_sets)
  score <- matrix(0, nrow = length(clusters), ncol = length(cats), dimnames = list(paste0("C", clusters), cats))
  nhit <- matrix(0, nrow = length(clusters), ncol = length(cats), dimnames = list(paste0("C", clusters), cats))
  for (i in seq_along(clusters)) {
    cl <- clusters[i]
    path <- file.path(markers_dir, paste0("markers_cluster_", cl, ".csv"))
    if (!file.exists(path)) next
    df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
    names(df) <- tolower(names(df))
    df$gene <- as.character(df$gene)
    df$avg_log2fc <- safe_numeric(df$avg_log2fc)
    df <- df[nzchar(df$gene) & is.finite(df$avg_log2fc), , drop = FALSE]
    for (j in seq_along(cats)) {
      hits <- intersect(marker_sets[[cats[j]]], df$gene)
      nhit[i, j] <- length(hits)
      score[i, j] <- if (length(hits) > 0) max(df$avg_log2fc[df$gene %in% hits], na.rm = TRUE) else 0
    }
  }

  pal <- ezhu_div_pal(64)
  zlim <- range(score, finite = TRUE)
  if (!all(is.finite(zlim)) || diff(zlim) == 0) zlim <- c(0, 1)
  image(seq_len(ncol(score)), seq_len(nrow(score)), t(score), axes = FALSE, col = pal, zlim = zlim,
        main = "Canonical marker support (max avg_log2FC)", xlab = "", ylab = "")
  axis(1, at = seq_len(ncol(score)), labels = colnames(score), las = 2, cex.axis = 0.9)
  axis(2, at = seq_len(nrow(score)), labels = rownames(score), las = 2, cex.axis = 0.9)
  box()
  pts <- expand.grid(x = seq_len(ncol(score)), y = seq_len(nrow(score)))
  size_vals <- as.vector(nhit)
  cex <- ifelse(size_vals > 0, 0.8 + 0.8 * pmin(size_vals, 3), 0.2)
  points(pts$x, pts$y, pch = 21, bg = grDevices::adjustcolor("black", 0.15), col = "grey40", cex = cex, lwd = 0.5)
  
  # Add size legend in the right margin
  usr <- par("usr")
  x_legend <- usr[2] + 0.15 * (usr[2] - usr[1])
  y_legend <- usr[3] + 0.2 * (usr[4] - usr[3])
  legend(x = x_legend, y = y_legend, legend = c("1 hit", "2 hits", "3+ hits"), 
         pt.cex = c(1.6, 2.4, 3.2), pch = 21, bg = NA, box.col = NA, cex = 0.6, title = "Hits", xpd = TRUE, xjust = 0)

  # Add colorbar
  xleft <- usr[2] + 0.05 * (usr[2] - usr[1])
  xright <- usr[2] + 0.1 * (usr[2] - usr[1])
  y_seq <- seq(usr[3], usr[4], length.out = 65)
  for (idx_c in 1:64) {
    rect(xleft, y_seq[idx_c], xright, y_seq[idx_c+1], col = pal[idx_c], border = NA, xpd = TRUE)
  }
  text(xright, c(usr[3], usr[4]), labels = round(zlim, 1), pos = 4, cex = 0.6, xpd = TRUE)
  panel_label("G", x_shift_lines = label_shift)

  # H: annotation summary (repeat-run p-values)
  par(mar = c(6, 12, 3, 4)) # Increased left margin
  if (file.exists(ann_path)) {
    ann <- utils::read.delim(ann_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    ann$cluster <- as.character(ann$cluster)
    ann$repeat_run_pvalue <- safe_numeric(ann$repeat_run_pvalue)
    ann <- ann[ann$cluster %in% clusters & is.finite(ann$repeat_run_pvalue), , drop = FALSE]
    ann <- ann[order(ann$repeat_run_pvalue, decreasing = FALSE), , drop = FALSE]
    
    vals <- ann$repeat_run_pvalue
    score <- -log10(pmax(vals, 1e-300))
    # Use context labels with better truncation handling
    labels <- paste0("C", ann$cluster, " ", short_context_label(ann$proposed_cell_type, 25L))
    
    pal2 <- ezhu_seq_pal(64)
    if (!all(is.finite(score))) score[!is.finite(score)] <- 0
    rng <- range(score, finite = TRUE)
    if (!is.finite(rng[1]) || diff(rng) == 0) rng <- c(0, 1)
    
    # Enhanced color mapping: non-linear to show more detail in significant range
    idx <- as.integer(1 + 63 * ((score - rng[1]) / diff(rng)))
    idx <- pmax(1, pmin(64, idx))
    cols <- pal2[idx]
    
    # Increase xlim significantly to ensure p-value text fits on the right
    x_max <- max(score, na.rm=TRUE)
    b <- barplot(
      rev(score),
      names.arg = rev(labels),
      horiz = TRUE,
      las = 1,
      col = rev(cols),
      border = "white",
      lwd = 0.5,
      xlab = expression(-log[10](italic(p))~"(repeat run)"),
      main = "Top annotated contexts",
      cex.names = 0.8,
      xlim = c(0, x_max * 1.6) 
    )
    # Add actual p-value text with a subtle offset to avoid overlap
    p_labels <- format_p_compact(rev(vals))
    text(x = rev(score), y = b, labels = p_labels, pos = 4, cex = 0.75, font = 3, col = "black")
    abline(v = -log10(0.05), lty = 2, col = "grey60")
  } else {
    ezhu_placeholder_plot("Annotated contexts ranked by repeat-run p-value", c(paste0("Missing: ", ann_path)))
  }
  panel_label("H", x_shift_lines = label_shift)
}

ann_path <- Sys.getenv("EZHU_CLUSTER_ANNOTATIONS", unset = "results/annotation/cluster_annotations_filled.tsv")

ezhu_write_stage_metadata(
  "07_figures_figure4",
  params = list(
    formats = TARGET_FORMATS,
    inputs = list(
      markers_dir = "results/annotation/markers",
      annotation = ann_path
    )
  ),
  seed = seed
)

ezhu_save_plot(
  function() draw_figure4_markers(ann_path),
  "figure4_top_markers_text",
  height = 11
)
