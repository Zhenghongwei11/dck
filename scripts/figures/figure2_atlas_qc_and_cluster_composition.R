source("scripts/figures/_setup.R")
source("scripts/figures/_common.R")

draw_figure2_atlas <- function(qc_df, retention_df, cluster_df, cell_scores, umap_df = NULL, qc_cells = NULL) {
  layout(matrix(1:8, nrow = 4, byrow = TRUE))

  # Global ordering: rank samples by retained fraction (descending) for Panels A/B/E/H.
  retention_df$pct_retained <- safe_numeric(retention_df$pct_retained)
  retention_df$raw_cells <- safe_numeric(retention_df$raw_cells)
  retention_df$retained_cells <- safe_numeric(retention_df$retained_cells)
  retention_df$group <- as.character(retention_df$group)
  o <- order(ifelse(is.finite(retention_df$pct_retained), retention_df$pct_retained, -Inf), decreasing = TRUE)
  retention_df <- retention_df[o, , drop = FALSE]

  # Cluster palette (used across Panels D/F/G when UMAP anchors are available).
  cluster_pal <- NULL
  if (!is.null(umap_df) && nrow(umap_df) > 0 && "cluster" %in% names(umap_df)) {
    cluster_pal <- make_cluster_palette(umap_df$cluster)
  }

  # A: raw vs retained
  par(mar = c(6, 8, 3, 4))
  short <- short_sample_label(retention_df$sample)
  barplot(
    rbind(retention_df$raw_cells, retention_df$retained_cells),
    beside = TRUE,
    names.arg = sanitize_ascii(short),
    las = 2,
    cex.names = 0.85,
    col = c("grey80", "grey30"),
    main = "Atlas input vs retained cells",
    ylab = ""
  )
  mtext("Cells", side = 2, line = 4.5, cex = 0.85)
  opar <- par(xpd = FALSE)
  legend("topleft", inset = 0.01, legend = c("Raw", "Retained"), fill = c("grey80", "grey30"), bty = "n", cex = 0.85)
  par(opar)
  panel_label("A")

  # B: retention percentage
  par(mar = c(6, 8, 3, 4))
  labs <- short_sample_label(retention_df$sample)
  y_pos <- seq_len(nrow(retention_df))
  cols <- ifelse(retention_df$group == "diabetes", "#b2182b", "grey60")
  plot(
    retention_df$pct_retained * 100,
    y_pos,
    yaxt = "n",
    pch = 16,
    col = cols,
    xlab = "Retained (%)",
    ylab = "",
    main = "Retention per sample",
    xlim = c(0, 100)
  )
  axis(2, at = y_pos, labels = labs, las = 2, cex.axis = 0.85)
  abline(v = 20, lty = 2, col = "grey60")
  legend("bottomright", legend = c("control", "diabetes"), pch = 16, col = c("grey60", "#b2182b"), bty = "n", cex = 0.8)
  # Place cutoff label at the bottom to avoid hitting sample names
  mtext("Cutoff: 20%", side = 1, line = -1.2, at = 20, adj = -0.1, cex = 0.65, col = "grey40")
  panel_label("B")

  # C: atlas UMAP (by group)
  par(mar = c(4, 8, 3, 4))
  if (!is.null(umap_df) && nrow(umap_df) > 0) {
    col <- ifelse(umap_df$group == "diabetes", "#b2182b", "grey60")
    plot(
      umap_df$UMAP_1,
      umap_df$UMAP_2,
      pch = 16,
      cex = 0.35,
      col = grDevices::adjustcolor(col, alpha.f = 0.35),
      xlab = "UMAP_1",
      ylab = "UMAP_2",
      main = "Atlas UMAP (by group)"
    )
    legend("topright", legend = c("control", "diabetes"), pch = 16, col = c("grey60", "#b2182b"), bty = "n", cex = 0.85)
    n_control <- sum(umap_df$group == "control", na.rm = TRUE)
    n_diabetes <- sum(umap_df$group == "diabetes", na.rm = TRUE)
    usr <- par("usr")
    text(
      usr[2] - 0.02 * (usr[2] - usr[1]),
      usr[3] + 0.02 * (usr[4] - usr[3]),
      labels = paste0("N=", nrow(umap_df), " (control=", n_control, "; diabetes=", n_diabetes, ")"),
      adj = c(1, 0),
      cex = 0.7,
      col = "grey40"
    )
  } else {
    ezhu_placeholder_plot("Atlas UMAP (by group)", c("Missing UMAP anchor table.", "Run: export UMAP anchors from cloud."))
  }
  panel_label("C")

  # D: atlas UMAP (by cluster)
  par(mar = c(4, 8, 3, 4))
  if (!is.null(umap_df) && nrow(umap_df) > 0 && "cluster" %in% names(umap_df)) {
    cl <- as.character(umap_df$cluster)
    cl_levels <- sort(unique(cl))
    pal <- if (is.null(cluster_pal)) grDevices::hcl.colors(length(cl_levels), palette = "Dark 3") else unname(cluster_pal[cl_levels])
    col <- pal[match(cl, cl_levels)]
    plot(
      umap_df$UMAP_1,
      umap_df$UMAP_2,
      pch = 16,
      cex = 0.35,
      col = grDevices::adjustcolor(col, alpha.f = 0.35),
      xlab = "UMAP_1",
      ylab = "UMAP_2",
      main = "Atlas UMAP (by cluster)"
    )
    # Legend: top clusters by cell count (avoid misleading "lowest cluster IDs" selection).
    cl_count <- sort(table(cl), decreasing = TRUE)
    legend_clusters <- names(cl_count)[seq_len(min(8, length(cl_count)))]
    legend_labs <- paste0("C", legend_clusters)
    if ("cluster_label" %in% names(umap_df)) {
      m <- tapply(as.character(umap_df$cluster_label), umap_df$cluster, function(x) {
        x <- x[nzchar(x)]
        if (length(x) == 0) return(NA_character_)
        x[1]
      })
      mapped <- unname(m[legend_clusters])
      mapped[is.na(mapped) | !nzchar(mapped)] <- legend_labs[is.na(mapped) | !nzchar(mapped)]
      legend_labs <- sanitize_ascii(mapped)
    }
    legend("topright", legend = legend_labs, col = pal[match(legend_clusters, cl_levels)], pch = 16, bty = "n", cex = 0.65)
  } else {
    ezhu_placeholder_plot("Atlas UMAP (by cluster)", c("Missing cluster labels in UMAP anchor table."))
  }
  par(mar = c(4, 4, 2, 1))
  panel_label("D")

  # E: cell-level QC distribution (nFeature_RNA) using anchor tables
  par(mar = c(8, 8, 3, 4))
  if (!is.null(qc_cells) && nrow(qc_cells) > 0) {
    qc_cells$sample <- factor(qc_cells$sample, levels = unique(retention_df$sample))
    labs <- short_sample_label(levels(qc_cells$sample))
    sample_group <- setNames(as.character(retention_df$group), as.character(retention_df$sample))
    samp_levels <- levels(qc_cells$sample)
    cols <- ifelse(sample_group[samp_levels] == "diabetes", "#b2182b", "grey70")
    boxplot(
      qc_cells$nFeature_RNA ~ qc_cells$sample,
      las = 2,
      outline = FALSE,
      col = cols,
      main = "QC: detected genes per cell (nFeature_RNA)",
      ylab = "nFeature_RNA",
      xlab = "",
      xaxt = "n"
    )
    axis(1, at = seq_along(labs), labels = sanitize_ascii(labs), las = 2, cex.axis = 0.75)
  } else {
    ezhu_placeholder_plot("QC distributions (nFeature_RNA)", c("Missing QC anchor table.", "Run: make anchors-qc"))
  }
  panel_label("E")

  # F: cluster composition per sample (top clusters)
  par(mar = c(6, 8, 3, 4))
  if (!is.null(umap_df) && nrow(umap_df) > 0 && all(c("sample", "cluster") %in% names(umap_df))) {
    umap_df$cluster <- as.character(umap_df$cluster)
    tab_sc <- table(sample = umap_df$sample, cluster = umap_df$cluster)
    tot <- colSums(tab_sc)
    top <- names(sort(tot, decreasing = TRUE))[seq_len(min(8, length(tot)))]
    tab_sc <- tab_sc[, top, drop = FALSE]
    # Order samples: control first, then diabetes (if available), otherwise keep input order.
    sample_levels <- unique(retention_df$sample)
    sample_levels <- sample_levels[sample_levels %in% rownames(tab_sc)]
    if (length(sample_levels) > 0) tab_sc <- tab_sc[sample_levels, , drop = FALSE]

    pal_sc <- if (is.null(cluster_pal)) grDevices::hcl.colors(ncol(tab_sc), palette = "Dark 3") else unname(cluster_pal[colnames(tab_sc)])
    barplot(
      t(tab_sc),
      beside = FALSE,
      names.arg = short_sample_label(rownames(tab_sc)),
      las = 2,
      col = pal_sc,
      border = NA,
      main = "Cluster composition per sample (top clusters)",
      ylab = "Cells"
    )
    legend_labs <- paste0("C", colnames(tab_sc))
    if ("cluster_label" %in% names(umap_df)) {
      mapped <- umap_df$cluster_label[match(colnames(tab_sc), umap_df$cluster)]
      mapped <- ifelse(nzchar(mapped), mapped, legend_labs)
      legend_labs <- sanitize_ascii(mapped)
    }
    legend("topright", legend = legend_labs, fill = pal_sc, bty = "n", cex = 0.65, ncol = 1, xpd = TRUE, inset = c(-0.05, 0))
  } else {
    top <- aggregate(list(n_cells = cluster_df$n_cells), by = list(cluster = cluster_df$cluster), FUN = sum)
    top <- top[order(top$n_cells, decreasing = TRUE), , drop = FALSE]
    top <- utils::head(top$cluster, 12)
    sub <- cluster_df[cluster_df$cluster %in% top, , drop = FALSE]
    sub$cluster <- factor(sub$cluster, levels = top)
    mat <- xtabs(n_cells ~ cluster + group, data = sub)
    mat <- as.matrix(mat)
    colnames(mat) <- sanitize_ascii(colnames(mat))
    barplot(
      t(mat),
      beside = TRUE,
      names.arg = paste0("C", rownames(mat)),
      las = 2,
      col = c("grey30", "#b2182b"),
      border = NA,
      main = "Cluster composition by group (top clusters)",
      ylab = "Cells"
    )
  }
  panel_label("F")

  # G: cluster composition by group (%), computed from the same cell set when available
  par(mar = c(6, 8, 3, 4))
  if (!is.null(umap_df) && nrow(umap_df) > 0 && all(c("group", "cluster") %in% names(umap_df))) {
    umap_df$cluster <- as.character(umap_df$cluster)
    tab_g <- table(cluster = umap_df$cluster, group = umap_df$group)
    # Keep same top clusters as Panel F when possible
    if (exists("top", inherits = FALSE) && length(top) > 0) {
      tab_g <- tab_g[rownames(tab_g) %in% top, , drop = FALSE]
      tab_g <- tab_g[top[top %in% rownames(tab_g)], , drop = FALSE]
    }
    tab_g <- as.matrix(tab_g)
    col_order <- intersect(c("control", "diabetes"), colnames(tab_g))
    col_order <- c(col_order, setdiff(colnames(tab_g), col_order))
    tab_g <- tab_g[, col_order, drop = FALSE]
    pct <- tab_g / pmax(rowSums(tab_g), 1)
    barplot(
      t(pct) * 100,
      beside = TRUE,
      names.arg = paste0("C", rownames(pct)),
      las = 2,
      col = c("grey30", "#b2182b")[seq_len(ncol(pct))],
      border = NA,
      main = "Cluster composition by group (%)",
      ylab = "% cells"
    )
    abline(h = 0, col = "black", lwd = 1)
  } else {
    pct <- mat / pmax(rowSums(mat), 1)
    barplot(
      t(pct) * 100,
      beside = TRUE,
      names.arg = paste0("C", rownames(mat)),
      las = 2,
      col = c("grey30", "#b2182b"),
      border = NA,
      main = "Cluster composition by group (%)",
      ylab = "% cells"
    )
    abline(h = 0, col = "black", lwd = 1)
  }
  panel_label("G")

  # H: cell-level QC distribution (percent.mt), flagged at 20% mt as a common QC threshold
  par(mar = c(8, 8, 3, 4))
  if (!is.null(qc_cells) && nrow(qc_cells) > 0) {
    mt <- qc_cells$percent.mt
    if (all(!is.finite(mt)) || all(mt == 0, na.rm = TRUE)) {
      qc_cells$sample <- factor(qc_cells$sample, levels = unique(retention_df$sample))
      labs <- short_sample_label(levels(qc_cells$sample))
      sample_group <- setNames(as.character(retention_df$group), as.character(retention_df$sample))
      samp_levels <- levels(qc_cells$sample)
      cols <- ifelse(sample_group[samp_levels] == "diabetes", "#b2182b", "grey70")
      boxplot(
        qc_cells$nCount_RNA ~ qc_cells$sample,
        las = 2,
        outline = FALSE,
        col = cols,
        main = "QC: UMI counts per cell (nCount_RNA)",
        ylab = "nCount_RNA",
        xlab = "",
        xaxt = "n"
      )
      axis(1, at = seq_along(labs), labels = sanitize_ascii(labs), las = 2, cex.axis = 0.75)
      mtext("percent.mt unavailable (all-zero); run: make anchors-qc", side = 3, line = 0.2, cex = 0.7, col = "grey40")
    } else {
      qc_cells$sample <- factor(qc_cells$sample, levels = unique(retention_df$sample))
      labs <- short_sample_label(levels(qc_cells$sample))
      sample_group <- setNames(as.character(retention_df$group), as.character(retention_df$sample))
      samp_levels <- levels(qc_cells$sample)
      cols <- ifelse(sample_group[samp_levels] == "diabetes", "#b2182b", "grey70")
      boxplot(
        qc_cells$percent.mt ~ qc_cells$sample,
        las = 2,
        outline = FALSE,
        col = cols,
        main = "QC: mitochondrial fraction (percent.mt)",
        ylab = "percent.mt",
        xlab = "",
        xaxt = "n"
      )
      axis(1, at = seq_along(labs), labels = sanitize_ascii(labs), las = 2, cex.axis = 0.75)
      if (max(mt, na.rm = TRUE) > 20) {
        abline(h = 20, lty = 2, col = "grey60")
        usr <- par("usr")
        text(usr[2], 20, labels = "Ref: 20%", adj = c(1, -0.5), cex = 0.65, col = "grey40")
      }
    }
  } else {
    ezhu_placeholder_plot("QC: mitochondrial fraction (percent.mt)", c("Missing QC anchor table.", "Run: make anchors-qc"))
  }
  panel_label("H")
}

run_figure2 <- function() {
  write_manifest_snapshot_anchor()
  qc_df <- read_singlecell_qc()
  cell_scores_repeat <- read_cell_scores("results/scpagwas/repeat_run/cell_scores.tsv")
  retained_counts <- table(cell_scores_repeat$sample)
  retention_df <- write_sample_retention_anchor(qc_df, retained_counts)
  cluster_df <- read_cluster_by_group()
  qc_cells <- NULL
  if (file.exists("results/figures/gse131882_qc_cells.tsv")) {
    qc_cells <- read_qc_cells_anchor("results/figures/gse131882_qc_cells.tsv")
    qc_cells <- qc_cells[qc_cells$cell %in% cell_scores_repeat$cell, , drop = FALSE]
  }
  umap_df <- NULL
  if (file.exists("results/figures/gse131882_umap_cells.tsv")) {
    umap_df0 <- read_umap_cells("results/figures/gse131882_umap_cells.tsv")
    umap_df <- umap_df0[umap_df0$cell %in% cell_scores_repeat$cell, , drop = FALSE]
  }
  # Optional: map cluster IDs to annotated context labels for legends.
  ann_path <- Sys.getenv("EZHU_CLUSTER_ANNOTATIONS", unset = "results/annotation/cluster_annotations_filled.tsv")
  ann <- NULL
  if (!is.null(umap_df) && file.exists(ann_path)) {
    ann <- utils::read.delim(ann_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    ann$cluster <- as.character(ann$cluster)
    ann$proposed_cell_type <- as.character(ann$proposed_cell_type)
    ann$proposed_state <- as.character(ann$proposed_state)
    ctx <- paste0(ann$proposed_cell_type, " / ", ann$proposed_state)
    ctx <- trimws(ctx)
    ctx <- short_context_label(ctx, 35L)
    map <- setNames(paste0("C", ann$cluster, " ", ctx), ann$cluster)
    umap_df$cluster <- as.character(umap_df$cluster)
    umap_df$cluster_label <- map[umap_df$cluster]
    umap_df$cluster_label[!nzchar(umap_df$cluster_label)] <- paste0("C", umap_df$cluster[!nzchar(umap_df$cluster_label)])
  }

  ezhu_write_stage_metadata(
    "07_figures_figure2",
    params = list(
      formats = TARGET_FORMATS,
      inputs = list(
        qc = "results/singlecell/gse131882_counts_qc.tsv",
        cell_scores_repeat = "results/scpagwas/repeat_run/cell_scores.tsv",
        cluster_by_group = "results/scpagwas/repeat_run/cluster_by_group_counts.tsv",
        umap_cells = "results/figures/gse131882_umap_cells.tsv",
        manifest_snapshot = "results/figures/manifest_snapshot.tsv",
        annotation = ann_path
      )
    ),
    seed = seed
  )

  ezhu_save_plot(
    function() draw_figure2_atlas(qc_df, retention_df, cluster_df, cell_scores_repeat, umap_df = umap_df, qc_cells = qc_cells),
    "figure2_atlas_qc_and_cluster_composition",
    height = 11
  )
}

tryCatch(
  run_figure2(),
  error = function(e) {
    ezhu_write_stage_metadata(
      "07_figures_figure2",
      params = list(error = conditionMessage(e)),
      seed = seed
    )
    ezhu_save_plot(
      function() ezhu_placeholder_plot("Figure 2 (atlas QC) (failed)", c(conditionMessage(e))),
      "figure2_atlas_qc_and_cluster_composition",
      height = 11
    )
  }
)