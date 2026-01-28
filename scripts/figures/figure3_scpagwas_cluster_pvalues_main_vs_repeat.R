source("scripts/figures/_setup.R")
source("scripts/figures/_common.R")

plot_top_pvalues <- function(df, title, top_n = 12L, cell_counts = NULL) {
  df <- df[order(df$pvalue, decreasing = FALSE), , drop = FALSE]
  df <- utils::head(df, top_n)
  if (nrow(df) == 0) {
    plot.new()
    title(main = title)
    return(invisible(NULL))
  }
  scores <- -log10(pmax(df$pvalue, 1e-300))
  labels <- paste0("C", df$celltype)
  if (!is.null(cell_counts) && length(cell_counts) > 0) {
    n_cells <- as.integer(cell_counts[as.character(df$celltype)])
    n_cells[is.na(n_cells)] <- 0L
    labels <- paste0(labels, " (n=", n_cells, ")")
  }
  pal <- ezhu_seq_pal(100)
  cols <- pal[cut(scores, breaks = 100)]

  opar <- par(mar = c(4, 8, 3, 4))
  on.exit(par(opar), add = TRUE)

  y_pos <- seq_along(scores)
  plot(NA, xlim = c(0, max(scores) * 1.15), ylim = c(0.5, length(scores) + 0.5),
       yaxt = "n", xlab = expression(-log[10](italic(p))), ylab = "", main = title, frame.plot = FALSE)
  abline(v = pretty(c(0, max(scores))), col = "#F0F0F0", lwd = 1)
  abline(v = -log10(0.05), lty = 2, col = "#C7C7C7")
  segments(0, y_pos, scores, y_pos, col = cols, lwd = 3)
  sig <- rep(FALSE, length(scores))
  if ("fdr" %in% names(df)) sig <- is.finite(df$fdr) & df$fdr <= 0.05
  points(scores, y_pos, pch = 21, bg = cols, col = ifelse(sig, "black", "#6E6E6E"), cex = 1.2, lwd = 0.8)
  axis(2, at = y_pos, labels = labels, las = 1, cex.axis = 0.9)
  if ("fdr" %in% names(df)) {
    legend("topleft", legend = c("FDR <= 0.05", "FDR > 0.05 / NA"), pch = 21, pt.bg = "grey80",
           col = c("black", "#6E6E6E"), bty = "n", cex = 0.75)
  }
  box()
  invisible(NULL)
}

read_annotation_table <- function(path) {
  if (!file.exists(path)) return(NULL)
  ann <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  ann$cluster <- as.character(ann$cluster)
  ann$proposed_cell_type <- as.character(ann$proposed_cell_type)
  ann$proposed_state <- as.character(ann$proposed_state)
  ann$main_run_pvalue <- safe_numeric(ann$main_run_pvalue)
  ann$repeat_run_pvalue <- safe_numeric(ann$repeat_run_pvalue)
  ann
}

context_labels_from_ann <- function(ann) {
  if (is.null(ann) || nrow(ann) == 0) return(character())
  label <- paste0(ann$proposed_cell_type, " / ", ann$proposed_state)
  short_context_label(label, 34)
}

draw_figure3_scpagwas <- function(main_df, repeat_df, ann_path, cell_scores, qc_df) {
  layout(matrix(1:8, nrow = 4, byrow = TRUE))
  cell_counts <- cluster_counts_from_cell_scores(cell_scores)

  # A: main top clusters
  plot_top_pvalues(main_df, "Main run: top clusters", cell_counts = cell_counts)
  panel_label("A")

  # B: repeat top clusters
  plot_top_pvalues(repeat_df, "Repeat run: top clusters", cell_counts = cell_counts)
  panel_label("B")

  # C: concordance
  par(mar = c(4, 8, 3, 4))
  merged <- merge(main_df, repeat_df, by = "celltype", suffixes = c("_main", "_repeat"))
  x <- -log10(pmax(merged$pvalue_main, 1e-300))
  y <- -log10(pmax(merged$pvalue_repeat, 1e-300))
  plot(x, y, pch = 16, col = "grey40", xlab = "-log10(p) main", ylab = "-log10(p) repeat",
       main = "Cluster-level concordance (main vs repeat)")
  abline(0, 1, lty = 2, col = "grey60")
  # Add Spearman Correlation
  cor_val <- stats::cor(x, y, method = "spearman")
  usr <- par("usr")
  text(usr[2] - 0.05*(usr[2]-usr[1]), usr[3] + 0.1*(usr[4]-usr[3]), 
       labels = paste0("Spearman r = ", round(cor_val, 3)), pos = 2, cex = 0.8, font = 2)
  label_idx <- which(pmax(x, y) >= 3)
  if (length(label_idx) > 0) {
    text(x[label_idx], y[label_idx], labels = paste0("C", merged$celltype[label_idx]), pos = 3, cex = 0.7)
  }
  panel_label("C")

  # D: marker-supported contexts (main vs repeat)
  par(mar = c(6, 12, 3, 4))
  ann <- read_annotation_table(ann_path)
  if (!is.null(ann) && nrow(ann) > 0) {
    ann$context <- paste0(ann$proposed_cell_type, " / ", ann$proposed_state)
    ann$context <- trimws(ann$context)
    ann <- ann[nzchar(ann$context) & is.finite(ann$main_run_pvalue) & is.finite(ann$repeat_run_pvalue), , drop = FALSE]
    ann <- ann[order(ann$repeat_run_pvalue, decreasing = FALSE), , drop = FALSE]
    ann <- utils::head(ann, 8)
    if (nrow(ann) > 0) {
      x_main <- -log10(pmax(ann$main_run_pvalue, 1e-300))
      x_rep <- -log10(pmax(ann$repeat_run_pvalue, 1e-300))
      y <- rev(seq_len(nrow(ann)))
      labels <- wrap_label(paste0("C", ann$cluster, "  ", ann$context), 30)
      xlim <- range(c(0, x_main, x_rep), finite = TRUE)
      plot(
        NA,
        xlim = c(0, max(xlim) * 1.05),
        ylim = c(0.5, length(y) + 0.5),
        yaxt = "n",
        xlab = expression(-log[10](italic(p))),
        ylab = "",
        main = "Marker-supported contexts (main vs repeat)",
        frame.plot = FALSE
      )
      axis(2, at = y, labels = labels, las = 1, cex.axis = 0.8)
      ezhu_light_grid(nx = NULL, ny = NA)
      segments(pmin(x_main, x_rep), y, pmax(x_main, x_rep), y, col = "grey70", lwd = 2)
      points(x_main, y, pch = 16, col = "grey40", cex = 0.9)
      points(x_rep, y, pch = 16, col = "#D55E00", cex = 0.9)
      legend("topleft", legend = c("main", "repeat"), pch = 16, col = c("grey40", "#D55E00"), bty = "n", cex = 0.8)
    } else {
      plot.new()
      title(main = "Marker-supported contexts (main vs repeat)")
      mtext("No annotated contexts found.", side = 3, line = -1.2, cex = 0.8, col = "grey40")
    }
  } else {
    plot.new()
    title(main = "Marker-supported contexts (main vs repeat)")
    mtext("Missing annotation table.", side = 3, line = -1.2, cex = 0.8, col = "grey40")
  }
  panel_label("D")

  # E: cell-level p-values
  par(mar = c(7, 8, 3, 4))
  p <- safe_numeric(cell_scores$Random_Correct_BG_p)
  p <- p[is.finite(p) & p >= 0 & p <= 1]
  if (length(p) == 0) {
    plot.new()
    title(main = "Cell-level p-values")
  } else {
    n_total <- length(p)
    n_p0 <- sum(p == 0, na.rm = TRUE)
    n_p1 <- sum(p == 1, na.rm = TRUE)
    n_sig <- sum(p > 0 & p <= 0.05, na.rm = TRUE)
    n_mid <- sum(p > 0.05 & p < 1, na.rm = TRUE)
    n_other <- n_total - n_p0 - n_p1 - n_sig - n_mid
    counts <- c(`p=0` = n_p0, `p<=0.05` = n_sig, `0.05<p<1` = n_mid, `p=1` = n_p1)
    pal <- okabe_ito_palette()
    cols <- c(pal["orange"], pal["grey"], pal["light_grey"], pal["blue"])
    barplot(counts, col = cols, border = NA, main = "Cell-level p-value categories", ylab = "Cells (n)",
            las = 2, cex.names = 0.85)
    usr <- par("usr")
    text(usr[1], usr[4], "p=0: underflow; p=1: null", adj = c(0, 1), cex = 0.75, col = "grey40")
    txt <- paste0("N=", n_total)
    mtext(txt, side = 3, line = 0.2, adj = 1, cex = 0.75, col = "grey30")
  }
  panel_label("E")

  # F: cell-level z-scores
  par(mar = c(4, 8, 3, 4))
  z <- safe_numeric(cell_scores$Random_Correct_BG_z)
  z <- z[is.finite(z)]
  if (length(z) < 2) {
    plot.new()
    title(main = "Cell-level z-scores")
  } else {
    z_cap <- pmax(pmin(z, 10), -10)
    d <- stats::density(z_cap, na.rm = TRUE)
    plot(d$x, d$y, type = "l", lwd = 2, col = "grey30", main = "Cell-level z-scores",
         xlab = "Random_Correct_BG_z (capped at \u00B110)", ylab = "density")
    polygon(d$x, d$y, col = grDevices::adjustcolor("#9ECAE1", 0.4), border = NA)
    lines(d$x, d$y, col = "grey30", lwd = 2)
    abline(v = 0, lty = 2, col = "#BDBDBD")
    # Clarification text inside plot
    text(8, max(d$y)*0.9, "Capped at \u00B110", cex = 0.7, font = 3, col = "grey40")
  }
  panel_label("F")

  # G: TRS vs gPAS
  par(mar = c(4, 8, 3, 4))
  tx <- cell_scores$scPagwas.TRS.Score
  ty <- cell_scores$scPagwas.gPAS.score
  plot(tx, ty, pch = 16, col = rgb(0, 0, 0, 0.12),
       xlab = "TRS score", ylab = "gPAS score", main = "Cell-level score concordance")
  # Add Spearman Correlation
  cor_val_g <- stats::cor(tx, ty, method = "spearman", use = "complete.obs")
  text(max(tx, na.rm = TRUE), max(ty, na.rm = TRUE), labels = paste0("Spearman r = ", round(cor_val_g, 3)), pos = 2, cex = 0.8, font = 2)
  panel_label("G")

  # H: TRS by sample
  par(mar = c(6, 8, 3, 4))
  cell_scores$sample_short <- short_sample_label(cell_scores$sample)
  boxplot(scPagwas.TRS.Score ~ sample_short, data = cell_scores, las = 2, col = "#D9D9D9",
          main = "TRS score by sample", ylab = "TRS score", xlab = "")
  panel_label("H")
}

run_figure3 <- function() {
  qc_df <- read_singlecell_qc()
  cell_scores_main <- read_cell_scores("results/scpagwas/main/cell_scores.tsv")
  main_df <- read_pvalue_table("results/scpagwas/main/merged_celltype_pvalue.csv")
  repeat_df <- read_pvalue_table("results/scpagwas/repeat_run/merged_celltype_pvalue.csv")
  ann_path <- Sys.getenv("EZHU_CLUSTER_ANNOTATIONS", unset = "results/annotation/cluster_annotations_filled.tsv")

  ezhu_write_stage_metadata(
    "07_figures_figure3",
    params = list(
      formats = TARGET_FORMATS,
      inputs = list(
        main_pvalues = "results/scpagwas/main/merged_celltype_pvalue.csv",
        repeat_pvalues = "results/scpagwas/repeat_run/merged_celltype_pvalue.csv",
        cell_scores_main = "results/scpagwas/main/cell_scores.tsv",
        annotation = ann_path
      )
    ),
    seed = seed
  )

  ezhu_save_plot(
    function() draw_figure3_scpagwas(main_df, repeat_df, ann_path, cell_scores_main, qc_df),
    "figure3_scpagwas_cluster_pvalues_main_vs_repeat",
    height = 11
  )
}

tryCatch(
  run_figure3(),
  error = function(e) {
    ezhu_write_stage_metadata(
      "07_figures_figure3",
      params = list(error = conditionMessage(e)),
      seed = seed
    )
    ezhu_save_plot(
      function() ezhu_placeholder_plot("Figure 3 (scPagwas) (failed)", c(conditionMessage(e))),
      "figure3_scpagwas_cluster_pvalues_main_vs_repeat",
      height = 11
    )
  }
)
