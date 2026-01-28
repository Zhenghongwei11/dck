source("scripts/figures/_setup.R")
source("scripts/figures/_common.R")

read_controls_table <- function(path) {
  if (!file.exists(path)) stop("Missing table: ", path, call. = FALSE)
  df <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- tolower(names(df))
  required <- c("cluster", "pvalue", "fdr", "gwas_source")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) stop("Missing columns in ", path, ": ", paste(missing, collapse = ", "), call. = FALSE)
  df$cluster <- as.character(df$cluster)
  df$pvalue <- suppressWarnings(as.numeric(df$pvalue))
  df$fdr <- suppressWarnings(as.numeric(df$fdr))
  df$gwas_source <- as.character(df$gwas_source)
  df <- df[nzchar(df$cluster) & is.finite(df$pvalue) & nzchar(df$gwas_source), , drop = FALSE]
  df
}

read_cluster_labels <- function(path = "results/annotation/cluster_annotations_filled.tsv") {
  if (!file.exists(path)) return(list(label = character(), type = character(), state = character()))
  ann <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  names(ann) <- tolower(names(ann))
  if (!all(c("cluster", "proposed_cell_type") %in% names(ann))) return(list(label = character(), type = character(), state = character()))
  ann$cluster <- as.character(ann$cluster)
  ann$proposed_cell_type <- as.character(ann$proposed_cell_type)
  ann$proposed_state <- as.character(ann$proposed_state %||% "")

  full <- paste0(
    short_context_label(ann$proposed_cell_type, 40L),
    ifelse(nzchar(ann$proposed_state), paste0(" / ", short_context_label(ann$proposed_state, 28L)), "")
  )
  list(
    label = stats::setNames(full, ann$cluster),
    type = stats::setNames(ann$proposed_cell_type, ann$cluster),
    state = stats::setNames(ann$proposed_state, ann$cluster)
  )
}

pretty_gwas_label <- function(id) {
  id <- as.character(id)
  if (identical(id, "Main")) return("DKD\n(GCST005881)")
  if (identical(id, "DKDGENIE2022_ALLVCTRL")) return("DKD meta\n(GENIE)")
  if (identical(id, "GCST008064")) return("CKD\n(CKDGen)")
  if (identical(id, "GCST008058")) return("eGFR\n(CKDGen)")
  id
}

draw_figure6_multi_gwas <- function() {
  top_path <- "results/scpagwas/gwas_controls/controls_top_contexts.tsv"
  full_path <- "results/scpagwas/gwas_controls/controls_summary.tsv"
  if (!file.exists(top_path) || !file.exists(full_path)) {
    ezhu_placeholder_plot(
      "Cross-GWAS robustness (missing tables)",
      c(
        paste0("Missing: ", top_path, " and/or ", full_path),
        "Run: Rscript --vanilla scripts/11_scpagwas_gwas_controls_summary.R"
      )
    )
    return(invisible(NULL))
  }

  top_df <- read_controls_table(top_path)
  full_df <- read_controls_table(full_path)

  # Define GWAS column order (human-readable).
  gwas_levels <- c("Main", "DKDGENIE2022_ALLVCTRL", "GCST008064", "GCST008058")
  gwas_levels <- gwas_levels[gwas_levels %in% unique(full_df$gwas_source)]
  if (!"Main" %in% gwas_levels) stop("controls tables missing 'Main' rows.", call. = FALSE)

  # Contexts: use primary significant clusters (already pre-filtered into controls_top_contexts.tsv).
  main_top <- top_df[top_df$gwas_source == "Main", , drop = FALSE]
  main_top <- main_top[order(main_top$pvalue, decreasing = FALSE), , drop = FALSE]
  clusters <- unique(as.character(main_top$cluster))

  # Attach cluster labels (marker-supported).
  label_map <- read_cluster_labels()
  cluster_labels <- unname(label_map$label[clusters])
  missing_idx <- which(is.na(cluster_labels) | !nzchar(cluster_labels))
  if (length(missing_idx) > 0) cluster_labels[missing_idx] <- paste0("Cluster ", clusters[missing_idx])

  short_override <- c(
    `12` = "C12 Podocyte (glomerular)",
    `3`  = "C3 Proximal tubule (S3-like)",
    `0`  = "C0 CD PC cell/Principal-CD",
    `7`  = "C7 TAL (TAL/Loop)",
    `13` = "C13 Injured epithelial (mixed)",
    `6`  = "C6 Distal nephron/CNT-like",
    `9`  = "C9 VCAM1+ adhesion/stress",
    `11` = "C11 CCN2/ECM-associated",
    `2`  = "C2 Distal convoluted tubule",
    `4`  = "C4 CD IC (prov.)/Acid-base",
    `5`  = "C5 Cluster 5",
    `15` = "C15 Pericyte/vascular mural"
  )

  fallback <- short_context_label(cluster_labels, 28L)
  y_labels <- ifelse(clusters %in% names(short_override), short_override[clusters], paste0("C", clusters, "  ", fallback))
  y_labels <- wrap_label(y_labels, width = 30)

  # Build matrices for plotting.
  get_val <- function(df, cl, gwas, col) {
    sub <- df[df$cluster == cl & df$gwas_source == gwas, , drop = FALSE]
    if (nrow(sub) == 0) return(NA_real_)
    suppressWarnings(as.numeric(sub[[col]][1]))
  }

  logp_mat <- matrix(NA_real_, nrow = length(clusters), ncol = length(gwas_levels), dimnames = list(clusters, gwas_levels))
  fdr_mat <- matrix(NA_real_, nrow = length(clusters), ncol = length(gwas_levels), dimnames = list(clusters, gwas_levels))
  for (i in seq_along(clusters)) {
    cl <- clusters[i]
    for (j in seq_along(gwas_levels)) {
      gw <- gwas_levels[j]
      p <- get_val(full_df, cl, gw, "pvalue")
      f <- get_val(full_df, cl, gw, "fdr")
      logp_mat[i, j] <- ifelse(is.finite(p), -log10(pmax(p, 1e-300)), NA_real_)
      fdr_mat[i, j] <- f
    }
  }

  # Cap the dynamic range for visual comparability.
  cap <- 8
  logp_cap <- pmin(logp_mat, cap)

  # Persist an audit-friendly anchor matrix for Figure 6.
  out_anchor <- data.frame(
    cluster = clusters,
    context_label = cluster_labels,
    stringsAsFactors = FALSE
  )
  for (gw in gwas_levels) {
    out_anchor[[paste0("p_", gw)]] <- suppressWarnings(as.numeric(logp_mat[, gw]))
  }
  utils::write.table(out_anchor, file = "results/figures/figure6_multi_gwas_logp_matrix.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE))

  # Panel A: dot plot (color = -log10(p), outline = FDR).
  par(mar = c(7, 12, 3.5, 5))
  x_labels <- vapply(gwas_levels, pretty_gwas_label, character(1))
  n_x <- length(gwas_levels)
  n_y <- length(clusters)

  plot(NA, xlim = c(0.5, n_x + 0.5), ylim = c(0.5, n_y + 0.5), xaxt = "n", yaxt = "n",
       xlab = "", ylab = "", main = "Cross-GWAS robustness (top DKD-localized contexts)", frame.plot = FALSE)
  axis(1, at = seq_len(n_x), labels = x_labels, las = 2, cex.axis = 0.8, lwd = 0, lwd.ticks = 0)
  axis(2, at = seq_len(n_y), labels = rev(y_labels), las = 1, cex.axis = 0.75, lwd = 0, lwd.ticks = 0)

  pal <- ezhu_seq_pal(64)
  zlim <- c(0, cap)
  for (i in seq_len(n_y)) {
    y <- n_y - i + 1
    for (j in seq_len(n_x)) {
      v <- logp_cap[i, j]
      if (!is.finite(v)) v <- 0
      idx <- as.integer(1 + 63 * (v - zlim[1]) / diff(zlim))
      idx <- pmax(1, pmin(64, idx))
      f <- fdr_mat[i, j]
      border_col <- if (is.finite(f) && f <= 0.05) "black" else "grey80"
      points(j, y, pch = 21, bg = pal[idx], col = border_col, cex = 1.8, lwd = 0.8)
    }
  }

  # Colorbar
  usr <- par("usr")
  xleft <- usr[2] + 0.15
  xright <- usr[2] + 0.35
  y_seq <- seq(usr[3], usr[4], length.out = 65)
  for (k in 1:64) {
    rect(xleft, y_seq[k], xright, y_seq[k + 1], col = pal[k], border = NA, xpd = TRUE)
  }
  text(xright, c(usr[3], usr[4]), labels = c(0, cap), pos = 4, cex = 0.7, xpd = TRUE)
  mtext(expression(-log[10](italic(p))~"(capped)"), side = 4, line = 3.2, cex = 0.75)
  mtext("Outline: FDR <= 0.05", side = 1, line = 5.2, cex = 0.7, col = "grey40")
  panel_label("A")

  # Panel B: rank concordance vs DKD (Main), lollipop style.
  par(mar = c(5, 6, 2.5, 2))
  all_clusters <- sort(unique(full_df$cluster))
  main_all <- full_df[full_df$gwas_source == "Main", , drop = FALSE]
  main_all <- main_all[match(all_clusters, main_all$cluster), , drop = FALSE]
  main_logp <- -log10(pmax(main_all$pvalue, 1e-300))

  other <- setdiff(gwas_levels, "Main")
  rho <- numeric(length(other))
  names(rho) <- other
  for (k in seq_along(other)) {
    gw <- other[k]
    sub <- full_df[full_df$gwas_source == gw, , drop = FALSE]
    sub <- sub[match(all_clusters, sub$cluster), , drop = FALSE]
    logp <- -log10(pmax(sub$pvalue, 1e-300))
    rho[k] <- suppressWarnings(stats::cor(main_logp, logp, method = "spearman", use = "complete.obs"))
  }

  pal_ok <- okabe_ito_palette()
  cols <- c(pal_ok["orange"], pal_ok["grey"], pal_ok["sky"])[seq_along(other)]
  y <- seq_along(other)
  plot(rho, y,
       xlim = c(-1, 1),
       ylim = c(0.5, length(other) + 0.5),
       yaxt = "n",
       xlab = "Spearman rho (cluster -log10(p))",
       ylab = "",
       main = "Rank concordance vs DKD (Main)",
       frame.plot = FALSE)
  axis(2, at = y, labels = vapply(other, pretty_gwas_label, character(1)), las = 1, cex.axis = 0.85, lwd = 0, lwd.ticks = 0)
  axis(1, lwd = 0, lwd.ticks = 0)
  abline(v = 0, lty = 2, col = "grey60")
  segments(0, y, rho, y, col = "grey70", lwd = 2)
  points(rho, y, pch = 21, bg = cols, col = "grey30", cex = 1.6, lwd = 0.6)
  text(x = rho, y = y, labels = sprintf("%.2f", rho), pos = ifelse(rho >= 0, 4, 2), cex = 0.8, font = 2)
  panel_label("B")

  # Panel C/D: scatter concordance against DKD (Main)
  scatter_concordance <- function(other_id, title, letter) {
    sub <- full_df[full_df$gwas_source == other_id, , drop = FALSE]
    sub <- sub[match(all_clusters, sub$cluster), , drop = FALSE]
    other_logp <- -log10(pmax(sub$pvalue, 1e-300))

    pal_ok <- okabe_ito_palette()
    sig_main <- is.finite(main_all$fdr) & main_all$fdr <= 0.05
    point_bg <- ifelse(sig_main, pal_ok["orange"], pal_ok["light_grey"])
    point_col <- ifelse(sig_main, pal_ok["orange"], "grey70")

    lim <- range(c(main_logp, other_logp), finite = TRUE)
    if (!all(is.finite(lim)) || diff(lim) == 0) lim <- c(0, 1)
    lim <- c(0, max(lim) * 1.05)

    plot(
      main_logp, other_logp,
      xlim = lim, ylim = lim,
      pch = 21, bg = point_bg, col = point_col, cex = 1.0,
      xlab = "DKD (Main) -log10(p)",
      ylab = paste0(gsub("\n", " ", pretty_gwas_label(other_id)), " -log10(p)"),
      main = title,
      frame.plot = FALSE
    )
    abline(0, 1, lty = 2, col = "grey70")
    rho_val <- suppressWarnings(stats::cor(main_logp, other_logp, method = "spearman", use = "complete.obs"))
    text(x = lim[2] * 0.05, y = lim[2] * 0.92, labels = sprintf("rho=%.2f", rho_val), adj = 0, cex = 0.85)
    panel_label(letter)
  }

  par(mar = c(5, 6, 2.5, 2))
  scatter_concordance("DKDGENIE2022_ALLVCTRL", "DKD vs DKD meta", "C")

  par(mar = c(5, 6, 2.5, 2))
  scatter_concordance("GCST008058", "DKD vs eGFR", "D")
}

ezhu_write_stage_metadata(
  "07_figures_figure6",
  params = list(
    formats = TARGET_FORMATS,
    inputs = list(
      controls_top_contexts = "results/scpagwas/gwas_controls/controls_top_contexts.tsv",
      controls_summary = "results/scpagwas/gwas_controls/controls_summary.tsv",
      cluster_annotations = "results/annotation/cluster_annotations_filled.tsv"
    ),
    outputs = list(anchor_matrix = "results/figures/figure6_multi_gwas_logp_matrix.tsv")
  ),
  seed = seed
)

ezhu_save_plot(
  function() draw_figure6_multi_gwas(),
  "figure6_scpagwas_multi_gwas_robustness",
  width = 11,
  height = 7.2
)
