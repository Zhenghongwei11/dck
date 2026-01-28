source("scripts/figures/_setup.R")
source("scripts/figures/_common.R")

draw_figure5_mr <- function(mr_path = "results/causal/mr_validation.tsv", mr_meta_path = "results/metadata/06_causal_inference_mr__metadata.json") {
  if (!file.exists(mr_path)) {
    ezhu_placeholder_plot(
      "Causal follow-up (MR) (missing input)",
      c(paste0("Missing: ", mr_path), "Run: make causal")
    )
    return(invisible(NULL))
  }
  df <- utils::read.delim(mr_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- tolower(names(df))

  meta_note <- NULL
  if (file.exists(mr_meta_path) && requireNamespace("jsonlite", quietly = TRUE)) {
    meta <- try(jsonlite::fromJSON(mr_meta_path), silent = TRUE)
    if (!inherits(meta, "try-error") && !is.null(meta$params)) {
      outcomes <- meta$params$outcomes %||% character()
      p1 <- meta$params$p1 %||% NA_real_
      clump <- meta$params$clump %||% NA
      note_parts <- character()
      if (length(outcomes) > 0) note_parts <- c(note_parts, paste0("Outcomes: ", paste(outcomes, collapse = ", ")))
      if (is.finite(p1)) note_parts <- c(note_parts, paste0("p1=", format(p1, scientific = TRUE)))
      if (!is.na(clump)) note_parts <- c(note_parts, paste0("clump=", tolower(as.character(clump))))
      if (length(note_parts) > 0) meta_note <- paste(note_parts, collapse = "  |  ")
    }
  }

  df$ivw_b <- safe_numeric(df$ivw_b)
  df$ivw_se <- safe_numeric(df$ivw_se)
  df$ivw_pval <- safe_numeric(df$ivw_pval)
  df$ivw_fdr <- safe_numeric(df$ivw_fdr)
  df$wm_b <- safe_numeric(df$wm_b)
  df$het_ivw_q_pval <- safe_numeric(df$het_ivw_q_pval)
  df$egger_intercept_pval <- safe_numeric(df$egger_intercept_pval)

  df <- df[is.finite(df$ivw_pval) & is.finite(df$ivw_b) & is.finite(df$ivw_se), , drop = FALSE]
  if (nrow(df) == 0) {
    ezhu_placeholder_plot("Causal follow-up (MR)", c("No valid IVW rows found in the MR table."))
    return(invisible(NULL))
  }

  df$gene <- ifelse(nzchar(df$exposure_gene_symbol_ensembl), df$exposure_gene_symbol_ensembl, df$exposure_ensg)
  df$outcome_short <- sub(" \\|\\| id:.*$", "", as.character(df$outcome_trait))
  df$ci_low <- df$ivw_b - 1.96 * df$ivw_se
  df$ci_high <- df$ivw_b + 1.96 * df$ivw_se

  df <- df[order(df$ivw_fdr, df$ivw_pval, decreasing = FALSE, na.last = TRUE), , drop = FALSE]
  hits <- df[is.finite(df$ivw_fdr) & df$ivw_fdr <= 0.05, , drop = FALSE]
  if (nrow(hits) == 0) hits <- utils::head(df, 12)
  hits <- utils::head(hits, 12)

  fig_hits <- hits[, c("gene", "exposure_ensg", "outcome_id", "outcome_short", "ivw_nsnp", "ivw_b", "ivw_se", "ci_low", "ci_high", "ivw_pval", "ivw_fdr")]
  utils::write.table(fig_hits, file = "results/figures/mr_top_hits.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

  layout(matrix(1:8, nrow = 4, byrow = TRUE))
  pal <- okabe_ito_palette()

  # A: Top MR signals (Bubble Lollipop)
  par(mar = c(6, 18, 3.6, 4)) # Increased left margin
  use_fdr <- is.finite(fig_hits$ivw_fdr) & fig_hits$ivw_fdr >= 0 & fig_hits$ivw_fdr <= 1
  if (any(use_fdr)) {
    score <- -log10(pmax(fig_hits$ivw_fdr, 1e-300))
    xlab <- expression(-log[10](IVW~FDR))
    cutoff <- -log10(0.05)
    note <- "Size ~ nSNPs; X-axis ~ -log10(FDR)"
  } else {
    score <- -log10(pmax(fig_hits$ivw_pval, 1e-300))
    xlab <- expression(-log[10](IVW~p))
    cutoff <- -log10(0.05)
    note <- "Size ~ nSNPs; X-axis ~ -log10(p)"
  }
  
  labels <- paste0(fig_hits$gene, "\n(", fig_hits$outcome_short, ")")
  y_pos <- rev(seq_along(score))
  
  # Scale point size by nSNPs
  nsnp <- safe_numeric(fig_hits$ivw_nsnp)
  nsnp[!is.finite(nsnp)] <- 1
  cex_vals <- 1.2 + 1.5 * (nsnp / max(nsnp, na.rm=TRUE))
  
  plot(NA, xlim=c(0, max(score)*1.15), ylim=c(0.5, length(y_pos)+0.5), yaxt="n", xlab=xlab,
       ylab="", main="Top MR signals (IVW)", frame.plot=FALSE)
  axis(2, at=y_pos, labels=sanitize_ascii(labels), las=1, cex.axis=0.7)
  ezhu_light_grid(nx=NULL, ny=NA)
  
  segments(0, y_pos, rev(score), y_pos, col="grey60", lwd=1.5)
  points(rev(score), y_pos, pch=21, bg=pal["blue"], col="white", cex=rev(cex_vals), lwd=0.5)
  
  abline(v = cutoff, lty = 2, col = "grey70")
  mtext(note, side = 3, line = 0.45, cex = 0.65, col = "grey40")
  panel_label("A")

  # B
  par(mar = c(6, 14, 3, 4)) # Increased left margin
  y <- seq_len(nrow(fig_hits))
  plot(fig_hits$ivw_b, y, xlim = range(c(fig_hits$ci_low, fig_hits$ci_high), finite = TRUE), yaxt = "n",
       ylab = "", xlab = "IVW effect (b) with 95% CI", main = "IVW forest plot (screening)")
  axis(2, at = y, labels = sanitize_ascii(labels), las = 2, cex.axis = 0.7)
  abline(v = 0, lty = 2, col = "grey60")
  segments(fig_hits$ci_low, y, fig_hits$ci_high, y, col = "grey30", lwd = 2)
  points(fig_hits$ivw_b, y, pch = 16, col = pal["blue"])
  if (any(is.finite(fig_hits$ivw_fdr))) {
    x_text <- max(fig_hits$ci_high, fig_hits$ivw_b, finite = TRUE)
    x_text <- x_text + 0.03 * diff(range(c(fig_hits$ci_low, fig_hits$ci_high), finite = TRUE))
    txt <- paste0("FDR=", format_fdr_compact(fig_hits$ivw_fdr))
    text(x_text, y, labels = txt, adj = 0, cex = 0.65, col = "grey30", xpd = TRUE)
  }
  mtext("Effect scale follows the outcome GWAS.", side = 3, line = -1.2, cex = 0.7, col = "grey40")
  if (!is.null(meta_note)) {
    # Wrap text if extremely long
    meta_note <- paste(strwrap(meta_note, width = 90), collapse = "\n")
    mtext(meta_note, side = 1, line = 4.5, cex = 0.65, col = "grey40")
  }
  panel_label("B")

  # C
  par(mar = c(4, 14, 3, 4))
  ok <- is.finite(hits$wm_b)
  if (sum(ok) < 2) {
    plot.new()
    title(main = "Direction concordance (IVW vs WM)")
    text(0.05, 0.8, "Not enough WM estimates in the current top hits.", adj = c(0, 1), cex = 0.9)
  } else {
    plot(hits$ivw_b[ok], hits$wm_b[ok], pch = 16, col = "grey30", xlab = "IVW b", ylab = "WM b",
         main = "Direction concordance (IVW vs WM)")
    abline(0, 1, lty = 2, col = "grey60")
  }
  panel_label("C")

  # D
  par(mar = c(4, 14, 3, 4))
  hv <- df$het_ivw_q_pval
  hv <- hv[is.finite(hv) & hv >= 0 & hv <= 1]
  if (length(hv) < 2) {
    plot.new()
    title(main = "Heterogeneity (Q) p-values")
    text(0.05, 0.8, "Not available for most rows in the current table.", adj = c(0, 1), cex = 0.9)
  } else {
    hist(hv, breaks = 30, col = "#D9D9D9", border = "white", main = "Heterogeneity (Q) p-values", xlab = "het_ivw_q_pval")
    abline(v = 0.05, lty = 2, col = pal["orange"])
  }
  panel_label("D")

  # E
  par(mar = c(4, 14, 3, 4))
  ev <- df$egger_intercept_pval
  ev <- ev[is.finite(ev) & ev >= 0 & ev <= 1]
  if (length(ev) < 2) {
    plot.new()
    title(main = "Egger intercept p-values")
    text(0.05, 0.8, "Not available for most rows in the current table.", adj = c(0, 1), cex = 0.9)
  } else {
    hist(ev, breaks = 30, col = "#D9D9D9", border = "white", main = "Egger intercept p-values", xlab = "egger_intercept_pval")
    abline(v = 0.05, lty = 2, col = pal["orange"])
  }
  panel_label("E")

  # F
  par(mar = c(8, 14, 3, 9)) # More bottom margin for wrapped x labels, right for legend
  u_genes <- unique(fig_hits$gene)
  u_out <- unique(fig_hits$outcome_short)
  mat <- matrix(0, nrow = length(u_genes), ncol = length(u_out), dimnames = list(u_genes, u_out))
  for (i in seq_len(nrow(fig_hits))) {
    g <- fig_hits$gene[i]
    o <- fig_hits$outcome_short[i]
    mat[g, o] <- sign(fig_hits$ivw_b[i]) * (-log10(pmax(fig_hits$ivw_pval[i], 1e-300)))
  }
  pal_div <- ezhu_div_pal(64)
  zlim <- c(-max(abs(mat)), max(abs(mat)))
  if (zlim[1] == zlim[2]) zlim <- c(-1, 1)
  image(seq_along(u_out), seq_along(u_genes), t(mat), axes = FALSE, col = pal_div, zlim = zlim, xlab = "", ylab = "", main = "Signed -log10(p) (IVW)")
  axis(1, at = seq_along(u_out), labels = wrap_label(u_out, 14), las = 2, cex.axis = 0.55)
  axis(2, at = seq_along(u_genes), labels = u_genes, las = 2, cex.axis = 0.8)
  box()
  
  # Add direction legend in the right margin
  usr <- par("usr")
  x_legend <- usr[2] + 0.25 * (usr[2] - usr[1])
  y_legend <- usr[3] + 0.2 * (usr[4] - usr[3])
  legend(x = x_legend, y = y_legend, legend = c("Positive", "Negative"), bty = "n", cex = 0.6,
         text.col = c(pal["orange"], pal["blue"]), title = "Effect", xpd = TRUE, xjust = 0)

  # Add colorbar
  xleft <- usr[2] + 0.05 * (usr[2] - usr[1])
  xright <- usr[2] + 0.1 * (usr[2] - usr[1])
  y_seq <- seq(usr[3], usr[4], length.out = 65)
  for (idx_f in 1:64) {
    rect(xleft, y_seq[idx_f], xright, y_seq[idx_f+1], col = pal_div[idx_f], border = NA, xpd = TRUE)
  }
  text(xright, c(usr[3], (usr[3]+usr[4])/2, usr[4]), labels = round(c(zlim[1], 0, zlim[2]), 1), pos = 4, cex = 0.6, xpd = TRUE)
  panel_label("F")

  # G: number of instruments (nsnp) for top hits
  par(mar = c(6, 18, 3, 4))
  nsnp <- safe_numeric(fig_hits$ivw_nsnp)
  nsnp[!is.finite(nsnp)] <- 0
  barplot(rev(nsnp), names.arg = rev(sanitize_ascii(labels)), horiz = TRUE, las = 1, col = "#BDBDBD", border = NA,
          xlab = "IVW instruments (n SNPs)", main = "Instrument counts", cex.names = 0.7)
  panel_label("G")

  # H: pleiotropy vs heterogeneity snapshot
  par(mar = c(4, 14, 3, 4))
  xh <- -log10(pmax(df$het_ivw_q_pval, 1e-300))
  yh <- -log10(pmax(df$egger_intercept_pval, 1e-300))
  okh <- is.finite(xh) & is.finite(yh)
  if (sum(okh) < 2) {
    plot.new()
    title(main = "Pleiotropy vs heterogeneity")
  } else {
    colh <- ifelse(df$ivw_b >= 0, pal["orange"], pal["blue"])
    sz <- safe_numeric(df$ivw_nsnp)
    sz[!is.finite(sz)] <- 1
    cexh <- 0.7 + 0.15 * pmin(sz, 10)
    plot(xh[okh], yh[okh], pch = 16, col = grDevices::adjustcolor(colh[okh], 0.6), cex = cexh[okh],
         xlab = expression(-log[10](Q~pvalue)), ylab = expression(-log[10](Egger~intercept~pvalue)),
         main = "Pleiotropy vs heterogeneity")
    abline(v = -log10(0.05), lty = 2, col = "grey70")
    abline(h = -log10(0.05), lty = 2, col = "grey70")
  }
  panel_label("H")
}

ezhu_write_stage_metadata(
  "07_figures_figure5",
  params = list(
    formats = TARGET_FORMATS,
    inputs = list(
      mr_validation = "results/causal/mr_validation.tsv",
      mr_metadata = "results/metadata/06_causal_inference_mr__metadata.json"
    )
  ),
  seed = seed
)

ezhu_save_plot(
  function() draw_figure5_mr(),
  "figure5_mr_top_hits",
  height = 11
)
