source("scripts/figures/_setup.R")
source("scripts/figures/_common.R")

find_latest_ccc_dir <- function(root = "results/ccc") {
  if (!dir.exists(root)) return("")
  dirs <- list.dirs(root, full.names = TRUE, recursive = FALSE)
  dirs <- dirs[dir.exists(dirs)]
  if (length(dirs) == 0) return("")
  has_tables <- vapply(dirs, function(d) {
    file.exists(file.path(d, "cellchat_interactions.tsv")) &&
      file.exists(file.path(d, "cellchat_group_counts.tsv"))
  }, logical(1))
  dirs <- dirs[has_tables]
  if (length(dirs) == 0) return("")

  # Prefer the newest *interaction table* timestamp, not the directory mtime,
  # because some cloud sync/copy operations can preserve directory mtimes.
  stamp <- vapply(dirs, function(d) {
    f <- file.path(d, "cellchat_interactions.tsv")
    info <- suppressWarnings(file.info(f))
    as.numeric(info$mtime %||% 0)
  }, numeric(1))
  dirs[order(stamp, decreasing = TRUE)][1]
}

draw_figure7_ccc <- function() {
  ccc_dir <- Sys.getenv("EZHU_CCC_DIR", unset = "")
  if (!nzchar(ccc_dir)) ccc_dir <- find_latest_ccc_dir()
  if (!nzchar(ccc_dir)) {
    ezhu_placeholder_plot(
      "Cell-cell communication (CellChat) (missing input)",
      c("No results/ccc/<run_id>/ directory found.", "Run: make ccc")
    )
    return(invisible(NULL))
  }

  int_path <- file.path(ccc_dir, "cellchat_interactions.tsv")
  path_path <- file.path(ccc_dir, "cellchat_pathway_summary.tsv")
  grp_path <- file.path(ccc_dir, "cellchat_group_counts.tsv")
  if (!file.exists(int_path) || !file.exists(grp_path)) {
    ezhu_placeholder_plot(
      "Cell-cell communication (CellChat) (missing tables)",
      c(paste0("Dir: ", ccc_dir), paste0("Missing: ", paste(c(int_path, grp_path)[!file.exists(c(int_path, grp_path))], collapse = ", ")))
    )
    return(invisible(NULL))
  }

  comm <- utils::read.delim(int_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  comm$prob <- safe_numeric(comm$prob)
  comm$pval <- safe_numeric(comm$pval)
  comm <- comm[is.finite(comm$prob) & is.finite(comm$pval), , drop = FALSE]
  comm$fdr <- stats::p.adjust(comm$pval, method = "BH")
  comm <- comm[order(comm$fdr, -comm$prob), , drop = FALSE]

  # Persist derived table (audit-friendly, cheap to compute locally).
  out_with_fdr <- file.path(ccc_dir, "cellchat_interactions_with_fdr.tsv")
  try(utils::write.table(comm, file = out_with_fdr, sep = "\t", quote = FALSE, row.names = FALSE), silent = TRUE)

  # Deduplicate LR names for plotting (avoid repeated labels dominating the plot).
  lr_name <- as.character(comm$interaction_name_2 %||% comm$interaction_name %||% "")
  lr_name[!nzchar(lr_name)] <- paste0(as.character(comm$ligand), "_", as.character(comm$receptor))[!nzchar(lr_name)]
  comm$lr_name <- lr_name
  comm_lr <- aggregate(
    cbind(prob = comm$prob, pval = comm$pval, fdr = comm$fdr) ~ lr_name,
    data = comm,
    FUN = function(x) c(max = max(x, na.rm = TRUE), min = min(x, na.rm = TRUE))
  )
  if (is.matrix(comm_lr$prob)) {
    comm_lr$prob <- comm_lr$prob[, "max"]
  } else {
    comm_lr$prob <- vapply(comm_lr$prob, `[[`, numeric(1), "max")
  }
  if (is.matrix(comm_lr$pval)) {
    comm_lr$pval <- comm_lr$pval[, "min"]
  } else {
    comm_lr$pval <- vapply(comm_lr$pval, `[[`, numeric(1), "min")
  }
  if (is.matrix(comm_lr$fdr)) {
    comm_lr$fdr <- comm_lr$fdr[, "min"]
  } else {
    comm_lr$fdr <- vapply(comm_lr$fdr, `[[`, numeric(1), "min")
  }
  comm_lr <- comm_lr[order(comm_lr$fdr, -comm_lr$prob), , drop = FALSE]
  comm_top_lr <- utils::head(comm_lr, 12)

  comm_top <- utils::head(comm, 12)

  grps <- utils::read.delim(grp_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  grps$n_cells <- safe_numeric(grps$n_cells)
  grps <- grps[order(grps$n_cells, decreasing = TRUE), , drop = FALSE]
  grp_cluster <- sub("^C", "", as.character(grps$group))
  cl_pal <- make_cluster_palette(grp_cluster)
  grp_cols <- unname(cl_pal[grp_cluster])
  grp_cols[is.na(grp_cols)] <- "grey50"

  pathway <- NULL
  if (file.exists(path_path)) {
    pathway <- utils::read.delim(path_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    pathway$prob_sum <- safe_numeric(pathway$prob_sum)
    pathway <- pathway[is.finite(pathway$prob_sum), , drop = FALSE]
    pathway <- pathway[order(pathway$prob_sum, decreasing = TRUE), , drop = FALSE]
    if ("pathway_name" %in% names(pathway)) {
      pathway_sum <- aggregate(prob_sum ~ pathway_name, data = pathway, FUN = sum)
      pathway_sum <- pathway_sum[order(pathway_sum$prob_sum, decreasing = TRUE), , drop = FALSE]
      pathway <- pathway_sum
    }
  }

  layout(matrix(1:8, nrow = 4, byrow = TRUE))
  pal_ok <- okabe_ito_palette()

  # A: group sizes
  par(mar = c(6, 18, 3, 4))
  grp_labels <- short_context_label(grps$group, max_chars = 35L)
  barplot(rev(grps$n_cells), names.arg = rev(grp_labels), horiz = TRUE, las = 1, col = rev(grp_cols), border = NA,
          main = "CellChat input groups", xlab = "Cells (n)")
  panel_label("A")

  # B: aggregate communication heatmap (source x target; prob sum)
  par(mar = c(8, 12, 3, 4))
  mat <- xtabs(prob ~ source + target, data = comm)
  mat <- as.matrix(mat)
  mat[!is.finite(mat)] <- 0
  if (all(mat == 0)) mat[1, 1] <- 0
  pal <- ezhu_seq_pal(64)
  zlim <- range(mat, finite = TRUE)
  if (!all(is.finite(zlim)) || diff(zlim) == 0) zlim <- c(0, 1)
  image(seq_len(ncol(mat)), seq_len(nrow(mat)), t(mat), axes = FALSE, col = pal, zlim = zlim,
        xlab = "", ylab = "", main = "Aggregate communication (prob sum)")
  axis(1, at = seq_len(ncol(mat)), labels = short_context_label(colnames(mat), 30L), las = 2, cex.axis = 0.7)
  axis(2, at = seq_len(nrow(mat)), labels = short_context_label(rownames(mat), 30L), las = 2, cex.axis = 0.7)
  box()
  # Colorbar
  usr <- par("usr")
  xleft <- usr[2] + 0.03 * (usr[2] - usr[1])
  xright <- usr[2] + 0.07 * (usr[2] - usr[1])
  y_seq <- seq(usr[3], usr[4], length.out = 65)
  for (idx_b in 1:64) {
    rect(xleft, y_seq[idx_b], xright, y_seq[idx_b + 1], col = pal[idx_b], border = NA, xpd = TRUE)
  }
  text(xright, c(usr[3], usr[4]), labels = signif(zlim, 2), pos = 4, cex = 0.6, xpd = TRUE)
  mtext("prob_sum", side = 4, line = 2.4, cex = 0.75)
  panel_label("B")

  # C: top ligand-receptor interactions (Modern Dot Chart)
  par(mar = c(6, 12, 3, 4))
  lr <- sanitize_ascii(comm_top_lr$lr_name)
  lr <- ifelse(nchar(lr) > 45, paste0(substr(lr, 1, 42), "..."), lr)
  
  # Y-axis positions
  y_pos <- seq_along(lr)
  # X-axis is -log10(FDR)
  x_vals <- -log10(pmax(comm_top_lr$fdr, 1e-300))
  x_vals_cap <- pmin(x_vals, 50)
  
  # Point size based on interaction probability (strength)
  # Scale it from 1 to 3
  cex_vals <- 1 + 2 * (comm_top_lr$prob / max(comm_top_lr$prob, na.rm=TRUE))
  
  plot(NA, xlim = c(0, 60), ylim = c(0.5, length(lr) + 0.5),
       yaxt = "n", xlab = expression(-log[10](BH~FDR)), ylab = "", 
       main = "Top L-R interactions (Sig. & Strength)", frame.plot = FALSE)
  
  axis(2, at = y_pos, labels = rev(sanitize_ascii(lr)), las = 1, cex.axis = 0.8)
  ezhu_light_grid(nx = NULL, ny = NA)
  
  # Add segments for visual guidance (dot-chart style)
  segments(0, y_pos, rev(x_vals_cap), y_pos, col = "grey90", lty = 1, lwd = 1)
  
  # Plot points: Color is significance, Size is probability
  points(rev(x_vals_cap), y_pos, pch = 21, bg = pal_ok["blue"], col = "black", cex = rev(cex_vals))
  
  # Add legends
  mtext("Size ~ Prob; X-axis ~ -log10(FDR)", side = 1, line = 3.6, cex = 0.7, col = "grey40")
  panel_label("C")

  # D: top pathways (summed probability) - Lollipop Chart
  par(mar = c(6, 12, 3, 4))
  if (!is.null(pathway) && nrow(pathway) > 0 && "pathway_name" %in% names(pathway)) {
    pw <- sanitize_ascii(pathway$pathway_name)
    pw <- ifelse(nchar(pw) > 40, paste0(substr(pw, 1, 37), "..."), pw)
    pw_sum <- data.frame(pathway = pw, prob_sum = safe_numeric(pathway$prob_sum), stringsAsFactors = FALSE)
    pw_sum <- pw_sum[is.finite(pw_sum$prob_sum) & nzchar(pw_sum$pathway), , drop = FALSE]
    pw_sum <- pw_sum[order(pw_sum$prob_sum, decreasing = TRUE), , drop = FALSE]
    pw_sum <- utils::head(pw_sum, 12)
    
    y_pos <- rev(seq_len(nrow(pw_sum)))
    xlim <- c(0, max(pw_sum$prob_sum) * 1.1)
    
    plot(NA, xlim = xlim, ylim = c(0.5, length(y_pos) + 0.5), yaxt="n", xlab="Summed Probability", 
         ylab="", main="Top pathways (strength)", frame.plot=FALSE)
    axis(2, at=y_pos, labels=pw_sum$pathway, las=1, cex.axis=0.8)
    ezhu_light_grid(nx=NULL, ny=NA)
    
    # Lollipop stems
    segments(0, y_pos, pw_sum$prob_sum, y_pos, col="grey70", lwd=1.5)
    # Lollipop heads
    points(pw_sum$prob_sum, y_pos, pch=21, bg=pal_ok["blue"], col="white", cex=1.5, lwd=0.5)
    
  } else {
    plot.new()
    title(main = "Top pathways (summed probability)")
    mtext("Pathway annotations not available in anchor table.", side = 3, line = -1.2, cex = 0.8, col = "grey40")
  }
  panel_label("D")

  # E: prob vs -log10(p)
  par(mar = c(4, 12, 4, 4))
  f_vals <- -log10(pmax(comm$fdr, 1e-300))
  f_vals_cap <- pmin(f_vals, 50)
  plot(comm$prob, f_vals_cap, pch = 16, col = rgb(0, 0, 0, 0.12),
       xlab = "prob", ylab = expression(-log[10](BH~FDR)), main = "Interaction strength vs significance")
  mtext("BH-FDR across all interaction tests; -log10(FDR) capped at 50", side = 3, line = 0.5, cex = 0.7, col = "grey40")
  panel_label("E")

  # F: top sender->receiver edges (prob sum)
  par(mar = c(6, 12, 3, 10))
  edges <- aggregate(prob ~ source + target, data = comm, FUN = sum)
  pmin_tbl <- aggregate(pval ~ source + target, data = comm, FUN = min)
  edges <- merge(edges, pmin_tbl, by = c("source", "target"), all.x = TRUE)
  edges$prob <- safe_numeric(edges$prob)
  edges$pval <- safe_numeric(edges$pval)
  edges$fdr <- stats::p.adjust(edges$pval, method = "BH")
  edges <- edges[is.finite(edges$prob), , drop = FALSE]
  edges <- edges[order(edges$prob, decreasing = TRUE), , drop = FALSE]
  edges <- utils::head(edges, 10)
  src <- unique(edges$source)
  tgt <- unique(edges$target)
  ys <- setNames(seq_len(length(src)), src)
  yt <- setNames(seq_len(length(tgt)), tgt)
  plot(NA, xlim = c(0, 1), ylim = c(0.5, max(length(src), length(tgt)) + 0.5),
       xaxt = "n", yaxt = "n", xlab = "", ylab = "", main = "Top sender->receiver edges")
  axis(2, at = seq_along(src), labels = short_context_label(src, 30L), las = 1, cex.axis = 0.7)
  axis(4, at = seq_along(tgt), labels = short_context_label(tgt, 30L), las = 1, cex.axis = 0.75)
  mtext("source", side = 2, line = 5.0, cex = 0.8)
  mtext("target", side = 4, line = 6.0, cex = 0.8)
  segments(0.2, seq_along(src), 0.2, seq_along(src), col = "grey40", lwd = 2)
  segments(0.8, seq_along(tgt), 0.8, seq_along(tgt), col = "grey40", lwd = 2)

  w <- edges$prob / max(edges$prob)
  w[!is.finite(w)] <- 0
  lwd <- 1 + 4 * w
  sig <- -log10(pmax(edges$fdr, 1e-300))
  sig[!is.finite(sig)] <- 0
  pal_edge <- grDevices::colorRampPalette(c("#E0E0E0", "#9ECAE1", "#3182BD"))(64)
  cols <- pal_edge[pmax(1, pmin(64, as.integer(1 + 63 * (sig / max(sig + 1e-6)))))]

  for (i in seq_len(nrow(edges))) {
    arrows(0.22, ys[[edges$source[i]]], 0.78, yt[[edges$target[i]]],
           length = 0.08, angle = 20, col = cols[i], lwd = lwd[i])
  }
  mtext("Top 10 edges only; axes show involved groups.", side = 1, line = 4.3, cex = 0.7, col = "grey40")
  mtext("Edge width ~ prob_sum; color ~ -log10(BH FDR)", side = 1, line = 3.2, cex = 0.7, col = "grey40")
  panel_label("F")

  # G: outgoing signaling by source (Lollipop)
  par(mar = c(6, 12, 3, 4))
  out_sum <- rowSums(mat, na.rm = TRUE)
  out_df <- data.frame(group = names(out_sum), prob_sum = as.numeric(out_sum), stringsAsFactors = FALSE)
  out_df <- out_df[order(out_df$prob_sum, decreasing = TRUE), , drop = FALSE]
  out_df <- utils::head(out_df, 6)
  
  y_g <- rev(seq_len(nrow(out_df)))
  x_g <- out_df$prob_sum
  plot(NA, xlim=c(0, max(x_g)*1.15), ylim=c(0.5, length(y_g)+0.5), yaxt="n", xlab="Summed Probability",
       ylab="", main="Outgoing signaling (Sender)", frame.plot=FALSE)
  axis(2, at=y_g, labels=short_context_label(out_df$group, 30L), las=1, cex.axis=0.8)
  ezhu_light_grid(nx=NULL, ny=NA)
  segments(0, y_g, x_g, y_g, col="grey60", lwd=1.5)
  points(x_g, y_g, pch=21, bg=pal_ok["orange"], col="white", cex=1.6)
  panel_label("G")

  # H: incoming signaling by target (Lollipop)
  par(mar = c(6, 12, 3, 4))
  in_sum <- colSums(mat, na.rm = TRUE)
  in_df <- data.frame(group = names(in_sum), prob_sum = as.numeric(in_sum), stringsAsFactors = FALSE)
  in_df <- in_df[order(in_df$prob_sum, decreasing = TRUE), , drop = FALSE]
  in_df <- utils::head(in_df, 6)
  
  y_h <- rev(seq_len(nrow(in_df)))
  x_h <- in_df$prob_sum
  plot(NA, xlim=c(0, max(x_h)*1.15), ylim=c(0.5, length(y_h)+0.5), yaxt="n", xlab="Summed Probability",
       ylab="", main="Incoming signaling (Receiver)", frame.plot=FALSE)
  axis(2, at=y_h, labels=short_context_label(in_df$group, 30L), las=1, cex.axis=0.8)
  ezhu_light_grid(nx=NULL, ny=NA)
  segments(0, y_h, x_h, y_h, col="grey60", lwd=1.5)
  points(x_h, y_h, pch=21, bg=pal_ok["green"], col="white", cex=1.6)
  panel_label("H")
}

ezhu_write_stage_metadata(
  "07_figures_figure7",
  params = list(
    formats = TARGET_FORMATS,
    inputs = list(ccc_dir = Sys.getenv("EZHU_CCC_DIR", unset = "auto-latest"))
  ),
  seed = seed
)

ezhu_save_plot(
  function() draw_figure7_ccc(),
  "figure7_cellchat_ccc",
  height = 11
)
