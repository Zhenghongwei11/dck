#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

get_env_num <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  raw <- trimws(raw)
  if (!nzchar(raw)) return(default)
  val <- suppressWarnings(as.numeric(raw))
  if (!is.finite(val)) return(default)
  val
}

get_env_int <- function(name, default) {
  val <- get_env_num(name, default)
  as.integer(val)
}

cell_scores_path <- Sys.getenv("EZHU_SUBSAMPLING_CELL_SCORES", unset = "results/scpagwas/repeat_run/cell_scores.tsv")
subsample_frac <- get_env_num("EZHU_SUBSAMPLING_FRAC", 0.8)
n_iter <- get_env_int("EZHU_SUBSAMPLING_N", 200)
seed <- get_env_int("EZHU_SUBSAMPLING_SEED", 20251227)
skip_plots <- identical(Sys.getenv("EZHU_SUBSAMPLING_SKIP_PLOTS", unset = "0"), "1")

out_dir <- "results/scpagwas/robustness"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_stability <- file.path(out_dir, "bootstrap_cluster_stability.tsv")
out_corr <- file.path(out_dir, "bootstrap_replicate_correlations.tsv")

fig_pdf <- "plots/publication/figureS6_scpagwas_subsampling_stability.pdf"
fig_png <- "plots/publication/png/figureS6_scpagwas_subsampling_stability.png"
dir.create(dirname(fig_pdf), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(fig_png), recursive = TRUE, showWarnings = FALSE)

if (!file.exists(cell_scores_path)) {
  stop("Missing cell_scores input: ", cell_scores_path, call. = FALSE)
}

df <- data.table::fread(cell_scores_path, sep = "\t", data.table = FALSE, showProgress = FALSE)
trs_col <- intersect(c("scPagwas.TRS.Score1", "scPagwas.TRS.Score", "scPagwas.TRS.Mean"), colnames(df))
trs_col <- trs_col[1]
required <- c("sample", "cluster", trs_col)
missing <- setdiff(required, colnames(df))
if (length(missing) > 0) {
  stop("cell_scores.tsv missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

df <- df %>%
  mutate(
    sample = as.character(sample),
    cluster = as.character(cluster),
    trs = suppressWarnings(as.numeric(.data[[trs_col]]))
  ) %>%
  filter(!is.na(trs))

full_means <- df %>%
  group_by(cluster) %>%
  summarise(full_mean_trs = mean(trs, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(full_mean_trs)) %>%
  mutate(full_rank = rank(-full_mean_trs, ties.method = "first"))

clusters <- full_means$cluster
full_vec <- setNames(full_means$full_mean_trs, full_means$cluster)

set.seed(seed)

rank_mat <- matrix(NA_integer_, nrow = length(clusters), ncol = n_iter, dimnames = list(clusters, NULL))
top10_mat <- matrix(FALSE, nrow = length(clusters), ncol = n_iter, dimnames = list(clusters, NULL))

rep_cor <- data.frame(
  replicate = seq_len(n_iter),
  subsample_frac = rep(subsample_frac, n_iter),
  spearman_cor_mean_trs = rep(NA_real_, n_iter),
  ok = rep(FALSE, n_iter),
  stringsAsFactors = FALSE
)

split_by_sample <- split(seq_len(nrow(df)), df$sample)

for (i in seq_len(n_iter)) {
  idx <- unlist(lapply(split_by_sample, function(rows) {
    n <- length(rows)
    k <- max(1L, as.integer(floor(n * subsample_frac)))
    sample(rows, size = k, replace = FALSE)
  }), use.names = FALSE)

  sub <- df[idx, , drop = FALSE]
  means <- sub %>%
    group_by(cluster) %>%
    summarise(mean_trs = mean(trs, na.rm = TRUE), .groups = "drop")

  m <- setNames(means$mean_trs, means$cluster)
  # Ensure a value for every cluster; missing clusters get -Inf and rank last.
  m_full <- setNames(rep(-Inf, length(clusters)), clusters)
  m_full[names(m)] <- m

  ranks <- rank(-m_full, ties.method = "first")
  rank_mat[, i] <- as.integer(ranks[clusters])
  top10_mat[, i] <- ranks[clusters] <= 10

  # Spearman correlation vs full means (exclude missing clusters: those set to -Inf).
  ok_clusters <- names(m)
  ok <- length(ok_clusters) >= 2
  if (ok) {
    cor_val <- suppressWarnings(cor(m[ok_clusters], full_vec[ok_clusters], method = "spearman"))
    rep_cor$spearman_cor_mean_trs[i] <- cor_val
    rep_cor$ok[i] <- is.finite(cor_val)
  }
}

stability <- tibble::tibble(
  cluster = clusters,
  full_mean_trs = as.numeric(full_vec[clusters]),
  full_rank = as.integer(full_means$full_rank[match(clusters, full_means$cluster)]),
  bootstrap_top10_freq = rowMeans(top10_mat),
  bootstrap_rank_median = apply(rank_mat, 1, stats::median, na.rm = TRUE),
  bootstrap_rank_q25 = apply(rank_mat, 1, stats::quantile, probs = 0.25, na.rm = TRUE),
  bootstrap_rank_q75 = apply(rank_mat, 1, stats::quantile, probs = 0.75, na.rm = TRUE)
) %>%
  arrange(full_rank)

data.table::fwrite(stability, out_stability, sep = "\t", quote = FALSE, na = "NA")
data.table::fwrite(rep_cor, out_corr, sep = "\t", quote = FALSE, na = "NA")

  # Figure: left = Rank stability (Point + Error Bar); right = correlation histogram.
plot_one <- function() {
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  par(mfrow = c(1, 2), mar = c(8, 4, 3, 1))

  top_n <- min(20, nrow(stability))
  plot_df <- stability[seq_len(top_n), , drop = FALSE]
  
  # Color points by Top-10 frequency stability
  freq <- plot_df$bootstrap_top10_freq
  pal <- grDevices::colorRampPalette(c("grey80", "#d55e00"))(100)
  col_idx <- as.integer(freq * 99) + 1
  cols <- pal[col_idx]

  # Plot: Median Rank with IQR error bars
  # Note: Y-axis reversed because Rank 1 is "top"
  y_min <- min(plot_df$bootstrap_rank_q25, 1)
  y_max <- max(plot_df$bootstrap_rank_q75, 20)
  
  plot(NA, xlim = c(0.5, top_n + 0.5), ylim = c(y_max + 1, 0), # Reversed Y
       xaxt = "n", xlab = "", ylab = "Bootstrap Rank (Median \u00B1 IQR)",
       main = "Rank Stability (Subsampling)")
  
  axis(1, at = 1:top_n, labels = paste0("C", plot_df$cluster), las = 2, cex.axis = 0.8)
  abline(h = c(1, 5, 10, 15), col = "grey90", lty = 2)
  
  # Error bars (Q25 to Q75)
  segments(1:top_n, plot_df$bootstrap_rank_q25, 1:top_n, plot_df$bootstrap_rank_q75, col = "grey50", lwd = 1.5)
  
  # Median points
  points(1:top_n, plot_df$bootstrap_rank_median, pch = 21, bg = cols, col = "grey30", cex = 1.5, lwd = 0.8)
  
  legend("bottomright", legend = c("High Stability", "Low Stability"), 
         pt.bg = c("#d55e00", "grey80"), pch = 21, bty = "n", cex = 0.8, title = "Freq. in Top 10")

  vals <- rep_cor$spearman_cor_mean_trs[is.finite(rep_cor$spearman_cor_mean_trs)]
  hist(
    vals,
    breaks = 20,
    col = "#9ECAE1", # Nicer blue
    border = "white",
    main = "Correlation vs full ranking",
    xlab = "Spearman correlation (mean TRS)"
  )
  med <- stats::median(vals, na.rm = TRUE)
  abline(v = med, col = "#d55e00", lwd = 2, lty=2)
  mtext(sprintf("median = %.3f", med), side = 3, line = -1.5, adj = 0.05, cex = 0.85)
}
if (!isTRUE(skip_plots)) {
  pdf(fig_pdf, width = 10, height = 4.5, onefile = TRUE)
  plot_one()
  dev.off()

  png(fig_png, width = 2000, height = 900, res = 200)
  plot_one()
  dev.off()
}

message("Wrote:")
message("- ", out_stability)
message("- ", out_corr)
if (!isTRUE(skip_plots)) {
  message("- ", fig_pdf)
  message("- ", fig_png)
}
