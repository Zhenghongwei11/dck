ezhu_save_plot <- function(draw_fn, stem, width = 11, height = 7) {
  for (fmt in TARGET_FORMATS) {
    if (fmt == "pdf") {
      out_path <- file.path("plots/publication", paste0(stem, ".pdf"))
      message("Generating PDF: ", out_path)
      grDevices::pdf(out_path, width = width, height = height)
    } else if (fmt == "png") {
      out_path <- file.path("plots/publication/png", paste0(stem, ".png"))
      message("Generating PNG: ", out_path)
      grDevices::png(out_path, width = width, height = height, units = "in", res = 300)
    }

    tryCatch(draw_fn(), finally = grDevices::dev.off())
  }
}

parse_fig_list <- function(x) {
  x <- trimws(as.character(x %||% ""))
  if (!nzchar(x)) return(character())
  parts <- unlist(strsplit(x, "[,;\\s]+"))
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

panel_label <- function(letter, x_shift_lines = 0, line = 1.0, cex = 1.8) {
  letter <- as.character(letter)
  if (!nzchar(letter)) return(invisible(NULL))

  # Publication quality: Top-left alignment.
  # Automatically compensate for wide left margins to keep labels aligned vertically
  # across the column and closer to the outer edge (above text) for wide-margin plots.
  
  current_mar <- par("mar")[2]
  # We define a "standard" margin (e.g., 5 lines) as the baseline anchor.
  # If the margin is wider (e.g., 18), we shift the label left.
  target_mar <- 5 
  
  x_pos <- par("usr")[1] # Default: left edge of plot
  
  if (current_mar > target_mar || x_shift_lines != 0) {
    # Calculate width of one 'line' in user coordinates
    # Note: this conversion depends on the current plot scale
    line_width_user <- diff(grconvertX(c(0, 1), "lines", "user"))
    
    # Shift left by the excess margin lines
    shift_lines <- (current_mar - target_mar) + x_shift_lines
    x_pos <- x_pos - (shift_lines * line_width_user)
  }

  graphics::mtext(
    text = letter,
    side = 3,
    line = line,
    at = x_pos,
    adj = 0,
    cex = cex,
    font = 2,
    xpd = NA
  )
}

sanitize_ascii <- function(x) {
  if (!is.character(x)) x <- as.character(x)
  out <- suppressWarnings(iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "-"))
  out[is.na(out)] <- "-"
  out
}

okabe_ito_palette <- function() {
  c(
    blue = "#0072B2",
    orange = "#D55E00",
    green = "#009E73",
    sky = "#56B4E9",
    purple = "#CC79A7",
    yellow = "#F0E442",
    grey = "#4D4D4D",
    light_grey = "#D9D9D9"
  )
}

ezhu_seq_pal <- function(n = 64L) {
  grDevices::colorRampPalette(c("#F7FBFF", "#9ECAE1", "#3182BD"))(n)
}

ezhu_div_pal <- function(n = 64L) {
  grDevices::colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(n)
}

ezhu_light_grid <- function(nx = NULL, ny = NULL) {
  graphics::grid(nx = nx, ny = ny, col = "#F0F0F0", lty = 1)
}

short_sample_label <- function(x) {
  x <- as.character(x)
  x <- sub("\\.dgecounts\\.rds(\\.gz)?$", "", x)
  x <- sub("\\.rds(\\.gz)?$", "", x)
  x <- sub("^GSM[0-9]+_", "", x)
  x <- gsub("_", "-", x, fixed = TRUE)
  sanitize_ascii(x)
}

short_context_label <- function(x, max_chars = 40L) {
  x <- as.character(x)
  x <- gsub("[[:space:]]+/[[:space:]]+", "/", x)
  x <- gsub("[[:space:]]+", " ", x)
  x <- gsub("Collecting duct", "CD", x, fixed = TRUE)
  x <- gsub("Thick ascending limb", "TAL", x, fixed = TRUE)
  x <- gsub("loop of Henle", "Loop", x, fixed = TRUE)
  x <- gsub("intercalated", "IC", x, fixed = TRUE)
  x <- gsub("principal", "PC", x, fixed = TRUE)
  x <- gsub("Epithelial", "Epi", x, fixed = TRUE)
  x <- gsub("progenitor-like", "prog-like", x, fixed = TRUE)
  x <- gsub("provisional", "prov.", x, fixed = TRUE)
  x <- sanitize_ascii(x)
  ifelse(nchar(x) > max_chars, paste0(substr(x, 1, max_chars - 3), "..."), x)
}

wrap_label <- function(x, width = 18L) {
  x <- sanitize_ascii(x)
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

format_p_compact <- function(p) {
  p <- safe_numeric(p)
  ifelse(!is.finite(p), NA_character_, formatC(p, format = "e", digits = 2))
}

format_fdr_compact <- function(fdr) {
  fdr <- safe_numeric(fdr)
  ifelse(!is.finite(fdr), NA_character_, formatC(fdr, format = "f", digits = 3))
}

make_cluster_palette <- function(clusters) {
  clusters <- as.character(clusters)
  clusters <- clusters[nzchar(clusters)]
  lev <- unique(clusters)
  as_int <- suppressWarnings(as.integer(lev))
  if (all(is.finite(as_int))) {
    lev <- lev[order(as_int)]
  } else {
    lev <- sort(lev)
  }
  cols <- grDevices::hcl.colors(length(lev), palette = "Dark 3")
  setNames(cols, lev)
}

ezhu_placeholder_plot <- function(title, lines = character()) {
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  title(main = title)
  if (length(lines) > 0) {
    lines <- sanitize_ascii(lines)
    text(0.05, 0.95, paste(lines, collapse = "\n"), adj = c(0, 1), cex = 0.9)
  }
  invisible(NULL)
}

read_pvalue_table <- function(path) {
  if (!file.exists(path)) stop("Missing table: ", path, call. = FALSE)
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"pvalue" %in% names(df)) stop("Missing column 'pvalue' in: ", path, call. = FALSE)

  celltype_col <- NULL
  if ("celltype" %in% names(df)) celltype_col <- "celltype"
  if (is.null(celltype_col) && ncol(df) >= 2) celltype_col <- names(df)[2]
  if (is.null(celltype_col)) stop("Cannot locate celltype column in: ", path, call. = FALSE)

  df$celltype <- as.character(df[[celltype_col]])
  df$pvalue <- safe_numeric(df$pvalue)
  if ("fdr" %in% names(df)) df$fdr <- safe_numeric(df$fdr)
  df <- df[is.finite(df$pvalue) & nzchar(df$celltype), , drop = FALSE]
  df
}

read_singlecell_qc <- function(path = "results/singlecell/gse131882_counts_qc.tsv") {
  if (!file.exists(path)) stop("Missing table: ", path, call. = FALSE)
  df <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  df$sample <- as.character(df$sample)
  df$group <- as.character(df$group)
  df$n_genes <- safe_numeric(df$n_genes)
  df$n_cells <- safe_numeric(df$n_cells)
  if (mean(grepl("\\.rds$", df$sample)) < 0.9) df$sample <- paste0(df$sample, ".rds")
  df <- df[nzchar(df$sample) & nzchar(df$group) & is.finite(df$n_cells), , drop = FALSE]
  df
}

read_qc_cells_anchor <- function(path = "results/figures/gse131882_qc_cells.tsv") {
  if (!file.exists(path)) stop("Missing QC anchor table: ", path, call. = FALSE)
  df <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  need <- c("cell", "nCount_RNA", "nFeature_RNA", "percent.mt", "sample", "group")
  if (!all(need %in% names(df))) {
    stop("QC anchor table missing required columns: ", paste(setdiff(need, names(df)), collapse = ", "), call. = FALSE)
  }
  df$cell <- as.character(df$cell)
  df$sample <- as.character(df$sample)
  df$group <- as.character(df$group)
  df$nCount_RNA <- safe_numeric(df$nCount_RNA)
  df$nFeature_RNA <- safe_numeric(df$nFeature_RNA)
  df$percent.mt <- safe_numeric(df[["percent.mt"]])
  if ("seurat_clusters" %in% names(df)) df$seurat_clusters <- as.character(df$seurat_clusters)

  # Prefer a sample label consistent with scPagwas outputs (sample stub + ".rds").
  derived <- derive_sample_from_cell(df$cell)
  if (any(nzchar(derived))) {
    derived_has_rds <- mean(grepl("\\.rds$", derived)) > 0.9
    sample_has_rds <- mean(grepl("\\.rds$", df$sample)) > 0.9
    if (derived_has_rds && !sample_has_rds) df$sample <- derived
  }

  df <- df[nzchar(df$cell) & nzchar(df$sample) & nzchar(df$group) & is.finite(df$nCount_RNA) & is.finite(df$nFeature_RNA), , drop = FALSE]
  df
}

read_sample_index <- function(path = "results/singlecell/gse131882_sample_index.tsv") {
  if (!file.exists(path)) stop("Missing table: ", path, call. = FALSE)
  df <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  df$sample <- as.character(df$sample)
  df$group <- as.character(df$group)
  df$path <- as.character(df$path)
  df$size_bytes <- safe_numeric(df$size_bytes)
  df$md5 <- as.character(df$md5)
  df
}

derive_sample_from_cell <- function(cell) {
  sample_stub <- sub("_[^_]+$", "", cell)
  has_rds <- grepl("\\.rds$", sample_stub)
  ifelse(has_rds, sample_stub, paste0(sample_stub, ".rds"))
}

derive_group_from_sample <- function(sample) {
  is_control <- grepl("control", sample, fixed = TRUE)
  is_diabetes <- grepl("diabetes", sample, fixed = TRUE) | grepl("dn", sample, ignore.case = TRUE)
  ifelse(is_control, "control", ifelse(is_diabetes, "diabetes", NA_character_))
}

read_umap_cells <- function(path = "results/figures/gse131882_umap_cells.tsv") {
  if (!file.exists(path)) stop("Missing UMAP anchor table: ", path, call. = FALSE)
  df <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!all(c("cell", "UMAP_1", "UMAP_2") %in% names(df))) {
    stop("UMAP anchor table missing required columns: cell, UMAP_1, UMAP_2", call. = FALSE)
  }
  df$cell <- as.character(df$cell)
  df$UMAP_1 <- safe_numeric(df$UMAP_1)
  df$UMAP_2 <- safe_numeric(df$UMAP_2)
  if ("cluster" %in% names(df)) df$cluster <- as.character(df$cluster)

  df$sample <- derive_sample_from_cell(df$cell)
  if ("group" %in% names(df)) df$group <- as.character(df$group)
  if (!"group" %in% names(df) || all(!nzchar(df$group))) df$group <- derive_group_from_sample(df$sample)

  df <- df[nzchar(df$cell) & is.finite(df$UMAP_1) & is.finite(df$UMAP_2), , drop = FALSE]
  df
}

read_umap_features <- function(path = "results/figures/gse131882_umap_features.tsv") {
  if (!file.exists(path)) stop("Missing UMAP feature anchor table: ", path, call. = FALSE)
  df <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"cell" %in% names(df)) stop("UMAP feature table missing column: cell", call. = FALSE)
  df$cell <- as.character(df$cell)
  for (col in setdiff(names(df), "cell")) df[[col]] <- safe_numeric(df[[col]])
  df <- df[nzchar(df$cell), , drop = FALSE]
  df
}

read_cluster_by_group <- function(path = "results/scpagwas/repeat_run/cluster_by_group_counts.tsv") {
  if (!file.exists(path)) stop("Missing table: ", path, call. = FALSE)
  df <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  df$cluster <- as.character(df$cluster)
  df$group <- as.character(df$group)
  df$n_cells <- safe_numeric(df$n_cells)
  df <- df[nzchar(df$cluster) & nzchar(df$group) & is.finite(df$n_cells), , drop = FALSE]
  df
}

read_cell_scores <- function(path) {
  if (!file.exists(path)) {
    alts <- character()
    if (grepl("singlecell_results\\.csv$", path)) {
      alts <- c(
        sub("singlecell_results\\.csv$", "cell_scores.tsv", path),
        file.path(dirname(path), "cell_scores.tsv")
      )
    }
    alts <- unique(alts[file.exists(alts)])
    if (length(alts) > 0) {
      path <- alts[1]
    } else {
      stop("Missing score table: ", path, call. = FALSE)
    }
  }

  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("tsv", "txt")) {
    df <- utils::read.delim(path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }

  if (ncol(df) < 1) stop("Empty score table: ", path, call. = FALSE)
  if (!nzchar(names(df)[1])) names(df)[1] <- "cell"
  if (!"cell" %in% names(df)) df$cell <- as.character(df[[1]])
  df$cell <- as.character(df$cell)

  if ("sample" %in% names(df)) df$sample <- as.character(df$sample)
  if ("group" %in% names(df)) df$group <- as.character(df$group)
  if ("cluster" %in% names(df)) df$cluster <- as.character(df$cluster)

  if (!"sample" %in% names(df) || !nzchar(df$sample[1] %||% "")) {
    df$sample <- derive_sample_from_cell(df$cell)
  }
  if (!"group" %in% names(df) || !nzchar(df$group[1] %||% "")) {
    df$group <- derive_group_from_sample(df$sample)
  }
  if (mean(grepl("\\.rds$", df$sample)) < 0.9) {
    df$sample <- paste0(df$sample, ".rds")
  }

  df$scPagwas.TRS.Score <- safe_numeric(df[["scPagwas.TRS.Score"]] %||% df[["scPagwas.TRS.Score1"]] %||% NA_real_)
  df$scPagwas.gPAS.score <- safe_numeric(df[["scPagwas.gPAS.score"]] %||% NA_real_)
  df$Random_Correct_BG_p <- safe_numeric(df[["Random_Correct_BG_p"]] %||% NA_real_)
  df$Random_Correct_BG_adjp <- safe_numeric(df[["Random_Correct_BG_adjp"]] %||% NA_real_)
  df$Random_Correct_BG_z <- safe_numeric(df[["Random_Correct_BG_z"]] %||% NA_real_)
  df
}

cluster_counts_from_cell_scores <- function(cell_scores) {
  if (is.null(cell_scores) || !"cluster" %in% names(cell_scores)) return(integer())
  cl <- as.character(cell_scores$cluster)
  cl <- cl[nzchar(cl)]
  sort(table(cl), decreasing = TRUE)
}

write_sample_retention_anchor <- function(qc_df, retained_counts, out_path = "results/figures/gse131882_sample_retention.tsv") {
  qc_df$raw_cells <- qc_df$n_cells
  qc_df$retained_cells <- as.numeric(retained_counts[qc_df$sample])
  qc_df$retained_cells[is.na(qc_df$retained_cells)] <- 0
  qc_df$included_in_atlas <- qc_df$retained_cells > 0
  qc_df$pct_retained <- ifelse(qc_df$raw_cells > 0, qc_df$retained_cells / qc_df$raw_cells, NA_real_)
  out <- qc_df[, c("sample", "group", "raw_cells", "retained_cells", "pct_retained", "included_in_atlas", "n_genes")]
  utils::write.table(out, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  out
}

write_manifest_snapshot_anchor <- function(out_path = "results/figures/manifest_snapshot.tsv") {
  mf <- ezhu_read_manifest("data/manifest.tsv")
  notes <- lapply(as.character(mf$notes), ezhu_parse_notes_kv)
  build <- vapply(notes, function(x) x$build %||% "", character(1))
  ancestry <- vapply(notes, function(x) x$ancestry %||% "", character(1))

  out <- data.frame(
    id = as.character(mf$id),
    type = as.character(mf$type),
    name = as.character(mf$name),
    source = as.character(mf$source),
    accession = as.character(mf$accession),
    version_date = as.character(mf$version_date),
    local_path = as.character(mf$local_path),
    integrity = as.character(mf$integrity),
    build = build,
    ancestry = ancestry,
    stringsAsFactors = FALSE
  )
  utils::write.table(out, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  out
}
