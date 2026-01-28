args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
  repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
  if (dir.exists(repo_root)) setwd(repo_root)
}

source("scripts/utils/repro.R")
ezhu_set_repo_root()
ezhu_activate_renv()

seed_env <- Sys.getenv("EZHU_SEED", unset = "")
seed <- suppressWarnings(as.integer(seed_env))
if (is.na(seed)) seed <- 20251227L
set.seed(seed)

options(stringsAsFactors = FALSE, width = 120, warn = 0)

if (!requireNamespace("Seurat", quietly = TRUE)) stop("Missing package: Seurat", call. = FALSE)
if (!requireNamespace("reshape2", quietly = TRUE)) stop("Missing package: reshape2", call. = FALSE)

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

read_pvalue_table <- function(path) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- tolower(names(df))
  if (!"pvalue" %in% names(df)) stop("Missing column 'pvalue' in: ", path, call. = FALSE)

  celltype_col <- NULL
  if ("celltype" %in% names(df)) celltype_col <- "celltype"
  if (is.null(celltype_col) && ncol(df) >= 2) celltype_col <- names(df)[2]
  if (is.null(celltype_col)) stop("Cannot locate cluster/celltype column in: ", path, call. = FALSE)

  out <- data.frame(
    cluster = as.character(df[[celltype_col]]),
    pvalue = safe_numeric(df$pvalue),
    stringsAsFactors = FALSE
  )
  out <- out[nzchar(out$cluster) & is.finite(out$pvalue), , drop = FALSE]
  out
}

read_markers <- function(path, top_n = 100L) {
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- tolower(names(df))
  if (!"gene" %in% names(df)) stop("Marker table missing 'gene': ", path, call. = FALSE)
  if (!"avg_log2fc" %in% names(df)) stop("Marker table missing 'avg_log2FC': ", path, call. = FALSE)

  df$gene <- toupper(trimws(as.character(df$gene)))
  df$avg_log2fc <- safe_numeric(df$avg_log2fc)
  df <- df[nzchar(df$gene) & is.finite(df$avg_log2fc), , drop = FALSE]
  df <- df[order(df$avg_log2fc, decreasing = TRUE), , drop = FALSE]
  unique(utils::head(df$gene, top_n))
}

best_overlap <- function(rep_genes, primary_sets) {
  best <- list(
    primary_cluster = NA_character_,
    overlap_n = 0L,
    overlap_ratio = 0,
    jaccard = 0
  )
  rep_genes <- unique(rep_genes[nzchar(rep_genes) & !is.na(rep_genes)])
  if (length(rep_genes) == 0) return(best)

  for (cl in names(primary_sets)) {
    pg <- primary_sets[[cl]]
    pg <- unique(pg[nzchar(pg) & !is.na(pg)])
    if (length(pg) == 0) next
    inter <- length(intersect(rep_genes, pg))
    denom <- length(rep_genes) + length(pg) - inter
    jacc <- if (denom > 0) inter / denom else 0
    ratio <- inter / min(length(rep_genes), length(pg))
    if (ratio > best$overlap_ratio || (ratio == best$overlap_ratio && inter > best$overlap_n)) {
      best$primary_cluster <- cl
      best$overlap_n <- inter
      best$overlap_ratio <- ratio
      best$jaccard <- jacc
    }
  }
  best
}

output_dir <- Sys.getenv("EZHU_REPLICATION_DIR", unset = "results/replication")
ezhu_dir_create(output_dir)

primary_ann_path <- Sys.getenv("EZHU_PRIMARY_ANNOTATION_TSV", unset = "results/annotation/cluster_annotations_filled.tsv")
primary_ann <- NULL
if (file.exists(primary_ann_path)) {
  primary_ann <- utils::read.delim(primary_ann_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  primary_ann$cluster <- as.character(primary_ann$cluster)
  if (!all(c("cluster", "proposed_cell_type", "proposed_state") %in% names(primary_ann))) {
    warning("Primary annotation record missing expected columns; proceeding without labels: ", primary_ann_path, call. = FALSE)
    primary_ann <- NULL
  }
}

primary_main_pvalue_path <- Sys.getenv("EZHU_PRIMARY_MAIN_PVALUE", unset = "results/scpagwas/main/merged_celltype_pvalue.csv")
primary_repeat_pvalue_path <- Sys.getenv("EZHU_PRIMARY_REPEAT_PVALUE", unset = "results/scpagwas/repeat_run/merged_celltype_pvalue.csv")
if (!file.exists(primary_main_pvalue_path)) stop("Missing primary main pvalue table: ", primary_main_pvalue_path, call. = FALSE)
if (!file.exists(primary_repeat_pvalue_path)) stop("Missing primary repeat pvalue table: ", primary_repeat_pvalue_path, call. = FALSE)

primary_main <- read_pvalue_table(primary_main_pvalue_path)
primary_repeat <- read_pvalue_table(primary_repeat_pvalue_path)
primary_pvals <- merge(primary_main, primary_repeat, by = "cluster", suffixes = c("_main", "_repeat"))
if (nrow(primary_pvals) == 0) stop("No overlapping clusters between main and repeat pvalue tables.", call. = FALSE)
primary_pvals$score_pvalue <- pmin(primary_pvals$pvalue_main, primary_pvals$pvalue_repeat)
primary_pvals <- primary_pvals[order(primary_pvals$score_pvalue, decreasing = FALSE), , drop = FALSE]

primary_top <- suppressWarnings(as.integer(Sys.getenv("EZHU_PRIMARY_TOP_CLUSTERS", unset = "8")))
if (is.na(primary_top) || primary_top <= 0) primary_top <- 8L
primary_focus <- utils::head(primary_pvals, primary_top)
primary_clusters <- as.character(primary_focus$cluster)

primary_markers_dir <- Sys.getenv("EZHU_PRIMARY_MARKERS_DIR", unset = "results/annotation/markers")
primary_marker_top_n <- suppressWarnings(as.integer(Sys.getenv("EZHU_PRIMARY_MARKER_TOP_N", unset = "100")))
if (is.na(primary_marker_top_n) || primary_marker_top_n <= 0) primary_marker_top_n <- 100L

primary_clusters <- unique(primary_clusters[nzchar(primary_clusters)])
if (length(primary_clusters) == 0) stop("No primary clusters selected from scPagwas pvalue tables.", call. = FALSE)

primary_marker_sets <- stats::setNames(vector("list", length(primary_clusters)), primary_clusters)
for (i in seq_along(primary_clusters)) {
  cl <- primary_clusters[i]
  path <- file.path(primary_markers_dir, paste0("markers_cluster_", cl, ".csv"))
  if (!file.exists(path)) stop("Missing primary marker file: ", path, call. = FALSE)
  primary_marker_sets[[i]] <- read_markers(path, top_n = primary_marker_top_n)
}

rep_obj_path <- Sys.getenv("EZHU_REP_SINGLECELL_PROCESSED", unset = "data/processed/singlecell_gse195460_seurat_v4.rds")
if (!file.exists(rep_obj_path)) stop("Missing replication Seurat object: ", rep_obj_path, call. = FALSE)

rep_pvalue_path <- Sys.getenv("EZHU_REP_SCPAGWAS_PVALUE", unset = "")
if (!nzchar(rep_pvalue_path)) {
  candidates <- Sys.glob("results/scpagwas/gse195460/*Merged_celltype_pvalue.csv")
  if (length(candidates) == 0) stop("No replication pvalue table found under results/scpagwas/gse195460/.", call. = FALSE)
  info <- file.info(candidates)
  candidates <- rownames(info)[order(info$mtime, decreasing = TRUE)]
  rep_pvalue_path <- candidates[1]
}
if (!file.exists(rep_pvalue_path)) stop("Missing replication pvalue table: ", rep_pvalue_path, call. = FALSE)

rep_counts_path <- Sys.getenv("EZHU_REP_CLUSTER_COUNTS", unset = "results/scpagwas/gse195460/cluster_by_group_counts.tsv")
rep_counts <- NULL
if (file.exists(rep_counts_path)) {
  rep_counts <- utils::read.delim(rep_counts_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  rep_counts$cluster <- as.character(rep_counts$cluster)
  rep_counts$group <- as.character(rep_counts$group)
  rep_counts$n_cells <- safe_numeric(rep_counts$n_cells)
  rep_counts <- rep_counts[nzchar(rep_counts$cluster) & nzchar(rep_counts$group) & is.finite(rep_counts$n_cells), , drop = FALSE]
}

message("Reading replication Seurat object: ", rep_obj_path)
obj <- readRDS(rep_obj_path)
if (!inherits(obj, "Seurat")) stop("Replication object is not a Seurat object: ", rep_obj_path, call. = FALSE)

assay <- Sys.getenv("EZHU_REP_ASSAY", unset = "RNA")
if (!assay %in% Seurat::Assays(obj)) stop("Assay not present in replication Seurat object: ", assay, call. = FALSE)

idents_field <- Sys.getenv("EZHU_REP_IDENTS_FIELD", unset = "seurat_clusters")
if (idents_field %in% colnames(obj[[]])) {
  Seurat::Idents(obj) <- as.factor(obj[[idents_field, drop = TRUE]])
} else {
  stop("Missing idents field in replication meta.data: ", idents_field, call. = FALSE)
}

rep_marker_top_n <- suppressWarnings(as.integer(Sys.getenv("EZHU_REP_MARKER_TOP_N", unset = "100")))
if (is.na(rep_marker_top_n) || rep_marker_top_n <= 0) rep_marker_top_n <- 100L
rep_marker_min_pct <- suppressWarnings(as.numeric(Sys.getenv("EZHU_REP_MARKER_MIN_PCT", unset = "0.1")))
if (is.na(rep_marker_min_pct) || rep_marker_min_pct < 0) rep_marker_min_pct <- 0.1
rep_marker_logfc <- suppressWarnings(as.numeric(Sys.getenv("EZHU_REP_MARKER_LOGFC", unset = "0.25")))
if (is.na(rep_marker_logfc) || rep_marker_logfc < 0) rep_marker_logfc <- 0.25

min_overlap <- suppressWarnings(as.integer(Sys.getenv("EZHU_REP_MIN_MARKER_OVERLAP", unset = "3")))
if (is.na(min_overlap) || min_overlap < 1) min_overlap <- 3L

message("Computing replication markers (FindAllMarkers; assay=", assay, ", min.pct=", rep_marker_min_pct, ", logfc=", rep_marker_logfc, ")")
markers <- Seurat::FindAllMarkers(
  object = obj,
  assay = assay,
  only.pos = TRUE,
  logfc.threshold = rep_marker_logfc,
  min.pct = rep_marker_min_pct,
  test.use = "wilcox"
)
names(markers) <- tolower(names(markers))
if (!"cluster" %in% names(markers) || !"gene" %in% names(markers) || !"avg_log2fc" %in% names(markers)) {
  stop("Unexpected FindAllMarkers output; require columns: cluster, gene, avg_log2FC.", call. = FALSE)
}
markers$cluster <- as.character(markers$cluster)
markers$gene <- toupper(trimws(as.character(markers$gene)))
markers$avg_log2fc <- safe_numeric(markers$avg_log2fc)
markers <- markers[nzchar(markers$cluster) & nzchar(markers$gene) & is.finite(markers$avg_log2fc), , drop = FALSE]

rep_marker_sets <- split(markers, markers$cluster)
rep_marker_sets <- lapply(rep_marker_sets, function(df) {
  df <- df[order(df$avg_log2fc, decreasing = TRUE), , drop = FALSE]
  unique(utils::head(df$gene, rep_marker_top_n))
})

rep_pvalues <- read_pvalue_table(rep_pvalue_path)

mapping_rows <- lapply(names(rep_marker_sets), function(rep_cluster) {
  rep_genes <- rep_marker_sets[[rep_cluster]]
  best <- best_overlap(rep_genes, primary_marker_sets)

  mapped <- as.character(best$primary_cluster %||% "")
  if (!nzchar(mapped) || best$overlap_n < min_overlap) {
    mapped <- NA_character_
  }

  primary_row <- NULL
  if (!is.null(primary_ann) && !is.na(mapped)) {
    primary_row <- primary_ann[primary_ann$cluster == mapped, , drop = FALSE]
    if (nrow(primary_row) == 0) primary_row <- NULL
  }

  pv_row <- primary_pvals[primary_pvals$cluster == mapped, , drop = FALSE]
  primary_main_p <- if (nrow(pv_row) > 0) pv_row$pvalue_main[1] else NA_real_
  primary_repeat_p <- if (nrow(pv_row) > 0) pv_row$pvalue_repeat[1] else NA_real_

  rep_row <- rep_pvalues[rep_pvalues$cluster == rep_cluster, , drop = FALSE]
  rep_p <- if (nrow(rep_row) > 0) rep_row$pvalue[1] else NA_real_

  data.frame(
    rep_cluster = rep_cluster,
    rep_pvalue = rep_p,
    mapped_primary_cluster = mapped,
    mapped_cell_type = if (!is.null(primary_row)) primary_row$proposed_cell_type[1] %||% NA_character_ else NA_character_,
    mapped_state = if (!is.null(primary_row)) primary_row$proposed_state[1] %||% NA_character_ else NA_character_,
    marker_overlap_n = best$overlap_n,
    marker_overlap_ratio = best$overlap_ratio,
    marker_jaccard = best$jaccard,
    primary_main_pvalue = safe_numeric(primary_main_p),
    primary_repeat_pvalue = safe_numeric(primary_repeat_p),
    stringsAsFactors = FALSE
  )
})

cluster_map <- do.call(rbind, mapping_rows)
cluster_map$rep_neglog10p <- -log10(pmax(cluster_map$rep_pvalue, 1e-300))
cluster_map$rep_neglog10p[!is.finite(cluster_map$rep_neglog10p)] <- NA_real_

if (!is.null(rep_counts) && nrow(rep_counts) > 0) {
  totals <- aggregate(list(rep_n_cells_total = rep_counts$n_cells), by = list(rep_cluster = rep_counts$cluster), FUN = sum)
  cluster_map <- merge(cluster_map, totals, by = "rep_cluster", all.x = TRUE)

  wide <- reshape2::dcast(rep_counts, cluster ~ group, value.var = "n_cells", fun.aggregate = sum)
  names(wide)[1] <- "rep_cluster"
  cluster_map <- merge(cluster_map, wide, by = "rep_cluster", all.x = TRUE)
}

cluster_map <- cluster_map[order(cluster_map$rep_pvalue, decreasing = FALSE, na.last = TRUE), , drop = FALSE]

cluster_map_path <- file.path(output_dir, "gse195460_cluster_to_primary_context.tsv")
utils::write.table(cluster_map, file = cluster_map_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote: ", cluster_map_path)

context_summary <- cluster_map[!is.na(cluster_map$mapped_primary_cluster), , drop = FALSE]
context_summary <- split(context_summary, context_summary$mapped_primary_cluster)
context_summary <- lapply(context_summary, function(df) {
  df <- df[order(df$rep_pvalue, decreasing = FALSE, na.last = TRUE), , drop = FALSE]
  top <- df[1, , drop = FALSE]
  data.frame(
    primary_cluster = top$mapped_primary_cluster,
    primary_cell_type = top$mapped_cell_type,
    primary_state = top$mapped_state,
    primary_main_pvalue = top$primary_main_pvalue,
    primary_repeat_pvalue = top$primary_repeat_pvalue,
    rep_best_cluster = top$rep_cluster,
    rep_best_pvalue = top$rep_pvalue,
    rep_best_neglog10p = top$rep_neglog10p,
    rep_marker_overlap_n = top$marker_overlap_n,
    rep_marker_overlap_ratio = top$marker_overlap_ratio,
    stringsAsFactors = FALSE
  )
})
context_summary <- do.call(rbind, context_summary)
context_summary <- context_summary[order(context_summary$rep_best_pvalue, decreasing = FALSE, na.last = TRUE), , drop = FALSE]

context_path <- file.path(output_dir, "gse195460_concordance_context_summary.tsv")
utils::write.table(context_summary, file = context_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote: ", context_path)

ezhu_write_stage_metadata(
  "03c_replication_concordance",
  params = list(
    replication = list(
      tag = "gse195460",
      seurat_rds = rep_obj_path,
      pvalue_table = rep_pvalue_path,
      counts_table = if (file.exists(rep_counts_path)) rep_counts_path else NA_character_,
      assay = assay,
      idents_field = idents_field,
      marker_top_n = rep_marker_top_n,
      marker_min_pct = rep_marker_min_pct,
      marker_logfc_threshold = rep_marker_logfc,
      min_marker_overlap = min_overlap
    ),
    primary = list(
      annotation_tsv = primary_ann_path,
      markers_dir = primary_markers_dir,
      marker_top_n = primary_marker_top_n
    ),
    outputs = list(
      cluster_mapping_tsv = cluster_map_path,
      context_summary_tsv = context_path
    )
  ),
  seed = seed
)
