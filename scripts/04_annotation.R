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

options(stringsAsFactors = FALSE, width = 120)

if (!requireNamespace("Seurat", quietly = TRUE)) stop("Missing package: Seurat", call. = FALSE)
suppressPackageStartupMessages(library(Seurat))

parse_csv_env <- function(env_value) {
  value <- trimws(env_value)
  if (!nzchar(value)) return(character())
  parts <- unlist(strsplit(value, "[,;]"))
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

annotation_dir <- Sys.getenv("EZHU_ANNOTATION_DIR", unset = "results/annotation")
ezhu_dir_create(annotation_dir)
markers_dir <- Sys.getenv("EZHU_MARKERS_DIR", unset = file.path(annotation_dir, "markers"))
ezhu_dir_create(markers_dir)

scpagwas_rds <- Sys.getenv("EZHU_SCPAGWAS_RDS", unset = "data/processed/singlecell_gse131882_scpagwas_gse131882_gcst90435703.rds")
if (!file.exists(scpagwas_rds)) stop("Missing scPagwas Seurat object: ", scpagwas_rds, call. = FALSE)

message("Reading scPagwas object: ", scpagwas_rds)
obj <- readRDS(scpagwas_rds)
if (!inherits(obj, "Seurat")) stop("Expected a Seurat object: ", scpagwas_rds, call. = FALSE)

idents_field <- Sys.getenv("EZHU_IDENTS_FIELD", unset = "seurat_clusters")
if (idents_field %in% colnames(obj[[]])) {
  Seurat::Idents(obj) <- as.factor(obj[[idents_field, drop = TRUE]])
}

clusters_env <- Sys.getenv("EZHU_ANNOT_CLUSTERS", unset = "")
clusters <- parse_csv_env(clusters_env)
if (length(clusters) == 0) {
  clusters <- as.character(sort(unique(Seurat::Idents(obj))))
}

min_pct <- suppressWarnings(as.numeric(Sys.getenv("EZHU_MARKERS_MIN_PCT", unset = "0.25")))
if (is.na(min_pct) || min_pct < 0) min_pct <- 0.25

logfc_threshold <- suppressWarnings(as.numeric(Sys.getenv("EZHU_MARKERS_LOGFC_THRESHOLD", unset = "0.25")))
if (is.na(logfc_threshold)) logfc_threshold <- 0.25

only_pos <- !identical(Sys.getenv("EZHU_MARKERS_ONLY_POS", unset = "1"), "0")
test_use <- Sys.getenv("EZHU_MARKERS_TEST_USE", unset = "wilcox")
test_use <- trimws(test_use)
if (!nzchar(test_use)) test_use <- "wilcox"

message("Annotating clusters: ", paste(clusters, collapse = ", "))

message(
  "FindAllMarkers: only.pos=", only_pos,
  " min.pct=", min_pct,
  " logfc.threshold=", logfc_threshold,
  " test.use=", test_use
)
markers <- Seurat::FindAllMarkers(
  object = obj,
  only.pos = only_pos,
  min.pct = min_pct,
  logfc.threshold = logfc_threshold,
  test.use = test_use
)

names(markers) <- tolower(names(markers))
if (!"cluster" %in% names(markers)) stop("FindAllMarkers output missing column: cluster", call. = FALSE)
if (!"gene" %in% names(markers)) {
  markers$gene <- rownames(markers)
}
if ("avg_logfc" %in% names(markers) && !"avg_log2fc" %in% names(markers)) {
  markers$avg_log2fc <- markers$avg_logfc
}

markers$cluster <- as.character(markers$cluster)
markers$gene <- as.character(markers$gene)

if (length(clusters) > 0) {
  markers <- markers[markers$cluster %in% clusters, , drop = FALSE]
}

required_cols <- c("p_val", "avg_log2fc", "pct.1", "pct.2", "p_val_adj", "cluster", "gene")
missing_cols <- setdiff(required_cols, names(markers))
if (length(missing_cols) > 0) {
  stop("Marker table missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

detailed <- markers[, required_cols, drop = FALSE]
colnames(detailed)[colnames(detailed) == "avg_log2fc"] <- "avg_log2FC"

utils::write.csv(detailed, file = file.path(annotation_dir, "detailed_markers_top_risk.csv"), row.names = FALSE, quote = TRUE)

# Export per-cluster marker CSVs expected by downstream modules.
for (cl in sort(unique(as.character(detailed$cluster)))) {
  df <- detailed[detailed$cluster == cl, , drop = FALSE]
  df <- df[order(df$p_val_adj, df$p_val, decreasing = FALSE, na.last = TRUE), , drop = FALSE]
  utils::write.csv(
    df[, c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "cluster", "gene")],
    file = file.path(markers_dir, paste0("markers_cluster_", cl, ".csv")),
    row.names = FALSE,
    quote = TRUE
  )
}

top30 <- NULL
for (cl in unique(detailed$cluster)) {
  subset <- detailed[detailed$cluster == cl, , drop = FALSE]
  subset <- subset[order(subset$p_val_adj, subset$p_val, decreasing = FALSE, na.last = TRUE), , drop = FALSE]
  subset <- subset[seq_len(min(30L, nrow(subset))), , drop = FALSE]
  top30 <- rbind(top30, subset)
}
utils::write.csv(
  top30[, c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")],
  file = file.path(annotation_dir, "top30_markers_summary.csv"),
  row.names = FALSE,
  quote = TRUE
)

labeled <- detailed
labeled$label <- ""
utils::write.csv(
  labeled[, c("cluster", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "label")],
  file = file.path(annotation_dir, "top_risk_markers_labeled.csv"),
  row.names = FALSE,
  quote = TRUE
)

ezhu_write_stage_metadata(
  "04_annotation",
  params = list(
    scpagwas_rds = scpagwas_rds,
    idents_field = idents_field,
    annotation_dir = annotation_dir,
    clusters = clusters,
    markers = list(min_pct = min_pct, logfc_threshold = logfc_threshold),
    outputs = list(
      detailed_markers = file.path(annotation_dir, "detailed_markers_top_risk.csv"),
      markers_dir = markers_dir,
      top30_markers = file.path(annotation_dir, "top30_markers_summary.csv"),
      labeled_markers = file.path(annotation_dir, "top_risk_markers_labeled.csv")
    )
  ),
  seed = seed
)
