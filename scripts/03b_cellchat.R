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
if (is.na(seed)) seed <- 20260102L
set.seed(seed)

seurat_rds <- Sys.getenv("EZHU_SEURAT_RDS", unset = "")
if (!nzchar(seurat_rds)) {
  stop(
    "Missing EZHU_SEURAT_RDS.\n",
    "Example: EZHU_SEURAT_RDS=/path/to/singlecell_gse131882_seurat.rds",
    call. = FALSE
  )
}
if (!file.exists(seurat_rds)) stop("Missing Seurat RDS: ", seurat_rds, call. = FALSE)

group_col <- Sys.getenv("EZHU_CCC_GROUP_COL", unset = "cell_context")
min_cells <- suppressWarnings(as.integer(Sys.getenv("EZHU_CCC_MIN_CELLS", unset = "50")))
if (!is.finite(min_cells) || is.na(min_cells) || min_cells < 10) min_cells <- 50L
max_cells <- suppressWarnings(as.integer(Sys.getenv("EZHU_CCC_MAX_CELLS", unset = "20000")))
if (!is.finite(max_cells) || is.na(max_cells) || max_cells < 0) max_cells <- 20000L
run_id <- Sys.getenv("EZHU_CCC_RUN_ID", unset = format(Sys.time(), "%Y%m%d_%H%M%S"))

parse_csv <- function(x) {
  x <- trimws(as.character(x %||% ""))
  if (!nzchar(x)) return(character())
  parts <- unlist(strsplit(x, "[,;]"))
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

group_mode <- trimws(Sys.getenv("EZHU_CCC_GROUP_MODE", unset = ""))
only_annotated <- identical(Sys.getenv("EZHU_CCC_ONLY_ANNOTATED", unset = ""), "1") ||
  identical(tolower(group_mode), "annotated_only")
allowed_groups <- parse_csv(Sys.getenv("EZHU_CCC_ALLOWED_GROUPS", unset = ""))

ezhu_dir_create("results/ccc")
out_dir <- file.path("results/ccc", run_id)
ezhu_dir_create(out_dir)

ezhu_write_stage_metadata(
  "03b_cellchat",
  params = list(
    seurat_rds = seurat_rds,
    group_col = group_col,
    min_cells = min_cells,
    max_cells = max_cells,
    run_id = run_id,
    group_mode = group_mode,
    only_annotated = only_annotated,
    allowed_groups = allowed_groups
  ),
  seed = seed
)

require_or_stop <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Missing R package: ", pkg, "\n",
      "Install it in the active environment before running CellChat.",
      call. = FALSE
    )
  }
}

require_or_stop("Seurat")
require_or_stop("SeuratObject")
require_or_stop("CellChat")

suppressPackageStartupMessages({
  # CellChat (and some downstream code) expects Seurat generics like DefaultAssay()
  # to be on the search path, not just loaded via requireNamespace().
  library(Seurat)
  library(SeuratObject)
})

get_seurat_expr_matrix <- function(obj, prefer = c("data", "counts")) {
  assay <- tryCatch(DefaultAssay(obj), error = function(e) NULL)
  if (!nzchar(assay %||% "")) assay <- "RNA"

  prefer <- unique(as.character(prefer))
  for (layer in prefer) {
    m <- tryCatch(
      SeuratObject::GetAssayData(obj = obj, assay = assay, layer = layer),
      error = function(e) NULL
    )
    if (!is.null(m)) return(m)

    m <- tryCatch(
      SeuratObject::GetAssayData(obj = obj, assay = assay, slot = layer),
      error = function(e) NULL
    )
    if (!is.null(m)) return(m)
  }

  stop(
    "Could not extract an expression matrix from the Seurat object.\n",
    "Tried layers/slots: ", paste(prefer, collapse = ", "),
    call. = FALSE
  )
}

message("Loading Seurat: ", seurat_rds)
obj <- readRDS(seurat_rds)
if (!inherits(obj, "Seurat")) stop("Input is not a Seurat object: ", seurat_rds, call. = FALSE)

meta0 <- obj[[]]
if (!identical(rownames(meta0), colnames(obj))) {
  stop(
    "Seurat meta.data rownames do not match cell names in the object.\n",
    "This indicates an invalid Seurat object for downstream tools like CellChat.",
    call. = FALSE
  )
}

cluster_col <- Sys.getenv("EZHU_CLUSTER_COL", unset = "seurat_clusters")
if (!cluster_col %in% colnames(meta0)) {
  stop("Seurat meta.data missing cluster column: ", cluster_col, call. = FALSE)
}
cluster_ids <- as.character(obj[[cluster_col, drop = TRUE]])
cluster_ids[is.na(cluster_ids)] <- "NA"
cluster_ids <- trimws(cluster_ids)
cell_names <- rownames(meta0)

# Resolve a grouping column for CellChat.
# Default is a "cell_context" label, but we must ensure every cell has a non-empty group label.
#
# For "annotated_only" runs we only keep cells mapped by the cluster annotation file, regardless
# of any existing cell_context-like columns.
groups_existing <- NULL
if (group_col %in% colnames(meta0)) groups_existing <- as.character(meta0[[group_col]])
groups_existing <- as.character(groups_existing)
groups_existing[is.na(groups_existing)] <- ""
groups_existing <- trimws(groups_existing)
existing_has_any <- length(groups_existing) > 0 && any(nzchar(groups_existing))

mapping_path <- Sys.getenv(
  "EZHU_CLUSTER_ANNOTATIONS",
  unset = "results/annotation/cluster_annotations_filled.tsv"
)
map <- NULL
lookup <- NULL
if (file.exists(mapping_path)) {
  map <- utils::read.delim(mapping_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  names(map) <- tolower(names(map))
  if (all(c("cluster", "proposed_cell_type") %in% names(map))) {
    map$cluster <- trimws(as.character(map$cluster))
    map$proposed_cell_type <- trimws(as.character(map$proposed_cell_type))
    lookup <- stats::setNames(map$proposed_cell_type, map$cluster)
  }
}

use_cluster_annotations <- identical(tolower(group_mode), "from_cluster_annotations") ||
  only_annotated ||
  length(allowed_groups) > 0 ||
  !existing_has_any

if (use_cluster_annotations) {
  if (is.null(lookup)) {
    stop(
      "CellChat grouping requires a valid cluster annotation file, but none was found or parsed.\n",
      "Expected a TSV with columns: cluster, proposed_cell_type.\n",
      "Path: ", mapping_path,
      call. = FALSE
    )
  }
  groups <- as.character(lookup[cluster_ids])
} else {
  groups <- groups_existing
}

groups <- as.character(groups)
groups[is.na(groups)] <- ""
groups <- trimws(groups)

# Track which cells have a non-empty annotation label from cluster annotations.
has_annotation <- nzchar(groups)

# Fallback: for any cells still missing labels, use the cluster ID as a stable technical group.
# In annotated-only mode we do not create fallback groups.
if (!only_annotated) {
  missing <- !nzchar(groups)
  groups[missing] <- paste0("C", cluster_ids[missing])
}

obj[[group_col]] <- groups

# Group counts before any cell filtering (for auditability).
counts_raw <- sort(table(obj[[group_col, drop = TRUE]]), decreasing = TRUE)
counts_df_raw <- data.frame(group = names(counts_raw), n_cells = as.integer(counts_raw), stringsAsFactors = FALSE)
utils::write.table(
  counts_df_raw,
  file = file.path(out_dir, "cellchat_group_counts_raw.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Optional: keep only annotated contexts (drop fallback C* groups).
keep_mask <- rep(TRUE, length(cell_names))
if (only_annotated) {
  keep_mask <- keep_mask & has_annotation
}
if (length(allowed_groups) > 0) {
  keep_mask <- keep_mask & (groups %in% allowed_groups)
}
if (!all(keep_mask)) {
  keep_cells <- cell_names[keep_mask]
  obj <- subset(obj, cells = keep_cells)
}

counts_by_group <- sort(table(obj[[group_col, drop = TRUE]]), decreasing = TRUE)
counts_df <- data.frame(group = names(counts_by_group), n_cells = as.integer(counts_by_group), stringsAsFactors = FALSE)
utils::write.table(
  counts_df,
  file = file.path(out_dir, "cellchat_group_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

eligible_groups <- counts_df$group[counts_df$n_cells >= min_cells]
meta <- obj[[]]
meta_group <- as.character(meta[[group_col]])
meta_group[is.na(meta_group)] <- ""
meta_group <- trimws(meta_group)
keep_cells <- rownames(meta)[meta_group %in% eligible_groups]
obj <- subset(obj, cells = keep_cells)

if (max_cells > 0 && ncol(obj) > max_cells) {
  set.seed(seed)
  cells <- colnames(obj)
  sampled <- sample(cells, size = max_cells, replace = FALSE)
  obj <- subset(obj, cells = sampled)
}

if (ncol(obj) == 0) {
  stop(
    "No cells remain after applying group labeling and min-cells filtering.\n",
    "Try lowering EZHU_CCC_MIN_CELLS or inspect the group counts table written to:\n",
    file.path(out_dir, "cellchat_group_counts.tsv"),
    call. = FALSE
  )
}

message("CellChat input cells: ", ncol(obj), " | groups: ", length(unique(obj[[group_col, drop = TRUE]])))

db <- NULL
if (exists("CellChatDB.human", envir = asNamespace("CellChat"), inherits = FALSE)) {
  db <- get("CellChatDB.human", envir = asNamespace("CellChat"), inherits = FALSE)
} else {
  tryCatch(
    {
      data("CellChatDB.human", package = "CellChat", envir = environment())
      db <- get("CellChatDB.human", envir = environment(), inherits = FALSE)
    },
    error = function(e) NULL
  )
}
if (is.null(db)) {
  stop(
    "Could not load CellChatDB.human from the CellChat package.\n",
    "Check that the installed CellChat version includes the built-in human database.",
    call. = FALSE
  )
}

meta <- obj[[]]
if (!identical(rownames(meta), colnames(obj))) {
  stop("Seurat meta.data rownames do not match cell names after filtering.", call. = FALSE)
}

# IMPORTANT: pass a matrix input (not a Seurat object) to avoid Seurat v5 incompatibilities
# inside CellChat/SeuratObject APIs (e.g., GetAssayData(slot=...) being defunct).
expr <- get_seurat_expr_matrix(obj, prefer = c("data", "counts"))
if (!identical(colnames(expr), rownames(meta))) {
  expr <- expr[, rownames(meta), drop = FALSE]
}

cellchat <- CellChat::createCellChat(object = expr, meta = meta, group.by = group_col)
cellchat@DB <- db

cellchat <- CellChat::subsetData(cellchat)
cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)
cellchat <- CellChat::computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- CellChat::filterCommunication(cellchat, min.cells = min_cells)
cellchat <- CellChat::computeCommunProbPathway(cellchat)
cellchat <- CellChat::aggregateNet(cellchat)

ezhu_dir_create("data/processed/ccc")
saveRDS(cellchat, file = file.path("data/processed/ccc", paste0("cellchat_", run_id, ".rds")))

comm <- CellChat::subsetCommunication(cellchat)
comm$source <- as.character(comm$source)
comm$target <- as.character(comm$target)
comm$ligand <- as.character(comm$ligand)
comm$receptor <- as.character(comm$receptor)
comm$pathway_name <- as.character(comm$pathway_name)

comm$prob <- suppressWarnings(as.numeric(comm$prob))
comm$pval <- suppressWarnings(as.numeric(comm$pval))
comm <- comm[is.finite(comm$prob) & is.finite(comm$pval), , drop = FALSE]
comm <- comm[order(comm$pval, -comm$prob), , drop = FALSE]

utils::write.table(
  comm,
  file = file.path(out_dir, "cellchat_interactions.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

pathway_df <- comm
if (nrow(pathway_df) > 0 && "pathway_name" %in% names(pathway_df)) {
  pathway_df <- pathway_df[nzchar(pathway_df$pathway_name), , drop = FALSE]
  pathway_summary <- stats::aggregate(
    prob ~ pathway_name + source + target,
    data = pathway_df,
    FUN = sum
  )
  names(pathway_summary)[names(pathway_summary) == "prob"] <- "prob_sum"
  pathway_summary <- pathway_summary[order(pathway_summary$prob_sum, decreasing = TRUE), , drop = FALSE]
  utils::write.table(
    pathway_summary,
    file = file.path(out_dir, "cellchat_pathway_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

session_path <- file.path(out_dir, "sessionInfo.txt")
writeLines(capture.output(sessionInfo()), con = session_path, useBytes = TRUE)

message("Wrote: ", file.path(out_dir, "cellchat_interactions.tsv"))
message("Wrote: ", file.path(out_dir, "cellchat_pathway_summary.tsv"))
message("Wrote: ", file.path(out_dir, "cellchat_group_counts.tsv"))
