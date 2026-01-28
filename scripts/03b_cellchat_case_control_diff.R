#!/usr/bin/env Rscript
# scripts/03b_cellchat_case_control_diff.R
#
# Purpose:
# - Differential CellChat analysis between two groups (Control vs DKD/diabetes).
#
# Outputs:
# - results/ccc_diff/<run_id>/
#   - interactions_control.tsv, interactions_case.tsv
#   - interactions_long.tsv (all interactions with a group column)
#   - diff_interactions.tsv (joined control vs case with Δprob/log2FC)
#   - cellchat_merged.rds
# - plots/publication/
#   - figureS_cellchat_diff_interactions_summary_<run_id>.pdf
#   - figureS_cellchat_diff_ranknet_<run_id>.pdf

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
if (is.na(seed)) seed <- 20260117L
set.seed(seed)

options(stringsAsFactors = FALSE, width = 120)

require_or_stop <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Missing R package: ", pkg, call. = FALSE)
  }
}

require_or_stop("Seurat")
require_or_stop("SeuratObject")
require_or_stop("CellChat")
require_or_stop("ggplot2")

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(CellChat)
  library(ggplot2)
})

parse_csv <- function(x) {
  x <- trimws(as.character(x %||% ""))
  if (!nzchar(x)) return(character())
  parts <- unlist(strsplit(x, "[,;]"))
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

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

  stop("Could not extract an expression matrix from Seurat.", call. = FALSE)
}

load_cellchat_db_human <- function() {
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
  if (is.null(db)) stop("Could not load CellChatDB.human.", call. = FALSE)
  db
}

seurat_path <- Sys.getenv("EZHU_SEURAT_RDS", unset = "data/processed/singlecell_gse131882_seurat.rds")
if (!file.exists(seurat_path)) stop("Missing Seurat RDS: ", seurat_path, call. = FALSE)

case_control_col <- Sys.getenv("EZHU_CASE_CONTROL_COL", unset = "group")
group_col <- Sys.getenv("EZHU_CCC_GROUP_COL", unset = "cell_context")
cluster_col <- Sys.getenv("EZHU_CLUSTER_COL", unset = "seurat_clusters")
min_cells <- suppressWarnings(as.integer(Sys.getenv("EZHU_CCC_MIN_CELLS", unset = "50")))
if (!is.finite(min_cells) || is.na(min_cells) || min_cells < 10) min_cells <- 50L

run_id <- Sys.getenv("EZHU_CCC_DIFF_RUN_ID", unset = format(Sys.time(), "%Y%m%d_%H%M%S"))
out_dir <- file.path("results/ccc_diff", run_id)
ezhu_dir_create(out_dir)
ezhu_dir_create("plots/publication")

message("Loading Seurat: ", seurat_path)
obj <- readRDS(seurat_path)
if (!inherits(obj, "Seurat")) stop("Input is not a Seurat object: ", seurat_path, call. = FALSE)

meta0 <- obj[[]]
if (!identical(rownames(meta0), colnames(obj))) {
  stop("Seurat meta.data rownames do not match cell names in the object.", call. = FALSE)
}
if (!case_control_col %in% colnames(meta0)) {
  stop("Seurat meta.data missing case/control column: ", case_control_col, call. = FALSE)
}
if (!cluster_col %in% colnames(meta0)) {
  stop("Seurat meta.data missing cluster column: ", cluster_col, call. = FALSE)
}

# Resolve per-cell group labels for CellChat (cell contexts).
groups_existing <- NULL
if (group_col %in% colnames(meta0)) groups_existing <- as.character(meta0[[group_col]])
groups_existing <- as.character(groups_existing)
groups_existing[is.na(groups_existing)] <- ""
groups_existing <- trimws(groups_existing)
existing_has_any <- length(groups_existing) > 0 && any(nzchar(groups_existing))

cluster_ids <- trimws(as.character(meta0[[cluster_col]]))
cluster_ids[is.na(cluster_ids)] <- "NA"

mapping_path <- Sys.getenv("EZHU_CLUSTER_ANNOTATIONS", unset = "results/annotation/cluster_annotations_filled.tsv")
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

use_cluster_annotations <- !existing_has_any
groups <- NULL
if (use_cluster_annotations) {
  if (is.null(lookup)) {
    message("No usable '", group_col, "' column and no cluster annotation file; falling back to cluster IDs.")
    groups <- rep("", length(cluster_ids))
  } else {
    groups <- as.character(lookup[cluster_ids])
  }
} else {
  groups <- groups_existing
}

groups <- as.character(groups)
groups[is.na(groups)] <- ""
groups <- trimws(groups)
missing <- !nzchar(groups)
groups[missing] <- paste0("C", cluster_ids[missing])
obj[[group_col]] <- groups

# Identify the two case/control labels to compare.
cc <- meta0[[case_control_col]]
cc <- as.character(cc)
cc[is.na(cc)] <- ""
cc <- trimws(cc)
cc_lc <- tolower(cc)
unique_groups <- sort(unique(cc_lc[nzchar(cc_lc)]))
if (length(unique_groups) < 2) {
  stop("Need >=2 groups in ", case_control_col, " to run differential CellChat.", call. = FALSE)
}

control_labels <- tolower(parse_csv(Sys.getenv("EZHU_CONTROL_LABELS", unset = "control")))
case_labels <- tolower(parse_csv(Sys.getenv("EZHU_CASE_LABELS", unset = "diabetes,dn,case,disease")))

pick_first <- function(candidates) {
  hit <- intersect(candidates, unique_groups)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

control_key <- pick_first(control_labels)
case_key <- pick_first(case_labels)
if (!nzchar(control_key) || !nzchar(case_key) || identical(control_key, case_key)) {
  # Fallback: if exactly 2 unique groups, take them as (control, case) in alphabetical order.
  if (length(unique_groups) == 2) {
    control_key <- unique_groups[1]
    case_key <- unique_groups[2]
    message("Auto-selected groups: control=", control_key, " case=", case_key)
  } else {
    stop(
      "Could not unambiguously pick control/case groups from '", case_control_col, "'.\n",
      "Observed groups: ", paste(unique_groups, collapse = ", "), "\n",
      "Set EZHU_CONTROL_LABELS and EZHU_CASE_LABELS explicitly.",
      call. = FALSE
    )
  }
}

keep_mask <- cc_lc %in% c(control_key, case_key)
keep_cells <- rownames(meta0)[keep_mask]
obj <- subset(obj, cells = keep_cells)

meta <- obj[[]]
meta_cc <- tolower(trimws(as.character(meta[[case_control_col]])))
meta_cc[is.na(meta_cc)] <- ""

# Filter to groups with sufficient cells in BOTH conditions.
counts_control <- table(meta[[group_col]][meta_cc == control_key])
counts_case <- table(meta[[group_col]][meta_cc == case_key])
eligible <- intersect(names(counts_control[counts_control >= min_cells]), names(counts_case[counts_case >= min_cells]))
eligible <- sort(unique(as.character(eligible)))
if (length(eligible) < 2) {
  stop("Too few eligible cell groups after min_cells filtering (min_cells=", min_cells, ").", call. = FALSE)
}

keep2 <- meta[[group_col]] %in% eligible
obj <- subset(obj, cells = rownames(meta)[keep2])
meta <- obj[[]]
meta_cc <- tolower(trimws(as.character(meta[[case_control_col]])))

write_counts <- function(x, path) {
  df <- data.frame(group = names(x), n_cells = as.integer(x), stringsAsFactors = FALSE)
  utils::write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE)
}
write_counts(sort(table(meta[[group_col]][meta_cc == control_key]), decreasing = TRUE), file.path(out_dir, "group_counts_control.tsv"))
write_counts(sort(table(meta[[group_col]][meta_cc == case_key]), decreasing = TRUE), file.path(out_dir, "group_counts_case.tsv"))

message("Running CellChat on control=", control_key, " case=", case_key, " (groups=", length(eligible), ")")
db <- load_cellchat_db_human()

run_cellchat <- function(sub_obj) {
  sub_meta <- sub_obj[[]]
  expr <- get_seurat_expr_matrix(sub_obj, prefer = c("data", "counts"))
  if (!identical(colnames(expr), rownames(sub_meta))) expr <- expr[, rownames(sub_meta), drop = FALSE]
  cc <- CellChat::createCellChat(object = expr, meta = sub_meta, group.by = group_col)
  cc@DB <- db
  cc <- CellChat::subsetData(cc)
  cc <- CellChat::identifyOverExpressedGenes(cc)
  cc <- CellChat::identifyOverExpressedInteractions(cc)
  cc <- CellChat::computeCommunProb(cc, raw.use = TRUE)
  cc <- CellChat::filterCommunication(cc, min.cells = min_cells)
  cc <- CellChat::computeCommunProbPathway(cc)
  cc <- CellChat::aggregateNet(cc)
  cc
}

obj_control <- subset(obj, cells = rownames(meta)[meta_cc == control_key])
obj_case <- subset(obj, cells = rownames(meta)[meta_cc == case_key])

cc_control <- run_cellchat(obj_control)
cc_case <- run_cellchat(obj_case)

cellchat <- CellChat::mergeCellChat(list(control = cc_control, case = cc_case), add.names = c("control", "case"))
saveRDS(cellchat, file = file.path(out_dir, "cellchat_merged.rds"))

comm_control <- CellChat::subsetCommunication(cc_control)
comm_case <- CellChat::subsetCommunication(cc_case)

standardize_comm <- function(df) {
  df$source <- as.character(df$source)
  df$target <- as.character(df$target)
  df$ligand <- as.character(df$ligand)
  df$receptor <- as.character(df$receptor)
  if ("interaction_name" %in% names(df)) df$interaction_name <- as.character(df$interaction_name)
  if ("pathway_name" %in% names(df)) df$pathway_name <- as.character(df$pathway_name)
  df$prob <- suppressWarnings(as.numeric(df$prob))
  df$pval <- suppressWarnings(as.numeric(df$pval))
  df <- df[is.finite(df$prob) & is.finite(df$pval), , drop = FALSE]
  df
}
comm_control <- standardize_comm(comm_control)
comm_case <- standardize_comm(comm_case)

utils::write.table(comm_control, file = file.path(out_dir, "interactions_control.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(comm_case, file = file.path(out_dir, "interactions_case.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(comm_control, file = file.path(out_dir, paste0("interactions_", control_key, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
utils::write.table(comm_case, file = file.path(out_dir, paste0("interactions_", case_key, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

comm_control$group <- "control"
comm_case$group <- "case"
comm_long <- rbind(comm_control, comm_case)
utils::write.table(comm_long, file = file.path(out_dir, "interactions_long.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Build a joined diff table at the interaction level.
join_key <- intersect(
  c("source", "target", "ligand", "receptor", "interaction_name", "pathway_name"),
  intersect(names(comm_control), names(comm_case))
)
if (length(join_key) < 4) {
  join_key <- c("source", "target", "ligand", "receptor")
}

cc1 <- comm_control[, unique(c(join_key, "prob", "pval")), drop = FALSE]
cc2 <- comm_case[, unique(c(join_key, "prob", "pval")), drop = FALSE]
names(cc1)[names(cc1) == "prob"] <- "prob_control"
names(cc1)[names(cc1) == "pval"] <- "pval_control"
names(cc2)[names(cc2) == "prob"] <- "prob_case"
names(cc2)[names(cc2) == "pval"] <- "pval_case"

diff_tbl <- merge(cc1, cc2, by = join_key, all = TRUE, sort = FALSE)
diff_tbl$prob_control <- suppressWarnings(as.numeric(diff_tbl$prob_control))
diff_tbl$prob_case <- suppressWarnings(as.numeric(diff_tbl$prob_case))
diff_tbl$prob_diff <- diff_tbl$prob_case - diff_tbl$prob_control
diff_tbl$prob_log2fc <- log2((diff_tbl$prob_case + 1e-12) / (diff_tbl$prob_control + 1e-12))

utils::write.table(diff_tbl, file = file.path(out_dir, "diff_interactions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Figures (comparison summary + rankNet)
gg1 <- CellChat::compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2))
gg2 <- CellChat::compareInteractions(cellchat, show.legend = FALSE, group = c(1, 2), measure = "weight")

pdf(file.path("plots/publication", paste0("figureS_cellchat_diff_interactions_summary_", run_id, ".pdf")), width = 6, height = 4)
print(gg1 + gg2)
dev.off()

pdf(file.path("plots/publication", paste0("figureS_cellchat_diff_ranknet_", run_id, ".pdf")), width = 8, height = 7)
CellChat::rankNet(cellchat, mode = "comparison", stacked = TRUE, do.stat = TRUE)
dev.off()

ezhu_write_stage_metadata(
  "03b_cellchat_case_control_diff",
  params = list(
    seurat_rds = seurat_path,
    case_control_col = case_control_col,
    group_col = group_col,
    cluster_col = cluster_col,
    control_key = control_key,
    case_key = case_key,
    min_cells = min_cells,
    run_id = run_id,
    outputs = list(
      diff_interactions = file.path(out_dir, "diff_interactions.tsv"),
      interactions_control = file.path(out_dir, "interactions_control.tsv"),
      interactions_case = file.path(out_dir, "interactions_case.tsv"),
      figures = list(
        paste0("plots/publication/figureS_cellchat_diff_interactions_summary_", run_id, ".pdf"),
        paste0("plots/publication/figureS_cellchat_diff_ranknet_", run_id, ".pdf")
      )
    )
  ),
  seed = seed
)
