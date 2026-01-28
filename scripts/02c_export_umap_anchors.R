source("scripts/utils/repro.R")
ezhu_set_repo_root()
ezhu_activate_renv()

seed_env <- Sys.getenv("EZHU_SEED", unset = "")
seed <- suppressWarnings(as.integer(seed_env))
if (is.na(seed)) seed <- 20251227L
if (!nzchar(seed_env)) Sys.setenv(EZHU_SEED = as.character(seed))
set.seed(seed)

options(stringsAsFactors = FALSE, width = 120)

if (!requireNamespace("Seurat", quietly = TRUE)) stop("Missing package: Seurat", call. = FALSE)
if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("Missing package: SeuratObject", call. = FALSE)
if (!requireNamespace("Matrix", quietly = TRUE)) stop("Missing package: Matrix", call. = FALSE)
suppressPackageStartupMessages(library(Seurat))

safe_read_seurat <- function(path) {
  if (!file.exists(path)) stop("Missing Seurat RDS: ", path, call. = FALSE)
  readRDS(path)
}

get_assay_matrix <- function(obj, assay = "RNA", layer = "data") {
  ga <- NULL
  if (exists("GetAssayData", where = asNamespace("SeuratObject"), inherits = FALSE)) {
    ga <- get("GetAssayData", envir = asNamespace("SeuratObject"))
  }
  if (is.null(ga)) stop("Cannot locate SeuratObject::GetAssayData()", call. = FALSE)

  out <- try(ga(obj[[assay]], slot = layer), silent = TRUE) # Seurat v4
  if (!inherits(out, "try-error")) return(out)
  out <- try(ga(obj[[assay]], layer = layer), silent = TRUE) # Seurat v5+
  if (!inherits(out, "try-error")) return(out)
  NULL
}

list_assay_layers <- function(assay_obj) {
  layers <- character()
  if (requireNamespace("SeuratObject", quietly = TRUE) &&
      exists("Layers", where = asNamespace("SeuratObject"), inherits = FALSE)) {
    layers <- tryCatch(SeuratObject::Layers(assay_obj), error = function(e) character())
  }
  if (length(layers) == 0L) {
    layers <- tryCatch(names(assay_obj@layers), error = function(e) character())
  }
  unique(as.character(layers))
}

get_assay_layer_matrix <- function(obj, assay, layer_name) {
  ga <- NULL
  if (exists("GetAssayData", where = asNamespace("SeuratObject"), inherits = FALSE)) {
    ga <- get("GetAssayData", envir = asNamespace("SeuratObject"))
  }
  if (is.null(ga)) return(NULL)
  out <- try(ga(obj[[assay]], layer = layer_name), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  out
}

safe_fetch_expression <- function(obj, assay, genes, cells) {
  Seurat::DefaultAssay(obj) <- assay
  genes <- genes[genes %in% rownames(obj[[assay]])]
  if (length(genes) == 0) return(data.frame(cell = cells, stringsAsFactors = FALSE))

  fetched <- tryCatch(
    {
      df <- Seurat::FetchData(obj, vars = genes, cells = cells, slot = "data")
      df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
      df$cell <- rownames(df)
      df
    },
    error = function(e) NULL
  )
  if (!is.null(fetched)) return(fetched)

  layers <- list_assay_layers(obj[[assay]])
  data_layers <- layers[layers == "data" | grepl("^data\\.", layers)]
  if (length(data_layers) > 0) {
    out <- matrix(NA_real_, nrow = length(cells), ncol = length(genes), dimnames = list(cells, genes))
    for (layer_name in data_layers) {
      mat <- get_assay_layer_matrix(obj, assay = assay, layer_name = layer_name)
      if (is.null(mat)) next
      layer_cells <- intersect(colnames(mat), cells)
      if (length(layer_cells) == 0) next
      layer_genes <- intersect(genes, rownames(mat))
      if (length(layer_genes) == 0) next
      sub <- mat[layer_genes, layer_cells, drop = FALSE]
      sub <- as.matrix(sub)
      out[layer_cells, layer_genes] <- t(sub)
    }
    df <- data.frame(cell = cells, out, check.names = FALSE, stringsAsFactors = FALSE)
    return(df)
  }

  counts_layers <- layers[layers == "counts" | grepl("^counts\\.", layers)]
  if (length(counts_layers) > 0) {
    out <- matrix(NA_real_, nrow = length(cells), ncol = length(genes), dimnames = list(cells, genes))
    for (layer_name in counts_layers) {
      mat <- get_assay_layer_matrix(obj, assay = assay, layer_name = layer_name)
      if (is.null(mat)) next
      layer_cells <- intersect(colnames(mat), cells)
      if (length(layer_cells) == 0) next
      layer_genes <- intersect(genes, rownames(mat))
      if (length(layer_genes) == 0) next
      lib <- Matrix::colSums(mat[, layer_cells, drop = FALSE])
      lib[lib <= 0] <- NA_real_
      sub <- mat[layer_genes, layer_cells, drop = FALSE]
      sub <- as.matrix(sub)
      norm <- t(t(sub) / as.numeric(lib)) * 10000
      norm <- log1p(norm)
      out[layer_cells, layer_genes] <- t(norm)
    }
    df <- data.frame(cell = cells, out, check.names = FALSE, stringsAsFactors = FALSE)
    return(df)
  }

  mat <- get_assay_matrix(obj, assay = assay, layer = "data")
  if (!is.null(mat)) {
    sub <- mat[genes, cells, drop = FALSE]
    sub <- Matrix::Matrix(sub, sparse = FALSE)
    df <- data.frame(cell = colnames(sub), t(sub), check.names = FALSE, stringsAsFactors = FALSE)
    return(df)
  }

  counts <- get_assay_matrix(obj, assay = assay, layer = "counts")
  if (is.null(counts)) {
    warning("Cannot export feature expression: neither data nor counts layer available.", call. = FALSE)
    return(data.frame(cell = cells, stringsAsFactors = FALSE))
  }
  sub <- counts[genes, cells, drop = FALSE]
  lib <- Matrix::colSums(counts[, cells, drop = FALSE])
  lib[lib <= 0] <- NA_real_
  # LogNormalize(scale.factor=10000) on the fly for selected genes.
  norm <- t(t(sub) / lib) * 10000
  norm <- log1p(norm)
  norm <- Matrix::Matrix(norm, sparse = FALSE)
  df <- data.frame(cell = colnames(norm), t(norm), check.names = FALSE, stringsAsFactors = FALSE)
  df
}

export_umap_anchors <- function(
  seurat_rds,
  assay = "RNA",
  reduction = "umap",
  features = c("TNFRSF12A", "SLC12A1", "AQP2", "UMOD", "SLC4A1", "PECAM1"),
  out_cells = "results/figures/gse131882_umap_cells.tsv",
  out_features = "results/figures/gse131882_umap_features.tsv"
) {
  obj <- safe_read_seurat(seurat_rds)
  if (!inherits(obj, "Seurat")) stop("Expected a Seurat object: ", seurat_rds, call. = FALSE)
  if (!reduction %in% names(obj@reductions)) stop("Missing reduction '", reduction, "' in Seurat object: ", seurat_rds, call. = FALSE)
  if (!assay %in% Seurat::Assays(obj)) stop("Assay not present in Seurat object: ", assay, call. = FALSE)

  emb <- Seurat::Embeddings(obj, reduction = reduction)
  if (ncol(emb) < 2) stop("UMAP embedding has <2 dimensions.", call. = FALSE)

  md <- obj@meta.data
  md$cell <- rownames(md)
  emb_df <- data.frame(
    cell = rownames(emb),
    UMAP_1 = as.numeric(emb[, 1]),
    UMAP_2 = as.numeric(emb[, 2]),
    stringsAsFactors = FALSE
  )

  out <- merge(emb_df, md, by = "cell", all.x = TRUE, sort = FALSE)
  if (!"seurat_clusters" %in% names(out)) {
    out$seurat_clusters <- as.character(Seurat::Idents(obj))[match(out$cell, names(Seurat::Idents(obj)))]
  }
  out$cluster <- as.character(out$seurat_clusters)
  out$sample <- as.character(out$sample %||% NA_character_)
  out$group <- as.character(out$group %||% NA_character_)

  cells_out <- out[, c("cell", "UMAP_1", "UMAP_2", "cluster", "sample", "group"), drop = FALSE]
  cells_out <- cells_out[nzchar(cells_out$cell) & is.finite(cells_out$UMAP_1) & is.finite(cells_out$UMAP_2), , drop = FALSE]

  present <- features[features %in% rownames(obj[[assay]])]
  missing <- setdiff(features, present)
  if (length(missing) > 0) {
    warning("UMAP feature export: missing genes in assay rownames: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  feat_df <- data.frame(cell = cells_out$cell, stringsAsFactors = FALSE)
  if (length(present) > 0) {
    fetched <- safe_fetch_expression(obj, assay = assay, genes = present, cells = cells_out$cell)
    feat_df <- merge(feat_df, fetched, by = "cell", all.x = TRUE, sort = FALSE)
  }

  ezhu_dir_create(dirname(out_cells))
  utils::write.table(cells_out, file = out_cells, sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(feat_df, file = out_features, sep = "\t", quote = FALSE, row.names = FALSE)

  invisible(list(cells = cells_out, features = feat_df))
}

seurat_rds <- Sys.getenv("EZHU_SEURAT_RDS", unset = "data/processed/singlecell_gse131882_seurat.rds")
assay <- Sys.getenv("EZHU_UMAP_ASSAY", unset = "RNA")
reduction <- Sys.getenv("EZHU_UMAP_REDUCTION", unset = "umap")
features_env <- trimws(Sys.getenv("EZHU_UMAP_FEATURES", unset = "TNFRSF12A,SLC12A1,AQP2,UMOD,SLC4A1,PECAM1"))
features <- trimws(unlist(strsplit(features_env, "[,;\\s]+")))
features <- features[nzchar(features)]

ezhu_write_stage_metadata(
  "02c_export_umap_anchors",
  params = list(
    seurat_rds = seurat_rds,
    assay = assay,
    reduction = reduction,
    features = features,
    out_cells = "results/figures/gse131882_umap_cells.tsv",
    out_features = "results/figures/gse131882_umap_features.tsv"
  ),
  seed = seed
)

export_umap_anchors(seurat_rds, assay = assay, reduction = reduction, features = features)
