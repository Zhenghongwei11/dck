source("scripts/utils/repro.R")
ezhu_set_repo_root()
ezhu_activate_renv()

seed_env <- Sys.getenv("EZHU_SEED", unset = "")
seed <- suppressWarnings(as.integer(seed_env))
if (is.na(seed)) seed <- 20251227L
if (!nzchar(seed_env)) Sys.setenv(EZHU_SEED = as.character(seed))
set.seed(seed)

safe_read_seurat <- function(path) {
  if (!file.exists(path)) stop("Missing Seurat RDS: ", path, call. = FALSE)
  readRDS(path)
}

get_assay_counts <- function(obj, assay = "RNA") {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("Missing package: SeuratObject", call. = FALSE)
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Missing package: Matrix", call. = FALSE)

  assay_names <- character()
  try(assay_names <- names(obj@assays), silent = TRUE)
  if (length(assay_names) == 0L || !assay %in% assay_names) {
    stop("Assay not found in Seurat object: ", assay, call. = FALSE)
  }

  ga <- NULL
  if (exists("GetAssayData", where = asNamespace("SeuratObject"), inherits = FALSE)) {
    ga <- get("GetAssayData", envir = asNamespace("SeuratObject"))
  }
  if (is.null(ga)) stop("Cannot locate SeuratObject::GetAssayData()", call. = FALSE)

  out <- try(ga(obj[[assay]], slot = "counts"), silent = TRUE) # Seurat v4
  if (!inherits(out, "try-error")) return(out)
  out <- try(ga(obj[[assay]], layer = "counts"), silent = TRUE) # Seurat v5+
  if (!inherits(out, "try-error")) return(out)

  stop("Cannot extract assay counts from Seurat object (assay=", assay, ").", call. = FALSE)
}

compute_percent_mt_seurat <- function(obj, assay = "RNA") {
  # Seurat helper is convenient when rownames already use HGNC symbols (MT-*)
  if (!requireNamespace("Seurat", quietly = TRUE)) return(NULL)
  if (!exists("PercentageFeatureSet", where = asNamespace("Seurat"), inherits = FALSE)) return(NULL)

  fn <- get("PercentageFeatureSet", envir = asNamespace("Seurat"))
  out <- try(fn(obj, pattern = "^MT-", assay = assay), silent = TRUE)
  if (inherits(out, "try-error")) return(NULL)
  out <- suppressWarnings(as.numeric(out))
  if (length(out) == 0L || all(!is.finite(out))) return(NULL)
  out
}

compute_percent_mt <- function(counts, feature_ids, mapping_path = "data/references/singlecell/gse131882_gene_universe.tsv") {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Missing package: Matrix", call. = FALSE)

  feat_upper <- toupper(feature_ids)
  mito <- grepl("^MT-", feat_upper) | grepl("^MT\\.", feat_upper)

  if (!any(mito) && file.exists(mapping_path)) {
    map <- utils::read.delim(mapping_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    names(map) <- tolower(names(map))
    if (all(c("ensg", "gene_symbol") %in% names(map))) {
      sym_map <- setNames(as.character(map$gene_symbol), as.character(map$ensg))
      ensg_norm <- sub("\\..*$", "", feature_ids)
      sym <- sym_map[ensg_norm]
      mito <- !is.na(sym) & grepl("^MT-", toupper(sym))
    }
  }

  if (!any(mito)) return(rep(NA_real_, ncol(counts)))

  mt_counts <- Matrix::colSums(counts[mito, , drop = FALSE])
  total_counts <- Matrix::colSums(counts)
  pct <- ifelse(total_counts > 0, 100 * (mt_counts / total_counts), NA_real_)
  as.numeric(pct)
}

export_qc_anchors <- function(
  seurat_rds,
  out_cells = "results/figures/gse131882_qc_cells.tsv",
  out_sample = "results/figures/gse131882_qc_sample_summary.tsv"
) {
  obj <- safe_read_seurat(seurat_rds)
  counts <- tryCatch(get_assay_counts(obj, assay = "RNA"), error = function(e) NULL)

  md <- obj@meta.data
  md$cell <- rownames(md)

  derived_sample <- sub("_[^_]+$", "", md$cell)
  use_derived <- TRUE
  if ("sample" %in% names(md) && any(nzchar(as.character(md$sample)))) {
    current <- as.character(md$sample)
    # Use derived sample names if they carry explicit group/replicate tags.
    current_has_group <- mean(grepl("control|diabetes|dn", current, ignore.case = TRUE)) > 0.8
    derived_has_group <- mean(grepl("control|diabetes|dn", derived_sample, ignore.case = TRUE)) > 0.8
    if (current_has_group && !derived_has_group) use_derived <- FALSE
  }
  md$sample <- if (use_derived) derived_sample else as.character(md$sample)
  if (mean(grepl("\\.rds$", md$sample)) < 0.9) md$sample <- paste0(md$sample, ".rds")
  if (!"group" %in% names(md) || all(!nzchar(as.character(md$group)))) {
    md$group <- ifelse(
      grepl("control", md$sample, fixed = TRUE),
      "control",
      ifelse(grepl("diabetes", md$sample, fixed = TRUE) | grepl("dn", md$sample, ignore.case = TRUE), "diabetes", NA_character_)
    )
  }
  md$sample[is.na(md$sample) | !nzchar(md$sample)] <- "unknown.rds"
  md$group[is.na(md$group) | !nzchar(md$group)] <- "unknown"

  # Prefer Seurat-provided QC fields; only fall back to computing from counts if required.
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Missing package: Matrix", call. = FALSE)
  if (!"nCount_RNA" %in% names(md) || !all(is.finite(suppressWarnings(as.numeric(md$nCount_RNA))))) {
    if (is.null(counts)) {
      stop(
        "Cannot extract assay counts from Seurat object (assay=RNA), and meta.data lacks a usable nCount_RNA column.\n",
        "Fix options:\n",
        "  - Re-run preprocessing without stripping RNA counts.\n",
        "  - Or regenerate QC anchors from a Seurat object that retains RNA counts.\n",
        call. = FALSE
      )
    }
    md$nCount_RNA <- Matrix::colSums(counts)[md$cell]
  }
  if (!"nFeature_RNA" %in% names(md) || !all(is.finite(suppressWarnings(as.numeric(md$nFeature_RNA))))) {
    if (is.null(counts)) {
      stop(
        "Cannot extract assay counts from Seurat object (assay=RNA), and meta.data lacks a usable nFeature_RNA column.\n",
        "Fix options:\n",
        "  - Re-run preprocessing without stripping RNA counts.\n",
        "  - Or regenerate QC anchors from a Seurat object that retains RNA counts.\n",
        call. = FALSE
      )
    }
    md$nFeature_RNA <- Matrix::colSums(counts > 0)[md$cell]
  }

  # percent.mt: prefer meta.data (written during preprocess). Otherwise compute from the object.
  pct_mt <- NULL
  if ("percent.mt" %in% names(md)) {
    pct_mt <- suppressWarnings(as.numeric(md[["percent.mt"]]))
    if (length(pct_mt) != nrow(md) || all(!is.finite(pct_mt))) pct_mt <- NULL
  }
  if (is.null(pct_mt)) {
    pct_mt <- compute_percent_mt_seurat(obj, assay = "RNA")
  }
  if (is.null(pct_mt) || all(!is.finite(pct_mt)) || all(pct_mt == 0, na.rm = TRUE)) {
    if (!is.null(counts)) {
      pct_mt <- compute_percent_mt(counts, rownames(counts))
      names(pct_mt) <- colnames(counts)
      pct_mt <- pct_mt[md$cell]
    } else {
      pct_mt <- rep(NA_real_, nrow(md))
    }
  }
  md$percent.mt <- pct_mt

  if ("seurat_clusters" %in% names(md)) md$seurat_clusters <- as.character(md$seurat_clusters)

  out <- md[, c("cell", "nCount_RNA", "nFeature_RNA", "percent.mt", "sample", "group", "seurat_clusters")]
  out$nCount_RNA <- suppressWarnings(as.numeric(out$nCount_RNA))
  out$nFeature_RNA <- suppressWarnings(as.numeric(out$nFeature_RNA))
  out$percent.mt <- suppressWarnings(as.numeric(out$percent.mt))
  out$sample <- as.character(out$sample)
  out$group <- as.character(out$group)
  out$seurat_clusters <- as.character(out$seurat_clusters)

  ezhu_dir_create(dirname(out_cells))
  utils::write.table(out, file = out_cells, sep = "\t", quote = FALSE, row.names = FALSE)

  sm <- aggregate(
    cbind(nCount_RNA, nFeature_RNA, percent.mt) ~ sample + group,
    data = out,
    FUN = function(x) mean(x, na.rm = TRUE)
  )
  cell_count <- as.integer(table(out$sample)[sm$sample])
  sample_out <- data.frame(
    sample = as.character(sm$sample),
    group = as.character(sm$group),
    cell_count = cell_count,
    mean_nCount = as.numeric(sm$nCount_RNA),
    mean_nFeature = as.numeric(sm$nFeature_RNA),
    mean_percent_mt = as.numeric(sm$percent.mt),
    stringsAsFactors = FALSE
  )
  utils::write.table(sample_out, file = out_sample, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(list(cells = out, samples = sample_out))
}

seurat_rds <- Sys.getenv("EZHU_SEURAT_RDS", unset = "data/processed/singlecell_gse131882_seurat.rds")
ezhu_write_stage_metadata(
  "02d_export_qc_anchors",
  params = list(
    seurat_rds = seurat_rds,
    out_cells = "results/figures/gse131882_qc_cells.tsv",
    out_sample = "results/figures/gse131882_qc_sample_summary.tsv"
  ),
  seed = seed
)
export_qc_anchors(seurat_rds)
