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

# Matrix >= 1.7 removed the S4 virtual class "mMatrix". Seurat v4 / SeuratObject v4
# require this class to be defined at runtime to avoid validation errors
# during Graph object construction. Defining it early (before library calls)
# provides a compatibility shim.
if (!methods::isClass("mMatrix")) {
  tryCatch(
    methods::setClass("mMatrix", contains = "Matrix"),
    error = function(e) NULL
  )
}

seed_env <- Sys.getenv("EZHU_SEED", unset = "")
seed <- suppressWarnings(as.integer(seed_env))
if (is.na(seed)) seed <- 20251227L
set.seed(seed)

ezhu_write_stage_metadata("02_preprocess", params = list(), seed = seed)

options(stringsAsFactors = FALSE, width = 120)

if (!requireNamespace("Matrix", quietly = TRUE)) stop("Missing package: Matrix", call. = FALSE)
if (!requireNamespace("Seurat", quietly = TRUE)) stop("Missing package: Seurat", call. = FALSE)
if (!requireNamespace("SeuratObject", quietly = TRUE)) stop("Missing package: SeuratObject", call. = FALSE)

suppressPackageStartupMessages(library(Matrix))

suppressPackageStartupMessages(library(SeuratObject))
suppressPackageStartupMessages(library(Seurat))

matrix_version <- tryCatch(packageVersion("Matrix"), error = function(e) NULL)
seurat_version <- tryCatch(packageVersion("Seurat"), error = function(e) NULL)
seuratobject_version <- tryCatch(packageVersion("SeuratObject"), error = function(e) NULL)

if (!is.null(seuratobject_version) && seuratobject_version >= "5.0.0") {
  stop(
    "Incompatible SeuratObject version detected: ", as.character(seuratobject_version), "\n",
    "This pipeline requires SeuratObject < 5.0.0 (Seurat v4-era APIs). If you installed SeuratObject from CRAN, it likely drifted to v5.\n",
    "Fix (recommended): recreate a clean micromamba env from env/ezhu-r-seurat4.yml and run all steps via scripts/run_r_mamba.sh.\n",
    call. = FALSE
  )
}
if (!is.null(seurat_version) && seurat_version >= "5.0.0") {
  stop(
    "Incompatible Seurat version detected: ", as.character(seurat_version), "\n",
    "This pipeline requires Seurat < 5.0.0.\n",
    "Fix (recommended): recreate a clean micromamba env from env/ezhu-r-seurat4.yml and run all steps via scripts/run_r_mamba.sh.\n",
    call. = FALSE
  )
}
if (!is.null(matrix_version) && !is.null(seurat_version)) {
  # Matrix >= 1.7 removed the S4 virtual class "mMatrix". Seurat v4 / SeuratObject v4
  # can error when constructing Graph objects if Matrix is too new.
  if (seurat_version < "5.0.0" && matrix_version >= "1.7.0") {
    stop(
      "Incompatible Matrix version for Seurat v4.\n",
      "Detected Seurat=", as.character(seurat_version), " with Matrix=", as.character(matrix_version), ".\n",
      "Matrix >= 1.7 removes the S4 class 'mMatrix', which breaks Graph construction in Seurat v4.\n",
      "Recommended fix (micromamba/conda-forge):\n",
      "  micromamba install -y -n <ENV> -c conda-forge \"r-matrix<1.7\"\n",
      call. = FALSE
    )
  }
}

if (!is.null(seurat_version) && seurat_version < "5.0.0" &&
    !is.null(matrix_version) && matrix_version < "1.7.0") {
  m_ok <- methods::isClass("mMatrix") &&
    !is.null(methods::getClassDef("mMatrix")) &&
    methods::getClassDef("mMatrix")@package == "Matrix"
  if (!isTRUE(m_ok)) {
    stop(
      "Matrix is loaded but the S4 class 'mMatrix' is not defined by the Matrix package.\n",
      "This breaks Seurat v4 Graph validation (invalid class 'Graph': superclass 'mMatrix' not defined).\n",
      "Fix (recommended): recreate a clean micromamba env from env/ezhu-r-seurat4.yml (pins a known-good Matrix build).\n",
      call. = FALSE
    )
  }
}

ezhu_dir_create("data/processed")
singlecell_results_dir <- Sys.getenv("EZHU_SINGLECELL_RESULTS_DIR", unset = "results/singlecell")
ezhu_dir_create(singlecell_results_dir)

singlecell_tag <- Sys.getenv("EZHU_SINGLECELL_TAG", unset = "gse131882")
singlecell_tag <- tolower(trimws(singlecell_tag))
if (!nzchar(singlecell_tag)) singlecell_tag <- "gse131882"

input_dir <- Sys.getenv(
  "EZHU_SINGLECELL_INPUT_DIR",
  unset = Sys.getenv("EZHU_GSE131882_DIR", unset = "data/raw/singlecell/GSE131882/files")
)
files <- sort(
  list.files(
    input_dir,
    pattern = "\\.(rds(\\.gz)?|h5(\\.gz)?)$",
    full.names = TRUE,
    ignore.case = TRUE
  )
)
if (length(files) == 0) {
  stop("No .rds/.h5 files found under: ", input_dir, call. = FALSE)
}

has_h5 <- any(grepl("\\.h5(\\.gz)?$", files, ignore.case = TRUE))
if (has_h5 && !requireNamespace("hdf5r", quietly = TRUE)) {
  stop(
    "Missing package: hdf5r (required to read 10x .h5 matrices).\n",
    "Recommended fix (micromamba/conda-forge):\n",
    "  micromamba install -y -n <ENV> -c conda-forge r-hdf5r hdf5 zlib libaec\n",
    "Then rerun preprocess using scripts/run_r_mamba.sh.\n",
    call. = FALSE
  )
}

parse_csv_env <- function(env_value) {
  value <- trimws(env_value)
  if (!nzchar(value)) return(character())
  parts <- unlist(strsplit(value, "[,;]"))
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

extract_counts <- function(obj) {
  if (inherits(obj, "dgCMatrix") || inherits(obj, "dgTMatrix") || inherits(obj, "dgRMatrix")) {
    return(obj)
  }
  if (is.matrix(obj)) {
    return(Matrix::Matrix(obj, sparse = TRUE))
  }
  if (inherits(obj, "Matrix")) {
    return(obj)
  }

  if (!is.list(obj)) {
    stop("Unsupported object type for count extraction: ", paste(class(obj), collapse = ", "), call. = FALSE)
  }

  candidates <- list(
    c("umicount", "exon", "all"),
    c("umicount", "exon"),
    c("counts"),
    c("data"),
    c("expr"),
    c("expression")
  )

  for (path in candidates) {
    current <- obj
    ok <- TRUE
    for (key in path) {
      if (!is.list(current) || is.null(current[[key]])) {
        ok <- FALSE
        break
      }
      current <- current[[key]]
    }
    if (!ok) next
    if (inherits(current, "Matrix") || is.matrix(current)) {
      return(Matrix::Matrix(current, sparse = TRUE))
    }
  }

  stack <- list(obj)
  visited <- 0L
  max_nodes <- 5000L
  while (length(stack) > 0) {
    current <- stack[[1]]
    stack <- stack[-1]
    visited <- visited + 1L
    if (visited > max_nodes) break

    if (inherits(current, "Matrix") || is.matrix(current)) {
      current_mat <- Matrix::Matrix(current, sparse = TRUE)
      if (!is.null(dim(current_mat)) && length(dim(current_mat)) == 2L) {
        return(current_mat)
      }
    }

    if (is.list(current)) {
      for (child in current) {
        if (is.list(child) || inherits(child, "Matrix") || is.matrix(child)) {
          stack[[length(stack) + 1L]] <- child
        }
      }
    }
  }

  stop("Failed to locate a count matrix inside the RDS object (checked known paths and scanned nested lists).", call. = FALSE)
}

infer_group <- function(sample_name) {
  sample_lc <- tolower(sample_name)
  if (grepl("control", sample_lc, fixed = TRUE)) return("control")
  if (grepl("dn", sample_lc, fixed = TRUE) || grepl("dkd", sample_lc, fixed = TRUE)) return("diabetes")
  if (grepl("diabetes", sample_lc, fixed = TRUE) || grepl("disease", sample_lc, fixed = TRUE)) return("diabetes")
  NA_character_
}

normalize_sample_name <- function(path) {
  base <- basename(path)
  base <- sub("\\.rds(\\.gz)?$", "", base, ignore.case = TRUE)
  base <- sub("\\.h5(\\.gz)?$", "", base, ignore.case = TRUE)
  base <- sub("_filtered_feature_bc_matrix$", "", base, ignore.case = TRUE)
  base <- sub("filtered_feature_bc_matrix$", "", base, ignore.case = TRUE)
  base
}

decompress_h5_gz <- function(path) {
  out <- tempfile(fileext = ".h5")
  in_con <- gzfile(path, open = "rb")
  out_con <- file(out, open = "wb")
  on.exit({
    try(close(in_con), silent = TRUE)
    try(close(out_con), silent = TRUE)
  }, add = TRUE)

  repeat {
    chunk <- readBin(in_con, what = "raw", n = 1024 * 1024)
    if (length(chunk) == 0) break
    writeBin(chunk, out_con)
  }

  out
}

safe_read_rds <- function(path) {
  tryCatch(
    readRDS(path),
    error = function(e) {
      # Be robust to files that are gzip-compressed but lack a ".gz" extension.
      con <- file(path, open = "rb")
      magic <- try(readBin(con, what = "raw", n = 2), silent = TRUE)
      try(close(con), silent = TRUE)

      is_gzip <- !inherits(magic, "try-error") &&
        length(magic) == 2L &&
        as.integer(magic[[1]]) == 0x1f &&
        as.integer(magic[[2]]) == 0x8b

      if (!isTRUE(is_gzip)) stop(e)

      gz_con <- gzfile(path, open = "rb")
      on.exit(try(close(gz_con), silent = TRUE), add = TRUE)
      readRDS(gz_con)
    }
  )
}

sample_table <- data.frame(
  sample = vapply(files, normalize_sample_name, character(1)),
  group = vapply(vapply(files, normalize_sample_name, character(1)), infer_group, character(1)),
  file = files,
  stringsAsFactors = FALSE
)

exclude_samples <- parse_csv_env(Sys.getenv("EZHU_EXCLUDE_SAMPLES", unset = ""))
if (length(exclude_samples) > 0) {
  keep <- !sample_table$sample %in% exclude_samples
  dropped <- sample_table$sample[!keep]
  sample_table <- sample_table[keep, , drop = FALSE]
  message("Excluding ", length(dropped), " samples via EZHU_EXCLUDE_SAMPLES: ", paste(dropped, collapse = ", "))
}

message("Found ", nrow(sample_table), " sample files under ", input_dir)
if (nrow(sample_table) == 0) stop("No samples remain after exclusions.", call. = FALSE)

sample_objects <- vector("list", nrow(sample_table))
sample_qc <- data.frame(
  sample = sample_table$sample,
  group = sample_table$group,
  n_genes_raw = NA_integer_,
  n_cells_raw = NA_integer_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(sample_table))) {
  sample_name <- sample_table$sample[i]
  sample_file <- sample_table$file[i]
  message("Loading: ", sample_name)

  if (grepl("\\.h5(\\.gz)?$", sample_file, ignore.case = TRUE)) {
    h5_path <- sample_file
    if (grepl("\\.gz$", sample_file, ignore.case = TRUE)) {
      message("Decompressing: ", basename(sample_file))
      h5_path <- decompress_h5_gz(sample_file)
    }

    obj <- Seurat::Read10X_h5(h5_path)
    if (is.list(obj)) {
      if ("Gene Expression" %in% names(obj)) {
        obj <- obj[["Gene Expression"]]
      } else if ("RNA" %in% names(obj)) {
        obj <- obj[["RNA"]]
      } else {
        obj <- obj[[1]]
      }
    }
    counts <- Matrix::Matrix(obj, sparse = TRUE)
  } else {
    obj <- safe_read_rds(sample_file)
    counts <- extract_counts(obj)
  }

  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("Counts matrix must have rownames (genes) and colnames (cells): ", sample_file, call. = FALSE)
  }

	  ens_ids <- rownames(counts)
	  if (any(grepl("^ENSG", ens_ids))) {
      if (!requireNamespace("AnnotationDbi", quietly = TRUE) || !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        stop(
          "Missing Bioconductor annotation packages required to map Ensembl IDs to symbols.\n",
          "Install one of the following:\n",
          "  - micromamba (recommended): bioconductor-annotationdbi bioconductor-org.hs.eg.db\n",
          "  - R/BiocManager: BiocManager::install(c('AnnotationDbi','org.Hs.eg.db'))\n",
          call. = FALSE
        )
      }
	    message("Mapping Ensembl IDs to Symbols for: ", sample_name)
	    ens_keys <- sub("\\.[0-9]+$", "", ens_ids)
	    symbols <- tryCatch({
	       AnnotationDbi::mapIds(
	         org.Hs.eg.db::org.Hs.eg.db,
	         keys = ens_keys,
	         column = "SYMBOL",
	         keytype = "ENSEMBL",
	         multiVals = "first"
	       )
	    }, error = function(e) {
	       warning("Mapping failed: ", e$message)
	       return(NULL)
	    })
	    if (!is.null(symbols)) {
	       symbols[is.na(symbols)] <- ens_ids[is.na(symbols)]
	       rownames(counts) <- make.unique(as.character(symbols))
	    }
	  }

  colnames(counts) <- paste0(sample_name, "_", colnames(counts))

  sample_qc$n_genes_raw[i] <- nrow(counts)
  sample_qc$n_cells_raw[i] <- ncol(counts)

  seurat_obj <- Seurat::CreateSeuratObject(counts = counts, assay = "RNA", project = sample_name)
  seurat_obj$sample <- sample_name
  seurat_obj$group <- sample_table$group[i]
  sample_objects[[i]] <- seurat_obj
}

counts_qc <- data.frame(
  sample = as.character(sample_qc$sample),
  group = as.character(sample_qc$group),
  n_genes = as.integer(sample_qc$n_genes_raw),
  n_cells = as.integer(sample_qc$n_cells_raw),
  stringsAsFactors = FALSE
)
utils::write.table(
  counts_qc,
  file = file.path(singlecell_results_dir, paste0(singlecell_tag, "_counts_qc.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

merged <- sample_objects[[1]]
if (length(sample_objects) > 1) {
  for (j in 2:length(sample_objects)) {
    merged <- merge(merged, y = sample_objects[[j]])
  }
}

merged$group <- as.character(merged$group)
if (anyNA(merged$group)) {
  warning("Some cells have missing group labels; check sample naming conventions.", call. = FALSE)
}

assay <- Sys.getenv("EZHU_SC_ASSAY", unset = "RNA")
Seurat::DefaultAssay(merged) <- assay

mt_pattern <- Sys.getenv("EZHU_MT_PATTERN", unset = "^MT-")
merged[["percent.mt"]] <- Seurat::PercentageFeatureSet(merged, pattern = mt_pattern)

integration_method <- tolower(trimws(Sys.getenv("EZHU_INTEGRATION_METHOD", unset = "seurat")))
if (!nzchar(integration_method)) integration_method <- "seurat"
if (integration_method %in% c("off", "false", "0")) integration_method <- "none"
if (integration_method %in% c("merge", "none")) integration_method <- "merge"
if (!integration_method %in% c("seurat", "merge")) {
  stop(
    "Invalid EZHU_INTEGRATION_METHOD: ", integration_method, "\n",
    "Valid values: seurat, merge",
    call. = FALSE
  )
}

qc_min_features <- suppressWarnings(as.integer(Sys.getenv("EZHU_QC_MIN_FEATURES", unset = "200")))
if (is.na(qc_min_features)) qc_min_features <- 200L
qc_min_counts <- suppressWarnings(as.numeric(Sys.getenv("EZHU_QC_MIN_COUNTS", unset = "0")))
if (is.na(qc_min_counts)) qc_min_counts <- 0
qc_max_percent_mt <- suppressWarnings(as.numeric(Sys.getenv("EZHU_QC_MAX_PERCENT_MT", unset = "20")))
if (is.na(qc_max_percent_mt)) qc_max_percent_mt <- 20
qc_max_features <- suppressWarnings(as.integer(Sys.getenv("EZHU_QC_MAX_FEATURES", unset = "")))
has_max_features <- is.finite(qc_max_features) && !is.na(qc_max_features) && qc_max_features > 0

pre_cells <- ncol(merged)

subset_expr <- merged$nFeature_RNA >= qc_min_features &
  merged$nCount_RNA >= qc_min_counts &
  merged$percent.mt <= qc_max_percent_mt
if (has_max_features) subset_expr <- subset_expr & merged$nFeature_RNA <= qc_max_features

merged <- merged[, which(subset_expr)]
post_cells <- ncol(merged)

message("QC filter: kept ", post_cells, " / ", pre_cells, " cells")

pca_dims <- suppressWarnings(as.integer(Sys.getenv("EZHU_PCA_DIMS", unset = "30")))
if (is.na(pca_dims) || pca_dims <= 0) pca_dims <- 30L

resolution <- suppressWarnings(as.numeric(Sys.getenv("EZHU_CLUSTER_RESOLUTION", unset = "0.5")))
if (is.na(resolution) || resolution <= 0) resolution <- 0.5

dims_use <- seq_len(min(30L, pca_dims))

integration_min_cells <- suppressWarnings(as.integer(Sys.getenv("EZHU_INTEGRATION_MIN_CELLS", unset = "50")))
if (is.na(integration_min_cells) || integration_min_cells <= 0) integration_min_cells <- 50L

if (integration_method == "seurat" && length(sample_table$sample) > 1L) {
  min_cells_by_sample <- tryCatch(as.integer(min(table(merged$sample))), error = function(e) NA_integer_)
  if (is.finite(min_cells_by_sample) && min_cells_by_sample < integration_min_cells) {
    warning(
      "Integration disabled: too few cells in at least one sample after QC (min_cells=",
      min_cells_by_sample,
      ", required>=",
      integration_min_cells,
      "). Falling back to merge-only pipeline.",
      call. = FALSE
    )
    integration_method <- "merge"
  }
}

if (integration_method == "seurat" && length(sample_table$sample) > 1L) {
  message("Integration: Seurat anchors (split.by=sample)")

  obj_list <- Seurat::SplitObject(merged, split.by = "sample")
  obj_list <- lapply(obj_list, function(obj) {
    Seurat::DefaultAssay(obj) <- assay
    obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    obj
  })

  integration_nfeatures <- suppressWarnings(as.integer(Sys.getenv("EZHU_INTEGRATION_NFEATURES", unset = "2000")))
  if (is.na(integration_nfeatures) || integration_nfeatures <= 0) integration_nfeatures <- 2000L
  integration_nfeatures <- min(integration_nfeatures, nrow(merged))

  min_cells_by_sample <- min(vapply(obj_list, function(obj) as.integer(ncol(obj)), integer(1)), na.rm = TRUE)
  k_max <- max(1L, min_cells_by_sample - 1L)
  k_filter <- suppressWarnings(as.integer(Sys.getenv("EZHU_INTEGRATION_K_FILTER", unset = "200")))
  if (is.na(k_filter) || k_filter <= 0) k_filter <- 200L
  k_filter <- min(k_filter, k_max)
  k_anchor <- suppressWarnings(as.integer(Sys.getenv("EZHU_INTEGRATION_K_ANCHOR", unset = "5")))
  if (is.na(k_anchor) || k_anchor <= 0) k_anchor <- 5L
  k_anchor <- min(k_anchor, k_max)
  k_score <- suppressWarnings(as.integer(Sys.getenv("EZHU_INTEGRATION_K_SCORE", unset = "30")))
  if (is.na(k_score) || k_score <= 0) k_score <- 30L
  k_score <- min(k_score, k_max)

  k_weight <- suppressWarnings(as.integer(Sys.getenv("EZHU_INTEGRATION_K_WEIGHT", unset = "20")))
  if (is.na(k_weight) || k_weight <= 0) k_weight <- 20L
  k_weight <- min(k_weight, k_max)

  integration_features <- Seurat::SelectIntegrationFeatures(object.list = obj_list, nfeatures = integration_nfeatures)
  anchors <- Seurat::FindIntegrationAnchors(
    object.list = obj_list,
    anchor.features = integration_features,
    dims = dims_use,
    k.filter = k_filter,
    k.anchor = k_anchor,
    k.score = k_score
  )
  merged <- Seurat::IntegrateData(anchorset = anchors, dims = dims_use, k.weight = k_weight)

  Seurat::DefaultAssay(merged) <- "integrated"
  Seurat::VariableFeatures(merged) <- integration_features
  merged <- Seurat::ScaleData(merged, features = integration_features, verbose = FALSE)
  merged <- Seurat::RunPCA(merged, npcs = pca_dims, features = integration_features, verbose = FALSE)
  merged <- Seurat::FindNeighbors(merged, dims = dims_use, verbose = FALSE)
  merged <- Seurat::FindClusters(merged, resolution = resolution, verbose = FALSE)
  merged <- Seurat::RunUMAP(merged, dims = dims_use, verbose = FALSE)

  Seurat::DefaultAssay(merged) <- assay
} else {
  if (integration_method != "merge") {
    message("Integration: disabled (n_samples<=1) -> merge-only pipeline")
  } else {
    message("Integration: merge-only pipeline")
  }

  merged <- Seurat::NormalizeData(merged, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  merged <- Seurat::FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  merged <- Seurat::ScaleData(merged, features = Seurat::VariableFeatures(merged), verbose = FALSE)
  merged <- Seurat::RunPCA(merged, npcs = pca_dims, verbose = FALSE)
  merged <- Seurat::FindNeighbors(merged, dims = dims_use, verbose = FALSE)
  merged <- Seurat::FindClusters(merged, resolution = resolution, verbose = FALSE)
  merged <- Seurat::RunUMAP(merged, dims = dims_use, verbose = FALSE)
}

out_rds <- Sys.getenv("EZHU_SINGLECELL_PROCESSED", unset = "data/processed/singlecell_gse131882_seurat.rds")
saveRDS(merged, file = out_rds)

meta <- merged[[]]
meta$cluster <- as.character(meta$seurat_clusters)
meta$group <- as.character(meta$group)
meta$sample <- as.character(meta$sample)

cluster_by_condition <- aggregate(
  x = list(n_cells = rep(1L, nrow(meta))),
  by = list(cluster = meta$cluster, group = meta$group),
  FUN = sum
)
cluster_by_condition <- cluster_by_condition[order(cluster_by_condition$cluster, cluster_by_condition$group), , drop = FALSE]
utils::write.table(
  cluster_by_condition,
  file = file.path(singlecell_results_dir, paste0(singlecell_tag, "_cluster_by_condition.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

kept_by_sample <- aggregate(
  x = list(n_cells_kept = rep(1L, nrow(meta))),
  by = list(sample = meta$sample, group = meta$group),
  FUN = sum
)

qc_summary <- merge(sample_qc, kept_by_sample, by = c("sample", "group"), all.x = TRUE)
qc_summary$n_cells_kept[is.na(qc_summary$n_cells_kept)] <- 0L
qc_summary$pct_cells_kept <- ifelse(qc_summary$n_cells_raw > 0, qc_summary$n_cells_kept / qc_summary$n_cells_raw, NA_real_)
qc_summary$qc_min_features <- qc_min_features
qc_summary$qc_min_counts <- qc_min_counts
qc_summary$qc_max_percent_mt <- qc_max_percent_mt
qc_summary$qc_max_features <- if (has_max_features) qc_max_features else NA_integer_
qc_summary$cluster_resolution <- resolution
qc_summary$pca_dims <- pca_dims

utils::write.table(
  qc_summary,
  file = file.path(singlecell_results_dir, paste0(singlecell_tag, "_atlas_qc_summary.tsv")),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ezhu_write_stage_metadata(
  "02_preprocess",
  params = list(
    singlecell_tag = singlecell_tag,
    input_dir = input_dir,
    n_samples = nrow(sample_table),
    exclude_samples = if (length(exclude_samples) > 0) exclude_samples else NULL,
    qc = list(
      min_features = qc_min_features,
      min_counts = qc_min_counts,
      max_percent_mt = qc_max_percent_mt,
      max_features = if (has_max_features) qc_max_features else NA_integer_
    ),
    pre_cells = pre_cells,
    post_cells = post_cells,
    assay = assay,
    integration_method = integration_method,
    pca_dims = pca_dims,
    cluster_resolution = resolution,
    output_rds = out_rds,
    results_dir = singlecell_results_dir,
    outputs = list(
      counts_qc_tsv = file.path(singlecell_results_dir, paste0(singlecell_tag, "_counts_qc.tsv")),
      qc_summary_tsv = file.path(singlecell_results_dir, paste0(singlecell_tag, "_atlas_qc_summary.tsv")),
      cluster_by_condition_tsv = file.path(singlecell_results_dir, paste0(singlecell_tag, "_cluster_by_condition.tsv"))
    )
  ),
  seed = seed
)
