args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
  repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
  if (dir.exists(repo_root)) setwd(repo_root)
}

source("scripts/utils/repro.R")
source("scripts/utils/manifest.R")
ezhu_set_repo_root()
ezhu_activate_renv()

seed_env <- Sys.getenv("EZHU_SEED", unset = "")
seed <- suppressWarnings(as.integer(seed_env))
if (is.na(seed)) seed <- 20251227L
set.seed(seed)

options(
  stringsAsFactors = FALSE,
  width = 120,
  warn = 0
)

options(lifecycle_verbosity = "quiet")

if (!requireNamespace("scPagwas", quietly = TRUE)) stop("Missing package: scPagwas", call. = FALSE)
if (!requireNamespace("Seurat", quietly = TRUE)) stop("Missing package: Seurat", call. = FALSE)
if (!requireNamespace("bigreadr", quietly = TRUE)) stop("Missing package: bigreadr", call. = FALSE)

# Fail-fast: scPagwas pipeline is validated against the Seurat v4 stack.
# Mixed Seurat (v4) + SeuratObject (v5) environments can behave unpredictably.
require_seurat4 <- !identical(Sys.getenv("EZHU_REQUIRE_SEURAT4", unset = "1"), "0")
if (isTRUE(require_seurat4)) {
  so_ver <- tryCatch(as.character(utils::packageVersion("SeuratObject")), error = function(e) NA_character_)
  s_ver <- tryCatch(as.character(utils::packageVersion("Seurat")), error = function(e) NA_character_)
  m_ver <- tryCatch(as.character(utils::packageVersion("Matrix")), error = function(e) NA_character_)

  bad <- FALSE
  if (nzchar(so_ver) && !is.na(so_ver) && utils::compareVersion(so_ver, "5.0.0") >= 0) bad <- TRUE
  if (nzchar(s_ver) && !is.na(s_ver) && utils::compareVersion(s_ver, "5.0.0") >= 0) bad <- TRUE
  if (nzchar(m_ver) && !is.na(m_ver) && utils::compareVersion(m_ver, "1.7-0") >= 0) bad <- TRUE

  if (isTRUE(bad)) {
    stop(
      "Unsupported R package stack for scPagwas (expected Seurat/SeuratObject < 5; Matrix < 1.7).\n",
      "Detected:\n",
      "  Seurat=", s_ver, "\n",
      "  SeuratObject=", so_ver, "\n",
      "  Matrix=", m_ver, "\n",
      "Fix:\n",
      "  - Use micromamba env from env/ezhu-r-seurat4.yml and run via scripts/run_r_mamba.sh\n",
      "  - Or set EZHU_REQUIRE_SEURAT4=0 to bypass (not recommended).",
      call. = FALSE
    )
  }
}

parse_csv_env <- function(env_value) {
  value <- trimws(env_value)
  if (!nzchar(value)) return(character())
  parts <- unlist(strsplit(value, "[,;]"))
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

safe_get_pkg_object <- function(package, name) {
  ns <- asNamespace(package)
  if (exists(name, envir = ns, inherits = FALSE)) {
    return(get(name, envir = ns, inherits = FALSE))
  }

  pkg_env <- tryCatch(as.environment(paste0("package:", package)), error = function(e) NULL)
  if (!is.null(pkg_env) && exists(name, envir = pkg_env, inherits = FALSE)) {
    return(get(name, envir = pkg_env, inherits = FALSE))
  }

  data_env <- new.env(parent = emptyenv())
  loaded <- tryCatch(
    {
      utils::data(list = name, package = package, envir = data_env)
      TRUE
    },
    error = function(e) FALSE
  )
  if (isTRUE(loaded) && exists(name, envir = data_env, inherits = FALSE)) {
    return(get(name, envir = data_env, inherits = FALSE))
  }

  NULL
}

ezhu_dir_create("results/scpagwas")

manifest <- ezhu_read_manifest("data/manifest.tsv")
active <- function(df) {
  fields <- c("type", "name", "source", "accession", "url", "local_path")
  fields <- intersect(fields, names(df))
  if (length(fields) == 0) return(rep(FALSE, nrow(df)))
  apply(df[, fields, drop = FALSE], 1, function(row) any(!is.na(row) & nzchar(trimws(as.character(row)))))
}
active_rows <- active(manifest)
gwas_rows <- which(active_rows & tolower(manifest$type) == "gwas")
if (length(gwas_rows) == 0) stop("No GWAS rows found in data/manifest.tsv (type=gwas).", call. = FALSE)

resolve_manifest_gwas_row <- function(manifest, active_rows, gwas_rows) {
  select_id <- Sys.getenv("EZHU_GWAS_MANIFEST_ID", unset = Sys.getenv("EZHU_GWAS_ID", unset = ""))
  select_id <- trimws(as.character(select_id))
  if (nzchar(select_id)) {
    matched <- which(active_rows & tolower(manifest$type) == "gwas" & as.character(manifest$id) == select_id)
    if (length(matched) == 0) {
      stop("EZHU_GWAS_MANIFEST_ID did not match any active GWAS row in data/manifest.tsv: ", select_id, call. = FALSE)
    }
    if (length(matched) > 1) {
      stop("EZHU_GWAS_MANIFEST_ID matched multiple rows; manifest IDs must be unique: ", select_id, call. = FALSE)
    }
    return(manifest[matched, , drop = FALSE])
  }

  # Fallback: try to infer GWAS ID from harmonized file name pattern gwas_harmonized__<ID>.tsv.gz
  gwas_path_env <- Sys.getenv("EZHU_GWAS_HARMONIZED", unset = "data/processed/gwas_harmonized.tsv.gz")
  base <- basename(trimws(as.character(gwas_path_env)))
  inferred <- sub("^gwas_harmonized__", "", base)
  inferred <- sub("\\.tsv\\.gz$", "", inferred, ignore.case = TRUE)
  inferred <- trimws(inferred)
  if (nzchar(inferred) && !identical(inferred, base)) {
    matched <- which(active_rows & tolower(manifest$type) == "gwas" & as.character(manifest$id) == inferred)
    if (length(matched) == 1) {
      return(manifest[matched, , drop = FALSE])
    }
  }

  if (length(gwas_rows) > 1) {
    message("Multiple GWAS rows found; using the first one. Set EZHU_GWAS_MANIFEST_ID to choose a specific GWAS.")
  }
  manifest[gwas_rows[1], , drop = FALSE]
}

row <- resolve_manifest_gwas_row(manifest, active_rows, gwas_rows)
notes <- ezhu_parse_notes_kv(row$notes)
build <- tolower(trimws(notes$build %||% ""))
if (!build %in% c("hg37", "hg38")) stop("GWAS build must be declared as hg37 or hg38 in manifest notes.", call. = FALSE)

# Fail-fast: ensure the resolved manifest row is consistent with env/file naming.
select_id_env <- trimws(as.character(Sys.getenv("EZHU_GWAS_MANIFEST_ID", unset = Sys.getenv("EZHU_GWAS_ID", unset = ""))))
if (nzchar(select_id_env)) {
  if (!identical(as.character(row$id), select_id_env)) {
    stop(
      "GWAS selection mismatch:\n",
      "  EZHU_GWAS_MANIFEST_ID=", select_id_env, "\n",
      "  resolved manifest row id=", as.character(row$id), "\n",
      "This indicates a manifest lookup bug or a misconfigured environment.",
      call. = FALSE
    )
  }
}

gwas_path_env0 <- Sys.getenv("EZHU_GWAS_HARMONIZED", unset = "data/processed/gwas_harmonized.tsv.gz")
base0 <- basename(trimws(as.character(gwas_path_env0)))
inferred0 <- sub("^gwas_harmonized__", "", base0)
inferred0 <- sub("\\.tsv\\.gz$", "", inferred0, ignore.case = TRUE)
inferred0 <- trimws(inferred0)
if (nzchar(inferred0) && !identical(inferred0, base0)) {
  if (!identical(as.character(row$id), inferred0)) {
    stop(
      "GWAS selection mismatch:\n",
      "  resolved manifest row id=", as.character(row$id), "\n",
      "  EZHU_GWAS_HARMONIZED basename implies id=", inferred0, " (", base0, ")\n",
      "Fix: pass consistent EZHU_GWAS_MANIFEST_ID and EZHU_GWAS_HARMONIZED.",
      call. = FALSE
    )
  }
}

gwas_path <- Sys.getenv("EZHU_GWAS_HARMONIZED", unset = "data/processed/gwas_harmonized.tsv.gz")
if (!file.exists(gwas_path)) stop("Missing GWAS harmonized file: ", gwas_path, call. = FALSE)

singlecell_path <- Sys.getenv("EZHU_SINGLECELL_PROCESSED", unset = "data/processed/singlecell_gse131882_seurat.rds")
if (!file.exists(singlecell_path)) stop("Missing processed single-cell Seurat object: ", singlecell_path, call. = FALSE)

singlecell_tag <- Sys.getenv("EZHU_SINGLECELL_TAG", unset = "gse131882")
singlecell_tag <- tolower(trimws(singlecell_tag))
if (!nzchar(singlecell_tag)) singlecell_tag <- "gse131882"

default_prefix <- paste0(singlecell_tag, "_", tolower(as.character(row$id)))
output_prefix <- Sys.getenv("EZHU_SCPAGWAS_PREFIX", unset = default_prefix)
output_dir <- Sys.getenv("EZHU_SCPAGWAS_DIR", unset = "results/scpagwas")
ezhu_dir_create(output_dir)

n_cores <- suppressWarnings(as.integer(Sys.getenv("EZHU_N_CORES", unset = "1")))
if (is.na(n_cores) || n_cores <= 0) n_cores <- 1L

iters_singlecell <- suppressWarnings(as.integer(Sys.getenv("EZHU_SCPAGWAS_ITERS_SINGLECELL", unset = "50")))
if (is.na(iters_singlecell) || iters_singlecell < 0) iters_singlecell <- 50L

iters_celltype <- suppressWarnings(as.integer(Sys.getenv("EZHU_SCPAGWAS_ITERS_CELLTYPE", unset = "50")))
if (is.na(iters_celltype) || iters_celltype < 0) iters_celltype <- 50L

maf_filter <- suppressWarnings(as.numeric(Sys.getenv("EZHU_SCPAGWAS_MAF_FILTER", unset = "0.01")))
if (is.na(maf_filter) || maf_filter <= 0) maf_filter <- 0.01

assay <- Sys.getenv("EZHU_SCPAGWAS_ASSAY", unset = "RNA")

run_singlecell <- !identical(Sys.getenv("EZHU_SCPAGWAS_SINGLECELL", unset = "1"), "0")
run_celltype <- !identical(Sys.getenv("EZHU_SCPAGWAS_CELLTYPE", unset = "1"), "0")
if (!run_singlecell && !run_celltype) {
  stop("Both EZHU_SCPAGWAS_SINGLECELL and EZHU_SCPAGWAS_CELLTYPE were set to 0; nothing to run.", call. = FALSE)
}

message("Reading GWAS harmonized: ", gwas_path)
gwas_df <- bigreadr::fread2(gwas_path)

required <- c("chrom", "pos", "rsid", "se", "beta", "maf")
missing_cols <- setdiff(required, colnames(gwas_df))
if (length(missing_cols) > 0) {
  stop("GWAS harmonized file is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

normalize_gwas_chrom <- function(x) {
  x0 <- trimws(as.character(x))
  x0 <- sub("^chr", "", x0, ignore.case = TRUE)
  is_num <- grepl("^[0-9]+$", x0)
  if (any(is_num)) {
    x0[is_num] <- sub("^0+", "", x0[is_num])
    x0[is_num & !nzchar(x0)] <- NA_character_
  }
  x0
}

chrom0 <- as.character(gwas_df$chrom)
chrom1 <- normalize_gwas_chrom(chrom0)
if (any(is.na(chrom1))) {
  stop(
    "GWAS harmonized 'chrom' contains invalid/empty chromosome values after normalization.\n",
    "Hint: check for leading zeros (e.g., 01) or malformed chromosome strings in the raw GWAS.\n",
    "File: ", gwas_path,
    call. = FALSE
  )
}
if (!identical(chrom0, chrom1)) {
  message("NOTE: Normalized GWAS chromosome values (e.g., removed leading zeros).")
}
gwas_df$chrom <- chrom1

message("Reading single-cell object: ", singlecell_path)
single_obj <- readRDS(singlecell_path)
if (!inherits(single_obj, "Seurat")) stop("Processed single-cell object is not a Seurat object: ", singlecell_path, call. = FALSE)
if (!assay %in% Seurat::Assays(single_obj)) stop("Assay not present in Seurat object: ", assay, call. = FALSE)

# scPagwas internally relies on Seurat's DefaultAssay + VariableFeatures.
# Ensure these are aligned to the assay we're passing to scPagwas to avoid
# "no match for Pathway gene and VariableFeatures" failures.
Seurat::DefaultAssay(single_obj) <- assay

# Optional smoke test mode: subset cells to reach the scoring stages quickly.
smoke <- identical(Sys.getenv("EZHU_SCPAGWAS_SMOKE", unset = ""), "1")
smoke_n <- suppressWarnings(as.integer(Sys.getenv("EZHU_SCPAGWAS_SMOKE_N_CELLS", unset = "2000")))
if (isTRUE(smoke)) {
  if (!is.finite(smoke_n) || is.na(smoke_n) || smoke_n <= 0) smoke_n <- 2000L
  nc <- ncol(single_obj)
  if (nc > smoke_n) {
    set.seed(seed)
    keep_cells <- sample(colnames(single_obj), size = smoke_n, replace = FALSE)
    single_obj <- subset(single_obj, cells = keep_cells)
    message("SMOKE mode enabled: subset to ", ncol(single_obj), " cells (from ", nc, ").")
    # Also reduce iterations for faster execution.
    iters_singlecell <- min(iters_singlecell, 3L)
    iters_celltype <- min(iters_celltype, 3L)
    n_cores <- 1L
  }
}

exclude_samples <- parse_csv_env(Sys.getenv("EZHU_EXCLUDE_SAMPLES", unset = ""))
exclude_low_retention <- !identical(Sys.getenv("EZHU_EXCLUDE_LOW_RETENTION_SAMPLES", unset = "1"), "0")
min_sample_retention <- suppressWarnings(as.numeric(Sys.getenv("EZHU_MIN_SAMPLE_RETENTION", unset = "0.2")))
if (!is.finite(min_sample_retention) || min_sample_retention <= 0) min_sample_retention <- 0.2

if (length(exclude_samples) == 0 && isTRUE(exclude_low_retention)) {
  singlecell_results_dir <- Sys.getenv("EZHU_SINGLECELL_RESULTS_DIR", unset = "results/singlecell")
  qc_summary_path <- file.path(singlecell_results_dir, paste0(singlecell_tag, "_atlas_qc_summary.tsv"))
  if (file.exists(qc_summary_path)) {
    qc <- tryCatch(
      utils::read.delim(qc_summary_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE),
      error = function(e) NULL
    )
    if (!is.null(qc) && all(c("sample", "pct_cells_kept") %in% names(qc))) {
      qc$sample <- as.character(qc$sample)
      qc$pct_cells_kept <- suppressWarnings(as.numeric(qc$pct_cells_kept))
      low <- qc$sample[is.finite(qc$pct_cells_kept) & qc$pct_cells_kept < min_sample_retention]
      low <- unique(low[nzchar(low)])
      if (length(low) > 0) {
        exclude_samples <- low
        message(
          "Auto-excluding low-retention samples (pct_cells_kept < ",
          min_sample_retention,
          "): ",
          paste(exclude_samples, collapse = ", ")
        )
      }
    } else {
      warning("QC summary missing required columns (sample, pct_cells_kept): ", qc_summary_path, call. = FALSE)
    }
  } else {
    warning(
      "QC summary not found for auto-exclusion: ",
      qc_summary_path,
      " (set EZHU_EXCLUDE_LOW_RETENTION_SAMPLES=0 to silence).",
      call. = FALSE
    )
  }
}

if (length(exclude_samples) > 0) {
  if (!"sample" %in% colnames(single_obj[[]])) {
    stop("EZHU_EXCLUDE_SAMPLES was set, but the Seurat object has no 'sample' column in meta.data.", call. = FALSE)
  }
  keep <- !as.character(single_obj$sample) %in% exclude_samples
  dropped <- unique(as.character(single_obj$sample)[!keep])
  single_obj <- single_obj[, which(keep)]
  message("Excluded ", length(dropped), " samples via EZHU_EXCLUDE_SAMPLES: ", paste(dropped, collapse = ", "))
}

if (!"package:scPagwas" %in% search()) {
  suppressPackageStartupMessages(library(scPagwas))
}

idents_field <- Sys.getenv("EZHU_IDENTS_FIELD", unset = "seurat_clusters")
if (idents_field %in% colnames(single_obj[[]])) {
  Seurat::Idents(single_obj) <- as.factor(single_obj[[idents_field, drop = TRUE]])
} else {
  warning("Idents field not found in meta.data: ", idents_field, " (keeping existing Idents).\n", call. = FALSE)
}

block_name <- Sys.getenv(
  "EZHU_BLOCK_ANNOTATION_NAME",
  unset = if (build == "hg37") "block_annotation_hg37" else "block_annotation"
)
block_annotation <- safe_get_pkg_object("scPagwas", block_name)
if (is.null(block_annotation)) {
  stop(
    "Cannot load scPagwas block annotation '", block_name, "'. ",
    "Try setting EZHU_BLOCK_ANNOTATION_NAME or provide a file path via EZHU_BLOCK_ANNOTATION_PATH.",
    call. = FALSE
  )
}

chrom_ld_name <- Sys.getenv("EZHU_CHROM_LD_NAME", unset = "chrom_ld")
chrom_ld <- safe_get_pkg_object("scPagwas", chrom_ld_name)
if (is.null(chrom_ld)) {
  stop(
    "Cannot load scPagwas LD object '", chrom_ld_name, "'. ",
    "Try setting EZHU_CHROM_LD_NAME.",
    call. = FALSE
  )
}

pathway_name <- Sys.getenv("EZHU_PATHWAY_LIST_NAME", unset = "Genes_by_pathway_kegg")
pathway_list <- safe_get_pkg_object("scPagwas", pathway_name)
if (is.null(pathway_list)) {
  stop(
    "Cannot load scPagwas pathway list '", pathway_name, "'. ",
    "Try setting EZHU_PATHWAY_LIST_NAME.",
    call. = FALSE
  )
}

if (isTRUE(smoke)) {
  smoke_n_pathways <- suppressWarnings(as.integer(Sys.getenv("EZHU_SCPAGWAS_SMOKE_N_PATHWAYS", unset = "200")))
  if (!is.finite(smoke_n_pathways) || is.na(smoke_n_pathways) || smoke_n_pathways < 10) smoke_n_pathways <- 200L
  if (length(pathway_list) > smoke_n_pathways) {
    sizes <- vapply(pathway_list, function(x) length(unique(as.character(x))), integer(1))
    ord <- order(sizes, decreasing = TRUE, na.last = NA)
    ord <- ord[seq_len(min(smoke_n_pathways, length(ord)))]
    pathway_list <- pathway_list[ord]
    message(
      "SMOKE mode: subset pathway list to ", length(pathway_list),
      " largest pathways (EZHU_SCPAGWAS_SMOKE_N_PATHWAYS=", smoke_n_pathways, ")."
    )
  }
}

features <- rownames(single_obj[[assay]])
features_base <- sub("\\.[0-9]+$", "", features)

match_mode <- tolower(trimws(Sys.getenv("EZHU_PATHWAY_MATCH_MODE", unset = "legacy")))
min_pathway_genes <- suppressWarnings(as.integer(Sys.getenv("EZHU_PATHWAY_MIN_GENES", unset = "5")))
if (is.na(min_pathway_genes) || min_pathway_genes < 1) min_pathway_genes <- 5L

if (!match_mode %in% c("legacy", "prune_only", "feature_names")) {
  stop(
    "Invalid EZHU_PATHWAY_MATCH_MODE: ", match_mode, "\n",
    "Valid values: legacy, prune_only, feature_names",
    call. = FALSE
  )
}

pathway_overlap <- data.frame(
  pathway = names(pathway_list),
  n_input_genes = as.integer(vapply(pathway_list, length, integer(1))),
  n_matched_features = NA_integer_,
  stringsAsFactors = FALSE
)

if (match_mode != "legacy") {
  # Map pathway gene symbols onto the actual feature names present in the Seurat object.
  # This is robust to duplicated symbols made unique via ".1" suffixing.
  mapped <- list()
  kept_names <- character()
  for (nm in names(pathway_list)) {
    genes <- unique(as.character(pathway_list[[nm]]))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (length(genes) == 0) next

    if (match_mode == "feature_names") {
      sel <- features[features_base %in% genes]
    } else {
      # prune_only: keep original gene names, but drop any that have no match in the data
      sel <- genes[genes %in% features_base]
    }

    sel <- unique(sel)
    pathway_overlap$n_matched_features[pathway_overlap$pathway == nm] <- length(sel)

    if (length(sel) >= min_pathway_genes) {
      mapped[[nm]] <- sel
      kept_names <- c(kept_names, nm)
    }
  }
  pathway_list <- mapped
  pathway_overlap <- pathway_overlap[pathway_overlap$pathway %in% kept_names, , drop = FALSE]
}

pathway_genes <- unique(as.character(unlist(pathway_list, use.names = FALSE)))
pathway_genes_base <- unique(sub("\\.[0-9]+$", "", pathway_genes))
overlap <- sum(features %in% pathway_genes)
message("Pathway gene overlap: ", overlap, " / ", length(features), " features in assay ", assay)
if (overlap < 200L) {
  stop(
    "Too few pathway genes overlap with Seurat feature names (", overlap,
    "). Ensure assay rownames are gene symbols that match Genes_by_pathway_kegg (not Ensembl IDs with version suffix).",
    call. = FALSE
  )
}

# Guard: if VariableFeatures() is empty or doesn't overlap pathway genes, set it to all features.
vf <- tryCatch(Seurat::VariableFeatures(single_obj), error = function(e) character())
vf <- unique(trimws(as.character(vf)))
vf <- vf[nzchar(vf)]
vf_overlap <- if (length(vf) == 0) 0L else sum(vf %in% pathway_genes)
if (length(vf) == 0 || vf_overlap == 0L) {
  vf_target <- unique(intersect(features, pathway_genes))
  if (length(vf_target) == 0L) {
    # Try mapping pathway gene symbols to Seurat feature names with ".1" suffixes.
    vf_target <- unique(features[features_base %in% pathway_genes])
  }
  if (length(vf_target) == 0L) {
    vf_target <- features
  }

  Seurat::VariableFeatures(single_obj) <- vf_target
  message(
    "WARN: VariableFeatures() had no overlap with pathway genes (vf_n=", length(vf),
    ", vf_overlap=", vf_overlap, "); set VariableFeatures to ", length(vf_target), " features (pathway-aware fallback)."
  )
}

if (match_mode != "legacy") {
  overlap_path <- file.path(output_dir, paste0("pathway_overlap_", singlecell_tag, ".tsv"))
  utils::write.table(pathway_overlap, file = overlap_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote: ", overlap_path)
}

patch_pathway_annotation <- !identical(Sys.getenv("EZHU_SCPAGWAS_PATCH_PATHWAY_ANNOTATION", unset = "1"), "0")
if (patch_pathway_annotation) {
  ns <- asNamespace("scPagwas")
  if (exists("Pathway_annotation_input", envir = ns, inherits = FALSE)) {
    ezhu_patch_pathway_annotation_input <- function() {
      orig <- get("Pathway_annotation_input", envir = ns, inherits = FALSE)
      if (!is.function(orig)) return(invisible(FALSE))

      Pathway_annotation_input <- function(Pagwas, block_annotation, ...) {
        try_call <- function() orig(Pagwas = Pagwas, block_annotation = block_annotation, ...)

        tryCatch(
          try_call(),
          error = function(e) {
            msg <- conditionMessage(e)
            if (!grepl("no match for Pathway gene and VariableFeatures", msg, fixed = TRUE)) {
              stop(e)
            }

            pathway_genes <- unique(as.character(unlist(Pagwas$Pathway_list, use.names = FALSE)))
            pathway_genes <- pathway_genes[!is.na(pathway_genes) & nzchar(pathway_genes)]

            snp_df <- Pagwas$snp_gene_df
            if (is.null(snp_df) || !is.data.frame(snp_df) || nrow(snp_df) == 0 || !"label" %in% names(snp_df)) {
              stop(
                "scPagwas::Pathway_annotation_input failed: no SNP→gene mappings are available (snp_gene_df is empty).\n",
                "This usually indicates a chromosome/build mismatch between the GWAS and scPagwas block annotations.\n",
                "Quick checks:\n",
                "  - Are GWAS chromosomes formatted like 01/02 (should be 1/2)?\n",
                "  - Is the declared GWAS build correct (hg37 vs hg38)?\n",
                call. = FALSE
              )
            }

            snp_labels <- unique(trimws(as.character(snp_df$label)))
            snp_labels <- snp_labels[!is.na(snp_labels) & nzchar(snp_labels)]

            strip_suffix <- function(x) sub("\\.[0-9]+$", "", x)
            overlap_raw <- length(intersect(pathway_genes, snp_labels))
            overlap_base <- length(intersect(strip_suffix(pathway_genes), strip_suffix(snp_labels)))

            if (overlap_raw == 0L && overlap_base == 0L) {
              stop(
                "scPagwas::Pathway_annotation_input failed: SNP-mapped genes do not overlap pathway genes.\n",
                "Diagnostics:\n",
                "  pathway_genes=", length(pathway_genes), "\n",
                "  snp_gene_labels=", length(snp_labels), "\n",
                "  overlap_raw=", overlap_raw, " overlap_base=", overlap_base, "\n",
                "Examples:\n",
                "  pathway_genes[1:5]=", paste(utils::head(pathway_genes, 5), collapse = ","), "\n",
                "  snp_gene_df$label[1:5]=", paste(utils::head(snp_labels, 5), collapse = ","), "\n",
                "Hint: this is typically a gene-identifier mismatch (e.g., Ensembl vs HGNC) or an incorrect GWAS build.\n",
                call. = FALSE
              )
            }

            # Repair: ensure Pagwas$VariableFeatures overlaps the SNP-mapped gene labels.
            vf0 <- unique(trimws(as.character(Pagwas$VariableFeatures %||% character())))
            vf0 <- vf0[!is.na(vf0) & nzchar(vf0)]
            vf_retry <- unique(c(vf0, snp_labels, strip_suffix(snp_labels)))
            vf_retry <- vf_retry[!is.na(vf_retry) & nzchar(vf_retry)]
            Pagwas$VariableFeatures <- vf_retry

            message(
              "WARN: scPagwas::Pathway_annotation_input failed with 'no match'. ",
              "Retrying after forcing Pagwas$VariableFeatures to include SNP-mapped genes (n=", length(vf_retry), ")."
            )
            try_call()
          }
        )
      }

      unlockBinding("Pathway_annotation_input", ns)
      assign("Pathway_annotation_input", Pathway_annotation_input, envir = ns)
      lockBinding("Pathway_annotation_input", ns)
      invisible(TRUE)
    }

    patched_pa <- tryCatch(ezhu_patch_pathway_annotation_input(), error = function(e) FALSE)
    if (isTRUE(patched_pa)) message("Applied EZHU patch: Pathway_annotation_input (VariableFeatures fallback).")
  }
}

patch_perform_score <- !identical(Sys.getenv("EZHU_SCPAGWAS_PATCH_PERFORM_SCORE", unset = "1"), "0")
if (patch_perform_score && requireNamespace("Matrix", quietly = TRUE)) {
  # scPagwas_perform_score() (v2.0.0) builds `pathway_expr` from a list of per-pathway
  # vectors. When a pathway intersects the data with exactly 1 gene and
  # Pagwas$data_mat is a sparse Matrix, `Pagwas$data_mat[a, ]` returns a 1xN Matrix
  # instead of a length-N vector, which can yield malformed `pathway_expr` and a
  # downstream "non-conformable arrays" error at scoring time.
  #
  # This patch forces the 1-gene case to a numeric vector (without changing values).
  ezhu_patch_scpagwas_perform_score <- function() {
    ns <- asNamespace("scPagwas")
    if (!exists("scPagwas_perform_score", envir = ns, inherits = FALSE)) return(invisible(FALSE))

    scPagwas_perform_score <- function(Pagwas, remove_outlier = TRUE) {
      options(bigmemory.allow.dimnames = TRUE)

      Pathway_sclm_results <- Pagwas$Pathway_sclm_results

      pca_mat <- Pagwas$pca_scCell_mat
      cell_ids <- colnames(pca_mat)
      if (length(cell_ids) == 0 || !all(nzchar(cell_ids))) {
        stop("scPagwas_perform_score: pca_scCell_mat has no valid colnames (cell IDs).", call. = FALSE)
      }

      data_mat <- Pagwas$data_mat
      data_cols <- colnames(data_mat)
      if (!is.null(data_cols) && length(data_cols) > 0 && all(nzchar(data_cols))) {
        if (!all(cell_ids %in% data_cols)) {
          stop("scPagwas_perform_score: cell IDs in pca_scCell_mat are not all present in data_mat columns.", call. = FALSE)
        }
        data_mat <- data_mat[, cell_ids, drop = FALSE]
      } else {
        if (ncol(data_mat) != length(cell_ids)) {
          stop(
            "scPagwas_perform_score: data_mat column count (", ncol(data_mat),
            ") does not match pca_scCell_mat cell count (", length(cell_ids), ").",
            call. = FALSE
          )
        }
      }

      # Normalize Pathway_sclm_results orientation to: (cells x pathways)
      sclm_raw <- Pathway_sclm_results
      if (is.null(sclm_raw) || length(sclm_raw) == 0) {
        singlecell_env <- Sys.getenv("EZHU_SCPAGWAS_SINGLECELL", unset = "")
        hint <- if (identical(singlecell_env, "0")) {
          "\nHint: EZHU_SCPAGWAS_SINGLECELL=0 disables the single-cell layer; scPagwas still reached Stage 8.\nSet EZHU_SCPAGWAS_SINGLECELL=1 (or EZHU_SCPAGWAS_CONTROLS_SINGLECELL=1 in the runner) and rerun."
        } else {
          ""
        }
        stop("scPagwas_perform_score: Pathway_sclm_results is empty.", hint, call. = FALSE)
      }

      sclm_cells_by_pathways <- NULL
      if (!is.null(dim(sclm_raw)) && nrow(sclm_raw) == length(cell_ids)) {
        sclm_cells_by_pathways <- sclm_raw
      } else if (!is.null(dim(sclm_raw)) && ncol(sclm_raw) == length(cell_ids)) {
        # Likely (pathways x cells); transpose.
        if (!is.null(colnames(sclm_raw)) && all(cell_ids %in% colnames(sclm_raw))) {
          sclm_raw <- sclm_raw[, cell_ids, drop = FALSE]
        }
        sclm_cells_by_pathways <- t(sclm_raw)
      } else {
        # Fall back to assuming (cells x pathways), matching upstream scPagwas behavior.
        sclm_cells_by_pathways <- sclm_raw
      }

      Pathway_names <- colnames(sclm_cells_by_pathways)
      if (length(Pathway_names) == 0 || !all(nzchar(Pathway_names))) {
        stop("scPagwas_perform_score: Pathway_sclm_results has no valid pathway names.", call. = FALSE)
      }

      mean_by_cell <- function(m) {
        if (inherits(m, "Matrix")) return(Matrix::colMeans(m))
        if (is.matrix(m)) return(colMeans(m))
        if (is.data.frame(m)) return(colMeans(as.matrix(m)))
        biganalytics::apply(m, 2, mean)
      }

      pathway_list <- Pagwas$Pathway_list
      pathway_list_names <- names(pathway_list)
      name_overlap <- if (!is.null(pathway_list_names) && length(pathway_list_names) > 0) {
        sum(Pathway_names %in% pathway_list_names)
      } else {
        0L
      }

      use_positional_pathway_list <- FALSE
      if (name_overlap == 0L && length(pathway_list) == length(Pathway_names)) {
        use_positional_pathway_list <- TRUE
      }

      data_genes <- rownames(data_mat)
      data_genes_base <- sub("\\.[0-9]+$", "", data_genes)
      has_suffix_mapping <- any(!identical(data_genes, data_genes_base))
      base_to_rows <- if (isTRUE(has_suffix_mapping)) split(data_genes, data_genes_base) else NULL

      get_pathway_genes <- function(i, pa) {
        if (isTRUE(use_positional_pathway_list)) return(pathway_list[[i]])
        pathway_list[[pa]]
      }

      map_genes_to_rows <- function(genes) {
        genes <- unique(as.character(genes %||% character()))
        genes <- genes[!is.na(genes) & nzchar(genes)]
        if (length(genes) == 0) return(character())

        hit <- intersect(genes, data_genes)
        if (length(hit) > 0) return(hit)
        if (!isTRUE(has_suffix_mapping)) return(character())

        mapped <- unlist(base_to_rows[genes], use.names = FALSE)
        mapped <- unique(as.character(mapped))
        mapped <- mapped[!is.na(mapped) & nzchar(mapped)]
        intersect(mapped, data_genes)
      }

      pathway_expr_list <- lapply(seq_along(Pathway_names), function(i) {
        pa <- Pathway_names[[i]]
        genes <- map_genes_to_rows(get_pathway_genes(i, pa))
        if (length(genes) == 0) return(NULL)

        if (length(genes) == 1) {
          one <- data_mat[genes, , drop = FALSE]
          if (inherits(one, "Matrix")) return(as.numeric(Matrix::as.matrix(one)))
          return(as.numeric(one))
        }

        v <- mean_by_cell(data_mat[genes, , drop = FALSE])
        v <- suppressWarnings(as.numeric(v))
        if (length(v) != length(cell_ids)) return(NULL)
        v
      })

      keep <- !vapply(pathway_expr_list, is.null, logical(1))
      pathway_expr_list <- pathway_expr_list[keep]
      Pathway_names <- Pathway_names[keep]

      if (length(Pathway_names) == 0) {
        diag <- paste0(
          "scPagwas_perform_score: no pathways remain after expression filtering.\n",
          "Diagnostics:\n",
          "  pathways_in_sclm=", length(colnames(sclm_cells_by_pathways)), "\n",
          "  pathways_in_list=", length(pathway_list), "\n",
          "  pathway_name_overlap=", name_overlap, "\n",
          "  pathway_list_mode=", if (use_positional_pathway_list) "positional" else "by_name", "\n",
          "  data_mat_genes=", length(data_genes), "\n",
          "  has_suffix_mapping=", if (has_suffix_mapping) "yes" else "no", "\n"
        )
        stop(diag, call. = FALSE)
      }

      pathway_expr <- do.call(cbind, pathway_expr_list)
      rownames(pathway_expr) <- cell_ids
      colnames(pathway_expr) <- Pathway_names

      # Align to the set of pathways that actually exist in the PCA matrix.
      common_pathways <- intersect(Pathway_names, rownames(pca_mat))
      if (length(common_pathways) == 0) {
        stop("scPagwas_perform_score: no pathways overlap between Pathway_sclm_results and pca_scCell_mat.", call. = FALSE)
      }

      pca_sub <- pca_mat[common_pathways, , drop = FALSE]
      # Ensure PCA columns align to cell_ids ordering.
      if (!identical(colnames(pca_sub), cell_ids)) {
        if (all(cell_ids %in% colnames(pca_sub))) {
          pca_sub <- pca_sub[, cell_ids, drop = FALSE]
        } else {
          stop("scPagwas_perform_score: pca_scCell_mat columns do not match expected cell IDs.", call. = FALSE)
        }
      }

      expr_sub <- pathway_expr[, common_pathways, drop = FALSE]
      # Both matrices should be (cells x pathways) after transposing the PCA subset.
      pa_exp_mat <- t(pca_sub) * expr_sub
      rm(pathway_expr, pca_sub, expr_sub, pathway_expr_list)

      # Align Pathway_sclm_results rows to cell order when possible.
      sclm <- sclm_cells_by_pathways
      if (!is.null(rownames(sclm)) && length(rownames(sclm)) > 0 && all(nzchar(rownames(sclm)))) {
        if (all(cell_ids %in% rownames(sclm))) {
          sclm <- sclm[cell_ids, , drop = FALSE]
        } else if (nrow(sclm) != length(cell_ids)) {
          stop("scPagwas_perform_score: Pathway_sclm_results rownames do not match cells and row count differs.", call. = FALSE)
        }
      } else if (nrow(sclm) != length(cell_ids)) {
        stop("scPagwas_perform_score: Pathway_sclm_results row count does not match number of cells.", call. = FALSE)
      }

      sclm <- sclm[, common_pathways, drop = FALSE]
      Pagwas$Pathway_single_results <- data.matrix(sclm) * pa_exp_mat
      rownames(Pagwas$Pathway_single_results) <- colnames(Pagwas$pca_scCell_mat)

      message("* Get Pathways'rankPvalue for each celltypes!")
      cl <- unique(Pagwas$Celltype_anno$annotation)

      Pagwas$Pathway_single_results <- t(data.matrix(Pagwas$Pathway_single_results))

      Pathways_rankPvalue <- lapply(cl, function(ss) {
        print(ss)
        tt <- Pagwas$Celltype_anno$annotation == ss
        PathwayrankPvalue <- scGene_rankP(Pagwas$Pathway_single_results[, tt, drop = FALSE])
        PathwayrankPvalue$pValueHigh
      })

      Pagwas$scPathways_rankPvalue <- Reduce(function(dtf1, dtf2) cbind(dtf1, dtf2), Pathways_rankPvalue)
      Pagwas$scPathways_rankPvalue <- as.data.frame(Pagwas$scPathways_rankPvalue)
      colnames(Pagwas$scPathways_rankPvalue) <- cl
      rownames(Pagwas$scPathways_rankPvalue) <- rownames(Pagwas$Pathway_single_results)
      rm(Pathways_rankPvalue)

      message("* Get scPgwas score for each single cell")
      scPagwas.gPAS.score <- colSums(Pagwas$Pathway_single_results)
      gc()
      if (remove_outlier) {
        Pagwas$scPagwas.gPAS.score <- scPagwas_score_filter(scPagwas_score = scPagwas.gPAS.score)
      }
      names(Pagwas$scPagwas.gPAS.score) <- colnames(Pagwas$pca_scCell_mat)
      Pagwas
    }

    unlockBinding("scPagwas_perform_score", ns)
    assign("scPagwas_perform_score", scPagwas_perform_score, envir = ns)
    lockBinding("scPagwas_perform_score", ns)
    invisible(TRUE)
  }

  patched <- tryCatch(ezhu_patch_scpagwas_perform_score(), error = function(e) FALSE)
  if (isTRUE(patched)) message("Applied EZHU patch: scPagwas_perform_score (sparse 1-gene pathways).")
}

main_name <- if (exists("scPagwas_main2", envir = asNamespace("scPagwas"))) {
  "scPagwas_main2"
} else {
  "scPagwas_main"
}

main_fn <- if (identical(main_name, "scPagwas_main2")) {
  get("scPagwas_main2", envir = asNamespace("scPagwas"))
} else {
  get("scPagwas_main", envir = asNamespace("scPagwas"))
}

message("Running scPagwas (output.dirs=", output_dir, ", prefix=", output_prefix, ")")
pagwas_result <- main_fn(
  Pagwas = NULL,
  gwas_data = gwas_df,
  Single_data = single_obj,
  output.prefix = output_prefix,
  output.dirs = output_dir,
  block_annotation = block_annotation,
  assay = assay,
  Pathway_list = pathway_list,
  chrom_ld = chrom_ld,
  n.cores = n_cores,
  maf_filter = maf_filter,
  iters_celltype = iters_celltype,
  iters_singlecell = iters_singlecell,
  singlecell = run_singlecell,
  celltype = run_celltype,
  seurat_return = TRUE
)

ezhu_add_bh_fdr <- function(path) {
  if (!file.exists(path)) return(invisible(FALSE))

  df <- tryCatch(
    utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || !"pvalue" %in% names(df)) return(invisible(FALSE))

  p <- suppressWarnings(as.numeric(df$pvalue))
  ok <- is.finite(p)
  if (!any(ok)) return(invisible(FALSE))

  df$fdr <- NA_real_
  df$fdr[ok] <- stats::p.adjust(p[ok], method = "BH")

  # Standardize output: avoid an extra unnamed row-index column.
  utils::write.csv(df, file = path, row.names = FALSE, quote = TRUE)
  invisible(TRUE)
}

# Ensure scPagwas celltype-level output reports BH-FDR (across all tested contexts).
celltype_candidates <- c(
  file.path(output_dir, "merged_celltype_pvalue.csv"),
  file.path(output_dir, "Merged_celltype_pvalue.csv")
)
extra <- Sys.glob(file.path(output_dir, "*Merged_celltype_pvalue.csv"))
celltype_candidates <- unique(c(celltype_candidates, extra))
for (path in celltype_candidates) {
  patched <- tryCatch(ezhu_add_bh_fdr(path), error = function(e) FALSE)
  if (isTRUE(patched)) message("Wrote BH-FDR: ", path)
}

# Canonicalize the primary celltype-level output filename expected by downstream scripts.
canonical_pvalues <- file.path(output_dir, "merged_celltype_pvalue.csv")
if (!file.exists(canonical_pvalues)) {
  alt <- c(file.path(output_dir, "Merged_celltype_pvalue.csv"), extra)
  alt <- alt[file.exists(alt)]
  if (length(alt) > 0) {
    ok_copy <- tryCatch(file.copy(alt[1], canonical_pvalues, overwrite = TRUE), error = function(e) FALSE)
    if (isTRUE(ok_copy)) {
      message("Copied celltype pvalues to canonical path: ", canonical_pvalues, " (from ", alt[1], ")")
    }
  }
}

out_rds <- Sys.getenv(
  "EZHU_SCPAGWAS_RDS_OUT",
  unset = file.path("data/processed", paste0("singlecell_", singlecell_tag, "_scpagwas_", output_prefix, ".rds"))
)
saveRDS(pagwas_result, file = out_rds)

meta <- pagwas_result[[]]
meta$cell <- rownames(meta)
meta$cluster <- as.character(meta$seurat_clusters %||% Seurat::Idents(pagwas_result))
meta$group <- as.character(meta$group %||% NA_character_)
meta$sample <- as.character(meta$sample %||% NA_character_)

trs_cols <- grep("^scPagwas\\.TRS\\.Score[0-9]+$", colnames(meta), value = TRUE)
if (length(trs_cols) > 0) {
  trs_matrix <- as.matrix(meta[, trs_cols, drop = FALSE])
  meta$scPagwas.TRS.Mean <- rowMeans(trs_matrix, na.rm = TRUE)
  meta$scPagwas.TRS.SD <- apply(trs_matrix, 1, stats::sd, na.rm = TRUE)
}

score_cols <- intersect(
  c(
    "scPagwas.TRS.Score1",
    "scPagwas.TRS.Mean",
    "scPagwas.TRS.SD",
    "scPagwas.gPAS.score",
    "Random_Correct_BG_p",
    "Random_Correct_BG_adjp",
    "Random_Correct_BG_z"
  ),
  colnames(meta)
)

cell_scores <- meta[, c("cell", "sample", "group", "cluster", score_cols), drop = FALSE]
utils::write.table(
  cell_scores,
  file = file.path(output_dir, "cell_scores.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cluster_summary <- aggregate(
  x = list(n_cells = rep(1L, nrow(meta))),
  by = list(cluster = meta$cluster, group = meta$group),
  FUN = sum
)
cluster_summary <- cluster_summary[order(cluster_summary$cluster, cluster_summary$group), , drop = FALSE]
utils::write.table(
  cluster_summary,
  file = file.path(output_dir, "cluster_by_group_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ezhu_write_stage_metadata(
  "03_analysis",
  params = list(
    singlecell_tag = singlecell_tag,
    gwas_harmonized = gwas_path,
    gwas_manifest_id = as.character(row$id),
    build = build,
    singlecell_processed = singlecell_path,
    exclude_samples = if (length(exclude_samples) > 0) exclude_samples else NULL,
    assay = assay,
    idents_field = idents_field,
    scpagwas = list(
      prefix = output_prefix,
      output_dir = output_dir,
      maf_filter = maf_filter,
      iters_singlecell = iters_singlecell,
      iters_celltype = iters_celltype,
      singlecell = run_singlecell,
      celltype = run_celltype,
      n_cores = n_cores,
      main_function = main_name
    ),
    outputs = list(
      scpagwas_seurat_rds = out_rds,
      cell_scores_tsv = file.path(output_dir, "cell_scores.tsv"),
      cluster_by_group_counts_tsv = file.path(output_dir, "cluster_by_group_counts.tsv")
    )
  ),
  seed = seed
)
