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

options(stringsAsFactors = FALSE, width = 120)

seed_env <- Sys.getenv("EZHU_SEED", unset = "")
seed <- suppressWarnings(as.integer(seed_env))
if (is.na(seed)) seed <- 20251227L
set.seed(seed)

source("scripts/utils/manifest.R")

manifest_path <- Sys.getenv("EZHU_MANIFEST", unset = "data/manifest.tsv")
if (!file.exists(manifest_path)) stop("Missing manifest: ", manifest_path, call. = FALSE)

mf <- ezhu_read_manifest(manifest_path)
mf$id <- as.character(mf$id)

row <- mf[mf$id == "GSE131882", , drop = FALSE]
if (nrow(row) != 1) stop("Expected exactly one manifest row for id=GSE131882 in: ", manifest_path, call. = FALSE)

src_url <- trimws(as.character(row$url))
if (!nzchar(src_url)) stop("Manifest row GSE131882 missing url.", call. = FALSE)

raw_tar <- Sys.getenv("EZHU_GSE131882_TAR", unset = trimws(as.character(row$local_path)))
if (!nzchar(raw_tar)) raw_tar <- "data/raw/singlecell/GSE131882/GSE131882_RAW.tar"

raw_dir <- dirname(raw_tar)
ezhu_dir_create(raw_dir)

curl_bin <- Sys.getenv("CURL", unset = Sys.which("curl"))
if (!nzchar(curl_bin)) stop("curl not found on PATH.", call. = FALSE)

download_if_needed <- function(url, path) {
  if (file.exists(path) && file.info(path)$size > 0) return(invisible(path))
  ezhu_dir_create(dirname(path))
  cmd <- paste(
    shQuote(curl_bin),
    "-L",
    "--silent",
    "--show-error",
    "--fail",
    "--output", shQuote(path),
    shQuote(url)
  )
  status <- system(cmd)
  if (!identical(status, 0L) || !file.exists(path) || file.info(path)$size <= 0) {
    stop("Failed to download: ", url, call. = FALSE)
  }
  invisible(path)
}

download_if_needed(src_url, raw_tar)

extract_dir <- Sys.getenv("EZHU_GSE131882_EXTRACT_DIR", unset = file.path(raw_dir, "extracted_universe"))
extract_dir <- trimws(extract_dir)
if (!nzchar(extract_dir)) extract_dir <- file.path(raw_dir, "extracted_universe")
ezhu_dir_create(extract_dir)

safe_untar <- function(tar_path, exdir) {
  # If already extracted, skip.
  existing <- list.files(exdir, pattern = "\\.rds(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  if (length(existing) > 0) return(invisible(existing))

  utils::untar(tar_path, exdir = exdir)
  out <- list.files(exdir, pattern = "\\.rds(\\.gz)?$", recursive = TRUE, full.names = TRUE)
  if (length(out) == 0) stop("No .rds.gz files found after untar: ", tar_path, call. = FALSE)
  out
}

files <- safe_untar(raw_tar, extract_dir)
files <- sort(files)

find_matrix_like <- function(x, max_nodes = 2000L) {
  is_mat <- function(v) {
    if (is.matrix(v)) return(TRUE)
    if (inherits(v, "dgCMatrix")) return(TRUE)
    if (inherits(v, "Matrix")) return(TRUE)
    FALSE
  }

  # Common structure used in this dataset.
  if (is.list(x) && !is.null(x$umicount) && is.list(x$umicount) &&
      !is.null(x$umicount$exon) && is.list(x$umicount$exon) &&
      !is.null(x$umicount$exon$all) && is_mat(x$umicount$exon$all)) {
    return(x$umicount$exon$all)
  }

  # Breadth-first search through lists.
  queue <- list(x)
  nodes <- 0L
  while (length(queue) > 0 && nodes < max_nodes) {
    cur <- queue[[1]]
    queue <- queue[-1]
    nodes <- nodes + 1L

    if (is_mat(cur)) return(cur)
    if (!is.list(cur)) next

    kids <- cur
    if (length(kids) == 0) next
    queue <- c(queue, kids)
  }

  NULL
}

gene_universe <- character()
for (path in files) {
  obj <- NULL
  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    con <- gzfile(path, open = "rb")
    obj <- tryCatch(readRDS(con), error = function(e) NULL, finally = try(close(con), silent = TRUE))
  } else {
    obj <- tryCatch(readRDS(path), error = function(e) NULL)
  }
  if (is.null(obj)) next

  mat <- find_matrix_like(obj)
  rm(obj)
  if (is.null(mat)) {
    warning("Could not locate matrix-like object in: ", path, call. = FALSE)
    next
  }

  genes <- rownames(mat)
  rm(mat)
  if (!is.null(genes) && length(genes) > 0) {
    genes <- toupper(trimws(as.character(genes)))
    genes <- genes[nzchar(genes) & !is.na(genes)]
    gene_universe <- unique(c(gene_universe, genes))
  }

  gc(verbose = FALSE)
}

if (length(gene_universe) == 0) stop("Failed to extract any genes from GSE131882 raw objects.", call. = FALSE)
gene_universe <- sort(unique(gene_universe))

out_path <- Sys.getenv(
  "EZHU_GSE131882_GENE_UNIVERSE_OUT",
  unset = "data/references/singlecell/gse131882_gene_universe.tsv"
)
out_path <- trimws(out_path)
if (!nzchar(out_path)) stop("Invalid output path.", call. = FALSE)
ezhu_dir_create(dirname(out_path))

is_ensg <- mean(grepl("^ENSG[0-9]+$", gene_universe)) >= 0.90

out <- NULL
if (is_ensg) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop(
      "jsonlite is required to map ENSG IDs to gene symbols via Ensembl REST.\n",
      "On Ubuntu: `sudo apt-get install -y r-cran-jsonlite`.\n",
      "Or install in R: install.packages(\"jsonlite\").",
      call. = FALSE
    )
  }
  if (!nzchar(curl_bin)) stop("curl not found on PATH.", call. = FALSE)

  lookup_chunk <- function(ids) {
    body <- jsonlite::toJSON(list(ids = ids), auto_unbox = TRUE, null = "null")
    cmd <- paste(
      shQuote(curl_bin),
      "-sS",
      "--fail",
      "-X", "POST",
      "-H", shQuote("Content-Type: application/json"),
      "-H", shQuote("Accept: application/json"),
      "-d", shQuote(body),
      shQuote("https://rest.ensembl.org/lookup/id?content-type=application/json")
    )
    txt <- paste(system(cmd, intern = TRUE), collapse = "\n")
    if (!nzchar(txt)) return(list())
    jsonlite::fromJSON(txt, simplifyVector = FALSE)
  }

  chunk_size <- suppressWarnings(as.integer(Sys.getenv("EZHU_ENSEMBL_LOOKUP_CHUNK", unset = "500")))
  if (is.na(chunk_size) || chunk_size <= 0) chunk_size <- 500L

  ids <- gene_universe
  res <- vector("list", length(ids))
  names(res) <- ids

  idx <- seq_along(ids)
  chunks <- split(idx, ceiling(idx / chunk_size))
  for (c in chunks) {
    sub_ids <- ids[c]
    m <- lookup_chunk(sub_ids)
    if (!is.list(m) || is.null(names(m))) next
    for (k in names(m)) {
      rec <- m[[k]]
      sym <- ""
      if (is.list(rec) && !is.null(rec$display_name)) sym <- as.character(rec$display_name)
      sym <- toupper(trimws(sym))
      res[[k]] <- sym
    }
    Sys.sleep(0.05)
  }

  sym_vec <- unlist(res, use.names = TRUE)
  ensg_vec <- names(sym_vec)
  out <- data.frame(ensg = ensg_vec, gene_symbol = as.character(sym_vec), stringsAsFactors = FALSE)
  out <- out[nzchar(out$ensg), , drop = FALSE]
  out$gene_symbol <- toupper(trimws(out$gene_symbol))
  out <- out[nzchar(out$gene_symbol), , drop = FALSE]
  out <- out[!duplicated(out$ensg), , drop = FALSE]
  out <- out[order(out$gene_symbol, out$ensg), , drop = FALSE]
} else {
  out <- data.frame(gene_symbol = gene_universe, stringsAsFactors = FALSE)
}

utils::write.table(out, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)

ezhu_write_stage_metadata(
  "02b_export_gene_universe_gse131882",
  params = list(
    manifest = manifest_path,
    url = src_url,
    raw_tar = raw_tar,
    extract_dir = extract_dir,
    out_path = out_path,
    universe_rows = nrow(out),
    universe_is_ensg = is_ensg,
    n_unique_input_ids = length(gene_universe),
    files_n = length(files)
  ),
  seed = seed
)

message("Wrote: ", out_path, " (n_rows=", nrow(out), ")")
