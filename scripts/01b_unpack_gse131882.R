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

tar_path <- Sys.getenv("EZHU_GSE131882_TAR", unset = "data/raw/singlecell/GSE131882/GSE131882_RAW.tar")
out_dir <- Sys.getenv("EZHU_GSE131882_FILES_DIR", unset = "data/raw/singlecell/GSE131882/files")
index_tsv <- Sys.getenv("EZHU_GSE131882_INDEX_TSV", unset = "results/singlecell/gse131882_sample_index.tsv")

ezhu_dir_create(dirname(tar_path))
ezhu_dir_create(out_dir)
ezhu_dir_create("results/metadata")

if (!file.exists(tar_path)) {
  stop(
    "Missing GSE131882 tar archive: ", tar_path, "\n",
    "Run: make data (or download it per data/manifest.tsv), then rerun this script.",
    call. = FALSE
  )
}
if (!file.exists(index_tsv)) {
  stop(
    "Missing expected index TSV (used to map/verify extracted files): ", index_tsv, "\n",
    "This repo expects the 6 curated sample RDS files documented in that table.",
    call. = FALSE
  )
}

idx <- utils::read.delim(index_tsv, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
need <- c("sample", "path", "md5")
if (!all(need %in% names(idx))) {
  stop("Index TSV missing required columns: ", paste(setdiff(need, names(idx)), collapse = ", "), call. = FALSE)
}
idx$sample <- as.character(idx$sample)
idx$path <- as.character(idx$path)
idx$md5 <- tolower(trimws(as.character(idx$md5)))
idx <- idx[nzchar(idx$sample) & nzchar(idx$path) & nzchar(idx$md5), , drop = FALSE]
if (nrow(idx) == 0) stop("Index TSV contains no valid rows: ", index_tsv, call. = FALSE)

expected_rel <- idx$path
expected_files <- basename(expected_rel)
expected_md5 <- setNames(idx$md5, expected_files)

already_ok <- TRUE
for (f in expected_files) {
  dest <- file.path(out_dir, f)
  if (!file.exists(dest)) {
    already_ok <- FALSE
    break
  }
  md5 <- tolower(unname(tools::md5sum(dest)))
  if (!identical(md5, expected_md5[[f]])) {
    already_ok <- FALSE
    break
  }
}

if (!already_ok) {
  message("Unpacking GSE131882 raw tar: ", tar_path)
  tmp_dir <- tempfile("gse131882_unpack_")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  utils::untar(tar_path, exdir = tmp_dir)
  extracted <- list.files(tmp_dir, recursive = TRUE, full.names = TRUE)
  extracted <- extracted[!dir.exists(extracted)]
  extracted <- extracted[grepl("\\.rds(\\.gz)?$", extracted, ignore.case = TRUE)]
  if (length(extracted) == 0) {
    stop("No .rds/.rds.gz files found inside tar: ", tar_path, call. = FALSE)
  }

  expected_md5_values <- unique(as.character(expected_md5))

  decompress_gz_to_file <- function(path, out_path) {
    in_con <- gzfile(path, open = "rb")
    out_con <- file(out_path, open = "wb")
    on.exit({
      try(close(in_con), silent = TRUE)
      try(close(out_con), silent = TRUE)
    }, add = TRUE)

    repeat {
      chunk <- readBin(in_con, what = "raw", n = 1024 * 1024)
      if (length(chunk) == 0) break
      writeBin(chunk, out_con)
    }

    out_path
  }

  md5_hits <- list()
  for (path in extracted) {
    is_gz <- grepl("\\.gz$", path, ignore.case = TRUE)

    if (is_gz) {
      # Index TSV records MD5 for the *decompressed* .rds payload.
      tmp_rds <- tempfile("gse131882_rds_", tmpdir = tmp_dir, fileext = ".rds")
      decompress_gz_to_file(path, tmp_rds)
      md5 <- tolower(unname(tools::md5sum(tmp_rds)))
      if (md5 %in% expected_md5_values) {
        md5_hits[[md5]] <- unique(c(md5_hits[[md5]] %||% character(), tmp_rds))
      } else {
        unlink(tmp_rds)
      }
    } else {
      md5 <- tolower(unname(tools::md5sum(path)))
      if (md5 %in% expected_md5_values) {
        md5_hits[[md5]] <- unique(c(md5_hits[[md5]] %||% character(), path))
      }
    }
  }

  matched <- 0L
  for (f in expected_files) {
    target_md5 <- expected_md5[[f]]
    hits <- md5_hits[[target_md5]] %||% character()
    if (length(hits) == 0) {
      stop("Could not find expected file by md5 inside tar for: ", f, call. = FALSE)
    }
    if (length(hits) > 1) {
      warning("Multiple files match expected md5 for ", f, "; taking first match.", call. = FALSE)
    }
    src <- hits[1]
    dest <- file.path(out_dir, f)
    
    is_src_gz <- grepl("\\.gz$", src, ignore.case = TRUE)
    if (is_src_gz) {
      message("Decompressing matched file for ", f, "...")
      system(paste("gunzip -c", shQuote(src), ">", shQuote(dest)))
    } else {
      file.copy(src, dest, overwrite = TRUE)
    }
    matched <- matched + 1L
  }

  message("Unpacked and matched ", matched, " curated sample files into: ", out_dir)
}

# Final verification (hard fail if mismatch).
for (f in expected_files) {
  dest <- file.path(out_dir, f)
  if (!file.exists(dest)) stop("Missing expected sample file after unpack: ", dest, call. = FALSE)
  md5 <- tolower(unname(tools::md5sum(dest)))
  if (!identical(md5, expected_md5[[f]])) {
    stop("MD5 mismatch after unpack for ", dest, call. = FALSE)
  }
}

ezhu_write_stage_metadata(
  "01b_unpack_gse131882",
  params = list(
    tar_path = tar_path,
    out_dir = out_dir,
    index_tsv = index_tsv,
    n_samples = length(expected_files),
    outputs = list(sample_files_dir = out_dir)
  ),
  seed = seed
)

message("GSE131882 unpack OK.")
