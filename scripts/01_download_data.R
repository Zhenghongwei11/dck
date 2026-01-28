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

options(timeout = suppressWarnings(as.integer(Sys.getenv("EZHU_DOWNLOAD_TIMEOUT", unset = "3600"))))

seed_env <- Sys.getenv("EZHU_SEED", unset = "")
seed <- suppressWarnings(as.integer(seed_env))
if (is.na(seed)) seed <- 20251227L
set.seed(seed)

parse_csv_env <- function(env_value) {
  value <- trimws(as.character(env_value %||% ""))
  if (!nzchar(value)) return(character())
  parts <- unlist(strsplit(value, "[,;]"))
  parts <- trimws(parts)
  parts[nzchar(parts)]
}

download_ids <- parse_csv_env(Sys.getenv("EZHU_DOWNLOAD_IDS", unset = ""))
exclude_ids <- parse_csv_env(Sys.getenv("EZHU_DOWNLOAD_EXCLUDE_IDS", unset = ""))

manifest <- ezhu_read_manifest("data/manifest.tsv")
required <- ezhu_required_manifest_columns()
missing <- setdiff(required, names(manifest))
if (length(missing) > 0) {
  stop("Manifest is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

manifest$download_method <- tolower(trimws(manifest$download_method))

download_results <- data.frame(
  id = manifest$id,
  local_path = manifest$local_path,
  url = manifest$url,
  download_method = manifest$download_method,
  status = character(nrow(manifest)),
  message = character(nrow(manifest)),
  stringsAsFactors = FALSE
)

download_with_curl <- function(url, dest) {
  curl <- Sys.which("curl")
  if (!nzchar(curl)) return(FALSE)

  args <- c(
    "-L",
    "--fail",
    "--retry", "5",
    "--retry-delay", "5",
    "--connect-timeout", "30",
    "--continue-at", "-",
    "-o", dest,
    url
  )
  status <- suppressWarnings(system2(curl, args))
  isTRUE(identical(status, 0L))
}

for (i in seq_len(nrow(manifest))) {
  id <- trimws(as.character(manifest$id[i] %||% ""))
  local_path <- manifest$local_path[i]
  url <- manifest$url[i]
  method <- manifest$download_method[i]

  if (length(download_ids) > 0 && !id %in% download_ids) {
    download_results$status[i] <- "skipped"
    download_results$message[i] <- "not selected (EZHU_DOWNLOAD_IDS)"
    next
  }

  if (length(exclude_ids) > 0 && id %in% exclude_ids) {
    download_results$status[i] <- "skipped"
    download_results$message[i] <- "excluded (EZHU_DOWNLOAD_EXCLUDE_IDS)"
    next
  }

  if (is.na(local_path) || !nzchar(trimws(local_path))) {
    download_results$status[i] <- "skipped"
    download_results$message[i] <- "missing local_path"
    next
  }
  if (!ezhu_is_safe_relative_path(local_path)) {
    download_results$status[i] <- "failed"
    download_results$message[i] <- "unsafe local_path (must be relative, no '..')"
    next
  }

  if (!identical(method, "script")) {
    download_results$status[i] <- "skipped"
    download_results$message[i] <- paste0("download_method=", method %||% "NA")
    next
  }
  if (is.na(url) || !nzchar(trimws(url))) {
    download_results$status[i] <- "failed"
    download_results$message[i] <- "missing url for script download"
    next
  }

  ezhu_dir_create(dirname(local_path))

  integrity <- ezhu_parse_integrity(manifest$integrity[i])

  if (file.exists(local_path) && length(integrity) > 0) {
    verify <- ezhu_verify_integrity(local_path, integrity)
    if (isTRUE(verify$ok)) {
      download_results$status[i] <- "skipped"
      download_results$message[i] <- "already exists (integrity ok)"
      next
    }
    # Corrupt/partial file: keep it as a resumable .part when possible.
    part <- paste0(local_path, ".part")
    if (!file.exists(part)) {
      ok_move <- tryCatch(file.rename(local_path, part), warning = function(w) FALSE, error = function(e) FALSE)
      if (!isTRUE(ok_move)) {
        unlink(local_path, force = TRUE)
      }
    }
  } else if (file.exists(local_path)) {
    download_results$status[i] <- "skipped"
    download_results$message[i] <- "already exists"
    next
  }

  ok <- TRUE
  msg <- "downloaded"

  tryCatch(
    {
      part <- paste0(local_path, ".part")
      dest <- if (file.exists(part)) part else local_path

      downloaded <- download_with_curl(url, dest)
      if (!isTRUE(downloaded)) {
        utils::download.file(url = url, destfile = dest, mode = "wb", quiet = FALSE, method = "libcurl")
      }

      # Move completed download into place (atomic-ish).
      if (file.exists(part)) {
        file.rename(part, local_path)
      }
    },
    error = function(e) {
      ok <<- FALSE
      msg <<- conditionMessage(e)
    }
  )

  if (ok) {
    if (length(integrity) > 0) {
      verify <- ezhu_verify_integrity(local_path, integrity)
      if (!isTRUE(verify$ok)) {
        ok <- FALSE
        msg <- "integrity check failed"
      }
    }
  }

  download_results$status[i] <- if (ok) "ok" else "failed"
  download_results$message[i] <- msg
}

ezhu_dir_create("results/metadata")
out_path <- file.path("results/metadata", "01_download_data__log.tsv")
utils::write.table(download_results, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
ezhu_write_stage_metadata(
  "01_download_data",
  params = list(
    manifest_path = "data/manifest.tsv",
    total_rows = nrow(manifest),
    ok = sum(download_results$status == "ok", na.rm = TRUE),
    skipped = sum(download_results$status == "skipped", na.rm = TRUE),
    failed = sum(download_results$status == "failed", na.rm = TRUE),
    download_ids = download_ids,
    exclude_ids = exclude_ids,
    log_tsv = out_path
  ),
  seed = seed
)
