ezhu_set_repo_root <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  start_dir <- getwd()
  if (length(file_arg) == 1) {
    script_path <- sub("^--file=", "", file_arg)
    script_path <- tryCatch(normalizePath(script_path, winslash = "/", mustWork = TRUE), error = function(e) "")
    if (nzchar(script_path)) start_dir <- dirname(script_path)
  }

  cur <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)
  for (i in seq_len(20)) {
    if (file.exists(file.path(cur, "renv.lock")) || dir.exists(file.path(cur, ".git"))) {
      setwd(cur)
      return(invisible(cur))
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }

  invisible(getwd())
}

ezhu_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

ezhu_activate_renv <- function(project = getwd()) {
  if (identical(Sys.getenv("EZHU_DISABLE_RENV", unset = ""), "1")) {
    return(invisible(.libPaths()))
  }
  if (requireNamespace("renv", quietly = TRUE)) {
    tryCatch(renv::activate(project = project), error = function(e) NULL)
  }
  invisible(.libPaths())
}

ezhu_git_sha <- function() {
  if (!dir.exists(".git")) return(NA_character_)
  out <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = TRUE), error = function(e) character())
  sha <- trimws(out[1] %||% "")
  if (!nzchar(sha) || grepl("fatal:", sha, fixed = TRUE)) return(NA_character_)
  sha
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

ezhu_write_text <- function(path, lines) {
  ezhu_dir_create(dirname(path))
  writeLines(text = lines, con = path, useBytes = TRUE)
  invisible(path)
}

ezhu_write_json <- function(path, object) {
  ezhu_dir_create(dirname(path))
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    json <- jsonlite::toJSON(object, pretty = TRUE, auto_unbox = TRUE, null = "null")
    writeLines(text = json, con = path, useBytes = TRUE)
    return(invisible(path))
  }

  dput_path <- sub("\\.json$", ".dput", path)
  dput(object, file = dput_path)
  invisible(dput_path)
}

ezhu_write_stage_metadata <- function(stage, params = list(), seed = NA_integer_) {
  ezhu_dir_create("results/metadata")

  metadata <- list(
    stage = stage,
    timestamp_utc = format(as.POSIXct(Sys.time(), tz = "UTC"), "%Y-%m-%dT%H:%M:%SZ"),
    timestamp_local = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    seed = seed,
    git_sha = ezhu_git_sha(),
    r_version = R.version.string,
    platform = R.version$platform,
    libpaths = .libPaths(),
    locale = tryCatch(Sys.getlocale(), error = function(e) NA_character_),
    working_directory = normalizePath(getwd(), winslash = "/", mustWork = FALSE),
    args = commandArgs(trailingOnly = FALSE),
    params = params
  )

  ezhu_write_json(file.path("results/metadata", paste0(stage, "__metadata.json")), metadata)
  ezhu_write_text(
    file.path("results/metadata", paste0(stage, "__sessionInfo.txt")),
    capture.output(sessionInfo())
  )
  ezhu_write_text(file.path("results/metadata", paste0(stage, "__done.txt")), c(metadata$timestamp_local, metadata$git_sha))

  invisible(metadata)
}
