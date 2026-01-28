audit_check_figures <- function(level = "full") {
  if (!dir.exists("plots/publication")) {
    audit_warn("plots/publication/ not found; skipping figure checks.")
    return(invisible(TRUE))
  }

  files <- list.files("plots/publication", full.names = TRUE)
  files <- files[!grepl("\\.gitkeep$", files)]
  files <- files[!dir.exists(files)]
  files <- files[basename(files) != ".DS_Store"]
  if (length(files) == 0) return(invisible(TRUE))

  write_inventory <- function(paths) {
    out_dir <- file.path("docs", "audit_runs", "current")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    out_path <- file.path(out_dir, "figures_inventory.tsv")

    sha_cmd <- ""
    if (nzchar(Sys.which("shasum"))) sha_cmd <- "shasum"
    if (!nzchar(sha_cmd) && nzchar(Sys.which("sha256sum"))) sha_cmd <- "sha256sum"

    pdfinfo_cmd <- Sys.which("pdfinfo")

    rows <- lapply(paths, function(path) {
      info <- suppressWarnings(file.info(path))
      bytes <- if (!is.null(info$size) && is.finite(info$size)) as.integer(info$size) else NA_integer_
      mtime_epoch <- NA_integer_
      if (!is.null(info$mtime) && !is.na(info$mtime)) {
        mtime_epoch <- as.integer(as.POSIXct(info$mtime, tz = "UTC"))
      }

      sha256 <- ""
      if (nzchar(sha_cmd)) {
        args <- if (identical(sha_cmd, "shasum")) c("-a", "256", path) else c(path)
        out <- tryCatch(system2(sha_cmd, args, stdout = TRUE, stderr = TRUE), error = function(e) character())
        if (length(out) > 0) sha256 <- sub("\\s+.*$", "", out[1])
      }

      pages <- ""
      ext <- tolower(tools::file_ext(path))
      if (identical(ext, "pdf") && nzchar(pdfinfo_cmd)) {
        pi <- tryCatch(system2(pdfinfo_cmd, path, stdout = TRUE, stderr = TRUE), error = function(e) character())
        pages_line <- pi[grepl("^Pages:", pi)]
        if (length(pages_line) > 0) pages <- sub("^Pages:\\s*", "", pages_line[1])
      }

      c(
        path = path,
        bytes = as.character(bytes),
        sha256 = sha256,
        mtime_epoch = as.character(mtime_epoch),
        pages = pages
      )
    })

    header <- c("path", "bytes", "sha256", "mtime_epoch", "pages")
    lines <- c(paste(header, collapse = "\t"))
    for (r in rows) {
      lines <- c(lines, paste(r[header], collapse = "\t"))
    }
    writeLines(lines, out_path, useBytes = TRUE)
  }

  allowed_ext <- c("pdf", "tif", "tiff", "png")
  for (path in files) {
    ext <- tolower(tools::file_ext(path))
    if (!ext %in% allowed_ext) {
      audit_fail("Unexpected figure file type in plots/publication/: ", path)
    }
    if (identical(ext, "png")) {
      audit_warn("PNG found in plots/publication/ (consider PDF/TIFF for publication): ", path)
    }

    # Heuristic guardrail: avoid committing placeholder figures created on errors.
    # Placeholder PDFs in this repo are typically very small (a few KB) and contain
    # obvious failure strings in their text layer.
    if (identical(level, "full") && identical(ext, "pdf")) {
      info <- suppressWarnings(file.info(path))
      if (!is.null(info$size) && is.finite(info$size) && info$size > 0 && info$size < 8000) {
        audit_warn("Suspiciously small PDF (possible placeholder; inspect if this is expected): ", path)
      }
      if (nzchar(Sys.which("pdftotext"))) {
        txt <- tryCatch(system2("pdftotext", c("-q", path, "-"), stdout = TRUE, stderr = TRUE), error = function(e) character())
        txt <- paste(txt, collapse = "\n")
        if (grepl("failed", txt, ignore.case = TRUE) || grepl("missing", txt, ignore.case = TRUE)) {
          audit_fail("PDF appears to contain a placeholder error message: ", path)
        }
      }
    }
  }

  if (identical(level, "full")) {
    write_inventory(files)
  }

  if (identical(level, "ci")) {
    audit_info("CI mode: skipping figure index linkage checks.")
    return(invisible(TRUE))
  }

  index_path <- "docs/FIGURE_INDEX.md"
  index_lines <- audit_read_text(index_path)
  if (length(index_lines) == 0) {
    audit_fail("Figure index missing or empty: ", index_path)
  }

  for (path in files) {
    rel <- normalizePath(path, winslash = "/", mustWork = FALSE)
    rel <- sub(paste0("^", normalizePath(getwd(), winslash = "/", mustWork = FALSE), "/"), "", rel)
    if (!any(grepl(rel, index_lines, fixed = TRUE))) {
      audit_fail("Figure not referenced in docs/FIGURE_INDEX.md: ", rel)
    }
  }

  invisible(TRUE)
}
