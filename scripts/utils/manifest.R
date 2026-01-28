ezhu_read_manifest <- function(path = "data/manifest.tsv") {
  if (!file.exists(path)) {
    stop("Manifest not found: ", path, call. = FALSE)
  }

  df <- utils::read.delim(
    file = path,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    quote = "",
    comment.char = ""
  )

  names(df) <- tolower(gsub("[^a-zA-Z0-9]+", "_", names(df)))
  df
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

ezhu_required_manifest_columns <- function() {
  c(
    "id",
    "type",
    "name",
    "source",
    "accession",
    "url",
    "version_date",
    "download_method",
    "local_path",
    "integrity",
    "processing_scripts",
    "notes"
  )
}

ezhu_is_safe_relative_path <- function(path) {
  if (is.na(path) || !nzchar(path)) return(FALSE)
  if (grepl("^/", path)) return(FALSE)
  if (grepl("^[A-Za-z]:[\\\\/]", path)) return(FALSE)
  if (grepl("\\.\\.", path)) return(FALSE)
  TRUE
}

ezhu_parse_integrity <- function(integrity) {
  if (is.na(integrity) || !nzchar(trimws(integrity))) return(list())

  parts <- unlist(strsplit(integrity, "[;[:space:]]+"))
  parts <- parts[nzchar(parts)]

  out <- list()
  for (p in parts) {
    p <- trimws(p)
    if (!nzchar(p)) next

    if (grepl("=", p, fixed = TRUE)) {
      kv <- strsplit(p, "=", fixed = TRUE)[[1]]
      if (length(kv) < 2) next
      key <- kv[1]
      value <- paste(kv[-1], collapse = "=")
    } else if (grepl(":", p, fixed = TRUE)) {
      kv <- strsplit(p, ":", fixed = TRUE)[[1]]
      if (length(kv) < 2) next
      key <- kv[1]
      value <- paste(kv[-1], collapse = ":")
    } else {
      next
    }

    key <- tolower(trimws(key))
    value <- trimws(value)
    if (!nzchar(key) || !nzchar(value)) next
    out[[key]] <- value
  }

  out
}

ezhu_parse_notes_kv <- function(notes) {
  if (is.na(notes) || !nzchar(trimws(notes))) return(list())
  parts <- unlist(strsplit(notes, "[;\r\n]+"))
  parts <- trimws(parts)
  parts <- parts[nzchar(parts)]

  out <- list()
  for (p in parts) {
    p <- trimws(p)
    if (!nzchar(p)) next

    if (grepl("=", p, fixed = TRUE)) {
      kv <- strsplit(p, "=", fixed = TRUE)[[1]]
      if (length(kv) < 2) next
      key <- kv[1]
      value <- paste(kv[-1], collapse = "=")
    } else if (grepl(":", p, fixed = TRUE)) {
      kv <- strsplit(p, ":", fixed = TRUE)[[1]]
      if (length(kv) < 2) next
      key <- kv[1]
      value <- paste(kv[-1], collapse = ":")
    } else {
      next
    }

    key <- tolower(trimws(key))
    value <- trimws(value)
    if (!nzchar(key) || !nzchar(value)) next
    out[[key]] <- value
  }
  out
}

ezhu_checksum_sha256 <- function(path) {
  if (requireNamespace("digest", quietly = TRUE)) {
    return(digest::digest(file = path, algo = "sha256", serialize = FALSE))
  }

  if (nzchar(Sys.which("sha256sum"))) {
    out <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
    return(strsplit(out[1], "\\s+")[[1]][1])
  }
  if (nzchar(Sys.which("shasum"))) {
    out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE, stderr = TRUE)
    return(strsplit(out[1], "\\s+")[[1]][1])
  }
  if (nzchar(Sys.which("openssl"))) {
    out <- system2("openssl", c("dgst", "-sha256", shQuote(path)), stdout = TRUE, stderr = TRUE)
    return(sub("^.*=\\s*", "", out[1]))
  }

  stop("No SHA256 implementation available (digest/sha256sum/shasum/openssl).", call. = FALSE)
}

ezhu_checksum_md5 <- function(path) {
  as.character(tools::md5sum(path))
}

ezhu_verify_integrity <- function(path, integrity_spec) {
  if (length(integrity_spec) == 0) return(list(ok = TRUE, checks = list()))

  checks <- list()
  ok <- TRUE

  if (!is.null(integrity_spec$size)) {
    expected <- suppressWarnings(as.numeric(integrity_spec$size))
    actual <- file.info(path)$size
    pass <- is.finite(expected) && !is.na(actual) && actual == expected
    checks$size <- list(expected = expected, actual = actual, pass = pass)
    ok <- ok && pass
  }

  if (!is.null(integrity_spec$md5)) {
    expected <- tolower(integrity_spec$md5)
    actual <- tolower(ezhu_checksum_md5(path))
    pass <- identical(actual, expected)
    checks$md5 <- list(expected = expected, actual = actual, pass = pass)
    ok <- ok && pass
  }

  if (!is.null(integrity_spec$sha256)) {
    expected <- tolower(integrity_spec$sha256)
    actual <- tolower(ezhu_checksum_sha256(path))
    pass <- identical(actual, expected)
    checks$sha256 <- list(expected = expected, actual = actual, pass = pass)
    ok <- ok && pass
  }

  list(ok = ok, checks = checks)
}
