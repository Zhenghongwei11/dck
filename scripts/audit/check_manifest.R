audit_check_manifest <- function(manifest) {
  required <- ezhu_required_manifest_columns()
  missing <- setdiff(required, names(manifest))
  audit_assert(length(missing) == 0, "Manifest missing required columns: ", paste(missing, collapse = ", "))

  active <- audit_active_manifest_rows(manifest)
  if (!any(active)) {
    audit_warn("No active rows in manifest (data/manifest.tsv).")
    return(invisible(TRUE))
  }

  # The manifest tracks both on-disk artifacts and external sources (e.g., APIs).
  # Keep this list aligned with real values used in data/manifest.tsv.
  allowed_methods <- c("script", "manual", "gated", "api", "curated")

  for (i in which(active)) {
    row <- manifest[i, , drop = FALSE]

    id <- as.character(row$id)
    row_prefix <- paste0("Manifest row id=", id, ": ")

    type <- trimws(as.character(row$type))
    name <- trimws(as.character(row$name))
    source <- trimws(as.character(row$source))
    accession <- trimws(as.character(row$accession))
    url <- trimws(as.character(row$url))
    method <- tolower(trimws(as.character(row$download_method)))
    local_path <- trimws(as.character(row$local_path))

    audit_assert(nzchar(type), row_prefix, "missing type")
    audit_assert(nzchar(name), row_prefix, "missing name")
    audit_assert(nzchar(source), row_prefix, "missing source")
    audit_assert(nzchar(method), row_prefix, "missing download_method")
    audit_assert(method %in% allowed_methods, row_prefix, "invalid download_method (", method, "); allowed: ", paste(allowed_methods, collapse = ", "))
    if (!identical(method, "api")) {
      audit_assert(nzchar(local_path), row_prefix, "missing local_path")
      audit_assert(ezhu_is_safe_relative_path(local_path), row_prefix, "local_path must be a safe relative path (no absolute paths, no '..'): ", local_path)
    }

    if (identical(method, "api")) {
      audit_assert(nzchar(url), row_prefix, "download_method=api requires a url")
    } else if (identical(method, "script")) {
      audit_assert(nzchar(url), row_prefix, "download_method=script requires a url")
    } else {
      audit_assert(nzchar(url) || nzchar(accession), row_prefix, "manual/gated rows must provide at least url or accession")
    }

    integrity <- ezhu_parse_integrity(row$integrity)
    if (length(integrity) > 0) {
      known_keys <- c("sha256", "md5", "size")
      unknown <- setdiff(names(integrity), known_keys)
      audit_assert(length(unknown) == 0, row_prefix, "unknown integrity keys: ", paste(unknown, collapse = ", "))
    }

    notes_kv <- ezhu_parse_notes_kv(row$notes)
    if (tolower(type) %in% c("gwas", "eqtl", "scrna")) {
      audit_assert(nzchar(notes_kv$build %||% ""), row_prefix, "notes must include build=hg37|hg38 for type=", type)
    }
    if (tolower(type) %in% c("gwas", "chrom_ld", "ld")) {
      audit_assert(nzchar(notes_kv$ancestry %||% ""), row_prefix, "notes must include ancestry=... (e.g., EUR/EAS) for type=", type)
    }
  }

  invisible(TRUE)
}
