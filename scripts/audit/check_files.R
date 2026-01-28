audit_check_files <- function(manifest, level = "full") {
  active <- audit_active_manifest_rows(manifest)
  if (!any(active)) return(invisible(TRUE))

  if (identical(level, "ci")) {
    audit_info("CI mode: skipping local file presence and checksum checks.")
    return(invisible(TRUE))
  }

  for (i in which(active)) {
    row <- manifest[i, , drop = FALSE]
    id <- as.character(row$id)
    local_path <- trimws(as.character(row$local_path))
    if (!nzchar(local_path)) next

    method <- tolower(trimws(as.character(row$download_method)))
    if (!identical(method, "script")) {
      # Large upstream artifacts are intentionally not committed; for non-script rows
      # we warn on missing files but do not fail the audit.
      if (!file.exists(local_path)) {
        audit_warn("Manifest row id=", id, ": missing non-script local_path on disk (skipping): ", local_path)
      }
      next
    }

    row_prefix <- paste0("Manifest row id=", id, ": ")
    # Many "script" rows are large upstream inputs. For local audits we warn if they
    # are missing rather than failing the entire run (reproduction can be done on
    # cloud runners; see docs/ for guidance).
    if (!file.exists(local_path)) {
      audit_warn(row_prefix, "missing local_path on disk (skipping): ", local_path)
      next
    }

    if (file.info(local_path)$isdir) {
      if (nzchar(trimws(as.character(row$integrity)))) {
        audit_warn(row_prefix, "integrity specified for a directory; skipping integrity checks: ", local_path)
      }
      next
    }

    integrity <- ezhu_parse_integrity(row$integrity)
    if (length(integrity) == 0) next

    verify <- ezhu_verify_integrity(local_path, integrity)
    if (!isTRUE(verify$ok)) {
      is_upstream_raw <- grepl("^data/raw/", local_path)
      if (is_upstream_raw) {
        audit_warn(row_prefix, "integrity verification failed for upstream raw artifact (skipping): ", local_path)
      } else {
        audit_fail(row_prefix, "integrity verification failed for: ", local_path)
      }
    }
  }

  invisible(TRUE)
}
