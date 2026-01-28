audit_check_renv <- function() {
  audit_assert(file.exists("renv.lock"), "Missing renv.lock (required for reproducibility).")

  if (requireNamespace("jsonlite", quietly = TRUE)) {
    lock <- tryCatch(jsonlite::fromJSON("renv.lock"), error = function(e) NULL)
    audit_assert(!is.null(lock), "Failed to parse renv.lock as JSON.")
    audit_assert(!is.null(lock$R$Version), "renv.lock is missing R version metadata.")
  }

  invisible(TRUE)
}

