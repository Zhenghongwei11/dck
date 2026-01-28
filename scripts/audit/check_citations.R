audit_check_citations <- function() {
  if (!file.exists("docs/references.bib")) {
    audit_warn("docs/references.bib not found; skipping citations audit.")
    return(invisible(TRUE))
  }

  py <- Sys.which("python3")
  if (!nzchar(py)) {
    audit_warn("python3 not found; skipping citations audit.")
    return(invisible(TRUE))
  }

  script <- "scripts/audit/check_citations.py"
  if (!file.exists(script)) {
    audit_warn("Citations audit script not found; skipping: ", script)
    return(invisible(TRUE))
  }

  out_dir <- file.path("docs", "audit_runs", "current")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  report_path <- file.path(out_dir, "citations_report.md")

  args <- c(script, "--repo-root", ".", "--bib", "docs/references.bib", "--output", report_path)
  out <- suppressWarnings(system2(py, args, stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && !identical(status, 0L)) {
    audit_fail("Citations audit failed. See: ", report_path)
  }

  invisible(TRUE)
}

