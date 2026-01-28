audit_check_scripts <- function() {
  script_files <- list.files("scripts", pattern = "\\.(R|sh)$", recursive = TRUE, full.names = TRUE)
  # Audit scripts themselves must also obey the path contract.
  script_files <- script_files[!grepl("^scripts/audit/check_scripts\\.R$", script_files)]

  forbidden_patterns <- c(
    "/Users/",
    "C:\\\\",
    "~/"
  )

  for (path in script_files) {
    lines <- audit_read_text(path)
    if (length(lines) == 0) next

    for (pat in forbidden_patterns) {
      if (any(grepl(pat, lines, fixed = TRUE))) {
        audit_fail("Forbidden path pattern detected in ", path, ": ", pat)
      }
    }
  }

  invisible(TRUE)
}
