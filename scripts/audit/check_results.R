audit_check_results <- function(level = "full") {
  if (!dir.exists("results")) {
    audit_warn("results/ directory not found; skipping results checks.")
    return(invisible(TRUE))
  }

  files <- list.files("results", pattern = "\\.(csv|tsv)$", recursive = TRUE, full.names = TRUE)
  files <- files[!grepl("^results/metadata/", files)]

  if (length(files) == 0) {
    audit_warn("No results tables found under results/ (csv/tsv).")
    return(invisible(TRUE))
  }

  for (path in files) {
    sep <- if (grepl("\\.tsv$", path, ignore.case = TRUE)) "\t" else ","
    df <- tryCatch(
      utils::read.table(path, sep = sep, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, nrows = 10000),
      error = function(e) NULL
    )
    if (is.null(df)) {
      audit_warn("Failed to read results table (skipping): ", path)
      next
    }

    colnames_lower <- tolower(names(df))
    p_cols <- grep("^(p|pval|p_value|pvalue|p\\.value|padj|fdr|qvalue)$", colnames_lower, value = TRUE)
    if (length(p_cols) == 0) next

    for (col in p_cols) {
      values <- suppressWarnings(as.numeric(df[[col]]))
      values <- values[!is.na(values)]
      if (length(values) == 0) next

      if (any(values < 0 | values > 1)) {
        audit_fail("Invalid probability-like values in ", path, " column '", col, "' (expected in [0,1]).")
      }
    }
  }

  invisible(TRUE)
}

