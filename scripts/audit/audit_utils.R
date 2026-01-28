audit_info <- function(...) message("[INFO] ", paste0(..., collapse = ""))
audit_warn <- function(...) message("[WARN] ", paste0(..., collapse = ""))

audit_fail <- function(...) {
  msg <- paste0(..., collapse = "")
  stop("[FAIL] ", msg, call. = FALSE)
}

audit_assert <- function(condition, ...) {
  if (!isTRUE(condition)) audit_fail(...)
  invisible(TRUE)
}

audit_parse_level <- function(args) {
  level <- Sys.getenv("EZHU_AUDIT_LEVEL", unset = "")
  if (nzchar(level)) return(tolower(level))

  level_flag <- grep("^--level=", args, value = TRUE)
  if (length(level_flag) == 1) return(tolower(sub("^--level=", "", level_flag)))

  level_pos <- which(args == "--level")
  if (length(level_pos) == 1 && level_pos < length(args)) return(tolower(args[level_pos + 1]))

  "full"
}

audit_active_manifest_rows <- function(df) {
  fields <- c("type", "name", "source", "accession", "url", "local_path")
  fields <- intersect(fields, names(df))
  if (length(fields) == 0) return(rep(FALSE, nrow(df)))

  apply(df[, fields, drop = FALSE], 1, function(row) any(!is.na(row) & nzchar(trimws(as.character(row)))))
}

audit_read_text <- function(path) {
  if (!file.exists(path)) return(character())
  readLines(path, warn = FALSE, encoding = "UTF-8")
}

