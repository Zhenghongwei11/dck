audit_check_figure_index <- function(level = "full") {
  index_path <- "docs/FIGURE_INDEX.md"
  lines <- audit_read_text(index_path)
  if (length(lines) == 0) {
    audit_fail("Figure index missing or empty: ", index_path)
  }

  text <- paste(lines, collapse = "\n")
  # Split into sections by headings of the form "### Figure ..."
  section_starts <- gregexpr("(?m)^###\\s+Figure\\b", text, perl = TRUE)[[1]]
  if (length(section_starts) == 1 && section_starts[1] == -1) {
    audit_warn("No '### Figure ...' sections detected in ", index_path, "; skipping provenance checks.")
    return(invisible(TRUE))
  }

  # Build sections by slicing at heading offsets.
  starts <- as.integer(section_starts)
  starts <- c(starts, nchar(text) + 1L)
  sections <- vector("list", length(starts) - 1L)
  for (i in seq_len(length(sections))) {
    sections[[i]] <- substr(text, starts[i], starts[i + 1L] - 1L)
  }

  # Extract backticked paths in each section.
  extract_backticked <- function(s) {
    m <- gregexpr("`([^`]+)`", s, perl = TRUE)[[1]]
    if (length(m) == 1 && m[1] == -1) return(character())
    matches <- regmatches(s, list(m))[[1]]
    gsub("^`|`$", "", matches)
  }

  normalize_rel <- function(p) {
    p <- gsub("\\\\", "/", p, fixed = FALSE)
    p <- sub("^\\./", "", p)
    p
  }

  is_placeholder <- function(p) grepl("<[^>]+>", p)
  is_glob <- function(p) grepl("[*?\\[]", p)

  for (sec in sections) {
    entries <- normalize_rel(extract_backticked(sec))

    outputs <- entries[grepl("^plots/publication/", entries)]
    scripts <- entries[grepl("^scripts/", entries)]
    anchors <- entries[grepl("^results/", entries)]

    # If a section references publication outputs, require script provenance unless explicitly exempted.
    if (length(outputs) > 0) {
      exempt <- grepl("Not script-generated", sec, fixed = TRUE) ||
        grepl("prompt-based", sec, ignore.case = TRUE) ||
        grepl("cloud run", sec, ignore.case = TRUE)

      if (!exempt && length(scripts) == 0) {
        audit_fail("Missing script provenance in figure index section for output(s): ", paste(outputs, collapse = ", "))
      }

      # For numeric figures, anchor tables should exist; allow schematic/prompt-based exemptions.
      if (!exempt && length(anchors) == 0) {
        audit_warn("No anchor tables listed in figure index section for output(s): ", paste(outputs, collapse = ", "))
      }
    }

    if (identical(level, "ci")) next

    # Verify referenced paths exist (skip placeholders).
    for (p in unique(c(scripts, anchors, outputs))) {
      if (is_placeholder(p)) next
      if (is_glob(p)) {
        matches <- Sys.glob(p)
        if (length(matches) == 0) {
          audit_fail("Glob path referenced in docs/FIGURE_INDEX.md matched no files: ", p)
        }
        next
      }
      if (!file.exists(p)) {
        audit_fail("Path referenced in docs/FIGURE_INDEX.md does not exist on disk: ", p)
      }
    }
  }

  invisible(TRUE)
}
