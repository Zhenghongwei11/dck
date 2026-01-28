args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
  repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
  if (dir.exists(repo_root)) setwd(repo_root)
}

source("scripts/utils/repro.R")
ezhu_set_repo_root()
ezhu_activate_renv()

options(stringsAsFactors = FALSE, width = 120)

ezhu_read_opengwas_jwt <- function() {
  jwt <- Sys.getenv("OPENGWAS_JWT", unset = "")
  if (nzchar(jwt)) return(jwt)

  jwt_file <- Sys.getenv("EZHU_OPENGWAS_JWT_FILE", unset = "")
  jwt_file <- trimws(jwt_file)
  if (!nzchar(jwt_file)) {
    home <- Sys.getenv("HOME", unset = "")
    if (!nzchar(home)) return("")
    jwt_file <- file.path(home, ".config", "ezhu", "opengwas_jwt")
  }
  if (!file.exists(jwt_file)) return("")

  lines <- readLines(jwt_file, warn = FALSE)
  jwt <- trimws(paste(lines, collapse = ""))
  jwt
}

seed_env <- Sys.getenv("EZHU_SEED", unset = "")
seed <- suppressWarnings(as.integer(seed_env))
if (is.na(seed)) seed <- 20251227L
set.seed(seed)

mode <- Sys.getenv("EZHU_CAUSAL_MODE", unset = "candidates_only")
mode <- trimws(tolower(mode))
if (!mode %in% c("candidates_only", "opengwas_search", "mr", "smr")) {
  stop("EZHU_CAUSAL_MODE must be one of: candidates_only, opengwas_search, mr, smr. Got: ", dQuote(mode), call. = FALSE)
}

annotation_path <- Sys.getenv(
  "EZHU_CLUSTER_ANNOTATIONS",
  unset = "results/annotation/cluster_annotations_filled.tsv"
)
markers_dir <- Sys.getenv("EZHU_CLUSTER_MARKERS_DIR", unset = "results/annotation/markers")

top_contexts <- suppressWarnings(as.integer(Sys.getenv("EZHU_CAUSAL_TOP_CONTEXTS", unset = "3")))
if (is.na(top_contexts) || top_contexts <= 0) top_contexts <- 3L

top_markers <- suppressWarnings(as.integer(Sys.getenv("EZHU_CAUSAL_TOP_MARKERS_PER_CONTEXT", unset = "50")))
if (is.na(top_markers) || top_markers <= 0) top_markers <- 50L

if (!dir.exists(markers_dir)) stop("Missing markers directory: ", markers_dir, call. = FALSE)

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

read_pvalue_table <- function(path) {
  if (!file.exists(path)) stop("Missing pvalue table: ", path, call. = FALSE)
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- tolower(names(df))
  if (!"pvalue" %in% names(df)) stop("Missing column 'pvalue' in: ", path, call. = FALSE)

  celltype_col <- NULL
  if ("celltype" %in% names(df)) celltype_col <- "celltype"
  if (is.null(celltype_col) && ncol(df) >= 2) celltype_col <- names(df)[2]
  if (is.null(celltype_col)) stop("Cannot locate cluster/celltype column in: ", path, call. = FALSE)

  out <- data.frame(
    cluster = as.character(df[[celltype_col]]),
    pvalue = safe_numeric(df$pvalue),
    fdr = safe_numeric(df$fdr %||% NA_real_),
    stringsAsFactors = FALSE
  )
  out <- out[nzchar(out$cluster) & is.finite(out$pvalue), , drop = FALSE]
  out
}

primary_main_pvalue_path <- Sys.getenv("EZHU_PRIMARY_MAIN_PVALUE", unset = "results/scpagwas/main/merged_celltype_pvalue.csv")
primary_repeat_pvalue_path <- Sys.getenv("EZHU_PRIMARY_REPEAT_PVALUE", unset = "results/scpagwas/repeat_run/merged_celltype_pvalue.csv")
main_tbl <- read_pvalue_table(primary_main_pvalue_path)
repeat_tbl <- read_pvalue_table(primary_repeat_pvalue_path)
merged <- merge(main_tbl, repeat_tbl, by = "cluster", suffixes = c("_main", "_repeat"))
if (nrow(merged) == 0) stop("No overlapping clusters between primary main and repeat pvalue tables.", call. = FALSE)

merged$score_pvalue <- pmin(merged$pvalue_main, merged$pvalue_repeat)
merged <- merged[order(merged$score_pvalue, decreasing = FALSE), , drop = FALSE]
merged <- utils::head(merged, top_contexts)
if (nrow(merged) == 0) stop("No clusters remain after selecting top contexts.", call. = FALSE)

labels <- rep(NA_character_, nrow(merged))
if (file.exists(annotation_path)) {
  ann <- utils::read.delim(annotation_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  names(ann) <- tolower(names(ann))
  if (all(c("cluster", "proposed_cell_type") %in% names(ann))) {
    ann$cluster <- as.character(ann$cluster)
    ann$proposed_cell_type <- as.character(ann$proposed_cell_type)
    map <- setNames(ann$proposed_cell_type, ann$cluster)
    labels <- map[merged$cluster]
  }
}
labels[!nzchar(labels %||% "")] <- paste0("Cluster ", merged$cluster)

read_markers <- function(cluster_id, top_n, markers_dir) {
  path <- file.path(markers_dir, paste0("markers_cluster_", cluster_id, ".csv"))
  if (!file.exists(path)) stop("Missing marker file for cluster ", cluster_id, ": ", path, call. = FALSE)
  df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- tolower(names(df))
  if (!all(c("gene", "avg_log2fc") %in% names(df))) {
    stop("Marker file missing required columns gene/avg_log2FC: ", path, call. = FALSE)
  }

  df$gene <- toupper(trimws(df$gene))
  df$avg_log2fc <- suppressWarnings(as.numeric(df$avg_log2fc))
  df <- df[is.finite(df$avg_log2fc) & nzchar(df$gene), , drop = FALSE]
  df <- df[order(df$avg_log2fc, decreasing = TRUE), , drop = FALSE]
  df$marker_rank <- seq_len(nrow(df))
  df <- utils::head(df, top_n)
  df
}

rows <- list()
for (i in seq_len(nrow(merged))) {
  cluster_id <- merged$cluster[i]
  ctx_label <- labels[i]
  p_main <- merged$pvalue_main[i]
  p_rep <- merged$pvalue_repeat[i]

  mk <- read_markers(cluster_id, top_n = top_markers, markers_dir = markers_dir)
  mk$cluster <- cluster_id
  mk$context_label <- ctx_label
  mk$main_run_pvalue <- p_main
  mk$repeat_run_pvalue <- p_rep
  rows[[length(rows) + 1L]] <- mk[, c("gene", "cluster", "context_label", "marker_rank", "avg_log2fc", "main_run_pvalue", "repeat_run_pvalue")]
}

out <- do.call(rbind, rows)
colnames(out) <- c("gene_symbol", "cluster", "context_label", "marker_rank", "marker_avg_log2fc", "main_run_pvalue", "repeat_run_pvalue")

ezhu_dir_create("results/causal")
utils::write.table(
  out,
  file = "results/causal/candidate_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

ezhu_write_stage_metadata(
  "06_causal_inference",
  params = list(
    mode = mode,
    cluster_annotations = annotation_path,
    markers_dir = markers_dir,
    top_contexts = top_contexts,
    top_markers_per_context = top_markers
  ),
  seed = seed
)

message("Wrote: results/causal/candidate_genes.tsv (n=", nrow(out), ")")

if (mode == "candidates_only") {
  message("Causal inference mode is candidates_only; stopping after candidate list export.")
  quit(status = 0)
}

if (mode == "opengwas_search") {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required for opengwas_search mode.", call. = FALSE)
  }

  api_base <- Sys.getenv("EZHU_OPENGWAS_API", unset = "https://api.opengwas.io")
  q_terms <- Sys.getenv("EZHU_OPENGWAS_QUERY", unset = "eGFR")
  q_terms <- trimws(q_terms)
  if (!nzchar(q_terms)) q_terms <- "eGFR"

  ezhu_dir_create("results/causal/opengwas")

  url_encode <- function(x) utils::URLencode(x, reserved = TRUE)

  json_in <- Sys.getenv("EZHU_OPENGWAS_JSON_IN", unset = "")
  token <- ezhu_read_opengwas_jwt()
  curl_bin <- Sys.getenv("CURL", unset = Sys.which("curl"))

  fetch_json <- function(path, query) {
    if (!nzchar(token)) {
      stop(
      "OpenGWAS JWT is not set. ",
      "Provide OPENGWAS_JWT via env var, or write it to $HOME/.config/ezhu/opengwas_jwt, then rerun with EZHU_CAUSAL_MODE=opengwas_search.",
      call. = FALSE
    )
  }
    if (!nzchar(curl_bin)) stop("curl not found on PATH; required for opengwas_search.", call. = FALSE)

    url <- paste0(api_base, "/api/gwasinfo?q=", url_encode(query))
    out_file <- file.path("results/causal/opengwas", path)

    cmd <- paste(
      shQuote(curl_bin),
      "-sS",
      "-H", shQuote(paste0("Authorization: Bearer ", token)),
      shQuote(url)
    )
    txt <- paste(system(cmd, intern = TRUE), collapse = "\n")
    if (!nzchar(txt)) stop("OpenGWAS API returned an empty response for query: ", query, call. = FALSE)
    if (grepl("\"ERROR\"", txt, fixed = TRUE) || grepl("\"message\"", txt, fixed = TRUE)) {
      # Save raw response for audit/debug, then stop with a short message.
      writeLines(txt, out_file, useBytes = TRUE)
      stop(
        "OpenGWAS API returned an error for query: ",
        query,
        ". Raw response saved to: ",
        out_file,
        call. = FALSE
      )
    }

    writeLines(txt, out_file, useBytes = TRUE)
    jsonlite::fromJSON(txt, simplifyVector = TRUE)
  }

  gwasinfo <- NULL
  if (nzchar(json_in) && file.exists(json_in)) {
    txt <- paste(readLines(json_in, warn = FALSE), collapse = "\n")
    gwasinfo <- jsonlite::fromJSON(txt, simplifyVector = TRUE)
  } else {
    gwasinfo <- fetch_json("gwasinfo_query.json", q_terms)
  }

  # OpenGWAS currently returns a mapping { id -> record } for gwasinfo queries.
  # We export a compact, auditable table and a local filtered view (by trait substring).
  records_from_mapping <- function(x) {
    if (!is.list(x) || is.null(names(x))) return(data.frame())

    clean_chr <- function(v) {
      v <- as.character(v %||% "")
      v <- gsub("[\t\r\n]+", " ", v)
      trimws(v)
    }

    ids <- names(x)
    out <- vector("list", length(ids))
    for (i in seq_along(ids)) {
      rec <- x[[i]]
      if (!is.list(rec)) next

      out[[i]] <- data.frame(
        id = ids[[i]],
        trait = clean_chr(rec$trait),
        population = clean_chr(rec$population),
        sample_size = suppressWarnings(as.numeric(rec$sample_size %||% NA)),
        ncase = suppressWarnings(as.numeric(rec$ncase %||% NA)),
        ncontrol = suppressWarnings(as.numeric(rec$ncontrol %||% NA)),
        nsnp = suppressWarnings(as.numeric(rec$nsnp %||% NA)),
        year = suppressWarnings(as.integer(rec$year %||% NA)),
        build = clean_chr(rec$build),
        category = clean_chr(rec$category),
        subcategory = clean_chr(rec$subcategory),
        group_name = clean_chr(rec$group_name),
        consortium = clean_chr(rec$consortium),
        author = clean_chr(rec$author),
        pmid = suppressWarnings(as.integer(rec$pmid %||% NA)),
        doi = clean_chr(rec$doi),
        sex = clean_chr(rec$sex),
        unit = clean_chr(rec$unit),
        stringsAsFactors = FALSE
      )
    }

    out <- out[!vapply(out, is.null, logical(1))]
    if (length(out) == 0) return(data.frame())
    do.call(rbind, out)
  }

  tbl <- records_from_mapping(gwasinfo)
  if (nrow(tbl) == 0) {
    warning("Could not parse OpenGWAS response mapping; raw JSON saved under results/causal/opengwas/.", call. = FALSE)
  } else {
    tbl <- tbl[order(tbl$sample_size, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
    utils::write.table(
      tbl,
      file = "results/causal/opengwas/gwasinfo_query_full.tsv",
      sep = "\t",
      quote = TRUE,
      row.names = FALSE
    )
    message("Wrote: results/causal/opengwas/gwasinfo_query_full.tsv (n=", nrow(tbl), ")")

    q_pat <- Sys.getenv("EZHU_OPENGWAS_FILTER_REGEX", unset = "")
    if (!nzchar(q_pat)) q_pat <- q_terms
    q_pat <- trimws(q_pat)
    if (nzchar(q_pat)) {
      hit <- grepl(q_pat, tbl$trait, ignore.case = TRUE)
      hits <- tbl[hit, , drop = FALSE]
      utils::write.table(
        hits,
        file = "results/causal/opengwas/gwasinfo_query_filtered.tsv",
        sep = "\t",
        quote = TRUE,
        row.names = FALSE
      )
      message("Wrote: results/causal/opengwas/gwasinfo_query_filtered.tsv (n=", nrow(hits), "; filter=", dQuote(q_pat), ")")
    }
  }

  message("OpenGWAS search finished. Next: decide the outcome GWAS id(s) and register them in data/manifest.tsv.")
  quit(status = 0)
}

if (mode == "mr") {
  if (!requireNamespace("TwoSampleMR", quietly = TRUE) || !requireNamespace("ieugwasr", quietly = TRUE)) {
    stop(
      "MR mode requires TwoSampleMR + ieugwasr. Run `make mr-setup` first.\n",
      "Then set OPENGWAS_JWT in your shell and rerun with EZHU_CAUSAL_MODE=mr.",
      call. = FALSE
    )
  }

  token <- ezhu_read_opengwas_jwt()
  if (!nzchar(token)) {
    stop(
      "OpenGWAS JWT is not set. Provide OPENGWAS_JWT via env var, or write it to $HOME/.config/ezhu/opengwas_jwt, then rerun with EZHU_CAUSAL_MODE=mr.",
      call. = FALSE
    )
  }
  # Ensure downstream packages (ieugwasr/TwoSampleMR) can authenticate.
  Sys.setenv(OPENGWAS_JWT = token)

  # Choose outcome IDs (default: those registered in data/manifest.tsv under opengwas_id=...).
  if (!requireNamespace("jsonlite", quietly = TRUE)) stop("jsonlite is required.", call. = FALSE)
  source("scripts/utils/manifest.R")
  mf <- ezhu_read_manifest("data/manifest.tsv")
  mf$type <- tolower(trimws(as.character(mf$type)))
  notes <- lapply(as.character(mf$notes), ezhu_parse_notes_kv)
  opengwas_ids <- vapply(notes, function(x) x$opengwas_id %||% "", character(1))
  ok_out <- mf$type == "gwas" & nzchar(opengwas_ids)
  outcome_ids <- unique(opengwas_ids[ok_out])

  # Conservative default: just the two renal outcomes we pinned in the decision log.
  default_outcomes <- c("ebi-a-GCST90103634", "ebi-a-GCST006586")
  outcome_ids <- unique(c(default_outcomes, outcome_ids))

  max_genes <- suppressWarnings(as.integer(Sys.getenv("EZHU_MR_MAX_GENES", unset = "30")))
  if (is.na(max_genes) || max_genes <= 0) max_genes <- 30L

  mr_p1 <- suppressWarnings(as.numeric(Sys.getenv("EZHU_MR_P1", unset = "5e-8")))
  if (!is.finite(mr_p1) || mr_p1 <= 0 || mr_p1 >= 1) mr_p1 <- 5e-8

  mr_clump <- trimws(Sys.getenv("EZHU_MR_CLUMP", unset = "1"))
  mr_clump <- identical(mr_clump, "1") || identical(tolower(mr_clump), "true")

  # Map gene symbols -> ENSG for OpenGWAS eQTL exposures (eqtl-a-ENSG...).
  # We keep a small cache under data/references/ensembl/ to avoid repeated API calls.
  ensembl_cache_path <- Sys.getenv(
    "EZHU_ENSEMBL_SYMBOL_CACHE",
    unset = "data/references/ensembl/symbol_to_ensg.tsv"
  )

  ezhu_dir_create(dirname(ensembl_cache_path))

  cache <- data.frame(gene_symbol = character(), ensg = character(), stringsAsFactors = FALSE)
  if (file.exists(ensembl_cache_path)) {
    cache <- utils::read.delim(ensembl_cache_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    names(cache) <- tolower(names(cache))
    if (!all(c("gene_symbol", "ensg") %in% names(cache))) cache <- data.frame(gene_symbol = character(), ensg = character(), stringsAsFactors = FALSE)
    cache$gene_symbol <- toupper(trimws(cache$gene_symbol))
    cache$ensg <- trimws(cache$ensg)
  }

  curl_bin <- Sys.getenv("CURL", unset = Sys.which("curl"))
  if (!nzchar(curl_bin)) stop("curl not found; required for Ensembl symbol mapping.", call. = FALSE)

  ensembl_xref <- function(symbol) {
    symbol <- toupper(trimws(symbol))
    hit <- cache$ensg[cache$gene_symbol == symbol][1]
    if (!is.na(hit) && nzchar(hit)) return(hit)

    url <- paste0("https://rest.ensembl.org/xrefs/symbol/homo_sapiens/", utils::URLencode(symbol, reserved = TRUE), "?content-type=application/json")
    cmd <- paste(shQuote(curl_bin), "-sS", shQuote(url))
    txt <- paste(system(cmd, intern = TRUE), collapse = "\n")
    if (!nzchar(txt)) return(NA_character_)
    if (grepl("error", txt, ignore.case = TRUE) && grepl("message", txt, fixed = TRUE)) return(NA_character_)

    recs <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = TRUE), error = function(e) NULL)
    if (is.null(recs) || length(recs) == 0) return(NA_character_)

    if (is.data.frame(recs)) {
      if (!"id" %in% names(recs)) return(NA_character_)
      ids <- as.character(recs$id)
      ids <- ids[grepl("^ENSG", ids)]
      if (length(ids) == 0) return(NA_character_)
      return(ids[[1]])
    }

    # Fallback if parsed as list
    ids <- unlist(lapply(recs, function(x) as.character(x$id %||% "")), use.names = FALSE)
    ids <- ids[grepl("^ENSG", ids)]
    if (length(ids) == 0) return(NA_character_)
    ids[[1]]
  }

  genes <- unique(out$gene_symbol)
  genes <- genes[nzchar(genes)]

  # Prefer HGNC-like symbols; keep ENSG if already present.
  ensg <- vapply(genes, function(g) if (grepl("^ENSG", g)) g else ensembl_xref(g), character(1))
  map_df <- data.frame(gene_symbol = genes, ensg = ensg, stringsAsFactors = FALSE)
  map_df <- map_df[!is.na(map_df$ensg) & nzchar(map_df$ensg), , drop = FALSE]

  # Update cache (append new mappings).
  if (nrow(map_df) > 0) {
    add <- map_df[!map_df$gene_symbol %in% cache$gene_symbol, , drop = FALSE]
    if (nrow(add) > 0) {
      cache2 <- rbind(cache[, c("gene_symbol", "ensg"), drop = FALSE], add[, c("gene_symbol", "ensg"), drop = FALSE])
      cache2 <- cache2[order(cache2$gene_symbol), , drop = FALSE]
      utils::write.table(cache2, file = ensembl_cache_path, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }

  exposures <- unique(map_df$ensg)
  exposures <- exposures[seq_len(min(length(exposures), max_genes))]
  exposure_ids <- paste0("eqtl-a-", exposures)

  ezhu_dir_create("results/causal/mr")

  params <- list(
    mode = "mr",
    exposures_n = length(exposure_ids),
    outcomes = outcome_ids,
    outcome_n = length(outcome_ids),
    max_genes = max_genes,
    p1 = mr_p1,
    clump = mr_clump,
    exposure_prefix = "eqtl-a-",
    ensembl_cache = ensembl_cache_path
  )

  # 1) Extract instruments (OpenGWAS server-side clumping).
  instruments <- TwoSampleMR::extract_instruments(outcomes = exposure_ids, p1 = mr_p1, clump = mr_clump)
  if (is.null(instruments) || nrow(instruments) == 0) {
    stop("No instruments extracted for the selected exposures (n=", length(exposure_ids), ").", call. = FALSE)
  }

  # 2) Extract outcomes at instrument SNPs.
  outcome_dat <- TwoSampleMR::extract_outcome_data(snps = unique(instruments$SNP), outcomes = outcome_ids)
  if (is.null(outcome_dat) || nrow(outcome_dat) == 0) {
    stop("No outcome data returned for the selected outcomes (n=", length(outcome_ids), ").", call. = FALSE)
  }

  # 3) Harmonise and run MR.
  harmonised <- TwoSampleMR::harmonise_data(exposure_dat = instruments, outcome_dat = outcome_dat)
  mr_res <- TwoSampleMR::mr(harmonised)
  het <- tryCatch(TwoSampleMR::mr_heterogeneity(harmonised), error = function(e) NULL)
  pleio <- tryCatch(TwoSampleMR::mr_pleiotropy_test(harmonised), error = function(e) NULL)

  utils::write.table(instruments, file = "results/causal/mr/instruments.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(harmonised, file = "results/causal/mr/harmonised.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  utils::write.table(mr_res, file = "results/causal/mr/mr_results.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  if (!is.null(het)) utils::write.table(het, file = "results/causal/mr/heterogeneity.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
  if (!is.null(pleio)) utils::write.table(pleio, file = "results/causal/mr/pleiotropy.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

  # Also write compact “anchor tables” outside results/causal/mr/ so they can be committed
  # while keeping verbose intermediate outputs optional.
  mr_res$method <- as.character(mr_res$method)
  mr_res$outcome_id <- as.character(mr_res$id.outcome)
  mr_res$exposure_id <- as.character(mr_res$id.exposure)
  mr_res$exposure_ensg <- sub("^eqtl-a-", "", mr_res$exposure_id)
  mr_res$outcome_trait <- as.character(mr_res$outcome)

  mr_res_out <- mr_res[, c("exposure_ensg", "exposure_id", "outcome_id", "outcome_trait", "method", "nsnp", "b", "se", "pval")]
  utils::write.table(
    mr_res_out,
    file = "results/causal/mr_summary_all_methods.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  ivw <- mr_res_out[tolower(mr_res_out$method) == "inverse variance weighted", , drop = FALSE]
  if (nrow(ivw) > 0) {
    ivw$ivw_fdr <- stats::p.adjust(ivw$pval, method = "BH")
    ivw <- ivw[order(ivw$pval, decreasing = FALSE), , drop = FALSE]
  }
  utils::write.table(
    ivw,
    file = "results/causal/mr_summary_ivw.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # Build a publication-friendly validation table that combines:
  # - IVW/WM/Egger estimates
  # - heterogeneity (Q_pval) and Egger intercept p
  # - IVW FDR (BH) across all gene-outcome pairs
  #
  # This table is designed to be committed as a compact “anchor table”, while keeping
  # the verbose intermediate MR artifacts under results/causal/mr/ optional.
  method_key <- function(x) trimws(tolower(as.character(x)))
  mk_key <- function(exposure_ensg, outcome_id) paste0(exposure_ensg, "||", outcome_id)

  # Reverse-map ENSG -> input gene symbol (from candidate list), for readability.
  map_df$ensg <- as.character(map_df$ensg)
  map_df$gene_symbol <- as.character(map_df$gene_symbol)
  map_df <- map_df[nzchar(map_df$ensg) & nzchar(map_df$gene_symbol), , drop = FALSE]
  map_df <- map_df[order(map_df$ensg, map_df$gene_symbol), , drop = FALSE]
  ensg_to_symbol <- tapply(map_df$gene_symbol, map_df$ensg, function(v) v[[1]])

  base_keys <- unique(mr_res_out[, c("exposure_ensg", "exposure_id", "outcome_id", "outcome_trait")])
  base_keys$mr_pair_id <- mk_key(base_keys$exposure_ensg, base_keys$outcome_id)
  base_keys$exposure_gene_symbol_input <- unname(ensg_to_symbol[base_keys$exposure_ensg])
  base_keys$exposure_gene_symbol_input[is.na(base_keys$exposure_gene_symbol_input)] <- ""

  # Also fetch the current Ensembl display symbol for each ENSG (helps resolve aliases
  # such as histone gene synonyms).
  ensg_symbol_cache_path <- file.path(dirname(ensembl_cache_path), "ensg_to_symbol.tsv")
  ensg_cache <- data.frame(ensg = character(), gene_symbol_ensembl = character(), stringsAsFactors = FALSE)
  if (file.exists(ensg_symbol_cache_path)) {
    ensg_cache <- utils::read.delim(ensg_symbol_cache_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
    names(ensg_cache) <- tolower(names(ensg_cache))
    if (!all(c("ensg", "gene_symbol_ensembl") %in% names(ensg_cache))) {
      ensg_cache <- data.frame(ensg = character(), gene_symbol_ensembl = character(), stringsAsFactors = FALSE)
    } else {
      ensg_cache$ensg <- trimws(as.character(ensg_cache$ensg))
      ensg_cache$gene_symbol_ensembl <- trimws(as.character(ensg_cache$gene_symbol_ensembl))
    }
  }

  ensembl_lookup_symbol <- function(ensg_id) {
    ensg_id <- trimws(as.character(ensg_id))
    if (!nzchar(ensg_id)) return(NA_character_)
    hit <- ensg_cache$gene_symbol_ensembl[ensg_cache$ensg == ensg_id][1]
    if (!is.na(hit) && nzchar(hit)) return(hit)

    url <- paste0("https://rest.ensembl.org/lookup/id/", utils::URLencode(ensg_id, reserved = TRUE), "?content-type=application/json")
    cmd <- paste(shQuote(curl_bin), "-sS", shQuote(url))
    txt <- paste(system(cmd, intern = TRUE), collapse = "\n")
    if (!nzchar(txt)) return(NA_character_)
    rec <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = TRUE), error = function(e) NULL)
    if (is.null(rec) || !is.list(rec)) return(NA_character_)
    sym <- as.character(rec$display_name %||% "")
    sym <- trimws(sym)
    if (!nzchar(sym)) return(NA_character_)
    sym
  }

  uniq_ensg <- unique(base_keys$exposure_ensg)
  ensg_symbol <- vapply(uniq_ensg, ensembl_lookup_symbol, character(1))
  sym_df <- data.frame(ensg = uniq_ensg, gene_symbol_ensembl = ensg_symbol, stringsAsFactors = FALSE)
  sym_df <- sym_df[!is.na(sym_df$gene_symbol_ensembl) & nzchar(sym_df$gene_symbol_ensembl), , drop = FALSE]
  if (nrow(sym_df) > 0) {
    # Update cache (append new).
    add2 <- sym_df[!sym_df$ensg %in% ensg_cache$ensg, , drop = FALSE]
    if (nrow(add2) > 0) {
      cache2 <- rbind(ensg_cache[, c("ensg", "gene_symbol_ensembl"), drop = FALSE], add2[, c("ensg", "gene_symbol_ensembl"), drop = FALSE])
      cache2 <- cache2[order(cache2$ensg), , drop = FALSE]
      utils::write.table(cache2, file = ensg_symbol_cache_path, sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }

  base_keys$exposure_gene_symbol_ensembl <- ensg_symbol[match(base_keys$exposure_ensg, uniq_ensg)]
  base_keys$exposure_gene_symbol_ensembl[is.na(base_keys$exposure_gene_symbol_ensembl)] <- ""

  # Method-specific estimates.
  pick_method <- function(df, want_method) {
    m <- df[method_key(df$method) == method_key(want_method), , drop = FALSE]
    if (nrow(m) == 0) return(data.frame())
    out <- m[, c("exposure_ensg", "outcome_id", "nsnp", "b", "se", "pval"), drop = FALSE]
    out$mr_pair_id <- mk_key(out$exposure_ensg, out$outcome_id)
    out
  }

  ivw2 <- pick_method(mr_res_out, "Inverse variance weighted")
  if (nrow(ivw2) > 0) {
    ivw2$ivw_fdr <- stats::p.adjust(ivw2$pval, method = "BH")
    colnames(ivw2) <- c("exposure_ensg", "outcome_id", "ivw_nsnp", "ivw_b", "ivw_se", "ivw_pval", "mr_pair_id", "ivw_fdr")
  }

  wm <- pick_method(mr_res_out, "Weighted median")
  if (nrow(wm) > 0) colnames(wm) <- c("exposure_ensg", "outcome_id", "wm_nsnp", "wm_b", "wm_se", "wm_pval", "mr_pair_id")

  egger <- pick_method(mr_res_out, "MR Egger")
  if (nrow(egger) > 0) colnames(egger) <- c("exposure_ensg", "outcome_id", "egger_nsnp", "egger_b", "egger_se", "egger_pval", "mr_pair_id")

  # Heterogeneity and Egger intercept.
  het_ivw <- data.frame()
  het_egger <- data.frame()
  if (!is.null(het) && nrow(het) > 0) {
    het$method <- as.character(het$method)
    het$id.exposure <- as.character(het$id.exposure)
    het$id.outcome <- as.character(het$id.outcome)
    het$exposure_ensg <- sub("^eqtl-a-", "", het$id.exposure)
    het$mr_pair_id <- mk_key(het$exposure_ensg, het$id.outcome)

    het_ivw <- het[method_key(het$method) == method_key("Inverse variance weighted"), c("mr_pair_id", "Q_pval"), drop = FALSE]
    if (nrow(het_ivw) > 0) colnames(het_ivw) <- c("mr_pair_id", "het_ivw_q_pval")

    het_egger <- het[method_key(het$method) == method_key("MR Egger"), c("mr_pair_id", "Q_pval"), drop = FALSE]
    if (nrow(het_egger) > 0) colnames(het_egger) <- c("mr_pair_id", "het_egger_q_pval")
  }

  pleio2 <- data.frame()
  if (!is.null(pleio) && nrow(pleio) > 0) {
    pleio$id.exposure <- as.character(pleio$id.exposure)
    pleio$id.outcome <- as.character(pleio$id.outcome)
    pleio$exposure_ensg <- sub("^eqtl-a-", "", pleio$id.exposure)
    pleio$mr_pair_id <- mk_key(pleio$exposure_ensg, pleio$id.outcome)
    pleio2 <- pleio[, c("mr_pair_id", "egger_intercept", "se", "pval"), drop = FALSE]
    colnames(pleio2) <- c("mr_pair_id", "egger_intercept", "egger_intercept_se", "egger_intercept_pval")
  }

  validation <- base_keys
  if (nrow(ivw2) > 0) validation <- merge(validation, ivw2[, c("mr_pair_id", "ivw_nsnp", "ivw_b", "ivw_se", "ivw_pval", "ivw_fdr")], by = "mr_pair_id", all.x = TRUE)
  if (nrow(wm) > 0) validation <- merge(validation, wm[, c("mr_pair_id", "wm_nsnp", "wm_b", "wm_se", "wm_pval")], by = "mr_pair_id", all.x = TRUE)
  if (nrow(egger) > 0) validation <- merge(validation, egger[, c("mr_pair_id", "egger_nsnp", "egger_b", "egger_se", "egger_pval")], by = "mr_pair_id", all.x = TRUE)
  if (nrow(het_ivw) > 0) validation <- merge(validation, het_ivw, by = "mr_pair_id", all.x = TRUE)
  if (nrow(het_egger) > 0) validation <- merge(validation, het_egger, by = "mr_pair_id", all.x = TRUE)
  if (nrow(pleio2) > 0) validation <- merge(validation, pleio2, by = "mr_pair_id", all.x = TRUE)

  # Keep a clean column order and sort by IVW p-value.
  keep_cols <- c(
    "exposure_gene_symbol_ensembl",
    "exposure_gene_symbol_input",
    "exposure_ensg",
    "exposure_id",
    "outcome_id",
    "outcome_trait",
    "ivw_nsnp",
    "ivw_b",
    "ivw_se",
    "ivw_pval",
    "ivw_fdr",
    "wm_b",
    "wm_se",
    "wm_pval",
    "egger_b",
    "egger_se",
    "egger_pval",
    "het_ivw_q_pval",
    "het_egger_q_pval",
    "egger_intercept",
    "egger_intercept_se",
    "egger_intercept_pval"
  )
  keep_cols <- keep_cols[keep_cols %in% colnames(validation)]
  validation <- validation[, keep_cols, drop = FALSE]
  if ("ivw_pval" %in% colnames(validation)) {
    validation <- validation[order(validation$ivw_pval, decreasing = FALSE, na.last = TRUE), , drop = FALSE]
  }

  utils::write.table(
    validation,
    file = "results/causal/mr_validation.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  message("Wrote: results/causal/mr_validation.tsv (n=", nrow(validation), ")")

  ezhu_write_stage_metadata("06_causal_inference_mr", params = params, seed = seed)
  message("Wrote MR outputs under results/causal/mr/")
  quit(status = 0)
}

stop(
  "Causal inference mode '", mode, "' is scaffolded but not executed yet. ",
  "Provide an explicit eQTL resource + GWAS outcome definition and implement the chosen engine (MR or SMR) before running this mode.",
  call. = FALSE
)
