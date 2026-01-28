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

seed_env <- Sys.getenv("EZHU_SEED", unset = "")
seed <- suppressWarnings(as.integer(seed_env))
if (is.na(seed)) seed <- 20251227L
set.seed(seed)

dry_run <- identical(Sys.getenv("EZHU_COLOC_DRY_RUN", unset = "0"), "1")

if (!isTRUE(dry_run) && !requireNamespace("coloc", quietly = TRUE)) {
  stop("Missing package: coloc. Run: make coloc-setup", call. = FALSE)
}
if (!requireNamespace("dplyr", quietly = TRUE)) stop("Missing package: dplyr. Run: make coloc-setup", call. = FALSE)
if (!requireNamespace("ieugwasr", quietly = TRUE)) stop("Missing package: ieugwasr. Run: make coloc-setup", call. = FALSE)
if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Missing package: jsonlite. Run: make coloc-setup", call. = FALSE)

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

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

mr_validation_path <- Sys.getenv("EZHU_MR_VALIDATION", unset = "results/causal/mr_validation.tsv")
mr_instruments_path <- Sys.getenv("EZHU_MR_INSTRUMENTS", unset = "results/causal/mr/instruments.tsv")
region_mode <- trimws(tolower(Sys.getenv("EZHU_COLOC_REGION_MODE", unset = "auto")))
if (!region_mode %in% c("auto", "instrument", "gene")) region_mode <- "auto"

top_hits <- suppressWarnings(as.integer(Sys.getenv("EZHU_COLOC_TOP_HITS", unset = "5")))
if (is.na(top_hits) || top_hits <= 0) top_hits <- 5L

window_kb <- suppressWarnings(as.integer(Sys.getenv("EZHU_COLOC_WINDOW_KB", unset = "500")))
if (is.na(window_kb) || window_kb <= 0) window_kb <- 500L
window_bp <- as.integer(window_kb) * 1000L

min_overlap <- suppressWarnings(as.integer(Sys.getenv("EZHU_COLOC_MIN_OVERLAP", unset = "50")))
if (is.na(min_overlap) || min_overlap <= 0) min_overlap <- 50L

strict_no_results <- identical(Sys.getenv("EZHU_COLOC_STRICT", unset = "1"), "1")
if (isTRUE(dry_run)) strict_no_results <- FALSE

p1 <- safe_numeric(Sys.getenv("EZHU_COLOC_P1", unset = "1e-4"))
p2 <- safe_numeric(Sys.getenv("EZHU_COLOC_P2", unset = "1e-4"))
p12 <- safe_numeric(Sys.getenv("EZHU_COLOC_P12", unset = "1e-5"))
if (!is.finite(p1) || p1 <= 0) p1 <- 1e-4
if (!is.finite(p2) || p2 <= 0) p2 <- 1e-4
if (!is.finite(p12) || p12 <= 0) p12 <- 1e-5

use_vcf <- identical(Sys.getenv("EZHU_COLOC_USE_VCF", unset = "0"), "1")
if (isTRUE(use_vcf) && !requireNamespace("gwasvcf", quietly = TRUE)) {
  stop("Missing package: gwasvcf. Either run: make coloc-setup, or set EZHU_COLOC_USE_VCF=0 to use the API backend.", call. = FALSE)
}

if (!file.exists(mr_validation_path)) stop("Missing MR validation table: ", mr_validation_path, call. = FALSE)

mr <- utils::read.delim(mr_validation_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

names(mr) <- tolower(names(mr))
inst <- NULL
if (file.exists(mr_instruments_path)) {
  inst <- utils::read.delim(mr_instruments_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, quote = "")
  names(inst) <- tolower(names(inst))
}

required_mr <- c("exposure_id", "outcome_id")
if (!all(required_mr %in% names(mr))) stop("MR table missing columns: ", paste(setdiff(required_mr, names(mr)), collapse = ", "), call. = FALSE)
required_inst <- c("id.exposure", "chr.exposure", "pos.exposure", "snp", "pval.exposure")
has_instruments <- !is.null(inst) && all(required_inst %in% names(inst))

extract_ensg <- function(exposure_id, ensg_hint) {
  ensg_hint <- as.character(ensg_hint %||% "")
  if (nzchar(ensg_hint)) return(ensg_hint)
  exposure_id <- as.character(exposure_id %||% "")
  m <- regexpr("ENSG[0-9]+", exposure_id, perl = TRUE)
  if (m[1] > 0) {
    regmatches(exposure_id, m)
  } else {
    ""
  }
}

fetch_gene_grch37 <- function(ensg, curl_bin) {
  ensg <- trimws(as.character(ensg))
  if (!nzchar(ensg)) return(NULL)
  url <- paste0("https://grch37.rest.ensembl.org/lookup/id/", ensg, "?content-type=application/json")
  txt <- ezhu_http_get_text(url, curl_bin = curl_bin)
  if (!nzchar(txt)) return(NULL)

  parsed <- tryCatch(jsonlite::fromJSON(txt), error = function(e) NULL)
  if (is.null(parsed)) return(NULL)

  chr <- as.character(parsed$seq_region_name %||% "")
  start <- suppressWarnings(as.integer(parsed$start %||% NA))
  end <- suppressWarnings(as.integer(parsed$end %||% NA))
  strand <- suppressWarnings(as.integer(parsed$strand %||% NA))
  if (!nzchar(chr) || !is.finite(start) || !is.finite(end)) return(NULL)

  list(chr = chr, start = start, end = end, strand = strand %||% NA_integer_)
}

fetch_region_rsids_1000g <- function(chr, center, radius_bp, start, end, max_positions = 5000L, pop = "EUR") {
  chr <- trimws(as.character(chr))
  center <- as.integer(center)
  radius_bp <- as.integer(radius_bp)
  start <- as.integer(start)
  end <- as.integer(end)
  max_positions <- as.integer(max_positions)
  pop <- trimws(as.character(pop %||% "EUR"))
  if (!nzchar(chr) || !is.finite(center) || center <= 0) return(list(rsids = character(), maf = setNames(numeric(), character())))
  if (!is.finite(radius_bp) || radius_bp <= 0) return(list(rsids = character(), maf = setNames(numeric(), character())))
  if (!is.finite(start) || !is.finite(end) || start <= 0 || end <= 0 || end < start) return(list(rsids = character(), maf = setNames(numeric(), character())))
  if (!is.finite(max_positions) || max_positions <= 0) max_positions <- 5000L

  ensure_opengwas_token()

  chrpos_center <- paste0(chr, ":", center)
  vpos <- tryCatch(ieugwasr::variants_chrpos(chrpos_center, radius = radius_bp), error = function(e) NULL)
  if (is.null(vpos) || nrow(vpos) == 0 || !"pos" %in% names(vpos)) return(list(rsids = character(), maf = setNames(numeric(), character())))

  pos <- suppressWarnings(as.integer(vpos$pos))
  pos <- pos[is.finite(pos)]
  pos <- unique(pos[pos >= start & pos <= end])
  if (length(pos) == 0) return(list(rsids = character(), maf = setNames(numeric(), character())))

  # Sample a manageable set of variant positions, then map to 1000G rsIDs + MAF via afl2_chrpos.
  # This avoids Ensembl rsID mismatches while keeping the SNP set explicit and reproducible.
  if (length(pos) > max_positions) {
    set.seed(seed)
    pos <- sample(pos, max_positions)
  }
  chrpos <- paste0(chr, ":", pos)
  chunks <- split(chrpos, ceiling(seq_along(chrpos) / 64))

  rsids <- character()
  maf_map <- setNames(numeric(), character())
  for (i in seq_along(chunks)) {
    afl <- tryCatch(ieugwasr::afl2_chrpos(chunks[[i]], reference = "1000g"), error = function(e) NULL)
    if (is.null(afl) || nrow(afl) == 0 || !"id" %in% names(afl)) next

    ids <- trimws(as.character(afl$id))
    ids <- ids[nzchar(ids) & grepl("^rs[0-9]+$", ids)]
    if (length(ids) == 0) next

    rsids <- c(rsids, ids)

    af_col <- paste0("AF.", toupper(pop))
    if (af_col %in% names(afl)) {
      af <- safe_numeric(afl[[af_col]])
      maf <- pmin(af, 1 - af)
      ok <- is.finite(maf) & maf > 0 & maf < 0.5
      if (any(ok)) {
        maf_map[ids[ok]] <- maf[ok]
      }
    }
  }

  rsids <- unique(rsids)
  list(rsids = rsids, maf = maf_map)
}

mr$ivw_fdr <- safe_numeric(mr$ivw_fdr %||% NA_real_)
mr$ivw_pval <- safe_numeric(mr$ivw_pval %||% mr$ivw_pval %||% NA_real_)
mr$rank_score <- mr$ivw_fdr
mr$rank_score[!is.finite(mr$rank_score)] <- mr$ivw_pval[!is.finite(mr$rank_score)]
mr <- mr[is.finite(mr$rank_score), , drop = FALSE]
mr <- mr[order(mr$rank_score, decreasing = FALSE), , drop = FALSE]
mr <- utils::head(mr, top_hits)
if (nrow(mr) == 0) stop("No MR rows remain for colocalization (after filtering on ivw_fdr/ivw_pval).", call. = FALSE)

base_url <- Sys.getenv("EZHU_OPENGWAS_FILES_BASE", unset = "https://gwas.mrcieu.ac.uk/files")
curl_bin <- Sys.which(Sys.getenv("CURL", unset = "curl"))

ezhu_http_get_text <- function(url, curl_bin = "") {
  url <- trimws(as.character(url))
  if (!nzchar(url)) return("")

  if (nzchar(curl_bin)) {
    txt <- tryCatch(
      paste(system2(curl_bin, c("-sS", "-L", url), stdout = TRUE, stderr = FALSE), collapse = ""),
      error = function(e) ""
    )
    if (nzchar(txt)) return(txt)
  }

  old_timeout <- getOption("timeout")
  on.exit(options(timeout = old_timeout), add = TRUE)
  options(timeout = max(120L, suppressWarnings(as.integer(old_timeout)) %||% 120L))

  con <- tryCatch(base::url(url, open = "rb"), error = function(e) NULL)
  if (is.null(con)) return("")
  on.exit(try(close(con), silent = TRUE), add = TRUE)

  txt <- tryCatch(paste(readLines(con, warn = FALSE), collapse = ""), error = function(e) "")
  txt
}

opengwas_token <- ""
ensure_opengwas_token <- function() {
  if (nzchar(opengwas_token)) return(invisible(opengwas_token))
  token <- ezhu_read_opengwas_jwt()
  if (!nzchar(token)) {
    stop(
      "OpenGWAS JWT is not set. Write it to $HOME/.config/ezhu/opengwas_jwt (preferred) or set OPENGWAS_JWT, then rerun.",
      call. = FALSE
    )
  }
  opengwas_token <<- token
  Sys.setenv(OPENGWAS_JWT = token)
  invisible(token)
}

opengwas_sample_size <- function(id) {
  id <- trimws(as.character(id))
  if (!nzchar(id)) return(NA_real_)
  ensure_opengwas_token()

  info <- tryCatch(ieugwasr::gwasinfo(id), error = function(e) NULL)
  if (is.null(info) || nrow(info) == 0) return(NA_real_)
  names(info) <- tolower(names(info))

  for (key in c("sample_size", "samplesize", "n")) {
    if (key %in% names(info)) {
      n <- safe_numeric(info[[key]][1])
      if (is.finite(n) && n > 0) return(n)
    }
  }
  NA_real_
}

download_if_missing <- function(url, dest) {
  if (file.exists(dest)) return(invisible(dest))
  if (!nzchar(curl_bin)) stop("curl not found on PATH (required for VCF downloads).", call. = FALSE)
  ezhu_dir_create(dirname(dest))
  cmd <- c("-L", "-sS", "--fail", "--retry", "3", "--retry-delay", "2", "-o", dest, url)
  status <- system2(curl_bin, cmd)
  if (!identical(status, 0L) || !file.exists(dest)) stop("Failed to download: ", url, call. = FALSE)
  dest
}

opengwas_vcf_paths <- function(id) {
  id <- trimws(as.character(id))
  if (!nzchar(id)) stop("Empty OpenGWAS ID.", call. = FALSE)
  url <- paste0(base_url, "/", id, "/", id, ".vcf.gz")
  url_tbi <- paste0(url, ".tbi")
  local <- file.path("data/raw/opengwas/vcf", paste0(id, ".vcf.gz"))
  local_tbi <- paste0(local, ".tbi")
  list(url = url, url_tbi = url_tbi, path = local, path_tbi = local_tbi)
}

fetch_region <- function(vcf_path, chrom, start, end) {
  region <- paste0(chrom, ":", start, "-", end)
  v <- gwasvcf::query_gwas(vcf_path, chrompos = region)
  gr <- gwasvcf::vcf_to_granges(v)
  gr
}

extract_rsids <- function(gr) {
  if (is.null(gr) || length(gr) == 0) return(character())
  cand <- character()
  m <- tryCatch(S4Vectors::mcols(gr), error = function(e) NULL)
  if (!is.null(m)) {
    m_names <- names(m)
    for (key in c("rsid", "RSID", "SNP", "variant", "ID")) {
      if (key %in% m_names) cand <- c(cand, as.character(m[[key]]))
    }
  }
  cand <- c(cand, names(gr))
  cand <- unique(trimws(cand))
  cand <- cand[nzchar(cand)]
  cand <- cand[grepl("^rs[0-9]+$", cand)]
  cand
}

gr_to_table <- function(gr) {
  df <- as.data.frame(gr)
  if (!"seqnames" %in% names(df) || !"start" %in% names(df)) return(data.frame())

  # Common OpenGWAS VCF fields exposed by gwasvcf.
  es <- df$ES %||% df$es %||% NA_real_
  se <- df$SE %||% df$se %||% NA_real_
  af <- df$AF %||% df$af %||% NA_real_
  lp <- df$LP %||% df$lp %||% NA_real_
  ss <- df$SS %||% df$ss %||% NA_real_

  out <- data.frame(
    chromosome = as.character(df$seqnames),
    position = as.integer(df$start),
    beta = safe_numeric(es),
    se = safe_numeric(se),
    af = safe_numeric(af),
    lp = safe_numeric(lp),
    n = safe_numeric(ss),
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$beta) & is.finite(out$se) & is.finite(out$af) & out$af >= 0 & out$af <= 1, , drop = FALSE]
  if (nrow(out) == 0) return(out)

  out$chromosome <- sub("^chr", "", out$chromosome, ignore.case = TRUE)
  out$id <- paste0(out$chromosome, ":", out$position)
  out$maf <- pmin(out$af, 1 - out$af)
  out <- out[is.finite(out$maf) & out$maf > 0 & out$maf < 0.5, , drop = FALSE]

  # Remove multiallelic duplicates (multiple rows per id).
  dup <- ave(out$id, out$id, FUN = length)
  out <- out[dup == 1, , drop = FALSE]
  out
}

assoc_to_table <- function(df, default_n = NA_real_, maf_fallback = NULL) {
  if (is.null(df) || nrow(df) == 0) return(data.frame())
  names(df) <- tolower(names(df))

  get1 <- function(keys) {
    for (k in keys) {
      if (k %in% names(df)) return(df[[k]])
    }
    rep(NA, nrow(df))
  }

  chr <- as.character(get1(c("chr", "chromosome", "chrom")))
  pos <- safe_numeric(get1(c("position", "pos", "bp", "base_pair_location")))
  beta <- safe_numeric(get1(c("beta", "b", "es", "effect")))
  se <- safe_numeric(get1(c("se", "standard_error", "stderr")))
  eaf <- safe_numeric(get1(c("eaf", "effect_allele_frequency", "af")))
  n <- safe_numeric(get1(c("n", "samplesize", "ss")))
  if (all(!is.finite(n)) && is.finite(default_n)) n <- rep(default_n, length(beta))

  out <- data.frame(
    chromosome = sub("^chr", "", chr, ignore.case = TRUE),
    position = as.integer(pos),
    beta = beta,
    se = se,
    af = eaf,
    n = n,
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$beta) & is.finite(out$se) & is.finite(out$position) & nzchar(out$chromosome), , drop = FALSE]
  if (nrow(out) == 0) return(out)
  out$rsid <- as.character(get1(c("rsid", "snp", "variant", "id")))
  out$rsid <- trimws(out$rsid)
  out$rsid[!grepl("^rs[0-9]+$", out$rsid)] <- ""
  out$id <- paste0(out$chromosome, ":", out$position)

  out$maf <- ifelse(is.finite(out$af) & out$af >= 0 & out$af <= 1, pmin(out$af, 1 - out$af), NA_real_)
  if (!is.null(maf_fallback)) {
    maf_fallback <- maf_fallback[nzchar(names(maf_fallback))]
    if (length(maf_fallback) > 0) {
      use <- !is.finite(out$maf) | out$maf <= 0 | out$maf >= 0.5
      has_rsid <- nzchar(out$rsid) & out$rsid %in% names(maf_fallback)
      idx <- which(use & has_rsid)
      if (length(idx) > 0) {
        out$maf[idx] <- maf_fallback[out$rsid[idx]]
      }
    }
  }
  out <- out[is.finite(out$maf) & out$maf > 0 & out$maf < 0.5, , drop = FALSE]

  dup <- ave(out$id, out$id, FUN = length)
  out <- out[dup == 1, , drop = FALSE]
  out
}

run_one_coloc <- function(eqtl_tbl, gwas_tbl, p1, p2, p12) {
  if (isTRUE(dry_run)) return(NULL)
  merged <- merge(eqtl_tbl, gwas_tbl, by = "id", suffixes = c("_eqtl", "_gwas"))
  if (nrow(merged) < min_overlap) return(NULL)

  n1 <- suppressWarnings(stats::median(merged$n_eqtl, na.rm = TRUE))
  n2 <- suppressWarnings(stats::median(merged$n_gwas, na.rm = TRUE))
  if (!is.finite(n1) || !is.finite(n2)) return(NULL)

  eqtl_dataset <- list(
    beta = merged$beta_eqtl,
    varbeta = merged$se_eqtl^2,
    N = n1,
    MAF = merged$maf_eqtl,
    type = "quant",
    snp = merged$id
  )

  gwas_dataset <- list(
    beta = merged$beta_gwas,
    varbeta = merged$se_gwas^2,
    N = n2,
    MAF = merged$maf_gwas,
    type = "quant",
    snp = merged$id
  )

  coloc::coloc.abf(
    dataset1 = eqtl_dataset,
    dataset2 = gwas_dataset,
    p1 = p1,
    p2 = p2,
    p12 = p12
  )
}

ezhu_dir_create("results/causal/coloc")

attempts <- list()
summaries <- list()
for (i in seq_len(nrow(mr))) {
  exposure_id <- as.character(mr$exposure_id[i])
  outcome_id <- as.character(mr$outcome_id[i])
  gene <- as.character(mr$exposure_gene_symbol_input[i] %||% mr$exposure_gene_symbol_ensembl[i] %||% mr$exposure_ensg[i] %||% "")
  ensg <- extract_ensg(exposure_id, mr$exposure_ensg[i] %||% "")

  use_mode <- region_mode
  if (identical(use_mode, "auto")) use_mode <- if (isTRUE(has_instruments)) "instrument" else "gene"
  if (identical(use_mode, "instrument") && !isTRUE(has_instruments)) use_mode <- "gene"

  chr <- NA_character_
  pos <- NA_integer_
  rsid <- NA_character_
  region <- NULL

  if (identical(use_mode, "instrument")) {
    inst_sub <- inst[as.character(inst[["id.exposure"]]) == exposure_id, , drop = FALSE]
    inst_sub$pval.exposure <- safe_numeric(inst_sub$pval.exposure)
    inst_sub$chr.exposure <- safe_numeric(inst_sub$chr.exposure)
    inst_sub$pos.exposure <- safe_numeric(inst_sub$pos.exposure)
    inst_sub <- inst_sub[is.finite(inst_sub$pval.exposure) & is.finite(inst_sub$chr.exposure) & is.finite(inst_sub$pos.exposure), , drop = FALSE]
    if (nrow(inst_sub) == 0) use_mode <- "gene"
  }

  if (identical(use_mode, "instrument")) {
    inst_sub <- inst_sub[order(inst_sub$pval.exposure, decreasing = FALSE), , drop = FALSE]
    sentinel <- inst_sub[1, , drop = FALSE]
    chr <- as.character(as.integer(sentinel$chr.exposure[1]))
    pos <- as.integer(sentinel$pos.exposure[1])
    rsid <- as.character(sentinel$snp[1] %||% "")
    start <- max(1L, pos - window_bp)
    end <- pos + window_bp
    region <- list(chr = chr, start = start, end = end, center = pos, mode = "instrument")
  }

  if (identical(use_mode, "gene")) {
    gene_info <- fetch_gene_grch37(ensg, curl_bin = curl_bin)
    if (is.null(gene_info)) next
    chr <- as.character(gene_info$chr)
    start <- max(1L, as.integer(gene_info$start) - window_bp)
    end <- as.integer(gene_info$end) + window_bp
    pos <- as.integer(round((as.integer(gene_info$start) + as.integer(gene_info$end)) / 2))
    region <- list(chr = chr, start = start, end = end, center = pos, mode = "gene")
  }
  if (is.null(region)) next

  max_variants <- suppressWarnings(as.integer(Sys.getenv("EZHU_COLOC_MAX_VARIANTS", unset = "5000")))
  if (is.na(max_variants) || max_variants <= 0) max_variants <- 5000L

  # Default N values: prefer MR instruments (if available), otherwise query OpenGWAS metadata via API.
  default_n_exposure <- NA_real_
  default_n_outcome <- NA_real_
  if (isTRUE(has_instruments) && exists("inst_sub") && nrow(inst_sub) > 0) {
    if ("samplesize.exposure" %in% names(inst_sub)) default_n_exposure <- stats::median(safe_numeric(inst_sub$samplesize.exposure), na.rm = TRUE)
    if ("samplesize.outcome" %in% names(inst_sub)) default_n_outcome <- stats::median(safe_numeric(inst_sub$samplesize.outcome), na.rm = TRUE)
  }
  if (!is.finite(default_n_exposure)) default_n_exposure <- opengwas_sample_size(exposure_id)
  if (!is.finite(default_n_outcome)) default_n_outcome <- opengwas_sample_size(outcome_id)

  # API backend: build a stable rsID set by mapping OpenGWAS chrpos variants to 1000G rsIDs.
  # This avoids Ensembl rsID drift while keeping overlap sufficient for coloc.
  region_probe <- fetch_region_rsids_1000g(
    chr = region$chr,
    center = region$center,
    radius_bp = window_bp,
    start = region$start,
    end = region$end,
    max_positions = max_variants,
    pop = "EUR"
  )
  snps <- region_probe$rsids
  maf_fallback <- region_probe$maf
  if (length(snps) < min_overlap) {
    attempts[[length(attempts) + 1L]] <- data.frame(
      exposure_id = exposure_id,
      outcome_id = outcome_id,
      exposure_gene_symbol = gene,
      exposure_ensg = ensg,
      sentinel_rsid = rsid,
      region_hg37 = paste0(region$chr, ":", region$start, "-", region$end),
      region_center_chrpos = paste0(region$chr, ":", region$center),
      region_mode = region$mode,
      window_kb = window_kb,
      max_variants = max_variants,
      n_region_rsids = length(snps),
      n_variants_queried = 0L,
      n_assoc_exposure = 0L,
      n_assoc_outcome = 0L,
      n_eqtl_rows = 0L,
      n_gwas_rows = 0L,
      n_overlap = 0L,
      reason = paste0("region_rsids_lt_", min_overlap),
      stringsAsFactors = FALSE
    )
    next
  }

  chunk_size <- suppressWarnings(as.integer(Sys.getenv("EZHU_COLOC_API_CHUNK", unset = "250")))
  if (is.na(chunk_size) || chunk_size <= 0) chunk_size <- 250L

  opengwas_assoc_chunked <- function(id, variants) {
    out <- list()
    n <- length(variants)
    if (n == 0) return(NULL)
    chunks <- split(variants, ceiling(seq_len(n) / chunk_size))
    for (j in seq_along(chunks)) {
      v <- chunks[[j]]
      res <- NULL
      for (attempt in 1:4) {
        res <- tryCatch(
          ieugwasr::associations(variants = v, id = id),
          error = function(e) tryCatch(ieugwasr::associations(v, id), error = function(e2) NULL)
        )
        if (!is.null(res)) break
        Sys.sleep(2^attempt)
      }
      if (!is.null(res) && nrow(res) > 0) out[[length(out) + 1L]] <- res
    }
    if (length(out) == 0) return(NULL)
    dplyr::bind_rows(out)
  }

  assoc_exp <- opengwas_assoc_chunked(exposure_id, snps)
  assoc_out <- opengwas_assoc_chunked(outcome_id, snps)

  eqtl_tbl <- assoc_to_table(assoc_exp, default_n = default_n_exposure, maf_fallback = maf_fallback)
  gwas_tbl <- assoc_to_table(assoc_out, default_n = default_n_outcome, maf_fallback = maf_fallback)

  merged_inputs_all <- merge(eqtl_tbl, gwas_tbl, by = "id", suffixes = c("_eqtl", "_gwas"))
  overlap_n <- nrow(merged_inputs_all)
  reason <- ""
  if (nrow(eqtl_tbl) == 0) reason <- "no_exposure_associations_after_qc"
  if (nrow(gwas_tbl) == 0) reason <- if (nzchar(reason)) paste(reason, "no_outcome_associations_after_qc", sep = ";") else "no_outcome_associations_after_qc"
  if (!nzchar(reason) && overlap_n < min_overlap) reason <- paste0("overlap_lt_", min_overlap)
  if (is.null(assoc_exp) || nrow(assoc_exp) == 0) {
    reason <- if (nzchar(reason)) paste(reason, "exposure_api_no_rows", sep = ";") else "exposure_api_no_rows"
  } else if (nrow(assoc_exp) < min_overlap) {
    reason <- if (nzchar(reason)) paste(reason, paste0("exposure_api_rows_lt_", min_overlap), sep = ";") else paste0("exposure_api_rows_lt_", min_overlap)
  }
  if (is.null(assoc_out) || nrow(assoc_out) == 0) {
    reason <- if (nzchar(reason)) paste(reason, "outcome_api_no_rows", sep = ";") else "outcome_api_no_rows"
  } else if (nrow(assoc_out) < min_overlap) {
    reason <- if (nzchar(reason)) paste(reason, paste0("outcome_api_rows_lt_", min_overlap), sep = ";") else paste0("outcome_api_rows_lt_", min_overlap)
  }

  attempts[[length(attempts) + 1L]] <- data.frame(
    exposure_id = exposure_id,
    outcome_id = outcome_id,
    exposure_gene_symbol = gene,
    exposure_ensg = ensg,
    sentinel_rsid = rsid,
    region_hg37 = paste0(region$chr, ":", region$start, "-", region$end),
    region_center_chrpos = paste0(region$chr, ":", region$center),
    region_mode = region$mode,
    window_kb = window_kb,
    max_variants = max_variants,
    n_region_rsids = length(snps),
    n_variants_queried = length(snps),
    n_assoc_exposure = if (is.null(assoc_exp)) 0L else nrow(assoc_exp),
    n_assoc_outcome = if (is.null(assoc_out)) 0L else nrow(assoc_out),
    n_eqtl_rows = nrow(eqtl_tbl),
    n_gwas_rows = nrow(gwas_tbl),
    n_overlap = overlap_n,
    reason = reason,
    stringsAsFactors = FALSE
  )

  coloc_res <- run_one_coloc(eqtl_tbl, gwas_tbl, p1 = p1, p2 = p2, p12 = p12)
  if (is.null(coloc_res)) next

  merged_inputs <- merged_inputs_all
  merged_inputs$chromosome <- merged_inputs$chromosome_eqtl
  merged_inputs$position <- merged_inputs$position_eqtl
  merged_inputs <- merged_inputs[, c(
    "id", "chromosome", "position",
    "beta_eqtl", "se_eqtl", "maf_eqtl", "n_eqtl",
    "beta_gwas", "se_gwas", "maf_gwas", "n_gwas"
  )]

  slug <- paste0("chr", chr, "_", pos)
  locus_prefix <- paste0(gsub("[^A-Za-z0-9._-]+", "_", exposure_id), "__", gsub("[^A-Za-z0-9._-]+", "_", outcome_id), "__", slug)
  locus_path <- file.path("results/causal/coloc", paste0("locus_", locus_prefix, ".tsv.gz"))
  gz <- gzfile(locus_path, open = "wt")
  utils::write.table(merged_inputs, file = gz, sep = "\t", quote = FALSE, row.names = FALSE)
  close(gz)

  s <- as.data.frame(t(as.data.frame(coloc_res$summary)))
  s$exposure_id <- exposure_id
  s$outcome_id <- outcome_id
  s$exposure_gene_symbol <- gene
  s$exposure_ensg <- ensg
  s$sentinel_rsid <- rsid
  s$region_hg37 <- paste0(region$chr, ":", region$start, "-", region$end)
  s$region_center_chrpos <- paste0(region$chr, ":", region$center)
  s$region_mode <- region$mode
  s$window_kb <- window_kb
  s$n_snps_used <- nrow(merged_inputs)
  s$locus_table <- locus_path
  summaries[[length(summaries) + 1L]] <- s
}

attempts_path <- "results/causal/coloc/coloc_attempts.tsv"
attempts_tbl <- if (length(attempts) > 0) do.call(rbind, attempts) else data.frame()
if (nrow(attempts_tbl) > 0) {
  utils::write.table(attempts_tbl, file = attempts_path, sep = "\t", quote = FALSE, row.names = FALSE)
}

out_path <- "results/causal/coloc_summary.tsv"

if (length(summaries) == 0) {
  empty_cols <- c(
    "nsnps",
    "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf",
    "exposure_id", "outcome_id", "exposure_gene_symbol", "exposure_ensg",
    "sentinel_rsid", "region_hg37", "region_center_chrpos", "region_mode",
    "window_kb", "n_snps_used", "locus_table"
  )
  empty_summary <- as.data.frame(setNames(replicate(length(empty_cols), character(0), simplify = FALSE), empty_cols))
  utils::write.table(empty_summary, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)

  msg <- paste0(
    "No colocalization results produced (0 loci passed overlap/N checks). ",
    "Wrote: ", out_path,
    if (nrow(attempts_tbl) > 0) paste0(" and ", attempts_path) else "",
    ". Most common cause is insufficient SNP overlap after QC (see `n_overlap` and `reason` in the attempts table)."
  )
  ezhu_write_stage_metadata(
    "06b_colocalization",
    params = list(
      mr_validation = mr_validation_path,
      mr_instruments = mr_instruments_path,
      top_hits = top_hits,
      window_kb = window_kb,
      min_overlap = min_overlap,
      strict_no_results = strict_no_results,
      dry_run = dry_run,
      priors = list(p1 = p1, p2 = p2, p12 = p12),
      coloc_backend = if (isTRUE(use_vcf)) "vcf+api" else "api",
      opengwas_files_base = base_url,
      output_summary = out_path,
      output_attempts = attempts_path
    ),
    seed = seed
  )
  if (isTRUE(strict_no_results)) stop(msg, call. = FALSE)
  warning(msg, call. = FALSE)
  message("Wrote: ", out_path, " (n=0)")
  quit(status = 0)
}

summary_tbl <- do.call(rbind, summaries)
summary_tbl$PP.H4.abf <- safe_numeric(summary_tbl$PP.H4.abf %||% NA_real_)
summary_tbl <- summary_tbl[order(summary_tbl$PP.H4.abf, decreasing = TRUE), , drop = FALSE]

utils::write.table(summary_tbl, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)

ezhu_write_stage_metadata(
  "06b_colocalization",
  params = list(
    mr_validation = mr_validation_path,
    mr_instruments = mr_instruments_path,
    top_hits = top_hits,
    window_kb = window_kb,
    min_overlap = min_overlap,
    strict_no_results = strict_no_results,
    dry_run = dry_run,
    priors = list(p1 = p1, p2 = p2, p12 = p12),
    coloc_backend = if (isTRUE(use_vcf)) "vcf+api" else "api",
    opengwas_files_base = base_url,
    output_summary = out_path,
    output_attempts = attempts_path
  ),
  seed = seed
)

message("Wrote: ", out_path, " (n=", nrow(summary_tbl), ")")
