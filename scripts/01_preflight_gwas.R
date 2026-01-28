args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
if (length(file_arg) == 1) {
  script_path <- sub("^--file=", "", file_arg)
  repo_root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)
  if (dir.exists(repo_root)) setwd(repo_root)
}

source("scripts/utils/repro.R")
source("scripts/utils/manifest.R")
ezhu_set_repo_root()
ezhu_activate_renv()

seed_env <- Sys.getenv("EZHU_SEED", unset = "")
seed <- suppressWarnings(as.integer(seed_env))
if (is.na(seed)) seed <- 20251227L
set.seed(seed)

ezhu_dir_create("results/metadata")

max_variants_env <- Sys.getenv("EZHU_GWAS_PREFLIGHT_MAX_VARIANTS", unset = "")
max_variants <- suppressWarnings(as.integer(max_variants_env))
if (!is.finite(max_variants) || is.na(max_variants) || max_variants <= 0) max_variants <- NA_integer_

fill_missing_rsid <- identical(Sys.getenv("EZHU_PREFLIGHT_FILL_MISSING_RSID", unset = "0"), "1")

manifest <- ezhu_read_manifest("data/manifest.tsv")
active <- audit_active_manifest_rows <- function(df) {
  fields <- c("type", "name", "source", "accession", "url", "local_path")
  fields <- intersect(fields, names(df))
  if (length(fields) == 0) return(rep(FALSE, nrow(df)))
  apply(df[, fields, drop = FALSE], 1, function(row) any(!is.na(row) & nzchar(trimws(as.character(row)))))
}
active_rows <- active(manifest)
gwas_rows <- which(active_rows & tolower(manifest$type) == "gwas")

if (length(gwas_rows) == 0) {
  message("No GWAS rows found in data/manifest.tsv (type=gwas). Preflight skipped.")
  ezhu_write_stage_metadata("01_preflight_gwas", params = list(skipped = TRUE, reason = "no gwas rows"), seed = seed)
  quit(status = 0)
}

select_id <- Sys.getenv("EZHU_GWAS_MANIFEST_ID", unset = Sys.getenv("EZHU_GWAS_ID", unset = ""))
select_id <- trimws(as.character(select_id))
if (nzchar(select_id)) {
  matched <- which(active_rows & tolower(manifest$type) == "gwas" & as.character(manifest$id) == select_id)
  if (length(matched) == 0) {
    stop("EZHU_GWAS_MANIFEST_ID did not match any active GWAS row in data/manifest.tsv: ", select_id, call. = FALSE)
  }
  if (length(matched) > 1) {
    stop("EZHU_GWAS_MANIFEST_ID matched multiple rows; manifest IDs must be unique: ", select_id, call. = FALSE)
  }
  row <- manifest[matched, , drop = FALSE]
} else {
  if (length(gwas_rows) > 1) {
    message("Multiple GWAS rows found; using the first one. Set EZHU_GWAS_MANIFEST_ID to choose a specific GWAS.")
  }
  row <- manifest[gwas_rows[1], , drop = FALSE]
}
notes <- ezhu_parse_notes_kv(row$notes)
build <- tolower(trimws(notes$build %||% ""))
ancestry_raw <- notes$ancestry %||% notes$cestry %||% ""
ancestry <- toupper(trimws(ancestry_raw))

if (!build %in% c("hg37", "hg38")) {
  stop("GWAS preflight requires notes to include build=hg37|hg38 for the GWAS row.", call. = FALSE)
}
if (!nzchar(ancestry)) {
  stop("GWAS preflight requires notes to include ancestry=... (e.g., EUR/EAS) for the GWAS row.", call. = FALSE)
}

qc_path <- Sys.getenv("EZHU_GWAS_PREFLIGHT_QC_OUT", unset = "")
qc_path <- trimws(as.character(qc_path))
if (!nzchar(qc_path)) {
  suffix <- trimws(as.character(row$id %||% select_id))
  qc_path <- if (nzchar(suffix)) {
    file.path("results/metadata", paste0("01_preflight_gwas__qc__", suffix, ".tsv"))
  } else {
    file.path("results/metadata", "01_preflight_gwas__qc.tsv")
  }
}

write_preflight_qc <- function(path, kv) {
  ezhu_dir_create(dirname(path))
  keys <- names(kv)
  if (is.null(keys) || length(keys) == 0) return(invisible(FALSE))
  vals <- unname(unlist(kv))
  lines <- paste0(keys, "\t", vals)
  writeLines(lines, con = path, sep = "\n")
  invisible(TRUE)
}

gwas_path <- trimws(as.character(row$local_path))
if (!nzchar(gwas_path) || !file.exists(gwas_path)) {
  stop("GWAS local_path does not exist on disk: ", gwas_path, call. = FALSE)
}

message("Reading GWAS: ", gwas_path)
is_gz <- grepl("\\.gz$", gwas_path, ignore.case = TRUE)
use_shell_stream <- nzchar(Sys.which("bash")) &&
  nzchar(Sys.which("awk")) &&
  (isTRUE(!is_gz) || nzchar(Sys.which("gzip"))) &&
  isTRUE(file.info(gwas_path)$size > 50 * 1024^2)

if (isTRUE(use_shell_stream)) {
  header_con <- if (isTRUE(is_gz)) gzfile(gwas_path, open = "rt") else file(gwas_path, open = "rt")
  header_line <- readLines(header_con, n = 1)
  close(header_con)
  header_line <- trimws(header_line)
  has_tab <- grepl("\t", header_line, fixed = TRUE)
  header_cols <- if (isTRUE(has_tab)) {
    strsplit(header_line, "\t", fixed = TRUE)[[1]]
  } else {
    strsplit(header_line, "[[:space:]]+")[[1]]
  }
  header_lc <- tolower(header_cols)

  pick_idx <- function(candidates) {
    hit <- match(candidates, header_lc)
    hit <- hit[!is.na(hit)]
    if (length(hit) == 0) return(NA_integer_)
    hit[1]
  }

  chr_i <- pick_idx(c("chrom", "chr", "chromosome"))
  pos_i <- pick_idx(c("pos", "bp", "position", "base_pair_location", "pos_b37", "pos_b38"))
  rsid_i <- pick_idx(c("rsid", "rs_id", "snp", "rs", "marker", "markername", "variant_id"))
  se_i <- pick_idx(c("se", "stderr", "standard_error", "or_se", "se_or"))
  beta_i <- pick_idx(c("beta", "effect", "b"))
  or_i <- pick_idx(c("or", "odds_ratio"))
  af_i <- pick_idx(c("maf", "eaf", "effect_allele_frequency", "af", "freq1", "freqa1", "freq_a1", "eaf_a1"))

  fast_ok <- !any(is.na(c(chr_i, pos_i, rsid_i, se_i, af_i))) && (!is.na(beta_i) || !is.na(or_i))

  if (isTRUE(fast_ok)) {
    out_dir <- "data/processed"
    ezhu_dir_create(out_dir)
    out_path <- Sys.getenv("EZHU_GWAS_HARMONIZED_OUT", unset = "")
    out_path <- trimws(as.character(out_path))
    if (!nzchar(out_path)) {
      # Use the ID from the matched row, or fallback to the provided select_id
      suffix <- trimws(as.character(row$id %||% select_id))
      out_path <- if (nzchar(suffix)) {
        file.path(out_dir, paste0("gwas_harmonized__", suffix, ".tsv.gz"))
      } else {
        file.path(out_dir, "gwas_harmonized.tsv.gz")
      }
    }
    out_tmp <- paste0(out_path, ".tmp")

    beta_i_awk <- if (is.na(beta_i)) 0L else beta_i
    or_i_awk <- if (is.na(or_i)) 0L else or_i
    fs_awk <- if (isTRUE(has_tab)) "\\t" else "[[:space:]]+"

  awk_program <- paste(
      "BEGIN{OFS=\"\\t\";",
      "print \"chrom\",\"pos\",\"rsid\",\"se\",\"beta\",\"maf\";",
      "rows_in=0;rows_out=0;skip_missing=0;skip_rsid=0;skip_se=0;skip_af=0;skip_beta=0;fill_rsid=0}",
      "NR==1{next}",
      "{rows_in++;",
      "chr=$(chr_i); gsub(/^chr/ , \"\", chr); if (chr ~ /^[0-9]+$/) { sub(/^0+/, \"\", chr) }",
      "pos=$(pos_i);",
      "rsid=$(rsid_i);",
      "se_str=$(se_i); se=se_str+0;",
      "af_str=$(af_i); af=af_str+0;",
      "if(chr==\"\"||pos==\"\"||pos==\"NA\"){skip_missing++; next}",
      "if(rsid==\"\"||rsid==\"NA\"||rsid==\".\"){if(fill_missing_rsid==1){rsid=chr \":\" pos; fill_rsid++} else {skip_rsid++; next}}",
      "if(se_str==\"\"||se_str==\"NA\"||se<=0){skip_se++; next}",
      "if(af_str==\"\"||af_str==\"NA\"||af<=0||af>=1){skip_af++; next}",
      "maf=(af<0.5?af:1-af);",
      "if(maf<=0||maf>0.5){skip_af++; next}",
      if (beta_i_awk > 0) "beta_str=$(beta_i); if(beta_str==\"\"||beta_str==\"NA\"){skip_beta++; next}; beta=beta_str+0;" else "or_str=$(or_i); orv=or_str+0; if(or_str==\"\"||or_str==\"NA\"||orv<=0){skip_beta++; next}; beta=log(orv);",
      "rows_out++; print chr,pos,rsid,se,beta,maf;",
      if (!is.na(max_variants)) paste0("if(rows_out>=", max_variants, "){exit}"), "}",
      "END{print \"rows_in\\t\" rows_in > qc;",
      "print \"rows_out\\t\" rows_out >> qc;",
      "print \"skip_missing\\t\" skip_missing >> qc;",
      "print \"skip_rsid\\t\" skip_rsid >> qc;",
      "print \"fill_rsid\\t\" fill_rsid >> qc;",
      "print \"skip_se\\t\" skip_se >> qc;",
      "print \"skip_af\\t\" skip_af >> qc;",
      "print \"skip_beta\\t\" skip_beta >> qc}",
      sep = " "
    )

    stream_cmd <- if (isTRUE(is_gz)) {
      paste0("gzip -dc ", shQuote(gwas_path))
    } else {
      paste0("cat ", shQuote(gwas_path))
    }

    cmd <- paste0(
      "set -euo pipefail; ",
      stream_cmd, " | ",
      "awk -F'", fs_awk, "' -v OFS='\\t' ",
      "-v chr_i=", chr_i, " -v pos_i=", pos_i, " -v rsid_i=", rsid_i,
      " -v se_i=", se_i, " -v af_i=", af_i,
      " -v beta_i=", beta_i_awk, " -v or_i=", or_i_awk,
      " -v fill_missing_rsid=", if (isTRUE(fill_missing_rsid)) 1L else 0L,
      " -v qc=", shQuote(qc_path), " '", awk_program, "' | ",
      "gzip -c > ", shQuote(out_tmp), "; ",
      "mv ", shQuote(out_tmp), " ", shQuote(out_path)
    )

    message("Using shell streaming for GWAS harmonization (fast path).")
    system2("bash", c("-lc", cmd))
    if (!file.exists(out_path)) stop("Fast-path output not found: ", out_path, call. = FALSE)

    shell_qc <- NULL
    shell_qc <- tryCatch(
      {
        q <- utils::read.delim(qc_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
        as.list(setNames(as.numeric(q[[2]]), q[[1]]))
      },
      error = function(e) NULL
    )

    param_summary <- list(
      gwas_row_id = as.character(row$id),
      gwas_source = as.character(row$source),
      gwas_name = as.character(row$name),
      gwas_local_path = gwas_path,
      build = build,
      ancestry = ancestry,
      detected_columns = list(
        chrom = header_cols[chr_i],
        pos = header_cols[pos_i],
        rsid = header_cols[rsid_i],
        se = header_cols[se_i],
        beta = if (beta_i_awk > 0) header_cols[beta_i_awk] else NA_character_,
        or = if (or_i_awk > 0) header_cols[or_i_awk] else NA_character_,
        af = header_cols[af_i]
      ),
      qc = shell_qc %||% list(),
      harmonized_rows = as.numeric(shell_qc$rows_out %||% NA_real_),
      harmonized_output = out_path,
      preflight_mode = "shell_stream"
    )

    ezhu_write_stage_metadata("01_preflight_gwas", params = param_summary, seed = seed)
    message("GWAS preflight OK (fast path). Harmonized output: ", out_path)
    quit(status = 0)
  } else {
    message("Shell streaming not applicable (missing required columns); falling back to R read.")
  }
}

detect_sep <- function(path, is_gz) {
  con <- if (isTRUE(is_gz)) gzfile(path, open = "rt") else file(path, open = "rt")
  line <- readLines(con, n = 1)
  close(con)
  if (length(line) == 0) return("\t")
  if (grepl("\t", line, fixed = TRUE)) "\t" else ""
}
sep_detected <- detect_sep(gwas_path, is_gz)

gwas_con <- NULL
gwas <- tryCatch(
  {
    if (grepl("\\.gz$", gwas_path, ignore.case = TRUE)) {
      gwas_con <<- gzfile(gwas_path, open = "rt")
      utils::read.delim(
        gwas_con,
        sep = sep_detected,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    } else {
      utils::read.delim(
        gwas_path,
        sep = sep_detected,
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  },
  finally = {
    if (!is.null(gwas_con)) close(gwas_con)
  }
)

original_names <- names(gwas)
names(gwas) <- tolower(names(gwas))

pick_col <- function(candidates) {
  hit <- intersect(candidates, names(gwas))
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

col_chr <- pick_col(c("chrom", "chr", "chromosome"))
col_pos <- pick_col(c("pos", "bp", "position", "base_pair_location", "pos_b37", "pos_b38"))
col_rsid <- pick_col(c("rsid", "rs_id", "snp", "rs", "marker", "markername", "variant_id"))
col_beta <- pick_col(c("beta", "effect", "b"))
col_or <- pick_col(c("or", "odds_ratio"))
col_se <- pick_col(c("se", "stderr", "standard_error", "or_se", "se_or"))
col_maf <- pick_col(c("maf", "eaf", "effect_allele_frequency", "af", "freq1", "freqa1", "freq_a1", "eaf_a1"))

problems <- character()
if (is.na(col_chr)) problems <- c(problems, "missing chromosome column (expected one of: chrom/chr/chromosome)")
if (is.na(col_pos)) problems <- c(problems, "missing position column (expected one of: pos/bp/position/base_pair_location)")
if (is.na(col_rsid)) problems <- c(problems, "missing rsid column (expected one of: rsid/snp/rs/marker/markername/variant_id)")
if (is.na(col_se)) problems <- c(problems, "missing standard error column (expected one of: se/stderr/standard_error/or_se/se_or)")
if (is.na(col_maf)) problems <- c(problems, "missing allele frequency column (expected one of: maf/eaf/af/freq1/freqa1/effect_allele_frequency)")

if (is.na(col_beta) && is.na(col_or)) {
  problems <- c(problems, "missing effect column (need beta or OR)")
}

if (length(problems) > 0) {
  stop("GWAS format issues:\n- ", paste(problems, collapse = "\n- "), call. = FALSE)
}

chr <- as.character(gwas[[col_chr]])
chr <- sub("^chr", "", chr, ignore.case = TRUE)
is_numeric_chr <- grepl("^[0-9]+$", chr)
chr[is_numeric_chr] <- sub("^0+", "", chr[is_numeric_chr])
chr[is_numeric_chr & !nzchar(chr)] <- NA_character_
pos <- suppressWarnings(as.integer(gwas[[col_pos]]))
rsid <- as.character(gwas[[col_rsid]])
se <- suppressWarnings(as.numeric(gwas[[col_se]]))
af <- suppressWarnings(as.numeric(gwas[[col_maf]]))

if (is.na(col_beta)) {
  or <- suppressWarnings(as.numeric(gwas[[col_or]]))
  beta <- log(or)
} else {
  beta <- suppressWarnings(as.numeric(gwas[[col_beta]]))
}

if (any(!is.na(af) & (af <= 0 | af >= 1))) {
  stop("GWAS QC failed: allele frequency must be in (0, 1) (found values outside range).", call. = FALSE)
}
maf <- pmin(af, 1 - af)

missing_rsid <- is.na(rsid) | !nzchar(trimws(rsid)) | trimws(rsid) %in% c("NA", ".")
filled_rsid <- 0L
if (any(missing_rsid)) {
  if (isTRUE(fill_missing_rsid)) {
    rsid[missing_rsid] <- paste0(chr[missing_rsid], ":", pos[missing_rsid])
    filled_rsid <- sum(missing_rsid)
  } else {
    # Drop rows with missing rsID unless explicitly allowed to fill placeholders.
    rsid[missing_rsid] <- NA_character_
  }
}

qc <- list(
  rows = nrow(gwas),
  missing_chr = sum(is.na(chr) | !nzchar(chr)),
  missing_pos = sum(is.na(pos)),
  missing_rsid = sum(is.na(gwas[[col_rsid]]) | !nzchar(as.character(gwas[[col_rsid]]))),
  filled_rsid = filled_rsid,
  missing_se = sum(is.na(se)),
  missing_maf = sum(is.na(af)),
  se_nonpos = sum(!is.na(se) & se <= 0),
  maf_out_of_range = sum(!is.na(maf) & (maf <= 0 | maf > 0.5)),
  beta_missing = sum(is.na(beta))
)

if (qc$se_nonpos > 0) stop("GWAS QC failed: se must be > 0 for all rows (found ", qc$se_nonpos, " violations).", call. = FALSE)
if (qc$maf_out_of_range > 0) stop("GWAS QC failed: maf must be in (0, 0.5] (found ", qc$maf_out_of_range, " violations).", call. = FALSE)

harm <- data.frame(
  chrom = chr,
  pos = pos,
  rsid = rsid,
  se = se,
  beta = beta,
  maf = maf,
  stringsAsFactors = FALSE
)

harm <- harm[complete.cases(harm), , drop = FALSE]

out_dir <- "data/processed"
ezhu_dir_create(out_dir)
out_path <- Sys.getenv("EZHU_GWAS_HARMONIZED_OUT", unset = "")
out_path <- trimws(as.character(out_path))
if (!nzchar(out_path)) {
  suffix <- trimws(as.character(row$id))
  out_path <- if (nzchar(suffix)) {
    file.path(out_dir, paste0("gwas_harmonized__", suffix, ".tsv.gz"))
  } else {
    file.path(out_dir, "gwas_harmonized.tsv.gz")
  }
}
gz <- gzfile(out_path, open = "wt")
utils::write.table(harm, file = gz, sep = "\t", quote = FALSE, row.names = FALSE)
close(gz)

write_preflight_qc(
  qc_path,
  list(
    rows_in = qc$rows,
    rows_out = nrow(harm),
    missing_chr = qc$missing_chr,
    missing_pos = qc$missing_pos,
    missing_rsid = qc$missing_rsid,
    filled_rsid = qc$filled_rsid,
    missing_se = qc$missing_se,
    missing_maf = qc$missing_maf,
    se_nonpos = qc$se_nonpos,
    maf_out_of_range = qc$maf_out_of_range,
    beta_missing = qc$beta_missing
  )
)

param_summary <- list(
  gwas_row_id = as.character(row$id),
  gwas_source = as.character(row$source),
  gwas_name = as.character(row$name),
  gwas_local_path = gwas_path,
  build = build,
  ancestry = ancestry,
  detected_columns = list(
    chrom = col_chr, pos = col_pos, rsid = col_rsid, se = col_se, beta = col_beta, or = col_or, maf = col_maf
  ),
  qc = qc,
  harmonized_rows = nrow(harm),
  harmonized_output = out_path,
  qc_tsv = qc_path
)

ezhu_write_stage_metadata("01_preflight_gwas", params = param_summary, seed = seed)
if (isTRUE(qc$filled_rsid > 0)) {
  message(
    "NOTE: Filled ", qc$filled_rsid, " missing rsid values using chrom:pos placeholders.\n",
    "      If scPagwas later reports empty results, consider using a GWAS source with true rsIDs."
  )
}
message("GWAS preflight OK. Harmonized output: ", out_path)
